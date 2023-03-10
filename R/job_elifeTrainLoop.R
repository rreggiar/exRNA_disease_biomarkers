library(tidyverse)
library(glmnet)
library(pROC)

import::from(.directory = r_code.dir,
             .from = 'helper_analysisFunctions.R',
             calc_entropy, run.de.seq.split)

inverse.logit.calc <- function(y_hat) { exp(y_hat)/(1+exp(y_hat)) }

doMC::registerDoMC(cores = 24)

count_out_list <-
  setNames(
    list(
      elife_rmsk_quant$liver_healthy,
      elife_rmsk_quant$lung_healthy,
      elife_rmsk_quant$esoph_healthy,
      elife_rmsk_quant$colorectal_healthy,
      elife_rmsk_quant$stomach_healthy
    ),
    c('Liver',
      'Lung',
      'Esophagus',
      'Colorectal',
      'Stomach')
  )

imap(count_out_list, function(analysis_set, name) {
  message(name)
  
  calc_entropy(analysis_set$deseq$norm_counts.df,
               elife_meta_data,
               grouping_level = 'clade') %>% 
    mutate(clade = ifelse(clade %in% names(te_color_pal), clade, 'other')) %>%
    filter(!grepl('GENCODE', clade)) %>%
    group_by(sample, clade) %>%
    summarize(entropy = mean(entropy)) %>%
    pivot_wider(names_from = 'clade', values_from = 'entropy') %>%
    select_if( ~ !any(is.na(.))) %>%
    column_to_rownames('sample') %>% as.matrix() -> binom_model_input_data_entropy
  
  cbind(binom_model_input_data_entropy,
        t(analysis_set$deseq$norm_counts.df[!rownames(analysis_set$deseq$norm_counts.df) %in%
                                      reference_meta_data$ucsc_rmsk_insert_info.df$gene,])) ->
    binom_model_input_data_entropy_plus
  
  
  binom_model_input_data_gencode <-
    t(analysis_set$deseq$norm_counts.df[!rownames(analysis_set$deseq$norm_counts.df) %in%
                                  reference_meta_data$ucsc_rmsk_insert_info.df$gene,])
  
  binom_model_input_data_te <-
    t(analysis_set$deseq$norm_counts.df[rownames(analysis_set$deseq$norm_counts.df) %in%
                                  reference_meta_data$ucsc_rmsk_insert_info.df$gene,])
  
  binom_model_input_data_all <-
    t(analysis_set$deseq$norm_counts.df)
  
  binom_model_input_meta_data <-
    analysis_set$quant_meta %>% 
    mutate(diagnosis = ifelse(diagnosis == 'Healthy donor', 0, 1))
  
  ## calculate split for this sample
  set.seed(2022 - 03 - 08)
  trainSplitIndex <-
    caret::createDataPartition(
      factor(binom_model_input_meta_data$diagnosis,
             levels = c(0, 1)),
      p = .8,
      list = FALSE,
      times = 1
    )
  
  # extract metaData
  binom_elife_y_in <-
    binom_model_input_meta_data

  # extract train & test splits
  binom_elife_y_in_train <-
    binom_elife_y_in[trainSplitIndex, , drop = F] %>% pull(diagnosis)
  
  # create folds for this cancer
  fold_ids <-
    caret::createFolds(factor(binom_elife_y_in_train), 10, list = F)
  
  binom_elife_y_in_test <-
    binom_elife_y_in[-trainSplitIndex, , drop = F] %>% pull(diagnosis)
  # build train-specific DE names
  
  print(binom_elife_y_in[trainSplitIndex, , drop = F])
  
  binom_model_input_data_de_names <-
    run.de.seq.split(
      input_se = analysis_set,
      top_level = paste0(name, ' cancer'),
      split = binom_elife_y_in[trainSplitIndex, , drop = F]
    )
  
  # build train-specific DE feature sets
  binom_model_input_data_te_de <-
    binom_model_input_data_te[, colnames(binom_model_input_data_te) %in% 
                                binom_model_input_data_de_names]
  
  binom_model_input_data_gencode_de <-
    binom_model_input_data_gencode[, colnames(binom_model_input_data_gencode) %in% 
                                binom_model_input_data_de_names]
  
  binom_model_input_data_all_de <-
    binom_model_input_data_all[, colnames(binom_model_input_data_all) %in% 
                                     binom_model_input_data_de_names]
  
  ###
  
  feature_set_list <- setNames(
    list(
      binom_model_input_data_entropy,
      binom_model_input_data_entropy_plus,
      binom_model_input_data_gencode,
      binom_model_input_data_te,
      binom_model_input_data_all,
      binom_model_input_data_all_de,
      binom_model_input_data_te_de,
      binom_model_input_data_gencode_de
    ),
    c('entropy',
      'entropy_plus',
      'gencode',
      'te',
      'all',
      'all_de',
      'te_de',
      'gencode_de')
  )
  
  imap(feature_set_list, function(feature_set, feature_set_name) {
    message(feature_set_name)
    
    binom_elife_x_in <-
      feature_set[rownames(feature_set) %in% binom_model_input_meta_data$sample,
                       colSums(feature_set != 0) > 0]
    
    print(dim(binom_elife_x_in))
    print(dim(binom_model_input_meta_data))
    
    binom_elife_x_in <-
      feature_set[match(rownames(binom_elife_x_in),
                        binom_model_input_meta_data$sample),]
    
    binom_elife_x_in_train <-
      binom_elife_x_in[trainSplitIndex, , drop = F]
    
    binom_elife_x_in_test <-
      binom_elife_x_in[-trainSplitIndex, , drop = F]
    
    lapply(seq(0, 1, 0.1), function(alpha) {
      message(alpha)
      
      binom_elife_model <-
        glmnet::cv.glmnet(
          x = binom_elife_x_in_train,
          y = binom_elife_y_in_train,
          family = 'binomial',
          type.measure = 'deviance',
          nfolds = 10,
          alpha = alpha,
          foldid = fold_ids,
          parallel = TRUE,
          keep = TRUE
        )
      
      y_hat_matrix <-
        binom_elife_model$fit.preval[, which(binom_elife_model$lambda == binom_elife_model$lambda.min)]
      
      predictions <-
        ROCR::prediction(
          labels = binom_elife_y_in_train,
          predictions = y_hat_matrix,
          label.ordering = c(0, 1)
        )
      
      performances <-
        ROCR::performance(
          prediction.obj = predictions,
          measure = "tpr",
          x.measure = "fpr"
        )
      
      ROCR::performance(prediction.obj = predictions,
                        measure = "auc")@y.values -> auc
      
      cutoffs <-
        data.frame(
          cut = performances@alpha.values[[1]],
          fpr = performances@x.values[[1]],
          tpr = performances@y.values[[1]]
        )
      
      perf_list <-
        cutoffs[which.min(abs(cutoffs$fpr - (1 - 0.9))),]
      perf_list$cut <-  inverse.logit.calc(perf_list$cut)
      perf_list$auc <-  auc
      perf_list$alpha <- alpha
      perf_list$feature_set <- feature_set_name
      
      list(
        'perf_list' = perf_list,
        'model' = binom_elife_model,
        'test_x' = binom_elife_x_in_test,
        'test_y' = binom_elife_y_in_test,
        'train_x' = binom_elife_x_in_train,
        'train_y' = binom_elife_y_in_train
      )
      
    })
    
  })
  
}) -> final_alpha_new
