library(tidyverse)
library(glmnet)
library(pROC)

import::from(.directory = r_code.dir,
             .from = 'helper_analysisFunctions.R',
             calc_entropy)

doMC::registerDoMC(cores = 12)

count_out_list <-
  setNames(
    list(
      elife_rmsk_quant$liver_healthy$deseq$norm_counts.df,
      elife_rmsk_quant$lung_healthy$deseq$norm_counts.df,
      elife_rmsk_quant$esoph_healthy$deseq$norm_counts.df,
      elife_rmsk_quant$colorectal_healthy$deseq$norm_counts.df,
      elife_rmsk_quant$stomach_healthy$deseq$norm_counts.df
    ),
    c('liver',
      'lung',
      'esoph',
      'colorectal',
      'stomach')
  )

imap(count_out_list, function(norm_counts.df, name) {
  message(name)
  
  calc_entropy(norm_counts.df,
               elife_meta_data,
               grouping_level = 'clade') %>%
    mutate(clade = ifelse(clade %in% names(te_color_pal), clade, 'other')) %>%
    filter(!grepl('GENCODE', clade)) %>%
    group_by(sample, clade) %>%
    summarize(entropy = mean(entropy)) %>%
    pivot_wider(names_from = 'clade', values_from = 'entropy') %>%
    select_if(~ !any(is.na(.))) %>%
    column_to_rownames('sample') %>% as.matrix() -> binom_model_input_data_entropy
  
  cbind(binom_model_input_data_entropy,
        t(norm_counts.df[!rownames(norm_counts.df) %in%
                           reference_meta_data$ucsc_rmsk_insert_info.df$gene, ])) ->
    binom_model_input_data_entropy_plus
  
  
  binom_model_input_data_gencode <-
    t(norm_counts.df[!rownames(norm_counts.df) %in%
                       reference_meta_data$ucsc_rmsk_insert_info.df$gene, ])
  
  binom_model_input_data_te <-
    t(norm_counts.df[rownames(norm_counts.df) %in%
                       reference_meta_data$ucsc_rmsk_insert_info.df$gene, ])
  
  binom_model_input_data_all <-
    t(norm_counts.df)
  
  feature_set_list <- setNames(
    list(
      binom_model_input_data_entropy,
      binom_model_input_data_entropy_plus,
      binom_model_input_data_gencode,
      binom_model_input_data_te,
      binom_model_input_data_all
    ),
    c('entropy',
      'entropy_plus',
      'gencode',
      'te',
      'all')
  )
  
  imap(feature_set_list, function(feature_set, feature_set_name) {
    message(feature_set_name)
    
    binom_model_input_meta_data <-
      elife_meta_data[match(rownames(feature_set), elife_meta_data$sample), ] %>%
      mutate(diagnosis = ifelse(diagnosis == 'Healthy donor', 0, 1))
    
    set.seed(2022 - 03 - 08)
    trainSplitIndex <-
      caret::createDataPartition(
        factor(binom_model_input_meta_data$diagnosis,
               levels = c(0, 1)),
        p = .8,
        list = FALSE,
        times = 1
      )
    
    binom_elife_x_in <-
      feature_set[, colSums(feature_set != 0) > 0]
    
    binom_elife_y_in <-
      binom_model_input_meta_data
    
    binom_elife_x_in_train <-
      binom_elife_x_in[trainSplitIndex, , drop = F]
    
    binom_elife_y_in_train <-
      binom_elife_y_in[trainSplitIndex, , drop = F] %>% pull(diagnosis)
    
    binom_elife_x_in_test <-
      binom_elife_x_in[-trainSplitIndex, , drop = F]
    
    binom_elife_y_in_test <-
      binom_elife_y_in[-trainSplitIndex, , drop = F] %>% pull(diagnosis)
    
    set.seed(2022 - 03 - 08)
    fold_ids <-
      caret::createFolds(factor(binom_elife_y_in_train), 10, list = F)
    
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
        cutoffs[which.min(abs(cutoffs$fpr - (1 - 0.9))), ]
      perf_list$cut <-  inverse.logit.calc(perf_list$cut)
      perf_list$auc <-  auc
      perf_list$alpha <- alpha
      perf_list$feature_set <- feature_set_name
      
      list(
        'perf_list' = perf_list,
        'model' = binom_elife_model,
        'test_x' = binom_elife_x_in_test,
        'test_y' = binom_elife_y_in_test
      )
      
    })
    
  })
  
}) -> final_alpha
