save.manuscript.panel <- function(figure, 
				  figure.dir,
				  name, 
				  width = unit(89, 'mm'), 
				  height = unit(183, 'mm'), 
				  plot.in = last_plot()) {

	if(!dir.exists(paste0(figure.dir, '/','fig.',figure))) {
		dir.create(paste0(figure.dir, '/','fig.',figure))
	}
  
	  ggsave(paste0(figure.dir,
			'/',
			'fig.',
			figure,
			'/', 
			name, 
			'.',
			round(width, 3), 
			'x', 
			round(height, 3), 
			'.pdf'),
	       width = width, 
	       height = height, 
	       units = 'mm', 
	       plot = plot.in, device = cairo_pdf)
}

load.tximeta.object.list <- function(reference, output_data.dir = here::here('data/output_data')) {
	import::here(.from = "HDF5Array", loadHDF5SummarizedExperiment)
	import::here(.from = "purrr", imap)

	output_data_path <- file.path(output_data.dir, paste0(reference, '_quant'))
	message(output_data_path)

	list.files(output_data_path) -> tximeta_obj.list
	sub('_h5_se', '', tximeta_obj.list) -> tmp
	sub('_gene', '', tmp) -> tmp
	sub('_split', '', tmp) %>% unique() -> tximeta_sample.list
	names(tximeta_sample.list) <- tximeta_sample.list	

	print(tximeta_sample.list)

	tximeta_sample.list %>%
		imap(function(name, sample) {

			     txi <- loadHDF5SummarizedExperiment(dir=file.path(output_data_path, 
									       paste0(sample, '_h5_se')))
			     gxi <- loadHDF5SummarizedExperiment(dir=file.path(output_data_path, 
									       paste0(sample, '_gene_h5_se')))


			     if (reference == 'process.aware.salmon'){
				     gxi.split <- loadHDF5SummarizedExperiment(dir=file.path(output_data_path,
											     paste0(sample, '_gene_split_h5_se')))
				     return(list('txi' = txi, 'gxi' = gxi, 'gxi.split' = gxi.split)) 
			     } else {
				     return(list('txi' = txi, 'gxi' = gxi)) 
			     }
	       }
	)
}

parse.tximeta.quant.metadata <- function(se_list, qc_data.dir = here::here('data/output_data/rna_qc')) { 
  
  import::here(.from = "purrr", imap) 
  import::here(.from = "tibble", lst)
  import::here(.from = 'readr', read_tsv)
  
  salmon_quant_cols_of_interest <-
    c("frag_length_mean", "frag_length_mean", "frag_length_sd",
    "num_processed", "num_mapped", "num_decoy_fragments",
    "num_dovetail_fragments", "num_fragments_filtered_vm", "num_alignments_below_threshold_for_mapped_fragments_vm",
    "percent_mapped") 
  
  imap(se_list, function(data, project) { 
    
    message(project)

	  meta_obj <- data$gxi@metadata$quantInfo
	  names <- data$gxi$names

	  meta_obj[names(meta_obj) %in% salmon_quant_cols_of_interest] %>%
	    as.data.frame() -> salmon_meta_tmp

	  colnames(salmon_meta_tmp) <- paste0('salmon_', colnames(salmon_meta_tmp))

	  salmon_meta_tmp %>% mutate(sample = names) -> salmon_meta_tmp
	  
	  process_aware = F
	  if(grepl('process.aware', project)) {process_aware = T}

	  project <- sub('_salmon', '', project)
	  project <- sub('_ucsc.rmsk.salmon', '', project)
	  project <- sub('_process.aware.salmon', '', project)
	  
	  star_meta_tmp_path <- Sys.glob(file.path(qc_data.dir,
    					                     'star',
    					                     paste0('multiqc.star.', project, '.*_data'),
    					                     'multiqc_star.txt'))

	  star_meta_tmp <- read_tsv(star_meta_tmp_path)

	  colnames(star_meta_tmp) <- paste0('star_', colnames(star_meta_tmp))

	  star_meta_tmp %>%
	    dplyr::rename('sample' = star_Sample) %>%
	    mutate(sample = sub('_second_pass_out', '', sample)) -> star_meta_tmp

	  meta_out <- merge(star_meta_tmp, salmon_meta_tmp, by = 'sample')
	  
	  if (process_aware == T){

	    return(lst('txi' = data$txi, 
        	        'gxi' = data$gxi, 
        	        'gxi.split' = data$gxi.split,
        	        'quant_meta' = meta_out))
  
	  } else {
	    return(lst('txi' = data$txi,
        	        'gxi' = data$gxi,
        	        'quant_meta' = meta_out))
	  }
	  
	  }) -> return_lst 
  
  return_lst
}

subset.tximeta.se <- function(se, filter_list = qc_fails) { 
  
  import::here(.from = 'SummarizedExperiment', colData)
  import::here(.from = 'tibble', lst)
  
  process_aware = F
  if(grepl('process.aware', deparse(substitute(se)))) { process_aware = T }
  
  txi_return <- 
    SummarizedExperiment::subset(se$txi, 
                                 select = !colData(se$txi)$names %in% filter_list)
  
  gxi_return <- 
    SummarizedExperiment::subset(se$gxi, 
                                 select = !colData(se$gxi)$names %in% filter_list)
  
  quant_meta_return <- 
    se[['quant_meta']][! se[['quant_meta']]$sample %in% filter_list, ]
  
  lst('txi' = txi_return,
      'gxi' = gxi_return,
      'quant_meta' = quant_meta_return)
  
  if (process_aware == T){
    
    gxi_split_return <- 
      SummarizedExperiment::subset(se$gxi.split, 
                                   select = !colData(se$gxi)$names %in% filter_list)
    
    return(lst('txi' = txi_return, 
               'gxi' = gxi_return, 
               'gxi.split' = gxi_split_return,
               'quant_meta' = quant_meta_return))
    
  } else {
    return(lst('txi' = txi_return,
               'gxi' = gxi_return,
               'quant_meta' = quant_meta_return))
  }
  
}

build.analysis.set <- function(se_list, 
                               se_1,
                               se_2 = NULL,
                               analysis_set_2 = NULL,
                               sample_meta_in) { 
  
  import::here(.from = 'SummarizedExperiment', colData)
  
  process_aware = F
  if(grepl('process.aware', deparse(substitute(se_list)))) { process_aware = T }
  
  if(! "SparseSummarizedExperiment" %in% rownames(installed.packages())) { 
    devtools::install_github("PeteHaitch/SparseSummarizedExperiment")
  }
  
  if(is.null(analysis_set_2) & !is.null(se_2)) { 
    
    quant_meta_return <- 
      lapply(se_list[names(se_list) %in% c(se_1, se_2)], 
             function(se) { se$quant_meta }) %>% bind_rows()
    
    sample_meta_return <- 
      lapply(se_list[names(se_list) %in% c(se_1, se_2)], 
             function(se) { colData(se$gxi) %>% 
                             as.data.frame() %>% 
                             remove_rownames() }) %>% bind_rows()
    
    merge(quant_meta_return, 
          sample_meta_return, by.x = 'sample',by.y = 'names') %>% 
      merge(sample_meta_in %>% select(-condition),
            by.x = 'sample', by.y = 'patient') -> quant_meta_return
    
    gxi_return <-
      SparseSummarizedExperiment::cbind(se_list[[se_1]]$gxi,
                                        se_list[[se_2]]$gxi)
    
    txi_return <-
      SparseSummarizedExperiment::cbind(se_list[[se_1]]$txi,
                                        se_list[[se_2]]$txi)
    
    if (process_aware == T){
      
      gxi_split_return <- 
        SparseSummarizedExperiment::cbind(se_list[[se_1]]$gxi.split,
                                          se_list[[se_2]]$gxi.split)
      
      return(lst('txi' = txi_return, 
                 'gxi' = gxi_return, 
                 'gxi.split' = gxi_split_return,
                 'quant_meta' = quant_meta_return))
      
    } else {
      return(lst('gxi' = gxi_return,
                 'txi' = txi_return,
                 'quant_meta' = quant_meta_return))
    }
    
  } else {
    
    quant_meta_return <- 
      lapply(se_list[names(se_list) %in% c(se_1)], 
             function(se) { se$quant_meta }) %>% bind_rows()
    
    sample_meta_return <- 
      lapply(se_list[names(se_list) %in% c(se_1, se_2)], 
             function(se) { colData(se$gxi) %>% 
                 as.data.frame() %>% 
                 remove_rownames() }) %>% bind_rows()
    
    merge(quant_meta_return, 
          sample_meta_return, by.x = 'sample',by.y = 'names') %>% 
      merge(sample_meta_in %>% select(-condition),
            by.x = 'sample', by.y = 'patient') -> quant_meta_return
    
    quant_meta_return <- 
      rbind(quant_meta_return,
            analysis_set_2[['quant_meta']]) %>% 
      distinct()
    
    gxi_return <-
      SparseSummarizedExperiment::cbind(se_list[[se_1]]$gxi,
                                        analysis_set_2[['gxi']])
    txi_return <-
      SparseSummarizedExperiment::cbind(se_list[[se_1]]$txi,
                                        analysis_set_2[['txi']])
    
    if (process_aware == T){
      
      gxi_split_return <- 
        SparseSummarizedExperiment::cbind(se_list[[se_1]]$gxi.split,
                                          analysis_set_2[['gxi.split']])
      
      return(lst('txi' = txi_return, 
                 'gxi' = gxi_return, 
                 'gxi.split' = gxi_split_return,
                 'quant_meta' = quant_meta_return))
      
    } else {
      return(lst('gxi' = gxi_return,
                 'txi' = txi_return,
                 'quant_meta' = quant_meta_return))
    }
    
  }
  
}

run.de.seq <- function(type = 'gxi', base_level = 'ctrl', 
                       input_se = salmon_quant$analysis_set$ctrl_panc_covid, 
                       dds_formula = as.formula('~age+input_vol+condition'), 
                       ref_type = 'salmon',
                       reference_meta_in = reference_meta_data$total_tx_to_gene.df) {  
  
  import::here(.from = 'SummarizedExperiment', colData, assay, rowRanges, assays)
  import::here(.from = 'DESeq2', .all = T)
  import::here(.from = 'tibble', lst)
  
  if(! "apeglm" %in% rownames(installed.packages())) { 
    BiocManager::install('apeglm', update = FALSE)
  }
  
  scaled_quant_meta_for_de.df <- 
    input_se$quant_meta %>% 
    mutate(input_vol = as.numeric(input_vol),
           condition = relevel(as.factor(condition), ref = base_level)) %>% 
    mutate_if(is.numeric, ~scale(., center = T)) %>% 
    mutate_if(is.character, ~as.factor(.))
  
  print(levels(scaled_quant_meta_for_de.df$condition))
  
  if (ref_type == 'ucsc.rmsk') { 
    
    gencode_sex_filter.df <-
      rowRanges(input_se[[type]]) %>% as.data.frame() %>% 
      filter(seqnames %in% c('chrY'),
             grepl('ENSG', group_name))
    
  } else {
    gencode_sex_filter.df <-
      rowRanges(input_se[[type]]) %>% as.data.frame() %>% 
      filter(seqnames %in% c('chrY'))
    
  }
  
  count_matrix.df <- 
    assay(input_se[[type]], 'counts') %>% 
    as.data.frame() %>% 
    rownames_to_column('gene') %>% 
    mutate_if(is.numeric, round) %>% 
    column_to_rownames('gene')
  
  input_set_dds <- DESeqDataSetFromMatrix(countData = count_matrix.df, 
                                          colData = scaled_quant_meta_for_de.df, 
                                          design = dds_formula)
  
  input_set_dds <- estimateSizeFactors(input_set_dds)
  
  input_set_dds_norm_counts.df <- as.data.frame(counts(input_set_dds, normalized=T))
  
  condition_aware_age_cor_filter <- 
    input_set_dds_norm_counts.df %>% 
    rownames_to_column('ensg') %>% 
    gather(sample, count, -ensg) %>% 
    merge(input_se$quant_meta %>% select(sample, age, condition), by = 'sample') %>% 
    group_by(condition, ensg) %>% 
    summarize(age_cor = cor(count, age, method = 'pearson')) %>% 
    drop_na() %>% 
    filter(abs(age_cor) >= 0.70) %>% 
    spread(condition, age_cor)
  
  input_set_dds_vst_counts <-  as.data.frame(assay(vst(input_set_dds, blind = F)))
  
  de_filter <- rowSums(counts(input_set_dds, normalized=T) >= 10) >= ncol(input_set_dds_norm_counts.df)*0.75
  
  input_set_dds <- DESeq(input_set_dds[de_filter,])
  
  print(resultsNames(input_set_dds))
  
  tail(resultsNames(input_set_dds), n=2)
  
  coef_list <- tail(resultsNames(input_set_dds), n=2)
  
  if (ref_type != 'ucsc.rmsk') { 
    
    covid_vs_ctrl <- 
      lfcShrink(input_set_dds, coef = coef_list[[1]]) %>% 
      as.data.frame() %>% 
      rownames_to_column('ensg') %>% 
      filter(!ensg %in% c(condition_aware_age_cor_filter[, c('ctrl', 'covid')]$ensg,
                          gencode_sex_filter.df$gene_id)) %>%
      merge(rowRanges(input_se[[type]]) %>% as.data.frame(), by.x = 'ensg', by.y = 'gene_id') %>%
      merge(reference_meta_in %>% select(ensg, gene) %>% distinct(), by = 'ensg')
    
    panc_vs_ctrl <- 
      lfcShrink(input_set_dds, coef = coef_list[[2]]) %>% 
      as.data.frame() %>% 
      rownames_to_column('ensg') %>% 
      filter(!ensg %in% c(condition_aware_age_cor_filter[, c('ctrl', 'panc')]$ensg,
                          gencode_sex_filter.df$gene_id)) %>% 
      merge(rowRanges(input_se[[type]]) %>% as.data.frame(), by.x = 'ensg', by.y = 'gene_id') %>%
      merge(reference_meta_in %>% select(ensg, gene) %>% distinct(), by = 'ensg') 
  } else { 
    
    covid_vs_ctrl <- 
      lfcShrink(input_set_dds, coef = coef_list[[1]]) %>% 
      as.data.frame() %>% 
      rownames_to_column('ensg') %>% 
      filter(!ensg %in% c(condition_aware_age_cor_filter[, c('ctrl', 'covid')]$ensg,
                          gencode_sex_filter.df$gene_id)) %>%
      merge(reference_meta_in %>% select(ensg, gene) %>% distinct(), by = 'ensg')
    
    panc_vs_ctrl <- 
      lfcShrink(input_set_dds, coef = coef_list[[2]]) %>% 
      as.data.frame() %>% 
      rownames_to_column('ensg') %>% 
      filter(!ensg %in% c(condition_aware_age_cor_filter[, c('ctrl', 'panc')]$ensg,
                          gencode_sex_filter.df$gene_id)) %>%
      merge(reference_meta_in %>% select(ensg, gene) %>% distinct(), by = 'ensg')
    
  }
  
  lst('norm_counts.df' = input_set_dds_norm_counts.df,
      'vst_counts.df' = input_set_dds_vst_counts,
      'age_filter.df' = condition_aware_age_cor_filter,
      'sex_filter.df' = gencode_sex_filter.df,
      'dds_object' = input_set_dds,
      'panc_vs_ctrl' = panc_vs_ctrl,
      'covid_vs_ctrl' = covid_vs_ctrl,
      'scaled_quant_meta_for_de.df' = scaled_quant_meta_for_de.df)  
  
}

run.de.seq.individual <- function(type = 'gxi', base_level = 'ctrl', 
                                  top_level = 'panc',
                                   input_se = salmon_quant$analysis_set$ctrl_panc_covid, 
                                   dds_formula = as.formula('~age+input_vol+condition'), 
                                   ref_type = 'salmon',
                                   reference_meta_in = reference_meta_data$total_tx_to_gene.df) {  
  
  import::here(.from = 'SummarizedExperiment', colData, assay, rowRanges, assays)
  import::here(.from = 'DESeq2', .all = T)
  import::here(.from = 'tibble', lst)
  
  if(! "apeglm" %in% rownames(installed.packages())) { 
    BiocManager::install('apeglm', update = FALSE)
  }
  
  scaled_quant_meta_for_de.df <- 
    input_se$quant_meta %>% 
    mutate(input_vol = as.numeric(input_vol),
           condition = relevel(as.factor(condition), ref = base_level)) %>% 
    mutate_if(is.numeric, ~scale(., center = T)) %>% 
    mutate_if(is.character, ~as.factor(.))
  
  print(levels(scaled_quant_meta_for_de.df$condition))
  
  if (ref_type == 'ucsc.rmsk') { 
    
    gencode_sex_filter.df <-
      rowRanges(input_se[[type]]) %>% as.data.frame() %>% 
      filter(seqnames %in% c('chrY'),
             grepl('ENSG', group_name))
    
  } else if (ref_type == 'process.aware') {
    
    gencode_sex_filter.df <-
      rowRanges(input_se[['gxi']]) %>% as.data.frame() %>% 
      filter(seqnames %in% c('chrY'))
    
  } else {
    gencode_sex_filter.df <-
      rowRanges(input_se[[type]]) %>% as.data.frame() %>% 
      filter(seqnames %in% c('chrY'))
    
  }
  
  if (ref_type == 'process.aware') {
    
    count_matrix.df <- 
      assay(input_se[[type]], 'unspliced') %>% 
      as.data.frame() %>% 
      rownames_to_column('gene') %>% 
      mutate_if(is.numeric, round) %>% 
      column_to_rownames('gene')
    
  } else {
  
    count_matrix.df <- 
      assay(input_se[[type]], 'counts') %>% 
      as.data.frame() %>% 
      rownames_to_column('gene') %>% 
      mutate_if(is.numeric, round) %>% 
      column_to_rownames('gene')
    
  }
  
  input_set_dds <- DESeqDataSetFromMatrix(countData = count_matrix.df, 
                                          colData = scaled_quant_meta_for_de.df, 
                                          design = dds_formula)
  
  input_set_dds <- estimateSizeFactors(input_set_dds)
  
  input_set_dds_norm_counts.df <- as.data.frame(counts(input_set_dds, normalized=T))
  
  condition_aware_age_cor_filter <- 
    input_set_dds_norm_counts.df %>% 
    rownames_to_column('ensg') %>% 
    gather(sample, count, -ensg) %>% 
    merge(input_se$quant_meta %>% select(sample, age, condition), by = 'sample') %>% 
    group_by(condition, ensg) %>% 
    summarize(age_cor = cor(count, age, method = 'pearson')) %>% 
    drop_na() %>% 
    filter(abs(age_cor) >= 0.70) %>% 
    spread(condition, age_cor)
  
  input_set_dds_vst_counts <- as.data.frame(assay(vst(input_set_dds, blind = F)))
  
  de_filter <- rowSums(counts(input_set_dds, normalized=T) >= 10) >= ncol(input_set_dds_norm_counts.df)*0.75
  
  input_set_dds <- DESeq(input_set_dds[de_filter,])
  
  print(resultsNames(input_set_dds))
  
  tail(resultsNames(input_set_dds), n=1)
  
  coef_list <- tail(resultsNames(input_set_dds), n=1)
  
  print(coef_list[[1]])
  
  if (ref_type == 'salmon') { 
    
    message('gencode')
    
    de_out <- 
      lfcShrink(input_set_dds, coef = coef_list[[1]]) %>% 
      as.data.frame() %>% 
      # indiv comparisons return flipped values ?
      # mutate(log2FoldChange = -log2FoldChange) %>% 
      rownames_to_column('ensg') %>% 
      filter(!ensg %in% c(condition_aware_age_cor_filter[, c(base_level, top_level)]$ensg,
                          gencode_sex_filter.df$gene_id)) %>%
      merge(rowRanges(input_se[[type]]) %>% as.data.frame(), by.x = 'ensg', by.y = 'gene_id') %>%
      merge(reference_meta_in %>% select(ensg, gene) %>% distinct(), by = 'ensg')
    
  } else if (ref_type == 'process.aware') {
    
    message('process.aware')
    
    de_out <- 
      lfcShrink(input_set_dds, coef = coef_list[[1]]) %>% 
      as.data.frame() %>% 
      # indiv comparisons return flipped values ?
      # mutate(log2FoldChange = -log2FoldChange) %>%
      rownames_to_column('ensg') %>%
      # filter(!ensg %in% c(condition_aware_age_cor_filter[, c(base_level, top_level)]$ensg,
      #                     gencode_sex_filter.df$gene_id)) %>%
      merge(reference_meta_in %>% select(ensg, gene) %>% distinct(), by = 'ensg')
    
  } else { 
    
    message('ucsc.rmsk')
    
    de_out <- 
      lfcShrink(input_set_dds, coef = coef_list[[1]]) %>% 
      as.data.frame() %>% 
      # indiv comparisons return flipped values ?
      # mutate(log2FoldChange = -log2FoldChange) %>% 
      rownames_to_column('ensg') %>% 
      filter(!ensg %in% c(condition_aware_age_cor_filter[, c(base_level, top_level)]$ensg,
                          gencode_sex_filter.df$gene_id)) %>%
      merge(reference_meta_in %>% select(ensg, gene) %>% distinct(), by = 'ensg')
    
  }
  
  lst('norm_counts.df' = input_set_dds_norm_counts.df,
      'vst_counts.df' = input_set_dds_vst_counts,
      'age_filter.df' = condition_aware_age_cor_filter,
      'sex_filter.df' = gencode_sex_filter.df,
      'dds_object' = input_set_dds,
      'de_out' = de_out,
      'scaled_quant_meta_for_de.df' = scaled_quant_meta_for_de.df)  
  
}

run.pca <- function(input_de, 
                    identity_color_pal, 
                    plot_tag = 'A') {
  
  import::here(.from = 'tibble', lst)
  import::here(.from = 'stats', prcomp)
  import::here(.from = 'patchwork', .all = T)
  import::here(.from = 'ggplot2', labs)
  
  input_de$vst_counts.df[! rownames(input_de$vst_counts.df) %in% 
                           c(input_de$age_filter.df$ensg,input_de$sex_filter.df$gene_id), ] %>% 
    t() -> pca_tmp_in 
  
  pca_tmp_in_sd <- apply(pca_tmp_in, 2, sd)
  pca_tmp_in_zero_sd <- pca_tmp_in[,pca_tmp_in_sd!=0]
  non_zero_pca <- prcomp(pca_tmp_in_zero_sd, center = T, scale. = T, rank. = 50)
  
  summary(non_zero_pca)
  pcs.of.interest <- apply(non_zero_pca$x, 2, var)
  pcs.props <-  pcs.of.interest/sum(pcs.of.interest)
  cumsum(pcs.props)[c(1,2)]
  print(as.data.frame(pcs.props) %>% 
          rownames_to_column('pc') %>% 
          mutate(pc = as.numeric(sub('PC', '', pc))) %>% 
          ggplot(aes(pc, pcs.props)) + 
          geom_col())
  
  pca.out <- as.data.frame(non_zero_pca$x) %>% 
    rownames_to_column('sample') %>% 
    merge(input_de$scaled_quant_meta_for_de.df, by = 'sample')
  
  pca.out.summary <- 
    as.data.frame(summary(non_zero_pca)$importance) %>% 
    rownames_to_column('metric') %>% 
    gather(pc, value, -metric) %>% 
    spread(metric, value)
  
  pca.out %>% 
    ggplot(aes(PC1, PC2,
               color = condition)) +
    geom_point(size = rel(1), alpha = 0.8) + 
    # scale_shape_manual(values = c(21,25)) +
    xlab(paste('PC1', round(cumsum(pcs.props)[1], digits = 3), sep = ' ')) +
    ylab(paste('PC2', round(cumsum(pcs.props)[2] - cumsum(pcs.props)[1], digits = 3), sep = ' ')) +
    geom_vline(xintercept = 0, linetype = 'dotted', alpha = 0.3, size = 0.25) +
    geom_hline(yintercept = 0, linetype = 'dotted', alpha = 0.3, size = 0.25) + 
    scale_color_manual(values = identity_color_pal) + 
    # theme(axis.title.x = element_blank()) +
    xlim(-200,200) -> pca_1v2.plt
  
  pca.out %>% 
    ggplot(aes(PC2, PC3,
               color = condition)) +
    geom_point(size = rel(1), alpha = 0.8) + 
    # scale_shape_manual(values = c(21,25)) +
    xlab(paste('PC2', round(cumsum(pcs.props)[2] - cumsum(pcs.props)[1], digits = 3), sep = ' ')) +
    ylab(paste('PC3', round(cumsum(pcs.props)[3] - cumsum(pcs.props)[2], digits = 3), sep = ' ')) +
    geom_vline(xintercept = 0, linetype = 'dotted', alpha = 0.3, size = 0.25) +
    geom_hline(yintercept = 0, linetype = 'dotted', alpha = 0.3, size = 0.25) + 
    scale_color_manual(values = identity_color_pal) + 
    ylim(-200,200) -> pca_2v3.plt
  
  pca.out %>% 
    ggplot(aes(PC3, PC4,
               color = condition)) +
    geom_point(size = rel(1), alpha = 0.8) + 
    # scale_shape_manual(values = c(21,25)) +
    xlab(paste('PC3', round(cumsum(pcs.props)[3] - cumsum(pcs.props)[2], digits = 3), sep = ' ')) +
    ylab(paste('PC4', round(cumsum(pcs.props)[4] - cumsum(pcs.props)[3], digits = 3), sep = ' ')) +
    geom_vline(xintercept = 0, linetype = 'dotted', alpha = 0.3, size = 0.25) +
    geom_hline(yintercept = 0, linetype = 'dotted', alpha = 0.3, size = 0.25) + 
    scale_color_manual(values = identity_color_pal) + 
    ylim(-200,200) -> pca_3v4.plt
  
  lst('pca_out.df' = pca.out,
      'pca_1v2.plt' = pca_1v2.plt,
      'pca_2v3.plt' = pca_2v3.plt,
      'pca_3v4.plt' = pca_3v4.plt)
}


extract.meta.pca.correlates <- function(pca.out) {
  
  Filter(function(x) sd(x) !=0, pca.out %>% select(where(is.numeric))) %>% 
    cor(method = 'pearson') %>% 
    as.data.frame() %>% 
    select(-starts_with('PC')) %>% 
    rownames_to_column('pc') %>% 
    filter(grepl('^PC', pc)) %>% 
    gather(var, cor, -pc) %>% 
    mutate(pc = as.factor(as.numeric(str_remove(pc, 'PC')))) %>% 
    filter(pc %in% c(1,2,3)) -> pca_pearson_cor.df
  
  pca_pearson_cor.df
  
  
}

make.fgsea.ranks <- function(de.seq){
  
  import::here(.from = 'fgsea', fgsea)
  
  # build a ranked gene list in the appropriate format (named list)
  # using shrunken log2 Fold Change values from DESeqII
  de.seq <- 
    de.seq %>%
    mutate(rank = as.double(log2FoldChange),
           Gene = as.character(gene)) %>% 
    dplyr::select(Gene, rank) %>%
    mutate(rank = scale(rank)) 
  
  de.seq.rnk <- 
    de.seq$rank
  de.seq.rnk <- 
    setNames(de.seq$rank, de.seq$Gene)
  # run the GSEA implementation on the hallmark sets supplemented with custom sets
  message('fgsea')
  de.seq.fgsea.res <- 
    fgsea(pathways = msig.df, 
          stats = de.seq.rnk, eps = 0.0)
  
  message('padj filter and clean names')
  # convert to a tidy format, clean the names for display
  de.seq.fgsea.df <-
    as_tibble(de.seq.fgsea.res) %>%
    filter(padj < 0.05) %>%
    mutate(pathway = gsub('_', ' ', pathway))
  
  return(de.seq.fgsea.df)
  
}

generate.unsplice.data <- function(analysis_set) {
  
  import::here(SummarizedExperiment, assay)
  
  split.gxi <- analysis_set$gxi.split
  
  assay(split.gxi, 'unspliced') %>% 
    as.data.frame() %>% 
    rownames_to_column('gene') %>% 
    gather(sample, count, -gene) %>% 
    group_by(sample) %>% 
    summarize(unsplice_count = sum(count)) -> unsplice.df
  
  assay(split.gxi, 'spliced') %>% 
    as.data.frame() %>% 
    rownames_to_column('gene') %>% 
    gather(sample, count, -gene) %>% 
    group_by(sample) %>% 
    summarize(splice_count = sum(count)) -> splice.df
  
  merge(unsplice.df, splice.df, by = 'sample') %>% 
    mutate(unsplice_rate = unsplice_count / (unsplice_count + splice_count))
  
}