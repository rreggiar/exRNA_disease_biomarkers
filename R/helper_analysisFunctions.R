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
	       plot = plot.in)
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
											     paste0(sample, 
												    '_gene_split_h5_se')))
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

	  lst('txi' = data$txi,
	      'gxi' = data$gxi,
	      'quant_meta' = meta_out) 
	  
	  }) -> return_lst 
  
  return_lst
}

subset.tximeta.se <- function(se, filter_list = qc_fails) { 
  
  import::here(.from = 'SummarizedExperiment', colData)
  import::here(.from = 'tibble', lst)
  
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
  
}

build.analysis.set <- function(se_list, 
                               se_1,
                               se_2 = NULL,
                               analysis_set_2 = NULL) { 
  
  if(! "SparseSummarizedExperiment" %in% rownames(installed.packages())) { 
    devtools::install_github("PeteHaitch/SparseSummarizedExperiment")
  }
  
  if(is.null(analysis_set_2) & !is.null(se_2)) { 
    
    quant_meta_return <- 
      lapply(salmon_quant[names(salmon_quant) %in% c(se_1, se_2)], 
             function(se) { se$quant_meta }) %>% bind_rows()
    
    gxi_return <-
      SparseSummarizedExperiment::cbind(se_list[[se_1]]$gxi,
                                        se_list[[se_2]]$gxi)
    
    txi_return <-
      SparseSummarizedExperiment::cbind(se_list[[se_1]]$txi,
                                        se_list[[se_2]]$txi)
    
    lst('gxi' = gxi_return,
        'txi' = txi_return,
        'quant_meta' = quant_meta_return)
    
  } else {
    
    quant_meta_return <- 
      rbind(se_list[['quant_meta']][names(se_list[['quant_meta']]) %in% c(se_1)][[1]],
            analysis_set_2[['quant_meta']]) %>% 
      distinct()
    
    gxi_return <-
      SparseSummarizedExperiment::cbind(se_list[[se_1]]$gxi,
                                        analysis_set_2[['gxi']])
    txi_return <-
      SparseSummarizedExperiment::cbind(se_list[[se_1]]$txi,
                                        analysis_set_2[['txi']])
    
    lst('gxi' = gxi_return,
        'txi' = txi_return,
        'quant_meta' = quant_meta_return)
    
  }
  
}

