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

subset.tximeta.object <- function(tximeta.list, grep.string){

  if (grepl('\\|', grep.string)){

    subset.txi <- SummarizedExperiment::subset(tximeta.list$txi,
                                         select = !grepl(grep.string, colData(tximeta.list$txi)$names))

    subset.gxi <- SummarizedExperiment::subset(tximeta.list$gxi,
                                         select = !grepl(grep.string, colData(tximeta.list$gxi)$names))
  } else{

    subset.txi <- SummarizedExperiment::subset(tximeta.list$txi,
                                         select = grepl(grep.string, colData(tximeta.list$txi)$names))

    subset.gxi <- SummarizedExperiment::subset(tximeta.list$gxi,
                                         select = grepl(grep.string, colData(tximeta.list$gxi)$names))
  }

  return(list('txi' = subset.txi,
              'gxi' = subset.gxi))
}
