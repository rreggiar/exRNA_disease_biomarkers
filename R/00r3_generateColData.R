#!/user/bin/env Rscript

## rreggiar@ucsc.edu

## generate 'colData' for tximeta aggregation of data using the naming conventions established 
##+ for organizing RNAseq data

# data/
## input_data
### experiment_name/
#### sample_name -- condition.rep.input_vol* *(if needed)

#[tximeta vignette](https://bioconductor.org/packages/devel/bioc/vignettes/tximeta/inst/doc/tximeta.html#Analysis_starts_with_sample_table)
#[alevin process aware](https://combine-lab.github.io/alevin-tutorial/2020/alevin-velocity/)

suppressPackageStartupMessages({
library(tidyverse)
library(tximeta)
library(SummarizedExperiment)
library(HDF5Array)
})

parse.input <- function(data_path, output_name, gencode_ver) {

	files <- dir(file.path(data_path))
	print('files')
	print(files)
	output_files <- list.files(data_path, full.names = TRUE)
    	names(output_files) <- output.files

    	print('output.files:')
    	print(output_files)

	quit()
    
    # creates empty df to fill with sample info
    sample.info <- data.frame('files' = file.path(),
    			      'names' = character(),
                              'condition' = character(),
                              'rep' = numeric(),
                              'input_vol' = numeric())

    # iteratively add to the empyty df the parsed sample names
    # here with format 'condition.rep'
    lapply(output_files, function(x){
		   x.path <- Sys.glob(file.path(x,paste0("*v", gencode_ver, ".",output_name,"*"), "quant.sf"))[1]
		   x <- basename(x)
		   sample.split <- str_split(x, '[.]', n=4)[[1]]
		   condition <- sample.split[1]
		   rep <- sample.split[2]
		   input_vol <- as.numeric(paste0(sample.split[3], '.', sample.split[4]))

            temp.df <- data.frame('files' = x.path,
            			  'names' = x,
                                  'condition' = condition,
                                  'rep' = rep,
                                  'input_vol' = input_vol)

            sample.info <<- rbind(sample.info, temp.df)
    })

    head(sample.info)

    sample.info.df <<- 
        sample.info %>%
        bind_rows() %>%
        filter(! is.na(files))

    print(sample.info.df)
}

build.tximeta.obj <- function() {

	if (output_name == 'ucsc.rmsk.salmon') {
		metaSkip = TRUE
		print('import data with tximeta')
		txi <- tximeta::tximeta(coldata = sample.info.df, type = 'salmon', skipMeta=metaSkip)


		print('save tx h5')
		saveHDF5SummarizedExperiment(txi, dir=paste0(outpath,paste0(output_name, '_h5_se')) , replace=TRUE)

		gxi <- tximeta::tximeta(coldata = sample.info.df, 
			type = 'salmon', 
			skipMeta=metaSkip, 
			txOut=FALSE, 
			tx2gene=tx2gene)

		print('save gene h5')
		saveHDF5SummarizedExperiment(gxi, dir=paste0(outpath,paste0(output_name, '_gene_h5_se')), replace=TRUE)


	} else {
		metaSkip = FALSE
		print('import data with tximeta')
		txi <- tximeta::tximeta(coldata = sample.info.df, type = 'salmon', skipMeta=metaSkip)
		print('save tx h5')
		saveHDF5SummarizedExperiment(txi, dir=paste0(outpath,paste0(output_name, '_h5_se')) , replace=TRUE)

		gxi <- summarizeToGene(txi)

		print('save gene h5')
		saveHDF5SummarizedExperiment(gxi, dir=paste0(outpath,paste0(output_name, '_gene_h5_se')), replace=TRUE)

		if (!is.null(txome_tsv)){
			cg <- read.delim(txome_tsv, header = TRUE, as.is = TRUE)
		
			colnames(cg)[colnames(cg) == 'intron'] <- 'unspliced'
			split.gxi <- tximeta::splitSE(gxi, cg, assayName = 'counts')
		
			saveHDF5SummarizedExperiment(split.gxi, 
				dir=paste0(outpath,paste0(output_name, '_gene_split_h5_se')), 
				replace=TRUE)
		}

	}
}



main <- function() {

	paths.in <- scan(file=file("stdin", "r"), what="character", n=6)

	data_path <- paths.in[1] 
	txome_path <- paths.in[2]
	output_name <- paths.in[3]
	txome_tsv <- paths.in[4]
	outpath <- paths.in[5]
	gencode_ver <- paths.in[6]

	print('parse input dir')
	parse.input(data_path, output_name, gencode_ver)

	#print('load tximeta')
	#tximeta::loadLinkedTxome(txome_path)

	#print(paste0('load tx2gene for gencode v ', gencode_ver))
	#tx2gene_path <- file.path('/public/groupls/kimlab/genomes.annotation/gencode.', gencode_ver, paste0('gencode.v', gencode_ver, '.ucsc.rmsk.tx2gene.csv'))

	#tx2gene <- read_csv(tx2gene_path, col_names=F)
	
	#print('build tximeta object')
	#build.tximeta.obj()

}

main()
