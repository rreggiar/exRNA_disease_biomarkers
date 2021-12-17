#!/user/bin/env Rscript

## rreggiar@ucsc.edu

## generate 'colData' for tximeta aggregation of data using the naming conventions established 
##+ for organizing RNAseq data

# data/
## experiment_name/
### sample_name -- condition.rep.input_vol* *(if needed)

#[tximeta vignette](https://bioconductor.org/packages/devel/bioc/vignettes/tximeta/inst/doc/tximeta.html#Analysis_starts_with_sample_table)
#[alevin process aware](https://combine-lab.github.io/alevin-tutorial/2020/alevin-velocity/)

suppressPackageStartupMessages({
library(tidyverse)
library(Seurat)
library(tximeta)
library(SummarizedExperiment)
library(HDF5Array)
})

parse.input <- function(data.path, output.name){
    ## assembles metadata object from input data

    files <- dir(file.path(data.path))

    # names(files) <- files
    output.files <- list.files(data.path, full.names = TRUE)
    names(output.files) <- output.files

    #print('output.files:')
    #print(output.files)
    
    #output.files.samples <- list.files(output.files, full.names = TRUE)
    #names(output.files.samples) <- output.files.samples


    # creates empty df to fill with sample info
    sample.info <- data.frame('files' = file.path(),
    						  'names' = character(),
                              'condition' = character(),
                              'rep' = numeric())

    # iteratively add to the empyty df the parsed sample names
    # here with format 'condition.rep'
    lapply(output.files, function(x){

    		x.path <- Sys.glob(file.path(x,paste0("*v35.",output.name,"*"), "alevin", "quants_mat.gz"))[1]
    		x <- basename(x)
            sample.split <- str_split(x, '[.]', n=2)[[1]]
            condition <- sample.split[1]
            rep <- sample.split[2]
            temp.df <- data.frame('files' = x.path,
            					  'names' = x,
                                  'condition' = condition,
                                  'rep' = rep)
            sample.info <<- rbind(sample.info, temp.df)
    })

    head(sample.info)

    sample.info.df <<- 
        sample.info %>%
        bind_rows() %>%
        filter(! is.na(files)) #%>%
        # mutate(condition = relevel(as.factor(condition), ref = 'ctrl'))

    print(sample.info.df)
}


paths.in <- scan(file=file("stdin", "r"), what="character", n=4)
data.path <- paths.in[1]
txome.path <- paths.in[2]
output.name <- paths.in[3]
txome.tsv <- paths.in[4]


outpath <- file.path("/public/groups/kimlab/pittsburgh_pah_exoRNA/output.data/")


print('parse input dir')
parse.input(data.path, output.name)

print('load tximeta')
tximeta::loadLinkedTxome(txome.path)

tx.to.gene <- read_csv('/public/groups/kimlab/genomes.annotations/gencode.35/gencode.v35.ucsc.rmsk.tx.to.gene.csv', col_names=F)


print('import data with tximeta')
lapply(1:nrow(sample.info.df), function(index){

    print(sample.info.df[index, 1])

    txi <- 
        tximeta::tximeta(sample.info.df[index, 1],
                         type = 'alevin', useHub = F, skipMeta = T,
                         alevinArgs=list(filterBarcodes=TRUE))

    counts <- assay(txi, 'counts')
    # returns sparse matrix, maybe will work in Seurat maybe coerce to mat

    seurat <- Seurat::CreateSeuratObject(counts = counts, 
                                         min.cells = 0, min.features = 0,
                                         project = sample.info.df[index, 2])

    name <- sample.info.df[index, 2][[1]]

    return_list <- lst(seurat)

    print(name)

    names(return_list) <- name

    return_list

}) %>% purrr::flatten() -> sc_txi

saveRDS(sc_txi, file=paste0(outpath,paste0(output.name, '_sc.rds')))
