#!/user/bin/env Rscript

## RER
## dkim lab

##################
suppressMessages({
library(tidyverse)
library(DESeq2)
library(tximport)
library(argparser)
})
##################

argparse <- function(){
    parser <- arg_parser('Import tx level quantification from Salmon, produce normalized counts of tx and genes, perform gene-level DESeq across specified condition(s)')

    parser <- add_argument(parser, '--input', help='path to salmon output', 
                           type='character',
                           default='/public/groups/kimlab/exoRNA-biomarkers-panc/data/collected.patient.data/')
    parser <- add_argument(parser, '--tx2gene', help='path to tx2gene csv',
                           type='character', 
                           default='/public/groups/kimlab/genomes.annotations/gencode.34/gen.34.ucsc.rmsk.tx.2.gene.csv')
    parser <- add_argument(parser, '--output', help='path to output *dir*', 
                           type='character')
    parser <- add_argument(parser, '--name', help='name of model/experiment',
                           type='character')
    parser <- add_argument(parser, '--first', help='first condition for comparison',
                           type='character',
                           default='panc')
    parser <- add_argument(parser, '--second', help='second condition for comparison',
                           type='character',
                           default='ctrl')
    parser <- add_argument(parser, '--metadata', help='patient meta data',
                           type='character',
                           default=NULL)
    parser <- add_argument(parser, '--default_formula', help='feature used as intercept',
                           type='character',
                           default='~condition')
    parser <- add_argument(parser, '--counts_only', help='do you only want counts?',
                           type='store_true',
                           default=FALSE)
    parser <- add_argument(parser, '--generate_counts', help='do you want to save counts?',
                           type='store_true',
                           default=TRUE)

    argv <- parse_args(parser)

    # '<<-' assigns object to parent scope (which in this case is global)
    in.path <<- argv$input
    out.path <<- argv$output
    tx2gene.path <<- argv$tx2gene
    name <<- argv$name
    one <<- argv$first
    two <<- argv$second
    meta.file <<- argv$metadata
    counts_only <<- argv$counts_only
    generate_counts <<- argv$generate_counts
    default.formula <<- as.formula(argv$default_formula)

    cat(paste('running DESeq normalization and comparison on ',
                in.path, 
                name,
                "\n",
                sep="\n"))

    if (!is.na(meta.file)){
        print(paste('metadata file: ', meta.file))
    }
    if(counts_only == TRUE){
        print('only generating counts')
    }
    if(generate_counts == FALSE){
        print('not writing counts')
    }

    print('creating output dir')
    dir.create(file.path(out.path, name))
    out.path <<- paste0(out.path, name,'/')
    print(paste0('out: ', out.path))

}

parse.input <- function(salmon.out){
    ## assembles metadata object from input data

    files <- dir(file.path(in.path))

    names(files) <- files

    # creates empty df to fill with sample info
    sample.info <- data.frame('sample' = character(),
                              'condition' = character(),
                              'rep' = numeric(), 'input_vol' = numeric())

    # iteratively add to the empyty df the parsed sample names
    # here with format 'condition.rep'
    lapply(files, function(x){
            x.split <- str_split(x, '[.]', n = 3)[[1]]
            condition <- x.split[1]
            rep <- x.split[2]
            input.vol <- x.split[3]
            temp.df <- data.frame('sample' = x,
                                  'condition' = condition,
                                  'rep' = rep,
                                  'input_vol' = input.vol)
            sample.info <<- rbind(sample.info, temp.df)
    })

    # convert to data frame
    # merge in the rest of the metadata
    duplicates <- c('panc.10.2.5', 'panc.6.2.5', 'panc.7.2.5', 'panc.9.2.5',
                    'ctrl.3.2.5', 'ctrl.2.3.4', 'ctrl.1.3.0', 'ctrl.1.5.0',
                    'ctrl.1.2.0', 'ctrl.1.1.5' , 'ctrl.1.0.5')


    sample.info.df <- 
        sample.info %>%
        filter(! sample %in% duplicates) %>%
        mutate(condition = relevel(as.factor(condition), ref = 'ctrl'))

    # sample.info.df <- sample.info.df[c(1:8, 20:23), ]

    print(sample.info.df)


    if (!is.na(meta.file)){
        meta.check <<- TRUE
        metadata.file <- read_csv(meta.file)

        sample.info.df <-
            sample.info.df %>%
            merge(metadata.file, by.x = 'sample', by.y = 'patient') %>%
            select(-condition.y) %>%
            rename(condition.x = 'condition') %>%
            mutate(input_vol = as.numeric(input_vol),
                   batch = as.factor(batch), 
                   cluster = relevel(as.factor(cluster), ref = two)) %>%
            mutate_if(is.numeric, ~ scale(., center = TRUE)) %>%
            mutate_if(is.numeric, as.factor) 
    }

    print(sample.info.df)

    # get paths to salmon output
    salmon.files <- 
        file.path(in.path, sample.info.df$sample, salmon.out, 'quant.sf')

    names(salmon.files) <- sample.info.df$sample

    tx2gene <- read_csv(tx2gene.path, col_names = F)
    names(tx2gene) <- c("Name", "GeneID")

    print('importing as transcripts')

    txi.tx <<- tximport(salmon.files, 
                type="salmon", 
                tx2gene = tx2gene,
                txOut = T)

    print('importing as genes')

    txi.gene <<- tximport(salmon.files, 
                    type="salmon", 
                    tx2gene = tx2gene,
                    txOut = F)

    return(sample.info.df) 
}

run.de.seq <- function(salmon.out, dds, formula_dir){
    ## takes argument that names the correct salmon output dir to look into
    ## some redundancy due to how tximport handles genes v txs
    ## following the vignette pretty closely, not reinventing the wheel

    def.feature <- as.character(default.formula)[[2]]

    coef <- paste0(def.feature, '_', one, '_vs_', two) 

    keep <- rowSums(counts(dds)) >= 10
    dds.txi <- dds[keep, ]
    dds.txi <- DESeq(dds.txi)
    print(resultsNames(dds.txi))
    res <- results(dds.txi, name = coef) 
    
    res05 <- results(dds.txi, alpha = 0.05)

    res05[order(res05$padj),]

    summary(res05)

    dir.create(paste0(out.path,formula_dir))

    pdf(paste0(out.path,formula_dir,'/','res.lfc.', name, '.ma.pdf'))

    plotMA(res, main='res')
    
    dev.off()

    print('shrinking log2fc')

    resLFC <- lfcShrink(dds.txi, 
                        coef = coef,
                        type = 'apeglm')

    pdf(paste0(out.path,formula_dir,'/','res.lfc.shrink.', name, '.ma.pdf'))

    plotMA(resLFC, main='shrunken')

    dev.off()

    return(resLFC)

}

formula.builder <- function(metadata){
    ## in: metadata object containing sample names and features
    ## out: list of possible DESeq formulae based on features

    # these will always be constants, not features, just descriptive labels
    constants <- c('condition', 'sample', 'rep', 'cluster')
    # empty list to populate with potential features
    variables <- c()

    cols_of_interest <- colnames(metadata %>% select(-all_of(constants)))
    
    lapply(cols_of_interest, function(x){
               variables <<- append(variables, x)
                        }
    )
    
    # rerwip: I need to understand this better. got from S/O
    # https://stackoverflow.com/questions/40049313/generate-all-combinations-of-all-lengths-in-r-from-a-vector
    variables <- 
        do.call("c", 
            lapply(seq_along(variables), 
                   function(i) combn(variables, i, FUN = list)))
    # empty list to store the formula objects
    formulae <- c()
    def.feature <- as.character(default.formula)[[2]]
    # add combinations to the formula list, pasted with the constant ~condition
    # feature at the end
    lapply(variables, function(comb){
            formulae <<- append(formulae, 
                                as.formula(
                                           paste(
                                                 paste0('~', 
                                                        paste(comb, collapse = ' + ')), 
                                                 def.feature, sep = ' + ')))
            }
    )
    
    return(formulae)
}

generate.norm.counts <- function(metadata, formula, salmon.out){
    ## in: medata object, formula, salmon output location
    ## out: deSeq objects holding counts, etc for tx's and genes in a list

    formula_dir <- sub(' [+] ', '_', as.character(formula)[[2]]) 
    #dir.create(paste0(out.path,formula_dir))

    print(paste0('path: ', salmon.out))

    if (generate_counts == TRUE){
        dds.tx <- DESeqDataSetFromTximport(txi.tx, 
                                       colData = metadata, 
                                       design = formula)

        dds.tx <- estimateSizeFactors(dds.tx)


        dds.tx.df <- as.data.frame(counts(dds.tx, normalized=T))
    }

    dds.gene <- DESeqDataSetFromTximport(txi.gene, 
                                         colData = metadata, 
                                         design = formula)

    dds.gene <- estimateSizeFactors(dds.gene)

    dds.gene.df <- as.data.frame(counts(dds.gene, normalized=T))


    print('returning counts and dds objects')
    
    if (generate_counts == TRUE){

        return(list('dds_gene' = dds.gene,
                    'gene_counts' = dds.gene.df,
                    'tx_counts' = dds.tx.df))
    } else{
        return(list('dds_gene' = dds.gene))
    }
}

write.data <- function(dds_list, salmon.out, formula_dir){
    ## in: the final deSeq output stored in named list object, salmon output,
    ## and parsed formula name for an outpur subdirectory
    ## out: writes csv's of counts and DE to formula subdir in output dir

    if (generate_counts == TRUE){

        print('writing tx level counts ..')

        write.csv(dds_list$tx_counts,
            paste0(out.path,
                   formula_dir,
                   '/',
                   'tx.', 
                   salmon.out, 
                   '.',
                   name,
                   '.',
                   one,
                   '.v.',
                   two,
                   '.de-seq.counts.csv'))
    
        print('done')


        print('writing gene level counts ..')

        write.csv(dds_list$gene_counts,
            paste0(out.path,
                   formula_dir,
                   '/',
                   'gene.', 
                   salmon.out, 
                   '.',
                   name,
                   '.',
                   one,
                   '.v.',
                   two,
                   '.de-seq.counts.csv'))
    }


    #if (length(dds_list$deSeq) != 0){

    print('writing gene level DE')
    
    tryCatch({
    write.csv(as.data.frame(dds_list$deSeq), 
        paste0(out.path,
                    formula_dir,
                    '/',
                    salmon.out, 
                    '.',
                    name,
                    '.',
                    one,
                    '.v.',
                    two,
                    '.de-seq.csv'))
    })

    print(' all done')
} 

run.and.collect <- function(ref){
    ## in: salmon output dir name
    ## out: optmized model DESeq output
    # implements tryCatch() to run all of the possible formulas built from the
    # metadata, most of which will fail due to co-linearity. Calls functions to
    # construct DESeq data objects and write finalized output.

    ref <- ref
    # parse the input and merge with metadata
    metadata <- parse.input(ref)
    # build out the possible formulas based on combinations of features
    design.formulae <- formula.builder(metadata)
    # empty DF that will be populated with formula performance and data
    # rerwip: this should be a class probably, like python, should look into
    # what I think are S3/4 classes and using them in R
    formula.performance <-
            tibble('formula' = character(),
                'sig_genes' = numeric(),
                'data_list' = list())
    # apply the run + collect algo to the reversed formula list to ensure that 
    # the longest (most complex) formula goes first. This ordering was initially
    # important as I was cross comparing the formulas to one another. As of now,
    # I'm simply comparing each formula to the base case of ~condition
    lapply(rev(design.formulae), function(formula){
               # seems necc. in order to catch and move on from err'd iterations
               skip_var <- FALSE
               # rerwip: this returns an object, probably a more elegant
               # solution hidden in that fact
               tryCatch({
                   # generate the norm counts, return the the different objects
                   # to a list `data.list`
                   formula.dir <- sub(' [+] ', '_', as.character(formula)[[2]])
                   data.list <- generate.norm.counts(metadata, formula, ref)
                   # run DESeq on the current formula, append the DESeq obj to
                   # the data list
                   data.list <- append(data.list,
                                list('deSeq' = run.de.seq(ref, 
                                                          data.list$dds_gene,
                                                          formula.dir)))
                   print('comparing that design to minimal design')
                   # use LRT to compare the current formula to the minimal
                   # ~condition design. H0 is that there is no difference
                   # between the models, sig padj indicates a difference in DE
                   # is detected in the more complex formula
                   form.compare <- DESeq(data.list$dds_gene, 
                                         test='LRT', 
                                         reduced=default.formula)
                   # transform DESeq object to DF
                   form.compare.res <- as.data.frame(results(form.compare))
                   print('# significant genes in LRT comparison to ~condition')
                   # count how many genes are sig diff between models
                   # rerwip: this may not be the best way to determine the
                   # 'best' model....
                   sig_gene_len <- length(which(form.compare.res$padj < 0.05))
                   #if (sig_gene_len > 50){
                   #print('more than 50 genes differ from ~condition, printing')
                   write.data(data.list, ref, formula.dir)
                   #}

                   # store the formula, the sig gene count, and all associated
                   # data objects in a row of the output data frame
                   formula.performance <<- formula.performance %>% 
                       add_row(tibble('formula' = as.character(formula)[[2]],
                                      'sig_genes' = sig_gene_len))
                   # this updates our skip value...
                   # rerwip: this probably can be refined, but has something to
                   # do with tryCatch() returning an object (?)
               }, error = function(e) {skip_var <<- TRUE}
               )
          }
    )

    print('summary of formula performance when compared to ~condition')
    print(formula.performance)
    write_csv(formula.performance %>% select(formula, sig_genes), 
              paste0(out.path, 
                     ref,
                     '.',
                     one,
                     '.',
                     two,
                     '.',
                     'formula.performance.summary.csv'))

}

### MAIN ###

# instantiate argparse objects
argparse()
# instantiate the salmon output directory names to search for data
ref_list <- c('gencode.34.salmon.out','te.locus.gencode.34.salmon.out')

formula.dir <- sub(' [+] ', '_', as.character(default.formula)[[2]])

if (!is.na(meta.file)){
# apply run and collect to both output directory names in the vector
    lapply(X = ref_list, FUN=run.and.collect)
} else{
    if (!dir.exists(file.path(out.path, formula.dir))){
        lapply(ref_list, function(ref){
            metadata <- parse.input(ref)
            data.list <- generate.norm.counts(metadata, default.formula, ref)
            print('running deSeq on first formula')
            # run DESeq on the current formula, append the DESeq obj to
            # the data list

            if (counts_only == TRUE){
                write.data(data.list, ref, formula.dir)

            } else{
                data.list <- append(data.list,
                            list('deSeq' = run.de.seq(ref, 
                                                      data.list$dds_gene,
                                                      formula.dir)))
                write.data(data.list, ref, formula.dir)
            }
        })

    } #else {
}

