#!/user/bin/env Rscript


##################
suppressMessages({
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(tximport)
library(data.table)
library(readr)
})
##################

## $ echo 'paths.in' | Rscript tximport.and.normalize.R

paths.in <- scan(file=file("stdin", "r"), what="character")
paths <- paths.in[1]
tx2genepath <- '/public/groups/kimlab/genomes.annotations/gencode.32/gen.32.ucsc.rmsk.tx.2.gene.csv'
files <- dir(file.path(paths))
print(files)
names(files)<- files

sample.info <- data.frame('sample' = character(),
                          'condition' = character(),
                          'rep' = numeric())


sample.info.list <- lapply(files, function(x){
            x.split <- str_split(x, '[.]', n = 2)[[1]]
            condition <- x.split[1]
            rep <- x.split[2]
            temp.df <- data.frame('sample' = x,
                                  'condition' = condition,
                                  'rep' = rep)
            sample.info <- rbind(sample.info, temp.df)
            return(sample.info)
            return(temp.df)
})

print(sample.info.list)

sample.info.df <- 
    as.data.frame(bind_rows(sample.info.list)) #%>%
    #filter(grepl('patient', sample))

print(sample.info.df)

files <- 
    file.path(paths, sample.info.df$sample, 'te.locus.out.for.sleuth', 'quant.sf')


names(files) <- sample.info.df$sample

tx2gene <- read_csv(tx2genepath, col_names = F)
names(tx2gene) <- c("Name", "GeneID")

txi <- tximport(files, 
                type="salmon", 
                tx2gene = tx2gene,
                txOut = F)

print(head(txi$counts))

dds <- DESeqDataSetFromTximport(txi, colData = sample.info.df, 
                                design = ~condition)

dds <- estimateSizeFactors(dds)

write.csv(as.data.frame(counts(dds, normalized=T)),
          '/public/groups/kimlab/aale.kras/data/bulk.rna.seq/exotic/output/exo.te.aale.kras.v.ctrl.de-seq.counts.csv')
quit()
keep <- rowSums(counts(dds)) >= 10
dds.txi <- dds[keep, ]
dds.txi <- DESeq(dds.txi)
#dds.txi <- estimateSizeFactors(dds.txi)
#dds.txi <- estimateDispersions(dds.txi)
#dds.txi <- nbinomWaldTest(dds.txi, maxit=500)

resultsNames(dds.txi)

res <- results(dds.txi, name = 'condition_kras_vs_ctrl')
#res <- results(dds.txi, name = 'group')
#res <- results(dds.txi, name = 'diabetes')
resLFC <- lfcShrink(dds.txi, 
                coef = 'condition_kras_vs_ctrl', 
                #coef = 'group',
                #coef = 'diabetes',
                type = 'apeglm')

write.csv(as.data.frame(resLFC), 
          '/public/groups/kimlab/aale.kras/data/bulk.rna.seq/aale/output/gencode.te.locus.aale.kras.v.ctrl.de-seq.csv')

