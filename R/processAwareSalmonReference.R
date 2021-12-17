#!/user/bin/env Rscript

# rreggiar@ucsc.edu

# [guide to generating with salmon](https://combine-lab.github.io/alevin-tutorial/2020/alevin-velocity/)

suppressPackageStartupMessages({
    library(Biostrings)
    library(eisaR)
    library(GenomicFeatures)
    library(BSgenome)
})

paths.in <- paths.in <- scan(file=file("stdin", "r"), what="character", n=1)

if (length(paths.in) < 1) {
    cat('please provide all arguments: [1] input directory')
    quit(status = 1)
}
 
input.dir <- paths.in[1]
# input.dir <- '/public/groups/kimlab/genomes.annotations/gencode.35'
gencode.version <- substring(input.dir, first = nchar(input.dir)-1, last = nchar(input.dir))

cat('gtf')
gtf <- file.path(input.dir, paste0('gencode.v', gencode.version, '.annotation.gtf.gz'))
cat('grl')
grl <- eisaR::getFeatureRanges(
    gtf = gtf,
    featureType = c("spliced", "intron"), 
    intronType = "separate", 
    flankLength = 90L, 
    joinOverlappingIntrons = FALSE, 
    verbose = TRUE
)

# cat('genome')
# genome <- Biostrings::readDNAStringSet(
#     file.path(input.dir, 'GRCh38.primary_assembly.genome.fa.gz')
# )

# cat('seqs')
# names(genome) <- sapply(strsplit(names(genome), " "), .subset, 1)
# seqs <- GenomicFeatures::extractTranscriptSeqs(
#     x = genome, 
#     transcripts = grl
# )

# cat('fasta')
# Biostrings::writeXStringSet(
#     seqs, filepath = file.path(input.dir, paste0('gencode.v', gencode.version, '.annotation.expanded.fa'))
# )
# cat('gtf out')
# eisaR::exportToGtf(
#     grl, 
#     filepath = file.path(input.dir, paste0('gencode.v', gencode.version, '.annotation.expanded.gtf'))
# )

write.table(
    metadata(grl)$corrgene, 
    file = file.path(input.dir, paste0('gencode.v', gencode.version, '.annotation.expanded.tsv')),
    row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t"
)

df <- eisaR::getTx2Gene(
    grl, filepath = file.path(input.dir, paste0('gencode.v', gencode.version, '.annotation.expanded.tx.to.gene.tsv'))
)