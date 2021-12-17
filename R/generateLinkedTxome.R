#!/user/bin/env Rscript

# rreggiar@ucsc.edu

# [guide to generating with salmon](https://combine-lab.github.io/alevin-tutorial/2020/alevin-velocity/)

suppressPackageStartupMessages({
  library(tximeta)
  library(rjson)
})

paths.in <- paths.in <- scan(file=file("stdin", "r"), what="character", n=2)
index.dir <- "/public/groups/kimlab/indexes"

# if (length(paths.in) < 2) {
#   cat('please provide all arguments: [1] input directory [2] salmon version')
#   quit(status = 1)
# }

input.dir <- file.path(paths.in[1])
salmon.version <- paths.in[2]
gencode.version <- substring(input.dir, first = nchar(input.dir)-1, last = nchar(input.dir))


# salmon.version <- '1.3.0'
# gencode.version <- '35'
# file.path(index.dir, paste0("sel.align.gencode.v", gencode.version, ".process.aware.salmon.v",salmon.version,".sidx"))

linkedTxomeWrap <- function(indexType) {

	tximeta::makeLinkedTxome(
	    indexDir = file.path(index.dir, paste0("sel.align.gencode.v", 
	    	gencode.version, 
	    	'.',
	    	indexType,
	    	'.v',
	    	salmon.version,
	    	".sidx")),
	    source = "de-novo", genome = "GRCh38",
	    organism = "Homo sapiens", release = "Hg38",
	    fasta = file.path(input.dir, paste0('gencode.v', 
	    	gencode.version,
	    	'.', 
	    	indexType,
	    	'.fa')),
	    gtf = file.path(input.dir, paste0('gencode.v', 
	    	gencode.version,
	    	'.', 
	    	indexType,
	    	'.gtf')),
	    write = TRUE, 
	    jsonFile = file.path(input.dir, paste0('gencode.v', 
	    	gencode.version, 
	    	'.',
	    	indexType,
	    	'.json'))
	)
}


if (!file.exists(file.path(input.dir, 
	paste0('gencode.v', gencode.version, '.annotation.expanded.json')))){
	if(!file.exists(file.path(input.dir, 
		paste0('gencode.v', gencode.version, '.process.aware.json')))){

		print('process aware')
		indexType = 'process.aware.json'
		linkedTxomeWrap(indexType)
	}
}

if (!file.exists(file.path(input.dir, 
	paste0('gencode.v', gencode.version, '.salmon.json')))){

	print('salmon')

	indexType = 'salmon'
	linkedTxomeWrap(indexType)
}

if (!file.exists(file.path(input.dir, 
	paste0('gencode.v', gencode.version, '.ucsc.rmsk.salmon.json')))){

	print('ucsc rmsk salmon')

	indexType = 'ucsc.rmsk.salmon'
	linkedTxomeWrap(indexType)
}