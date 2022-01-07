summarizeToGene_RER <- function(object, varReduce=FALSE, ...) {

  tximeta::missingMetadata(object, summarize=TRUE)

  txomeInfo <- metadata(object)$txomeInfo
  txdb <- tximeta::getTxDb(txomeInfo)
  message("obtaining transcript-to-gene mapping from database")

  suppressMessages({
    tx2gene <- select(txdb, keys(txdb, "TXNAME"), "GENEID", "TXNAME")
  })

  # TODO what to do about warnings about out-of-bound ranges? pass along somewhere?
  suppressWarnings({
    g <- getRanges(txdb=txdb, txomeInfo=txomeInfo, type="gene", single.strand.only = FALSE)
  })

  # need to add seqinfo for GENCODE and RefSeq
  if (all(is.na(seqlengths(g)))) {
    seqinfo(g) <- seqinfo(object)
  }
  
  txi <- list(
    abundance=assays(object)[["abundance"]],
    counts=assays(object)[["counts"]],
    length=assays(object)[["length"]],
    countsFromAbundance="no"
  )

  txi.gene <- summarizeToGene(object=txi, tx2gene=tx2gene, varReduce=varReduce, ...)
  
  # put 'counts' in front to facilitate DESeqDataSet construction
  assays <- txi.gene[c("counts","abundance","length")]
  if (varReduce) {
    assays <- c(assays, txi.gene["variance"])
  } else if ("infRep1" %in% assayNames(object)) {
    infReps <- txi.gene$infReps
    assays <- c(assays, infReps)
  }

  # here do the same check/subset but with the gene-level
  # tximport assay matrices and gene ranges
  assays <- checkAssays2Txps(assays, g)

  # TODO give a warning here if there are genes in TxDb not in Salmon index?
  g <- g[rownames(assays[["counts"]])]

  # store txp IDS
  tx_ids <- CharacterList(split(tx2gene$TXNAME, tx2gene$GENEID))
  if (all(names(g) %in% names(tx_ids))) {
    tx_ids <- tx_ids[names(g)]
    mcols(g)$tx_ids <- tx_ids
  }
  
  # calculate duplication
  if ("hasDuplicate" %in% colnames(mcols(object))) {
    stopifnot(all(rownames(object) %in% tx2gene[,1]))
    t2g.o <- tx2gene[match(rownames(object),tx2gene[,1]),]
    has.dup.list <- LogicalList(split(mcols(object)$hasDuplicate, t2g.o$GENEID))
    mcols(g)$numDupObjects <- sum(has.dup.list)
  }
  
  metadata <- metadata(object)
  # stash countsFromAbundance value
  metadata$countsFromAbundance <- txi.gene$countsFromAbundance
  metadata$level <- "gene"
  
  gse <- SummarizedExperiment(assays=assays,
                              rowRanges=g,
                              colData=colData(object),
                              metadata=metadata)  
  gse
}
