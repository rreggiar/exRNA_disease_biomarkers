# build or load a TxDb/EnsDb for the dataset
getTxDb <- function(txomeInfo, useHub=TRUE) {

  import::here(.from = 'helper_bfcLocation.R', .directory = '/public/groups/kimlab/exRNA_disease_biomarkers/R/', .all = T)
  import::here(.from = 'S4Vectors', .all = T)
  import::here(.from = 'SummarizedExperiment', .all = T)
  import::here(.from = 'GenomicRanges', seqnames)
  import::here(.from = 'tximport', tximport, summarizeToGene)
  import::here(.from = 'AnnotationDbi', .all = T)
  import::here(.from = 'GenomicFeatures', .all = T)
  import::here(.from = 'BiocFileCache', .all = T)
  import::here(.from = 'GenomeInfoDb', .all = T)
  import::here(.from = 'tools', R_user_dir)

  # TODO what if there are multiple GTF files?
  stopifnot(length(txomeInfo$gtf) == 1)
  stopifnot(txomeInfo$gtf != "")
  txdbName <- basename(txomeInfo$gtf)
  bfcloc <- getBFCLoc()
  bfc <- BiocFileCache(bfcloc)
  # look up txdbName
  q <- bfcquery(bfc, txdbName)
  # then filter for equality with rname
  q <- q[q$rname==txdbName,]

  ### No TxDb was found in the BiocFilecache ###
  if (bfccount(q) == 0) {

    # Ensembl and GENCODE best case we can find database on AnnotationHub
    hubSources <- c("Ensembl","GENCODE")
    srcName <- txomeInfo$source
    hubWorked <- FALSE
    if (srcName %in% hubSources) {
      ensSrc <- srcName == "Ensembl"
      dbType <- if (ensSrc) "EnsDb" else "TxDb"
      if (useHub) {
        message(paste("useHub=TRUE: checking for", dbType, "via 'AnnotationHub'"))
        ah <- AnnotationHub()
        # get records
        records <- query(ah, c(srcName, txomeInfo$organism, txomeInfo$release))
        # confirm source, organism, dbType through metadata columns
        records <- records[records$dataprovider==srcName &
                           records$species==txomeInfo$organism &
                           records$rdataclass==dbType,]        
        if (ensSrc) {
          # Confirm release number through grep on the title
          # EnsDb record titles look like "Ensembl 123 EnsDb for Homo sapiens"
          records <- records[grepl(paste(srcName, txomeInfo$release, dbType), records$title),]
        } else {
          # Narrow records based on the genome coordinates
          # GENCODE record titles look like "TxDb for Gencode v123 on hg38 coordinates"
          coords <- genome2UCSC(txomeInfo$genome)
          records <- records[grepl(coords, records$title),]
        }
        if (length(records) == 1) {
          message(paste("found matching", dbType, "via 'AnnotationHub'"))
          hubWorked <- TRUE
          txdb <- ah[[names(records)]]
          bfcadd(bfc, rname=txdbName, fpath=dbfile(dbconn(txdb)))
        } else {
          message(paste("did not find matching", dbType, "via 'AnnotationHub'"))
        }
      }
      # if check on AnnotationHub failed (or wasn't attempted)
      if (!hubWorked) {
        # build db for Ensembl
        if (ensSrc) {
          message("building EnsDb with 'ensembldb' package")
          # split code based on whether linkedTxome (bc GTF filename may be modified)
          if (!txomeInfo$linkedTxome) {
            # TODO what about suppressing all these warnings
            suppressWarnings({
              savepath <- ensDbFromGtf(
                txomeInfo$gtf,
                outfile = bfcnew(bfc, rname=txdbName, ext=".sqlite")
              )
            })
          } else {
            message("NOTE: linkedTxome with source='Ensembl', ensembldb will be used to parse GTF.
this may produce errors if the GTF is not from Ensembl, or has been modified")
            # for linkedTxome, because the GTF filename may be modified
            # we manually provide organism, genomeVersion, and version
            suppressWarnings({
              savepath <- ensDbFromGtf(
                txomeInfo$gtf,
                outfile = bfcnew(bfc, rname=txdbName, ext=".sqlite"),
                organism = txomeInfo$organism,
                genomeVersion = txomeInfo$genome,
                version = txomeInfo$release
              )
            })
          }
          txdb <- EnsDb(savepath)
        }
      }
    }

    # two cases left:
    # 1) Neither Ensembl or GENCODE source
    # 2) GENCODE source but AHub didn't work
    if ((!srcName %in% hubSources) | (srcName == "GENCODE" & !hubWorked)) {
      message("building TxDb with 'GenomicFeatures' package")
      txdb <- makeTxDbFromGFF(txomeInfo$gtf)
      saveDb(
        txdb,
        file = bfcnew(bfc, rname=txdbName, ext=".sqlite")
      )
    }

  } else {
    ### Yes, TxDb was found in the BiocFilecache ###
    loadpath <- bfcrpath(bfc, txdbName)
    if (txomeInfo$source == "Ensembl") {
      message(paste("loading existing EnsDb created:",q$create_time[1]))
      txdb <- EnsDb(loadpath)
    } else {
      message(paste("loading existing TxDb created:",q$create_time[1]))
      txdb <- loadDb(loadpath)
    }
  }
  
  txdb
}


checkAssays2Txps <- function(assays, txps) {
  assay.nms <- rownames(assays[["counts"]])
  txps.missing <- !assay.nms %in% names(txps)
  if (!all(assay.nms %in% names(txps))) {

    # it's probably ok that the messaging here uses the term 'txps',
    # because it's unlikely that we'd have genes not present in the GTF
    # which nevertheless had txps in the GTF...

    if (all(!assay.nms %in% names(txps))) {
      stop("none of the transcripts in the quantification files are in the GTF")
    } else {

      if (sum(txps.missing) > 3) {
        example.missing <- paste0("Example missing txps: [",
                                  paste(head(assay.nms[txps.missing],3),collapse=", "),
                                  ", ...]")
      } else {
        example.missing <- paste0("Missing txps: [",
                                  paste(assay.nms[txps.missing],collapse=", "), "]")
      }

      # TODO what to do here, GTF is missing some txps in FASTA for Ensembl
      warning(paste0("
Warning: the annotation is missing some transcripts that were quantified.
", sum(txps.missing), " out of ", nrow(assays[["counts"]]),
" txps were missing from GTF/GFF but were in the indexed FASTA.
(This occurs sometimes with Ensembl txps on haplotype chromosomes.)
In order to build a ranged SummarizedExperiment, these txps were removed.
To keep these txps, and to skip adding ranges, use skipMeta=TRUE
", example.missing, "
"))

      # after warning, then subset
      for (nm in names(assays)) {
        assays[[nm]] <- assays[[nm]][!txps.missing,,drop=FALSE]
      }

    }
  }
  assays
}

getRanges_RER <- function(txdb=txdb, txomeInfo=txomeInfo, type=c("txp","exon","cds","gene")) {

  import::here(.from = 'helper_bfcLocation.R', .directory = '/public/groups/kimlab/exRNA_disease_biomarkers/R/', .all = T)
  import::here(.from = 'S4Vectors', .all = T)
  import::here(.from = 'SummarizedExperiment', .all = T)
  import::here(.from = 'GenomicRanges', seqnames)
  import::here(.from = 'tximport', tximport, summarizeToGene)
  import::here(.from = 'AnnotationDbi', .all = T)
  import::here(.from = 'GenomicFeatures', .all = T)
  import::here(.from = 'BiocFileCache', .all = T)
  import::here(.from = 'GenomeInfoDb', .all = T)
  import::here(.from = 'tools', R_user_dir)

  long <- c(txp="transcript",exon="exon",cds="CDS",gene="gene")
  stopifnot(length(txomeInfo$gtf) == 1)
  stopifnot(txomeInfo$gtf != "")

  # TODO the entry in the BiocFileCache assumes that the GTF/GFF file
  # has a distinctive naming structure... works for GENCODE/Ensembl/RefSeq 
  rngsName <- paste0(type,"Rngs-",basename(txomeInfo$gtf))
  
  bfcloc <- getBFCLoc()
  bfc <- BiocFileCache(bfcloc)
  q <- bfcquery(bfc, rngsName)
  if (bfccount(q) == 0) {
    # now generate ranges
    message(paste("generating",long[type],"ranges"))
    # TODO what to do about warnings about out-of-bound ranges? pass along somewhere?

    if (type == "txp") {
      ################
      ## txp ranges ##
      ################

      if (txomeInfo$source == "Ensembl") {
        suppressWarnings({
          rngs <- transcripts(txdb)
        })
      } else {
        suppressWarnings({
          rngs <- transcripts(txdb, columns=c("tx_id","gene_id","tx_name"))
        })
      }
      names(rngs) <- rngs$tx_name
      # dammit de novo transcript annotation will have
      # the transcript names as seqnames (seqid in the GFF3)
      if (tolower(txomeInfo$source) == "dammit") {
        names(rngs) <- seqnames(rngs)
      }
    } else if (type == "exon") {
      #################
      ## exon ranges ##
      #################

      # TODO suppress warnings about out-of-bound ranges for now... how to pass this on
      suppressWarnings({
        rngs <- exonsBy(txdb, by="tx", use.names=TRUE)
      })
    } else if (type == "cds") {
      #################
      ## CDS ranges ##
      #################

      # TODO suppress warnings about out-of-bound ranges for now... how to pass this on
      suppressWarnings({
        rngs <- cdsBy(txdb, by="tx", use.names=TRUE)
      })
    } else if (type == "gene") {
      #################
      ## gene ranges ##
      #################

      message('using RER version of genes()')

      # TODO suppress warnings about out-of-bound ranges for now... how to pass this on
      suppressWarnings({
        rngs <- genes(txdb, single.strand.genes.only=FALSE)
      })
    }
    savepath <- bfcnew(bfc, rngsName, ext=".rds")
    saveRDS(rngs, file=savepath)
  } else {
    loadpath <- bfcrpath(bfc, rngsName)
    message(paste("loading existing",long[type],"ranges created:",q$create_time[1]))
    rngs <- readRDS(loadpath)
  }
  rngs
}



summarizeToGene_RER <- function(object, varReduce=FALSE, ...) {

  #missingMetadata(object, summarize=TRUE)

  import::here(.from = 'helper_bfcLocation.R', .directory = '/public/groups/kimlab/exRNA_disease_biomarkers/R/', .all = T)
  import::here(.from = 'S4Vectors', .all = T)
  import::here(.from = 'SummarizedExperiment', .all = T)
  import::here(.from = 'AnnotationDbi', .all = T)
  import::here(.from = 'GenomicRanges', seqnames)
  import::here(.from = 'tximport', tximport, summarizeToGene)
  import::here(.from = 'GenomicFeatures', .all = T)
  import::here(.from = 'BiocFileCache', .all = T)
  import::here(.from = 'GenomeInfoDb', .all = T)
  import::here(.from = 'tools', R_user_dir)


  txomeInfo <- metadata(object)$txomeInfo
  txdb <- getTxDb(txomeInfo)
  message("obtaining transcript-to-gene mapping from database")

  suppressMessages({
    tx2gene <- select(txdb, keys(txdb, "TXNAME"), "GENEID", "TXNAME")
  })

  # TODO what to do about warnings about out-of-bound ranges? pass along somewhere?
  suppressWarnings({
    g <- getRanges_RER(txdb=txdb, txomeInfo=txomeInfo, type="gene")
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
  if ("infRep1" %in% assayNames(object)) {
    infReps <- list(assays(object)[grep("infRep", assayNames(object))])
    if (varReduce) {
      # split from per replicate list into per sample list
      # (this is what tximport expects for varReduce=TRUE)
      infReps <- list(splitInfReps(infReps[[1]]))
    }
    txi <- c(txi, infReps=infReps)
  } else {
    if (varReduce) stop("cannot calculate inferential variance without inferential replicates")
  }

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
