#!/bin/bash

# rreggiar@ucsc.edu
# 2021_06_01
# run alevin on single cell fastq with preferred arguments


scriptName=$(basename $0)

if [ $# -lt 3 ]; then
    echo "error: usage $scriptName sampleDir indexDir gencodeDir "
    echo "example $scriptName /public/groups/kimlab/aale.kras/data/single.cell.rna.seq/input/ctrl.fastq \
    /public/groups/kimlab/indexes \
    /public/groups/kimlab/genomes.annotations/gencode.35"
    exit 1
fi


function runAlevin(){

	indexDir=$1 #/public/groups/kimlab/indexes
	sampleDir=$2 #/public/groups/kimlab/aale.kras/data/single.cell.rna.seq/input/ctrl.fastq (*_R1_*, *_R2_*)
	gencodeDir=$3 #/public/groups/kimlab/genomes.annotations/gencode.35 -- has tx.to.gene.tsv
	outdir=$4 

	# for index (process_aware, TE, gencode) available
	for index in "$indexDir"/sel.align.gencode.v35.*; do

		# ex: alevin_ctrl.fastq_sel.align.gencode.v35.salmon.v1.3.0.sidx
		outPath="$outdir"/alevin_$(basename "${sampleDir}")_$(basename "${index}")

		# if alevin output quants don't exist already
		if [[ ! -f "$outPath"/alevin/quants_mat_cols.txt ]]; then

			# RERWIP -- hard coding thr possible tx2genes here, could do better
			echo "INDEX:""$index"
			if [[ $index != '/public/groups/kimlab/indexes/sel.align.gencode.v35.process.aware.salmon.v1.3.0.sidx' ]]; then
				tx2gene="${gencodeDir}"/gencode.v35.ucsc.rmsk.tx.to.gene.tsv
			else
				tx2gene="${gencodeDir}"/gencode.v35.annotation.expanded.tx.to.gene.tsv
			fi

			echo "TX TO GENE:""$tx2gene"

			# rec args
			salmon alevin -i "$index" \
				-l ISR \
				-1 `ls "${sampleDir}"/*_R1_*.gz` \
				-2 `ls "${sampleDir}"/*_R2_*.gz` \
				--chromium \
				-p 16 \
				--tgMap "$tx2gene" \
				-o "$outPath"
		fi

	done
}

### MAIN

sampleDir="$1"
indexDir="$2"
gencodeDir="$3"

#outdir="/public/groups/kimlab/aale-KRAS-G12-transformation/output.data"
outdir="$sampleDir"

runAlevin "${indexDir}" "${sampleDir}" "${gencodeDir}" "${outdir}"
