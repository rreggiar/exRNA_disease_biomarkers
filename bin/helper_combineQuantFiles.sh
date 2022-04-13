#!/bin/bash


#inputPath='/public/groups/kimlab/seqData/PANC_TCGA_TE-AWARE_COUNTS/roman-panc/TCGA_PAAD_ControlledAccess_V1-0_DATA_edu_ucsc_kim_lab/quant.sf'

inputPath=$1
inputQuant=$2
outputPath=$3
countField=$4

tempQuant="$inputPath"/"$inputQuant"

cut -f 1 <(gzip -cd "$tempQuant") > quant_tmp.txt

for quantFile in "$inputPath"/*; do

	thisQuant=$(basename -s .quant.sf.gz "$quantFile")

	echo "$thisQuant"

	paste <(cat quant_tmp.txt) \
		<(gzip -cd "$quantFile" \
		| cut -f$countField \
		| sed  "1s/.*/"$thisQuant"/") > tmp && mv tmp quant_tmp.txt

done

rm tmp

mv quant_tmp.txt "$outputPath"
