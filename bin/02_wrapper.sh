#!/bin/bash

# array of jobs to execute the entire project's rna quant and qc workload

projectDir=${1:-'${PWD/*}'} # /public/groups/kimlab/exRNA_disease_biomarkers/
twoPass=${2:-'2pass=F'}
edit=${3:-'edit=F'}

for inputDataPath in $projectDir/data/input_data/*; do
	
	starGenome='/public/groups/kimlab/genomes.annotations/hg38_star_2.7.9a/'

	# input args
	#sampleDir="$1"
	#starGenome="$2"
	#2pass="$3"
	#edit="$4"

	nohup ./02_iterStarRun.sh $inputDataPath \
		$starGenome \
		$twoPass \
		$edit \
		$projectDir \
		2>&1 > $projectDir/tmp/logs/$(basename $inputDataPath)_star_align_log.txt &

	#nohup ./02b_starQualityCheck.sh "${inputDataPath}" "/public/groups/kimlab/exRNA_disease_biomarkers/data/output_data/rna_qc/star/" \
	#	2>&1 > $projectDir/tmp/logs/$(basename $inputDataPath)_star_qc_log.txt &

done


