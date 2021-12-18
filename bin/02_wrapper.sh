#!/bin/bash

# array of jobs to execute the entire project's rna quant and qc workload

projectDir=$1 # /public/groups/kimlab/exRNA_disease_biomarkers/
twoPass="2pass=T"
edit="edit=T"

for inputDataPath in $projectDir/data/input_data/*; do
	
	starGenome='/public/groups/kimlab/genomes.annotations/HG.38/'

	# input args
	#sampleDir="$1"
	#starGenome="$2"
	#2pass="$3"
	#edit="$4"

	nohup ./02_iterStarRun.sh $inputDataPath \
		$starGenome \
		$twoPass \
		$edit \
		2>&1 > $projectDir/tmp/logs/$(basename $inputDataPath)_star_align_log.txt &

done


