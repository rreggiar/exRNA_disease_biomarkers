#!/bin/bash

# array of jobs to execute the entire project's rna quant and qc workload

projectDir=$1 # /public/groups/kimlab/exRNA_disease_biomarkers/
indexType=$2 # ucsc.rmsk.salmon
indexVer=$3 # 39
salmonVersion=`salmon -v | cut -d' ' -f2`
index=`ls -d /public/groups/kimlab/indexes/*.v$indexVer.$indexType.v$salmonVersion.sidx`

echo $index

for inputDataPath in $projectDir/data/input_data/*; do
	
	adapterType='exo'
	
	if [[ $inputDataPath =~ 'intra' ]]; then adapterType='intra'; fi

	nohup ./01_trimAndSalmon.sh $inputDataPath \
		$adapterType \
		$index \
		$projectDir/data/output_data/rna_qc/$indexType \
		2>&1 > $projectDir/tmp/logs/$(basename $inputDataPath)_"$indexType"_quant_log.txt &

done


