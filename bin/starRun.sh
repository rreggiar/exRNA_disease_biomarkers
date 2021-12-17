#!/bin/bash

# rreggiar@ucsc.edu
# conda -- aale.analysis.env

## !!install STAR via conda!!
## $ conda config --add channels conda-forge
## $ conda install -c bioconda STAR
## $ STAR --version # double check you're at the current version

# run salmon on trimmed fastq files to quantify abundance of transcripts against provided 
# txome index. Resulting quant.sf files are used in DESeq

scriptName=$(basename $0)
if [ $# -lt 5 ]; then
    echo "error: usage $scriptName inputDir 2pass={T/F} edit={T/F} starGenome dateStamp"
    echo "example $scriptName {ctrl,kras}.{1,2,3..} /path/to/{name.of.star.genome} 2pass=F edit=T dateTime"
    echo "note: 2pass mode requires access to the entire experiments first pass output, should not be run until that has been generated"
    exit 1
fi

dateStamp="$5"
#set -x
echo "script: $scriptName"
echo "STAR version:" $(STAR --version)
echo "time: $dateStamp"
#set +x 

inputDir="$1"
twoPassArg="$2"
editArg="$3"
starGenome="$4"
firstPassDir="$inputDir"/$(basename "$starGenome")_${dateStamp}_star_out
rnaEditDir="$inputDir"/rna_editing_$(basename "$starGenome")_${dateStamp}_star_out
secondPassDir=$firstPassDir/second_pass_out

#set -x
echo "input:" "$inputDir"
echo "firstPassDir:" "$firstPassDir"
echo "secondPassDir:" "$secondPassDir"
#set +x 

## activate correct env
#function condaCheck() {
	# source the conda script so this shell has access
#	source /public/groups/kimlab/.install_bin/anaconda3/etc/profile.d/conda.sh

#	reqEnv="aale.analysis.env"
#	env=$(basename "$CONDA_PREFIX")

#	if [[ env != reqEnv ]]; then
#		echo "switching from "$env" to "$reqEnv""
#		conda activate $reqEnv
#	else
#		echo ""$reqEnv" is active"
#	fi
#}

function runStarFirstPass() {
	inputDir="$1"
	firstPassDir="$2"
	starGenome="$3"

	if [ ! -f "$firstPassDir"/Aligned.out.sam ]; then

		mkdir "$firstPassDir"

		trim_fwd=`ls "$inputDir"/*_output_forward_paired.fq.gz`
		trim_rev=`ls "$inputDir"/*_output_reverse_paired.fq.gz`

        STAR --genomeDir "$starGenome" \
	        --readFilesIn <(gunzip -c "$trim_fwd") <(gunzip -c "$trim_rev") \
	        --runThreadN 8 \
	        --outMultimapperOrder Random \
	        --outFilterMultimapNmax 50 \
	        --outFileNamePrefix "$firstPassDir"/ 

    fi
}

function runStarSecondPass() {

	inputDir="$1"
	secondPassDir="$2"
	starGenome="$3"

	cd "$inputDir"
	cd ../

	# starMasterDir == "$inputDir"/..

	starMasterDir=$PWD

	for inputSample in $PWD/*; do

		cd "$inputSample"

		if [ ! -d "$secondPassDir" ]; then

			mkdir "$secondPassDir"

		  trim_fwd=`ls "$inputDir"/*_output_forward_paired.fq.gz`
		  trim_rev=`ls "$inputDir"/*_output_reverse_paired.fq.gz`

		  STAR --genomeDir "$starGenome" \
			  --readFilesIn <(gunzip -c "$trim_fwd") <(gunzip -c "$trim_rev") \
			  --runThreadN 8 \
			  --sjdbFileChrStartEnd "$starMasterDir"/*/*_star_out/SJ.out.tab \
			  --outFilterMultimapNmax 50 \
			  --outReadsUnmapped Fastx \
			  --outMultimapperOrder Random \
			  --outSAMtype BAM SortedByCoordinate \
			  --outFileNamePrefix "$secondPassDir"/
    fi

  done
}

function runStarForRNAEditing() {
	inputDir="$1"
	rnaEditDir="$2"
	starGenome="$3"

	if [ ! -f "$rnaEditDir"/*.bam ]; then

		mkdir "$rnaEditDir"

		trim_fwd=`ls "$inputDir"/*_output_forward_paired.fq.gz`
		trim_rev=`ls "$inputDir"/*_output_reverse_paired.fq.gz`

        STAR --genomeDir $starGenome \
	        --readFilesIn <(gunzip -c $trim1) <(gunzip -c $trim2) \
	        --outFilterMatchNminOverLread 0.95 \
	        --outSAMtype BAM SortedByCoordinate \
	        --outFileNamePrefix "$rnaEditDir"/

    fi   
}

#condaCheck

runStarFirstPass "${inputDir}" "${firstPassDir}" "${starGenome}"

if [[ "$twoPassArg" == '2pass=T' ]]; then
	runStarSecondPass "${inputDir}" "${secondPassDir}" "${starGenome}"
fi

if [[ "$editArg" == "edit=T" ]]; then
	runStarForRNAEditing "${inputDir}" "${rnaEditDir}" "${starGenome}"
fi




