#!/bin/bash

# rreggiar@ucsc.edu
# conda -- aale.analysis.env

## !!install Trimmomatic via conda!!
## $ conda config --add channels conda-forge
## $ conda install -c bioconda trimmomatic
## $ trimmomatic -version # double check you're at the current version

# run trimmomatic to remove adapter sequences from the ends of reads
# the resulting paired output files will be used for all downstream analysis

scriptName=$(basename $0)
if [ $# -lt 2 ]; then
    echo "error: usage $scriptName  sample.input.dir adapterChoice"
    echo "example: $scriptName /path/to/{ctrl,kras}.{1,2,3..} {exo, intra}"
    exit 1
    # rerwip: Make the outDir argument make sense for symlinked databases
fi

# this is a constant 
adapterPath='/public/groups/kimlab/genomes.annotations/adapters'
# which adapters we use will depend on the library preparation
# the main situations this is relevant to are the different intracellular
# and extracellular preparations. This allows us to specify either of those as
# a shortcut or provide a custom argument that will be checked against our
# adapter library
if [[ "$2" == 'intra' ]]; then
	adapterChoice='TruSeq3-PE.fa'
elif [[ "$2" == 'exo' ]]; then
	adapterChoice='Nextera_pe_w_illumina_pe_universal_adapters.fa'
else
	if [[ -f "$adapterPath"/"$2" ]]; then
			adapterChoice="$2"
	else
		echo "you have provided an adapter file that doesn't exist (yet?): "$2""
		exit 1
	fi
fi

## don't want dateStamp functionality here
# dateStamp="$3"
# set -x
# echo "script: $scriptName"
# echo "time: $dateStamp"
# set +x 

inputDir="$1"

#set -x
echo "input: $inputDir"
echo "trimmomatic:" $(trimmomatic -version)
#set +x 

## activate correct env
# this will only run if you happen to be in the wrong env
# I think using directories as envNames could make this moot
# OR allow us to rely on "$(basename $PWD)" as $reqENV which would make the check
# super portable
# rerwip -- move this to util script
# function condaCheck() {
# 	# source the conda script so this shell has access
# 	source /public/groups/kimlab/.install_bin/anaconda3/etc/profile.d/conda.sh

# 	reqEnv="aale.analysis.env"
# 	env=$(basename "$CONDA_PREFIX")

# 	if [[ env != reqEnv ]]; then
# 		echo "switching from "$env" to "$reqEnv""
# 		conda activate $reqEnv
# 	else
# 		echo ""$reqEnv" is active"
# 	fi
# }

## since trimming is deterministic and reproducible, I'm going to treat these as un-dated
## output files for now so they can be accessed by various tools/pipelines and not be needlessly
## reproduced

## currently, this runs assuming you want output sent to $PWD, as in we are iterating over
## input dir's. This is somthing I need to talk to aedavids about...
function runTrimmomatic() {
	inputDir="$1"
	adapterPath="$2"
	adapterChoice="$3"

	# check for existence of canonical output file name
	if [[ ! -f "$inputDir"/$(basename "$inputDir")_output_forward_paired.fq.gz ]]; then
		# look for second paired read, if that exists approach as paired library
		if [ -f "$inputDir"/*_R2_001.fastq.gz ]; then
			#statements
			# if the read files exist and are symlinks
			if [ -L "$inputDir"/*_R1_001.fastq.gz ]; then
				# assign the reads  -- prepping for move to $SCRATCH data directories
				read_1=$(readlink -f "$inputDir"/*R1_001.fastq.gz)
				read_2=$(readlink -f "$inputDir"/*R2_001.fastq.gz)
			else
				read_1="$inputDir"/*R1_001.fastq.gz
				read_2="$inputDir"/*R2_001.fastq.gz
			fi
			
			set -x
			echo "read 1:" "$read_1"
			echo "read 2:" "$read_2"
			echo "Trimming fastq"

			# run trimm with recommended arguments
			trimmomatic PE  -threads 8 \
				"$read_1" "$read_2" \
				"$inputDir"/$(basename "$inputDir")_output_forward_paired.fq.gz \
				"$inputDir"/$(basename "$inputDir")_output_forward_unpaired.fq.gz \
				"$inputDir"/$(basename "$inputDir")_output_reverse_paired.fq.gz \
				"$inputDir"/$(basename "$inputDir")_output_reverse_unpaired.fq.gz \
				ILLUMINACLIP:"$adapterPath"/"$adapterChoice":1:30:10:4:true \
				LEADING:3 \
				TRAILING:3 \
				SLIDINGWINDOW:4:15 \
				MINLEN:36
			
			# remove the unpaired reads
			rm "$inputDir"/$(basename "$inputDir")*unpaired.fq.gz

			exitStatus=$?
			if [ "$exitStatus" -ne 0 ]; then

			    echo ERROR trimmomatic "$inputDir" returned exit status "$exitStatus"
			    continue

			fi

		elif [ ! -f "$inputDir"/$(basename "$inputDir")_output_single_end.fq.gz ]; then
			#statements
			# if read-pair not found, use single end arguments

			if [ -L "$inputDir"/*_1.fastq.gz ]; then
				# assign the reads  -- prepping for move to $SCRATCH data directories
				read_1=$(readlink "$inputDir"/*_1.fastq.gz)
				#read_2=$(readlink "$inputDir"/*_R2_001.fastq.gz)
			else
				read_1="$inputDir"/*_1.fastq.gz
				#read_2="$inputDir"/*_2*.fastq.gz
			fi
			
			set -x
			echo "read 1:" "$read_1"
			# echo "read 2:" "$read_2"
			echo "Trimming fastq"

			# run trimm with recommended arguments
			trimmomatic SE  -threads 8 \
				"$read_1"  \
				"$inputDir"/$(basename "$inputDir")_output_single_end.fq.gz \
				ILLUMINACLIP:"$adapterPath"/"$adapterChoice":1:30:10:4:true \
				LEADING:3 \
				TRAILING:3 \
				SLIDINGWINDOW:4:15 \
				MINLEN:36


			exitStatus=$?
			if [ "$exitStatus" -ne 0 ]; then

			    echo ERROR trimmomatic "$inputDir" returned exit status "$exitStatus"
			    continue

			fi	
		fi

		set +x

	fi
}

# condaCheck

runTrimmomatic "${inputDir}" "${adapterPath}" "${adapterChoice}" # specify file path

