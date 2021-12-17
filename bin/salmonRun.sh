#!/bin/bash

# rreggiar@ucsc.edu
# conda -- aale.analysis.env

## !!install salmon via conda!!
## $ conda config --add channels conda-forge
## $ conda install -c bioconda salmon
## $ salmon --version # double check you're at the current version

# run salmon on trimmed fastq files to quantify abundance of transcripts against provided 
# txome index. Resulting quant.sf files are used in DESeq

scriptName=$(basename $0)
if [ $# -lt 3 ]; then
    echo "error: usage $scriptName inputDir salmonIndex dateStamp"
    echo "example $scriptName {ctrl,kras}.{1,2,3..} /path/to/{name.of.salmon.index} pipeline_dateStamp"
    exit 1
fi

dateStamp="$3"
# set -x
echo "cmd: "$@""
echo "salmon version:" $(salmon --version)
# echo "time: $dateStamp"
# set +x 

# echo $@ -- print everything from cmdline

inputDir="$1"
salmonIndex="$2"
outputDir="$(basename "$inputDir")"_$(basename "$salmonIndex")_"$dateStamp"

# set -x
echo "input: $inputDir"
echo "output: $outputDir"
# set +x 

function runSalmon() {
	# runs salmon on one sample and outputs to that directory
	inputDir="$1"
	salmonIndex="$2"
	outputDir="$3"
	outputPath="$inputDir"/"$outputDir"

  	outpathGrob=$(echo "$inputDir"/*"$(basename "$salmonIndex")"*)

  	echo "$outpathGrob"
	
	set -x

  if [[ ! -f "$(find "$outpathGrob" -name 'quant.sf')" ]]; then

		mkdir "$outputPath"

		if [ -f "$inputDir"/"$(basename "$inputDir")"_output_single_end.fq.gz ]; then

			trim_read="$inputDir"/"$(basename "$inputDir")"_output_single_end.fq.gz

			salmon quant \
				-i "$salmonIndex" \
				--libType A \
				-r "$trim_read" \
				-p 8 \
				--validateMappings \
				--gcBias \
				--seqBias \
			  	--recoverOrphans \
			  	--rangeFactorizationBins 4 \
				--output "$outputPath" 

			exitStatus=$?
			if [ $? -ne 0 ]; then

			    echo ERROR salmon "$inputDir" returned exit status "$exitStatus"
			    continue

			fi

			set +x

		else

			trim_fwd="$inputDir"/"$(basename "$inputDir")"_output_forward_paired.fq.gz
			trim_rev="$inputDir"/"$(basename "$inputDir")"_output_reverse_paired.fq.gz


			salmon quant \
				-i "$salmonIndex" \
				--libType A \
				-1 "$trim_fwd" \
				-2 "$trim_rev" \
				-p 8 \
				--validateMappings \
				--gcBias \
				--seqBias \
			    --recoverOrphans \
			    --rangeFactorizationBins 4 \
				--output "$outputPath" 

			exitStatus=$?
			if [ $? -ne 0 ]; then

			    echo ERROR salmon "$inputDir" returned exit status "$exitStatus"
			    continue

			fi
		fi

		set +x
		
		echo SKIP salmon "$inputDir" already run for "$salmonIndex"
   fi 
}

runSalmon "${inputDir}" "${salmonIndex}" "${outputDir}"

