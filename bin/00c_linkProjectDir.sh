#!/bin/bash

scriptName=$(basename $0)
if [ $# -lt 3 ]; then
  echo "error: usage $scriptName inputDir project_name input_meta_file(subset to relevant rows)"
    echo "example: "$scriptName" /public/groups/kimlab/seqData/2020-09-30_exoRNABiomarkersPancAndCovid/bioIvt_covid\\
	    exRNA_disease_biomarkers\\
	    /public/groups/kimlab/seqData/2020-09-30_exoRNABiomarkersPancAndCovid/bioIvt_covid/ROSTER.csv\\
	    /public/groups/kimlab/exRNA_disease_biomarkers/data/input_data/bioIvt_covid"
    exit 1
fi

inputDir="$1"
projName="$2"
metaData="$3"
destinationDir="$4"

sampleName="$(basename "$inputDir")"

# laborious way of checking and generating the output directories in scratch
# I'll fix this when I'm less tired -- rerwip
groupsProjDir='/public/groups/kimlab/'"$projName"

if [[ ! -d "$groupsProjDir" ]]; then

  mkdir "$groupsProjDir"

fi

if [[ ! -d "$groupsProjDir"/data ]]; then

  mkdir "$groupsProjDir"/data

fi


if [[ ! -d "$destinationDir" ]]; then

  mkdir "$destinationDir"

fi

while IFS=, read run name; do 

  echo "$run" -- "$name"
  sampleDir="$name"

  if [[ ! -d "$destinationDir"/"$sampleDir" ]]; then

   mkdir "$destinationDir"/"$sampleDir"

  fi
   
  ln -s $inputDir/"$run"*fastq.gz "$destinationDir"/"$sampleDir"/

done<"${metaData}"
