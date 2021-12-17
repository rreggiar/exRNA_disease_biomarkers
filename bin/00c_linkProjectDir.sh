#!/bin/bash

scriptName=$(basename $0)
if [ $# -lt 3 ]; then
  echo "error: usage $scriptName indputDir project_name input_meta_file(subset to relevant rows)"
    echo "example $scriptName"
    exit 1
fi

inputDir="$1"
projName="$2"
metaData="$3"

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

scratchDir="$groupsProjDir"/data/"$sampleName"

if [[ ! -d "$scratchDir" ]]; then

  mkdir "$scratchDir"

fi

while IFS=, read run name; do 

  echo "$run" -- "$name"
  sampleDir="$name"

  if [[ ! -d "$scratchDir"/"$sampleDir" ]]; then

   mkdir "$scratchDir"/"$sampleDir"

  fi
   
  ln -s $inputDir/"$run"*fastq.gz "$scratchDir"/"$sampleDir"/

done<"${metaData}"
