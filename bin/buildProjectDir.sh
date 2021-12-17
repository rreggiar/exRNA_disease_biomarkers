#!/bin/bash

scriptName=$(basename $0)
if [ $# -lt 4 ]; then
  echo "error: usage $scriptName control_effect_name treatment_effect_name project_name input_meta_file(subset to relevant rows (for now))"
    echo "example $scriptName"
    exit 1
fi

inputDir="$1"
ctrl="$2"
treat="$3"
projName="$4"
metaData="$5"
# instantiate counters at 1 to avoid a zero sample number
ctrlCount=1
treatCount=1

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

scratchDir="$groupsProjDir"'/data/'"$sampleName"

if [[ ! -d "$scratchDir" ]]; then

  mkdir "$scratchDir"

fi

while IFS=, read run treatment; do 

  if [[ "$treatment" =~ "$treat" ]]; then

   echo "$run"."$treat"."$treatCount"
   sampleDir="$run"."$treat"."$treatCount"

   if [[ ! -d "$scratchDir"/"$sampleDir" ]]; then

     mkdir "$scratchDir"/"$sampleDir"

   fi
     
   ln -s $inputDir/"$run"*fastq.gz "$scratchDir"/"$sampleDir"/

   ((treatCount++))

  elif [[ "$treatment" =~ "$ctrl" ]]; then

   echo "$run"."$ctrl"."$ctrlCount"
   sampleDir="$run"."$ctrl"."$ctrlCount"

   if [[ ! -d "$scratchDir"/"$sampleDir" ]]; then

     mkdir "$scratchDir"/"$sampleDir"

   fi

   ln -s $inputDir/"$run"*fastq.gz "$scratchDir"/"$sampleDir"

   ((ctrlCount++))
  fi

 done<"${metaData}"
