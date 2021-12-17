#!/bin/bash

# rreggiar@ucsc.edu
# runs quick (hopefully) qc on fastq files (pref. trimmed to ensure success of adapter trimming)
# use multiqc output to assess; will aggregate into a nice webpage

scriptName=$(basename $0)
if [ $# -lt 2 ]; then
    echo "error: usage "$scriptName" sampleDir outputDir"
    echo "example "$scriptName" /scratch/kimlab/projects/exoRNA-biomarkers-panc/data/a549_0.2MOI_24hr"
    exit 1
fi

sampleDir="$1"

outputDir="$2"

dateStamp="$(bash dateStamp.sh)"

cmdList="$@"

echo "cmd: "$scriptName" "$cmdList""

fastqc -t 8 "$sampleDir"/*/*.fq.gz #--outdir "$outputDir"

multiqc "$sampleDir" --filename multiqc."$(basename "$inputDir")"."${dateStamp}" -o "$outputDir"
