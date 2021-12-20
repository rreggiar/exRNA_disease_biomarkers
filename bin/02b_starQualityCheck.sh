#!/bin/bash

# rreggiar@ucsc.edu
# runs quick (hopefully) qc on fastq files (pref. trimmed to ensure success of adapter trimming)
# use multiqc output to assess; will aggregate into a nice webpage

scriptName=$(basename $0)
if [ $# -lt 2 ]; then
    echo "error: usage "$scriptName" sampleDir outputDir"
    echo "example "$scriptName" /public/groups/kimlab/projects/exoRNA-biomarkers-panc/data/a549_0.2MOI_24hr /public/groups/kimlab/projects/exoRNA/exoRNA-biomarkers-panc/data/output_data_qc "
    exit 1
fi

sampleDir="$1"

outputDir="$2"

dateStamp="$(bash helper_dateStamp.sh)"

cmdList="$@"

echo "cmd: "$scriptName" "$cmdList""

multiqc ""$sampleDir"/hg38_star_2.7.9a_*_star_out/second_pass_out/*" --filename multiqc.star."$(basename "$sampleDir")"."${dateStamp}" -o "$outputDir"
