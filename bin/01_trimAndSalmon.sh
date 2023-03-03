#!/bin/bash

# rreggiar@ucsc.edu


scriptName=$(basename $0)
if [ $# -lt 3 ]; then
    echo "error: usage $scriptName sampleDir adapterChoice salmonIndex "
    echo "example: $scriptName /scratch/kimlab/projects/exoRNA-biomarkers-panc/data/a549_0.2MOI_24hr \
    intra /public/groups/kimlab/indexes/sel.aln.gen.34.ucsc.rmsk.index.salmon.1.2.1 /public/groups/kimlab/projects/x/data/output_data/qc"
    exit 1
fi

sampleDir="$1"
adapterChoice="$2"
salmonIndex="$3"
qc_output="$4"

dateStamp="$(bash helper_dateStamp.sh)"

echo "cmd: "$@""

for inputDir in "$sampleDir"/*; do

	set -x 

	echo "$inputDir"

	./01a_trimmomaticRun.sh "$inputDir" "$adapterChoice"

	./01b_salmonRun.sh "$inputDir" "$salmonIndex" "$dateStamp"

	exitStatus=$?
	if [ $? -ne 0 ]; then

	    echo ERROR "$scriptName" "$inputDir" returned exit status "$exitStatus"
	    continue

	fi

	set +x

done

#./01c_fastqQualityCheck.sh "$sampleDir" "$qc_output"
