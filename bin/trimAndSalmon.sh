#!/bin/bash

# rreggiar@ucsc.edu


scriptName=$(basename $0)
if [ $# -lt 3 ]; then
    echo "error: usage $scriptName sampleDir adapterChoice salmonIndex "
    echo "example $scriptName /scratch/kimlab/projects/exoRNA-biomarkers-panc/data/a549_0.2MOI_24hr \
    intra /public/groups/kimlab/indexes/sel.aln.gen.34.ucsc.rmsk.index.salmon.1.2.1"
    exit 1
fi

sampleDir="$1"
adapterChoice="$2"
salmonIndex="$3"

dateStamp="$(bash dateStamp.sh)"

echo "cmd: "$@""

for inputDir in "$sampleDir"/*; do

	set -x 

	./trimmomaticRun.sh "$inputDir" "$adapterChoice"

	./salmonRun.sh "$inputDir" "$salmonIndex" "$dateStamp"

	exitStatus=$?
	if [ $? -ne 0 ]; then

	    echo ERROR "$scriptName" "$inputDir" returned exit status "$exitStatus"
	    continue

	fi

	set +x

done
