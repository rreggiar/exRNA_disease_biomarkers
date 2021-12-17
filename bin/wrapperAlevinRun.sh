#!/bin/bash

# rreggiar@ucsc.edu


scriptName=$(basename $0)
if [ $# -lt 3 ]; then
    echo "error: usage $scriptName sampleDir indexDir gencodeDir"
    echo "example $scriptName /scratch/kimlab/projects/project/data/single_cell \
     /public/groups/kimlab/indexes \ 
     /public/groups/kimlab/genomes.annotations/gencode.35"
    exit 1
fi

sampleDir="$1"
indexDir="$2"
gencodeDir="$3"

dateStamp="$(bash dateStamp.sh)"

echo "cmd: "$@""

for inputDir in "$sampleDir"/*; do

	set -x 

  ./runAlevin.sh "$inputDir" "$indexDir" "$gencodeDir"

	exitStatus=$?
	if [ $? -ne 0 ]; then

	    echo ERROR "$scriptName" "$inputDir" returned exit status "$exitStatus"
	    continue

	fi

	set +x

done
