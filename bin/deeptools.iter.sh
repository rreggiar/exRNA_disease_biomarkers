#!/bin/bash

inputdir='/public/groups/kimlab/exoRNA-biomarkers-panc/data/'
outdir='/public/groups/kimlab/exoRNA-biomarkers-panc/output.data/bw.files.for.uscs.gb/'

for subdir in $inputdir/*/*/star.out/pass.2; do

  cd $subdir

  sample=$(echo $subdir | cut -d'/' -f 9)

  echo $sample

  samtools index Aligned.sortedByCoord.out.bam

  bamCoverage -b Aligned.sortedByCoord.out.bam -o $outdir$sample.bw

done
