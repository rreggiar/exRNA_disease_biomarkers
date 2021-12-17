#!/bin/bash

inputdir='/public/groups/kimlab/exoRNA-biomarkers-panc/data/panc.kras.intra'
outdir='/public/groups/kimlab/exoRNA-biomarkers-panc/output.data/bam.files.for.editing'

for subdir in $inputdir/*/star.out/; do

  cd $subdir

  sample=$(echo $subdir | cut -d'/' -f 8)

  echo $sample

  samtools view -S -b Aligned.out.sam > $sample.out.bam
  samtools sort -m 8G $sample.out.bam -o $sample.sorted.bam
  samtools index $sample.sorted.bam
  
  ln -s $subdir/$sample.sorted.bam $outdir/.

done
