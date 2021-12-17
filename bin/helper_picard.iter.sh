#!/bin/bash

inputdir='/public/groups/kimlab/exoRNA-biomarkers-panc/data/'

for subdir in $inputdir/*; do
  for sample in $subdir/*; do
    cd $sample
    name=$(basename $PWD)
    echo $name
    cd star.out
    cd pass.2
    echo 'Picard Read Groups + Sort'
    picard AddOrReplaceReadGroups I= Aligned.sortedByCoord.out.bam \
      O= $name.rg.sorted.out.bam SO= coordinate \
      RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample
    echo 'Picard Mark Duplicates'
    picard MarkDuplicates I= $name.rg.sorted.out.bam O= $name.dedupped.bam \
      CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics
  done
done
