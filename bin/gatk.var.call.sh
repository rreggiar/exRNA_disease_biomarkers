#!/bin/bash

inputdir='/public/groups/kimlab/exoRNA-biomarkers-panc/data/'
ref='/public/groups/kimlab/genomes.annotations/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'

for subdir in $inputdir/*; do
  for sample in $subdir/*; do

    cd $sample
    name=$(basename $sample)
    echo $name
    cd star.out/pass.2

    if [ ! -f $name.split.bam ]; then
      echo 'Trim Reads'
      gatk3 -T SplitNCigarReads -R $ref -I $name.dedupped.bam -o $name.split.bam \
        -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
    fi

    if [ ! -f $name.out.vcf ]; then
      echo 'Haplotype Caller'
      gatk3 -T HaplotypeCaller -R $ref -I $name.split.bam -dontUseSoftClippedBases \
        -stand_call_conf 20.0 -o $name.out.vcf
    fi

    echo 'Filter Variants'
    gatk3 -T VariantFiltration -R $ref -V $name.out.vcf -window 35 -cluster 3 \
      -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" \
      -o $name.out.vcf

  done
done
