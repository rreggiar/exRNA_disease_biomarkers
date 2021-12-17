#!/bin/bash
# reqs:  
# * all are included in exo-rna conda environment*
# bedtools (bed12toBed6)
# ucsc-gtfToGenePred
# ucsc-genePredToBed

INPUT=$1

NAME=$(basename "$INPUT" .gtf)

grep -e $'\texon\t' "$INPUT" > "$NAME.exon.gtf"

gtfToGenePred "$NAME.exon.gtf" "$NAME.exon.genePred"

genePredToBed "$NAME.exon.genePred" "$NAME.exon.bed12"

bed12ToBed6 -i "$NAME.exon.bed12" > "$NAME.exon.bed6"

rm "$NAME.exon.genePred" "$NAME.exon.gtf"
