#!/usr/bin/bash

# wrapper script to generate alu loci for editing

input_counts_file=$1
echo ""
output_tmp_counts=$2
output_locus_bed=$3

echo "$input_counts_file" "$output_tmp_counts" | Rscript alu.loci.for.editing.R

until [ -f "$output_tmp_counts" ]

do
	sleep 5
done

echo "creating bed file"

cut -d"'" -f1 "$output_tmp_counts" | cut -d'=' -f2 | sed 's/:/\t/g ; s/-/\t/g ; s/_5//g' | sort -k1,1 -k2,2n > "$output_locus_bed"

rm "$output_tmp_counts"
