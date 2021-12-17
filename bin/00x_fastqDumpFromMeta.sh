#!/bin/bash

## metadata file is useful to identify samples and label after import
metadata_run_csv="$1"
# extract first column 
for sra in $(cut -d, -f1 "$metadata_run_csv" | tail -n +2); do

  echo "$sra"
  # --split-files: split paired-end fastq files into respective read files
  # --gzip: gzip compress upon fastq generation
  # --origfmt: bug catch for gzip problems
  fastq-dump --split-files --origfmt --gzip "$sra"

done
