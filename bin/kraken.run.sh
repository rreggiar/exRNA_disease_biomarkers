#!/bin/bash

output="/public/groups/kimlab/exoRNA-biomarkers-panc/output.data/kraken.out"
input="/public/groups/kimlab/exoRNA-biomarkers-panc/data/collected.patient.data"
db_path="/public/groups/kimlab/.install_bin/"


for db_type in kraken_standard; do

  for sample in "$input"/*; do

    echo $(basename "$sample")
    echo $(basename "$db_path")

    sample_name=$(basename "$sample")

    mate_1="$sample/Unmapped.out.mate1"
    mate_2="$sample/Unmapped.out.mate2"
    
    if [ ! -d "$output/$db_type" ]; then
      
      mkdir "$output/$db_type"

    fi

    output_path="$output/$db_type"

    kraken2 --db "$db_path/$db_type" \
      --threads 32 \
      --use-names \
      --output "$output_path/$sample_name" \
      --report "$output_path/$sample_name.report" \
      --paired "$mate_1" "$mate_2"

    done

done


  
