#!/bin/bash

input_dir="$1"
input_type="$2"

for salmon_dir in "$input_dir"/*/*/*."$input_type".*/logs; do

  mapping_rate=$(grep 'Mapping rate' $salmon_dir/salmon_quant.log | cut -d'=' -f2)

  output_line=$(echo $(echo $salmon_dir | cut -d'/' -f8,9) $mapping_rate | \
    sed 's/ /,/g;s/\//,/g')

  echo $output_line
done


