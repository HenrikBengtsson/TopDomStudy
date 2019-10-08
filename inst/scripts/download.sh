#!/usr/bin/env bash

url_path="https://www.ncbi.nlm.nih.gov/geo/download"
samples=(GSM2254215_ML1 GSM2254216_ML2 GSM2254218_PL1 GSM2254219_PL2)
types=(percentages validPairs assignments)
dest_path=hicData/GSE84920
mkdir -p "$dest_path"
for sample in "${samples[@]}"; do
  for type in "${types[@]}"; do
    file=$sample.$type.txt.gz
    echo "File: $file"
    if [[ ! -f "$dest_path/$file" ]]; then
      url="$url_path/?acc=${sample//_*}&format=file&file=$file"
      curl "$url" -o "$dest_path/$file"
    fi
  done
done
