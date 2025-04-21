#! /bin/bash

while IFS= read -r line; do
    snakemake --forceall output/${line}_annotations.tsv
done < config/samples.txt