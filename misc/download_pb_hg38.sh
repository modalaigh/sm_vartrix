#! /bin/bash

cd ../input
wget -O human_GRCh38_no_alt_analysis_set.tar.gz https://downloads.pacbcloud.com/public/reference-genomes/human_GRCh38_no_alt_analysis_set.tar.2023-12-04.gz
tar -xzf human_GRCh38_no_alt_analysis_set.tar.gz
mv human_GRCh38_no_alt_analysis_set/human_GRCh38_no_alt_analysis_set.fasta* .
rm -rf human_GRCh38_no_alt_analysis_set/ human_GRCh38_no_alt_analysis_set.tar.gz
