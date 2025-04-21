# Introduction
This repo is for a snakemake pipeline to run the 10x-developed tool [Vartrix](https://github.com/10XGenomics/vartrix) to identify cells with a mutation of interest in PacBio Kinnex scRNA-seq data.
The primary output is a .tsv file containing celltype annotations (mutant/non-mutant) which can be imported as metadata to a Seurat object for downstream analysis. 
Please be aware that due to the technical limitations of scRNA-seq data, a `non-mutant` cell may not necessarily be unmutatated; the coverage of this gene in some mutated cells may have been too sparse to identify the mutation. 

### Set up conda environment
```
conda create -n sm_vartrix -c bioconda snakemake vartrix samtools conda-forge::r-essentials bioconda::gatk4
```

### Required inputs
Please place the following SMRT Link Read_Segmentation_and_Single-Cell_Iso-Seq outputs in the `input/` folder:
- `{sample}_barcodes.tsv`
- `{sample}_scisoseq.mapped.bam`
- `{sample}_scisoseq.mapped.bam.bai`

Where `{sample}` is the ID of a particular sample.

Also required in the inputs folder:
- A suitable reference genome. To download the reference genome used in the SMRT Link pipeline, please execute `misc/download_pb_hg38.sh`
- A properly-formatted VCF with a variant of interest. See `input/mutation.vcf` of an example of an NPM1 mutation

### Snakemake configuration file

This file is located at `config/config.yaml` and has two entries which may need to be altered:
- `filter_region`: The region in which to subset the BAM file. For the NPM1 example, a region +/- 13 bp from the mutation was selected. This was chosen to constrain the BAM file to reads aligning within the exon containing the mutation
- `reference_genome`: The path to the reference genome used. Needs to be altered if using a different reference than the one specified above

## Running the pipeline
Activate the conda environment:
```
conda activate sm_vartrix
```
Run the snakemake pipeline:
```
snakemake output/{sample}_annotations.tsv
```
Replace {sample} with the sample ID that you want to run the pipeline for

## Overview of the pipeline:
<img align="right" height=450 src="/imgs/pipeline.png">

1. `Filter`

   BAM file is subset to only focus on the region surrounding the mutation
   This region is specified in the `config/config.yaml` file

3. `Index`

   An index of the subsetted BAM file is generated with `samtools` 

5. `Trim`

   The "-1" is trimmed from barcodes as Kinnex BAMs don't contain this in the CB BAM tag

7. `Vartrix`

   `Vartrix` is used to assign mutations to single-cells
  
9. `Annotate`

   Outputs barcodes with their intrepreted assignments from Vartrix (mutant/non-mutant)
