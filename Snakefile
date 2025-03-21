configfile: "config/config.yaml"

rule filter:
    input: 
        "input/{sample}_scisoseq.mapped.bam" # BAMs from sciso-seq output
    output: 
        "filtered_bams/{sample}_scisoseq.mapped_filtered.bam"
    params:
        region=config["filter_region"]
    log:
        "logs/samtools_filter/{sample}.log"
    shell: 
        "samtools view -b -h -o {output} {input} {params.region} 2> {log}"

rule index:
    input:
        "filtered_bams/{sample}_scisoseq.mapped_filtered.bam"
    output:
        "filtered_bams/{sample}_scisoseq.mapped_filtered.bam.bai"
    shell:
        "samtools index {input}"

rule trim:
    input:
        "input/{sample}_barcodes.tsv"
    output:
        "input/{sample}_barcodes_trimmed.tsv"
    shell:
        "sed 's/-1$//' {input} > {output}"

rule vartrix:
    input: 
        barcodes="input/{sample}_barcodes_trimmed.tsv",
        filtered_bams="filtered_bams/{sample}_scisoseq.mapped_filtered.bam",
        filtered_bams_index="filtered_bams/{sample}_scisoseq.mapped_filtered.bam.bai",
        vcf="input/mutation.vcf",
        genome=config['reference_genome'] # Downloaded from PacBio GitHub repository with download_pb_hg38.sh
    output: 
        "vartrix/{sample}_vartrix.mtx"
    log:
        "logs/vartrix/{sample}.log"
    shell: 
        """
        vartrix \
          -v {input.vcf} \
          -b {input.filtered_bams} \
          -f {input.genome} \
          -c {input.barcodes} \
          -o {output} \
          2> {log}
        """

rule annotate:
    input: 
        "vartrix/{sample}_vartrix.mtx",
        "input/{sample}_barcodes_trimmed.tsv"
    output: 
        "output/{sample}_annotations.tsv"
    script: 
        "scripts/annotato.R"