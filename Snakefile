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

rule index_1:
    input:
        "filtered_bams/{sample}_scisoseq.mapped_filtered.bam"
    output:
        "filtered_bams/{sample}_scisoseq.mapped_filtered.bam.bai"
    shell:
        "samtools index {input}"

rule left_align:
    input:
        inbam="filtered_bams/{sample}_scisoseq.mapped_filtered.bam",
        inbai="filtered_bams/{sample}_scisoseq.mapped_filtered.bam.bai",
        genome=config['reference_genome']
    output:
        outbam="leftaligned_bams/{sample}_scisoseq.mapped_leftaligned.bam"
    log:
        "logs/left_align/{sample}.log"
    shell:
        """
        gatk LeftAlignIndels \
          -R {input.genome} \
          -I {input.inbam} \
          -O {output.outbam} \
          --disable-tool-default-read-filters \
          2> {log}
        """

rule index_2:
    input:
        "leftaligned_bams/{sample}_scisoseq.mapped_leftaligned.bam"
    output:
        "leftaligned_bams/{sample}_scisoseq.mapped_leftaligned.bam.bai"
    shell:
        "samtools index {input}"

rule trim:
    input:
        "input/{sample}_barcodes.tsv"
    output:
        "input/{sample}_barcodes_trimmed.tsv"
    shell:
        "sed 's/-1$/-1/' {input} > {output}"

rule vartrix:
    input: 
        barcodes="input/{sample}_barcodes_trimmed.tsv",
        filtered_bams="leftaligned_bams/{sample}_scisoseq.mapped_leftaligned.bam",
        filtered_bams_index="leftaligned_bams/{sample}_scisoseq.mapped_leftaligned.bam.bai",
        vcf="input/mutation.vcf",
        genome=config['reference_genome'] # Downloaded from PacBio GitHub repository with download_pb_hg38.sh
    output: 
        "vartrix/{sample}_vartrix.mtx"
    log:
        "logs/vartrix/{sample}.log"
    shell: 
        """
        if [[ "{wildcards.sample}" == *10x* ]]; then
            vartrix \
              -v {input.vcf} \
              -b {input.filtered_bams} \
              -f {input.genome} \
              -c {input.barcodes} \
              -o {output} \
              --umi \
              &> {log}
        else
            vartrix \
              -v {input.vcf} \
              -b {input.filtered_bams} \
              -f {input.genome} \
              -c {input.barcodes} \
              -o {output} \
              &> {log}
        fi
        """

rule annotate:
    input: 
        "vartrix/{sample}_vartrix.mtx",
        "input/{sample}_barcodes_trimmed.tsv"
    output: 
        "output/{sample}_annotations.tsv"
    script: 
        "scripts/annotato.R"
