# NOTE: this file is NOT included by the Snakefile. It is superseded by
# rules/alignment.smk, whose bwa_mem_se/bwa_mem_pe rules already pipe
# `bwa mem | samtools sort` to produce the sorted BAM in one step. Kept
# here only for reference; fixed up so it's at least valid Snakemake syntax
# if anyone wants to split sorting into its own rule later.
rule sort_bam:
    input:
        DATA + "/aligned/{sample}.sam"
    output:
        DATA + "/aligned/{sample}.sorted.bam"
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        samtools sort -o {output} {input}
        """


rule index_bam:
    input:
        DATA + "/aligned/{sample}.sorted.bam"
    output:
        DATA + "/aligned/{sample}.sorted.bam.bai"
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        samtools index {input}
        """

