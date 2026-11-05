rule mark_duplicates:
    input:
        DATA + "/aligned/{sample}.sorted.bam"
    output:
        bam=DATA + "/aligned/{sample}.dedup.bam",
        metrics=DATA + "/aligned/{sample}.metrics.txt"
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        gatk MarkDuplicates \
        -I {input} \
        -O {output.bam} \
        -M {output.metrics}
        """

rule index_dedup:
    input:
        DATA + "/aligned/{sample}.dedup.bam"
    output:
        DATA + "/aligned/{sample}.dedup.bam.bai"
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        samtools index {input}
        """