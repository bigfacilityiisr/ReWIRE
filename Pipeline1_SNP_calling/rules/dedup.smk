rule mark_duplicates:
    input:
        "data/aligned/{sample}.sorted.bam"
    output:
        bam="data/aligned/{sample}.dedup.bam",
        metrics="data/aligned/{sample}.metrics.txt"
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
        "data/aligned/{sample}.dedup.bam"
    output:
        "data/aligned/{sample}.dedup.bam.bai"
    shell:
        """
        samtools index {input}
        """