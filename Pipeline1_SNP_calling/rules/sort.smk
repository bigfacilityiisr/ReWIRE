rule sort_bam:
    input:
        f"{DATA}/aligned/{sample}.sam"
    output:
        f"{DATA}/aligned/{sample}.sorted.bam"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools sort -o {output} {input}
        """


rule index_bam:
    input:
        f"{DATA}/aligned/{sample}.sorted.bam"
    output:
        f"{DATA}/aligned/{sample}.sorted.bam.bai"
      conda:
        "../envs/samtools.yaml"
      shell:
        """
        samtools index {input}
        """

