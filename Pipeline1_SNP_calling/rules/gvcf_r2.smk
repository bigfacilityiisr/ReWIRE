rule haplotypecaller_r2:
    input:
        bam="data/aligned/{sample}.bqsr1.bam",
        ref=config["reference"]
    output:
        "results/variants/{sample}.r2.g.vcf"
    shell:
        """
        gatk HaplotypeCaller \
        -R {input.ref} \
        -I {input.bam} \
        -O {output} \
        -ERC GVCF
        """