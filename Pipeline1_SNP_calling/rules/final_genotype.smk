rule final_haplotypecaller:
    input:
        bam="data/aligned/{sample}.bqsr3.bam",
        ref=config["reference"]
    output:
        "results/variants/{sample}.final.g.vcf"
    shell:
        """
        gatk HaplotypeCaller \
        -R {input.ref} \
        -I {input.bam} \
        -O {output} \
        -ERC GVCF
        """


rule final_combine:
    input:
        expand("results/variants/{sample}.final.g.vcf", sample=SAMPLES)
    output:
        "results/variants/final_cohort.g.vcf"
    params:
        ref=config["reference"],
        gvcfs=lambda wildcards, input: " ".join(f"-V {i}" for i in input)
    shell:
        """
        gatk CombineGVCFs \
        -R {params.ref} \
        {params.gvcfs} \
        -O {output}
        """

rule final_genotype:
    input:
        "results/variants/final_cohort.g.vcf"
    output:
        "results/variants/final_variants.vcf"
    params:
        ref=config["reference"]
    shell:
        """
        gatk GenotypeGVCFs \
        -R {params.ref} \
        -V cohort.g.vcf\
        -O final_all_sites.vcf \
        --include-non-variant-sites
        """
