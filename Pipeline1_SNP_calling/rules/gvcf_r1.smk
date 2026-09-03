rule haplotypecaller_r1:
    input:
        bam="data/aligned/{sample}.dedup.bam",
        ref=config["reference"]
    output:
        "results/variants/{sample}.r1.g.vcf"
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        mkdir -p results/variants

        gatk HaplotypeCaller \
        -R {input.ref} \
        -I {input.bam} \
        -O {output} \
        -ERC GVCF
        """


rule combine_r1:
    input:
        expand("results/variants/{sample}.r1.g.vcf", sample=SAMPLES)
    output:
        "results/variants/cohort.r1.g.vcf"
    params:
        ref=config["reference"],
        gvcfs=lambda wildcards, input: " ".join(["-V " + f for f in input])
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        gatk CombineGVCFs \
        -R {params.ref} \
        {params.gvcfs} \
        -O {output}
        """


rule genotype_r1:
    input:
        "results/variants/cohort.r1.g.vcf"
    output:
        "results/variants/known_variants.vcf"
    params:
        ref=config["reference"]
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        gatk GenotypeGVCFs \
        -R {params.ref} \
        -V {input} \
        -O {output}
        """
