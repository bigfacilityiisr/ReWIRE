rule final_haplotypecaller:
    input:
        bam=DATA + "/aligned/{sample}.bqsr3.bam",
        ref=config["reference"]
    output:
        "results/variants/{sample}.final.g.vcf"
    conda:
        "../envs/gatk.yaml"
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
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        gatk CombineGVCFs \
        -R {params.ref} \
        {params.gvcfs} \
        -O {output}
        """

rule final_genotype:
    # NOTE: previously this rule ignored {input}/{output} entirely and
    # referenced files ("cohort.g.vcf", "final_all_sites.vcf") that don't
    # exist anywhere else in the pipeline. Fixed to use the real
    # input/output and to write a bgzipped, auto-indexed VCF so downstream
    # rules (workflow/rules/snp_postprocessing.smk) can consume it directly.
    input:
        "results/variants/final_cohort.g.vcf"
    output:
        vcf="results/variants/final_variants.vcf.gz",
        tbi="results/variants/final_variants.vcf.gz.tbi"
    params:
        ref=config["reference"]
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        gatk GenotypeGVCFs \
        -R {params.ref} \
        -V {input} \
        -O {output.vcf}
        """
