rule bqsr2:
    input:
        bam=DATA + "/aligned/{sample}.bqsr1.bam",
        ref=config["reference"],
        known="results/variants/known_variants.vcf"
    output:
        bam=DATA + "/aligned/{sample}.bqsr2.bam",
        table=DATA + "/aligned/{sample}.recal2.table"
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        gatk BaseRecalibrator \
        -R {input.ref} \
        -I {input.bam} \
        --known-sites {input.known} \
        -O {output.table}

        gatk ApplyBQSR \
        -R {input.ref} \
        -I {input.bam} \
        --bqsr-recal-file {output.table} \
        -O {output.bam}
        """