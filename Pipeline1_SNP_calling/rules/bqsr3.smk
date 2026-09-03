rule bqsr3:
    input:
        bam="data/aligned/{sample}.bqsr2.bam",
        ref=config["reference"],
        known="results/variants/known_variants.vcf"
    output:
        bam="data/aligned/{sample}.bqsr3.bam",
        table="data/aligned/{sample}.recal3.table"
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
