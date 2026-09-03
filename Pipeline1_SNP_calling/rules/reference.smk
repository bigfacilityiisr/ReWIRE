rule bwa_index:
    input:
        REF
    output:
        REF + ".amb",
        REF + ".ann",
        REF + ".bwt",
        REF + ".pac",
        REF + ".sa"
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        bwa index {input}
        """

rule faidx:
    input:
        REF
    output:
        REF + ".fai"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools faidx {input}
        """

rule dict:
    input:
        REF
    output:
        REF.replace(".fna",".dict")
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        gatk CreateSequenceDictionary \
        -R {input} \
        -O {output}
        """
