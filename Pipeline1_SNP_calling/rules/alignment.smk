SE = "|".join(config["se_samples"])
PE = "|".join(config["pe_samples"])

DATA = config["data_dir"]
REF = config["reference"]

rule bwa_mem_pe:
    input:
        ref=REF,
        r1=f"{DATA}/trimmed/{{sample}}_R1.fastq.gz",
        r2=f"{DATA}/trimmed/{{sample}}_R2.fastq.gz"
    output:
        f"{DATA}/aligned/{{sample}}.sorted.bam"
    threads: 12
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        mkdir -p {DATA}/aligned/tmp
        bwa mem -t {threads} \
        -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:ILLUMINA" \
        {input.ref} {input.r1} {input.r2} | \
        samtools sort -@ {threads} -m 1G \
        -T {DATA}/aligned/tmp -o {output}
        """

rule bwa_mem_se:
    input:
        ref=REF,
        reads=f"{DATA}/trimmed/{{sample}}.fastq.gz"
    output:
        f"{DATA}/aligned/{{sample}}.sorted.bam"
    threads: 12
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        mkdir -p {DATA}/aligned/tmp
        bwa mem -t {threads} {input.ref} {input.reads} | \
        samtools sort -@ {threads} -m 1G -T {DATA}/aligned/tmp -o {output}
        """
