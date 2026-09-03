DATA = config["data_dir"]

SE = "|".join(config["se_samples"])
PE = "|".join(config["pe_samples"])


rule trim_se:
    input:
        f"{DATA}/fastq/{{sample}}.fastq.gz"
    output:
        f"{DATA}/trimmed/{{sample}}.fastq.gz"
    wildcard_constraints:
        sample=SE
    threads: 4
    shell:
        """
        mkdir -p {DATA}/trimmed
        fastp \
        -i {input} \
        -o {output} \
        -w {threads}
        """


rule trim_pe:
    input:
        r1=f"{DATA}/fastq/{{sample}}_R1.fastq.gz",
        r2=f"{DATA}/fastq/{{sample}}_R2.fastq.gz"
    output:
        r1=f"{DATA}/trimmed/{{sample}}_R1.fastq.gz",
        r2=f"{DATA}/trimmed/{{sample}}_R2.fastq.gz"
    wildcard_constraints:
        sample=PE
    threads: 4
    shell:
        """
        mkdir -p {DATA}/trimmed
        fastp \
        -i {input.r1} \
        -I {input.r2} \
        -o {output.r1} \
        -O {output.r2} \
        -w {threads}
        """
