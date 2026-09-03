PROJECT = config["project_dir"]
rule fastqc_raw:
    input:
        f"{PROJECT}/raw_fastq/{{sample}}.fastq"
    output:
        html=f"{PROJECT}/qc/raw/{{sample}}_fastqc.html",
        zip=f"{PROJECT}/qc/raw/{{sample}}_fastqc.zip"
    conda:
        "../envs/preprocessing.yaml"
    threads: 4
    shell:
        """
        mkdir -p {PROJECT}/qc/raw
        fastqc {input} -o {PROJECT}/qc/raw
        """

rule multiqc_raw:
    input:
        expand(f"{PROJECT}/qc/raw/{{sample}}_fastqc.zip", sample=config["samples"])
    output:
        f"{PROJECT}/qc/raw_multiqc_report.html"
    conda:
        "../envs/preprocessing.yaml"
    shell:
        """
        multiqc {PROJECT}/qc/raw \
        -o {PROJECT}/qc \
        --force \
        -n raw_multiqc_report.html
        """

rule fastqc_trimmed:
    input:
        f"{PROJECT}/trimmed/{{sample}}.trimmed.fastq"
    output:
        html=f"{PROJECT}/trimmed/{{sample}}.trimmed_fastqc.html",
        zip=f"{PROJECT}/trimmed/{{sample}}.trimmed_fastqc.zip"
    conda:
        "../envs/preprocessing.yaml"
    threads: 4
    shell:
        """
        mkdir -p {PROJECT}/trimmed
        fastqc {input} -o {PROJECT}/trimmed
        """

rule multiqc_trimmed:
    input:
        expand(f"{PROJECT}/trimmed/{{sample}}.trimmed_fastqc.zip", sample=config["samples"])
    output:
        f"{PROJECT}/qc/multiqc_report.html"
    conda:
        "../envs/preprocessing.yaml"
    shell:
        """
        mkdir -p {PROJECT}/qc
        multiqc {PROJECT}/trimmed \
        -o {PROJECT}/qc \
        --force \
        -n multiqc_report.html
        """
