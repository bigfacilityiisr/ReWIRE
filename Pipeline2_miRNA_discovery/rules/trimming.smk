PROJECT=config["project_dir"]
ADAPTER = config["adapter"]
rule cutadapt_trim:
    input:
        f"{PROJECT}/raw_fastq/{{sample}}.fastq",
    output:
        f"{PROJECT}/trimmed/{{sample}}.trimmed.fastq"
    conda:
        "../envs/preprocessing.yaml"
    threads: 8    
    shell:
        """
        mkdir -p {PROJECT}/trimmed
        cutadapt \
            -a {ADAPTER} \
            -m 17 \
            -M 35 \
            -q 20 \
            -o {output} \
            {input}
        """
