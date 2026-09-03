PROJECT = config["project_dir"]

rule fastq_to_fasta:
    input:
        f"{PROJECT}/trimmed/{{sample}}.trimmed.fastq"
    output:
        f"{PROJECT}/filtered/{{sample}}.clean.fa"
    conda:
        "../envs/preprocessing.yaml"
    shell:
        """
        seqtk seq -A {input} > {output}
        """

