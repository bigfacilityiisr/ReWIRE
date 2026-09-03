DATA = config["data_dir"]

FASTQC_DIR = f"{DATA}/fastqc"
SE = "|".join(config["se_samples"])

rule fastqc_se:
    input:
        f"{DATA}/fastq/{{sample}}.fastq.gz"
    output:
        html=f"{FASTQC_DIR}/{{sample}}_fastqc.html",
        zip=f"{FASTQC_DIR}/{{sample}}_fastqc.zip"
    threads: 4
    wildcard_constraints:
        sample=SE
    shell:
        """
        mkdir -p {FASTQC_DIR}
        fastqc -t 4 --memory 4096 {input} --outdir {FASTQC_DIR}
        """
