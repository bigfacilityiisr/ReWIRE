DATA = config["data_dir"]

FASTQC_DIR = f"{DATA}/fastqc"
PE = "|".join(config["pe_samples"])

rule fastqc_pe:
    input:
        R1=f"{DATA}/fastq/{{sample}}_R1.fastq.gz",
        R2=f"{DATA}/fastq/{{sample}}_R2.fastq.gz"
    output:
        R1_html=f"{FASTQC_DIR}/{{sample}}_R1_fastqc.html",
        R1_zip=f"{FASTQC_DIR}/{{sample}}_R1_fastqc.zip",
        R2_html=f"{FASTQC_DIR}/{{sample}}_R2_fastqc.html",
        R2_zip=f"{FASTQC_DIR}/{{sample}}_R2_fastqc.zip"
    wildcard_constraints:
        sample=PE
    shell:
        """
        mkdir -p {FASTQC_DIR}
        fastqc {input.R1} {input.R2} --outdir {FASTQC_DIR}
        """
