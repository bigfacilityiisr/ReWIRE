rule multiqc:
    input:
        expand(f"{DATA}/fastqc/{{sample}}_fastqc.zip", sample=SE_SAMPLES),
        expand(f"{DATA}/fastqc/{{sample}}_R1_fastqc.zip", sample=PE_SAMPLES),
        expand(f"{DATA}/fastqc/{{sample}}_R2_fastqc.zip", sample=PE_SAMPLES)
    output:
        f"{DATA}/multiqc/multiqc_report.html"
    shell:
        """
        mkdir -p {DATA}/multiqc
        multiqc {DATA}/fastqc -o {DATA}/multiqc
        """
