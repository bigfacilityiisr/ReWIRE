rule fasterq_se:
    input:
        "/media/user/DATADRIVE1/SNP_Project/sra/{sample}/{sample}.sra"
    output:
        f"{DATA}/fastq/{{sample}}.fastq.gz"
    threads: 24
    wildcard_constraints:
        sample="|".join(SE_SAMPLES)
    shell:
        """
        mkdir -p {DATA}/fastq
        fasterq-dump {input} -O {DATA}/fastq
        gzip {DATA}/fastq/{wildcards.sample}.fastq
        """

rule fasterq_pe:
    input:
        "/media/user/DATADRIVE1/SNP_Project/sra/{sample}/{sample}.sra"
    output:
        R1=f"{DATA}/fastq/{{sample}}_R1.fastq.gz",
        R2=f"{DATA}/fastq/{{sample}}_R2.fastq.gz"
    threads: 24
    wildcard_constraints:
        sample="|".join(PE_SAMPLES)
    shell:
        """
        mkdir -p {DATA}/fastq
        fasterq-dump {input} --split-files -O {DATA}/fastq
        gzip {DATA}/fastq/{wildcards.sample}_1.fastq
        gzip {DATA}/fastq/{wildcards.sample}_2.fastq
        mv {DATA}/fastq/{wildcards.sample}_1.fastq.gz {output.R1}
        mv {DATA}/fastq/{wildcards.sample}_2.fastq.gz {output.R2}
        """
