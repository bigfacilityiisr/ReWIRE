rule fasterq_download:
    input:
        sra = f"{PROJECT}/sra/{{sample}}.sra"
    output:
        fastq = f"{PROJECT}/raw_fastq/{{sample}}.fastq"
    shell:
        """
        mkdir -p {PROJECT}/raw_fastq
        mkdir -p {PROJECT}/tmp
        fasterq-dump {input.sra} \
            -O {PROJECT}/raw_fastq \
            --threads 1 \
            --temp {PROJECT}/tmp \
            --bufsize 100M
        """
