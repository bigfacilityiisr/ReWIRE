PROJECT = config["project_dir"]
rule build_rfam_index:
    input:
        fasta = config["rfam_fasta"]
    output:
        expand(f"{PROJECT}/rfam_index.{{ext}}.ebwt",
               ext=["1","2","3","4","rev.1","rev.2"])
    shell:
        """
        bowtie-build {input.fasta} {PROJECT}/rfam_index
        """
rule filter_non_miRNA:
    input:
        reads = f"{PROJECT}/filtered/{{sample}}.clean.fa",
        index = expand(f"{PROJECT}/rfam_index.{{ext}}.ebwt",
                       ext=["1","2","3","4","rev.1","rev.2"])
    output:
        filtered = f"{PROJECT}/filtered/{{sample}}.filtered.fa"
    conda:
        "../envs/mirdeepP2.yaml"
    threads: 4
    shell:
        """
        mkdir -p {PROJECT}/filtered
        # Map against RFAM
        bowtie -f -v 1 -p {threads} -S \
            -x {PROJECT}/rfam_index \
            {input.reads} > {PROJECT}/filtered/{wildcards.sample}_rfam.sam
        # Extract unmapped read IDs
        grep -v "^@" {PROJECT}/filtered/{wildcards.sample}_rfam.sam \
        | awk '$3=="*"' \
        | cut -f1 \
        | sort | uniq > {PROJECT}/filtered/{wildcards.sample}_keep_ids.txt
        # Extract sequences + FIX FASTA FORMAT
        seqkit grep -f {PROJECT}/filtered/{wildcards.sample}_keep_ids.txt {input.reads} \
        | awk '/^>/ {{print $1}} !/^>/ {{print}}' \
        | seqkit seq -w 0 \
        > {output.filtered}
        """




