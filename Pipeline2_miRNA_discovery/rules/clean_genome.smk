rule prepare_genome_and_index:
    input:
        fasta="data/genome/GCF_018446385.1_Zo_v1.1_genomic.fna"
    output:
        chr1="data/genome/genome_chr1.fa",
        chr1_clean="data/genome/genome_chr1_clean.fa",
        genome_clean="data/genome/GCF_full_clean.fa",
        index=expand("genome_index/genome.{i}.ebwt", i=[1,2,3,4]),
        index_rev=expand("genome_index/genome.rev.{i}.ebwt", i=[1,2])
    params:
        chr_id="NC_055986.1",
        index_prefix="genome_index/genome"
    shell:
        r"""
        mkdir -p genome_index
        mkdir -p data/genome
        awk '/^>/{flag=($0 ~ "{params.chr_id}")} flag' {input.fasta} > {output.chr1}
        awk '/^>/ {{print $1; next}} {{print}}' {output.chr1} > {output.chr1_clean}
        awk '/^>/ {{print $1; next}} {{print}}' {input.fasta} > {output.genome_clean}
        bowtie-build {output.genome_clean} {params.index_prefix}
        """
