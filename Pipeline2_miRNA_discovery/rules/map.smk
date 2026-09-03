rule genome_map:
    input:
        fasta=f"{PROJECT}/filtered/{{sample}}.clean.fa",   
        index=expand(f"{PROJECT}/genome_index/genome.{{ext}}.ebwt",
                     ext=["1","2","3","4","rev.1","rev.2"])
    output:
        arf=f"{PROJECT}/mapped/{{sample}}.arf",
        processed=f"{PROJECT}/mapped/{{sample}}_processed.fa"
    conda:
        "../envs/mirdeepP2.yaml"
    shell:
        """
        mkdir -p {PROJECT}/mapped
        mapper.pl {input.fasta} \
            -c -j -m -l 18 \
            -p {PROJECT}/genome_index/genome \
            -s {output.processed} \
            -t {output.arf}
        """
