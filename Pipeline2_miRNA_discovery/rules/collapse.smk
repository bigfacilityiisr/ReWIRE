rule collapse_reads:
    input:
        fa = f"{PROJECT}/filtered/{{sample}}.filtered.fa",
        index = f"{PROJECT}/genome_index/genome.1.ebwt"
    output:
        collapsed = f"{PROJECT}/filtered/{{sample}}.collapsed.fa",
        arf = f"{PROJECT}/filtered/{{sample}}.arf"
    conda:
        "../envs/mirdeepP2.yaml"
    threads: 4
    shell:
        """
        mkdir -p {PROJECT}/filtered
        mapper.pl {input.fa} \
            -c \
            -j \
            -m \
            -l 18 \
            -p {PROJECT}/genome_index/genome \
            -s {output.collapsed} \
            -t {output.arf} \
            -v
        """
