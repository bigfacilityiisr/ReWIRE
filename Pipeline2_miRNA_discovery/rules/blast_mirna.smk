rule blast_mirna:
    input:
        query="results/mature_miRNAs.fa",
        db="reference/miRBase_db.nsq"
    output:
        "results/mirna_vs_miRBase.txt"
    conda:
        "../envs/blast.yaml"
    shell:
        """
        blastn \
            -task blastn-short \
            -query {input.query} \
            -db reference/miRBase_db \
            -out {output} \
            -outfmt 6
        """
