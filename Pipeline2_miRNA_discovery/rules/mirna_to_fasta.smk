rule mirna_to_fasta:
    input:
        "/media/user/DATADRIVE1/miRNA_Project/mirdeep/non_redundant_miRNAs.csv"
    output:
        "results/mature_miRNAs.fa"
    conda:
        "../envs/blast.yaml"
    shell:
        """
        python scripts/csv_to_fasta.py {input} {output}
        """
