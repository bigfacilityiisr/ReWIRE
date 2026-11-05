rule download_sra:
    output:
        sra=f"{DATA}/sra/{{sample}}/{{sample}}.sra"
    conda:
        "../envs/sra.yaml"
    shell:
        """
        mkdir -p {DATA}/sra/{wildcards.sample}
        prefetch {wildcards.sample} --output-directory {DATA}/sra
        """

