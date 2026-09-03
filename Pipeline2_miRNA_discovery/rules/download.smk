PROJECT = config["project_dir"]
rule download_sra:
    output:
        sra = f"{PROJECT}/sra/{{sample}}.sra"
    shell:
        """
        mkdir -p {PROJECT}/sra

        # run prefetch (ignore crash)
        prefetch {wildcards.sample} \
            --output-directory {PROJECT}/sra \
            --max-size 100G || true
        # FIX: move file to expected location
        if [ -f {PROJECT}/sra/{wildcards.sample}/{wildcards.sample}.sra ]; then
            mv {PROJECT}/sra/{wildcards.sample}/{wildcards.sample}.sra {output.sra}
        fi
        """
