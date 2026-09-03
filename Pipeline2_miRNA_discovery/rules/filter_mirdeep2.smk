PROJECT = config["project_dir"]
rule filter_mirdeep2:
    input:
        report=f"{PROJECT}/mirdeep/{{sample}}_report.txt"
    output:
        filtered=f"{PROJECT}/mirdeep/{{sample}}_filtered.csv"
    conda:
        "../envs/mirdeepP2.yaml"
    shell:
        """
        python scripts/filter_mirdeep2.py {input.report} {output.filtered}
        """
