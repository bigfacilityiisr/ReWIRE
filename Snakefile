ReWIRE — top-level Snakefile
Ties Module I (SNP calling), Module II (miRNA discovery), and Module III (integration) into a single DAG.

configfile: "config/config.yaml"

include: "workflow/rules/snp_calling.smk"
include: "workflow/rules/mirna_discovery.smk"
include: "workflow/rules/integration.smk"

rule all:
    input:
	    expand(f"{DATA}/variants_r3/{{sample}}.g.vcf.gz.tbi", sample=SAMPLES),
        expand(f"{PROJECT}/mirdeep/{{sample}}_report.txt", sample=SAMPLES),
        expand(f"{PROJECT}/mirdeep/{{sample}}_filtered.csv", sample=SAMPLES)
rule integration_go_kegg_comparison:
    input:
        wt=rules.integration_target_interproscan.output.wt_tsv,
        mut=rules.integration_target_interproscan.output.mut_tsv,
    output:
        "results/integration/go_kegg_comparison.tsv"
    conda:
        "../envs/integration/r_clusterprofiler.yaml"
    script:
        "../scripts/integration/go_kegg_comparison.R"