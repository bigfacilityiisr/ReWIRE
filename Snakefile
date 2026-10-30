# ============================================================================
# ReWIRE — top-level Snakefile
# Ties Module I (SNP calling + post-processing), Module II (miRNA discovery), and Module III (integration / WT vs MUT rewiring case study) into a single DAG, driven by the one config/config.yaml.

# Run from the repo root:
#   snakemake --use-conda --cores <N> -n     # dry run first
#   snakemake --use-conda --cores <N> -p
#
# Pipeline1_SNP_calling/ and Pipeline2_miRNA_discovery/ can each also be run
# standalone (from inside their own directory, with their own Snakefile +
# config/config.yaml) for isolated testing — this file does not touch those.
#
# Manual checkpoint: Module III stops at WT_mature.fa/MUT_mature.fa. Submit
# both to psRNATarget (https://www.zhaolab.org/psRNATarget/) against
# config["integration"]["transcriptome_fasta"], and save the results as
#   results/integration/WT_targets.txt
#   results/integration/MUT_targets.txt
# before asking Snakemake for anything past that point (see the
# [MANUAL CHECKPOINT] comment in workflow/rules/integration.smk).
# ============================================================================

configfile: "config/config.yaml"

# ---- variables consumed by Pipeline1_SNP_calling/rules/*.smk --------------
REF = config["reference"]
DATA = config["data_dir"]
SE_SAMPLES = config["se_samples"]
PE_SAMPLES = config["pe_samples"]
SAMPLES = SE_SAMPLES + PE_SAMPLES

# ---- variables consumed by Pipeline2_miRNA_discovery/rules/*.smk ----------
PROJECT = config["project_dir"]

# ============================================================================
# Module I — SNP calling (Pipeline1_SNP_calling/rules)
# ============================================================================
include: "Pipeline1_SNP_calling/rules/download.smk"
include: "Pipeline1_SNP_calling/rules/fasterq.smk"
include: "Pipeline1_SNP_calling/rules/fastqc_se.smk"
include: "Pipeline1_SNP_calling/rules/fastqc_pe.smk"
include: "Pipeline1_SNP_calling/rules/multiqc.smk"
include: "Pipeline1_SNP_calling/rules/trimming.smk"
include: "Pipeline1_SNP_calling/rules/reference.smk"
include: "Pipeline1_SNP_calling/rules/alignment.smk"
include: "Pipeline1_SNP_calling/rules/dedup.smk"
include: "Pipeline1_SNP_calling/rules/gvcf_r1.smk"
include: "Pipeline1_SNP_calling/rules/bqsr.smk"
include: "Pipeline1_SNP_calling/rules/gvcf_r2.smk"
include: "Pipeline1_SNP_calling/rules/bqsr2.smk"
include: "Pipeline1_SNP_calling/rules/bqsr3.smk"
include: "Pipeline1_SNP_calling/rules/final_genotype.smk"

# ============================================================================
# Module I post-processing — filtering, annotation, HIGH-impact proteins
# (workflow/rules). Depends on rules.final_genotype from Module I above.
# ============================================================================
include: "workflow/rules/snp_postprocessing.smk"

# ============================================================================
# Module II — miRNA discovery (Pipeline2_miRNA_discovery/rules)
# ============================================================================
include: "Pipeline2_miRNA_discovery/rules/download.smk"
include: "Pipeline2_miRNA_discovery/rules/fasterq.smk"
include: "Pipeline2_miRNA_discovery/rules/qc.smk"
include: "Pipeline2_miRNA_discovery/rules/trimming.smk"
include: "Pipeline2_miRNA_discovery/rules/fasta.smk"
include: "Pipeline2_miRNA_discovery/rules/filter_ncRNA.smk"
include: "Pipeline2_miRNA_discovery/rules/clean_genome.smk"
include: "Pipeline2_miRNA_discovery/rules/collapse.smk"
# rules/map.smk is intentionally NOT included — see the comment in
# Pipeline2_miRNA_discovery/Snakefile for why (superseded by collapse.smk).
include: "Pipeline2_miRNA_discovery/rules/mirdeep_P2.smk"
include: "Pipeline2_miRNA_discovery/rules/filter_mirdeep_P2.smk"
include: "Pipeline2_miRNA_discovery/rules/merge_mirna.smk"
include: "Pipeline2_miRNA_discovery/rules/mirna_to_fasta.smk"
include: "Pipeline2_miRNA_discovery/rules/mirbase_db.smk"
include: "Pipeline2_miRNA_discovery/rules/blast_mirna.smk"

# ============================================================================
# Module III — SNP-precursor integration / WT vs MUT rewiring case study
# Depends on rules.snp_extract_pass (Module I post-processing) and
# rules.mirna_merge_predictions (Module II) above.
# ============================================================================
include: "workflow/rules/integration.smk"


rule all:
    input:
        # ---- Module I: final annotated SNP set + HIGH-impact protein domains
        rules.snp_variant_stats.output,
        rules.snp_chromosome_distribution.output,
        rules.snp_tstv.output,
        rules.snp_annotation_summary.output,
        rules.snp_high_impact_interproscan.output,
        # ---- Module II: merged non-redundant miRNAs + miRBase homology search
        rules.mirna_merge_predictions.output.non_redundant,
        rules.blast_mirna.output,
        # ---- Module III: everything automatable up to the psRNATarget
        # manual checkpoint (WT/MUT precursor structure, seed-region check,
        # WT/MUT mature FASTAs ready to submit to psRNATarget)
        rules.integration_rnafold.output,
        rules.integration_seed_region_check.output,
        rules.integration_build_mature_fastas.output,
        # NOTE: rules.integration_compare_target_repertoire and everything
        # after it need results/integration/WT_targets.txt and
        # MUT_targets.txt from psRNATarget first (see the top-of-file
        # comment). Once those exist, run e.g.:
        #   snakemake --use-conda --cores <N> results/integration/go_kegg_comparison.tsv
