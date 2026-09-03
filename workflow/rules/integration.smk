# Module III — SNP-precursor integration and the WT vs MUT rewiring case
INT_OUT = config["integration"]["output_dir"]

rule mirna_precursor_bed:
    input:
        rules.mirna_merge_predictions.output   # all_samples_merged.csv, from Module II
    output:
        fasta=f"{INT_OUT}/precursor_miRNAs.fa",
        bed=f"{INT_OUT}/precursor_miRNA.bed",
    shell:
        r"""
        awk -F',' 'NR>1 {{print ">"$1"\n"$16}}' {input} > {output.fasta}
        awk -F',' '
        NR>1 {{
          split($17,a,":"); chr=a[1]
          split(a[2],b,"\\.\\."); start=b[1]-1; end=b[2]
          print chr"\t"start"\t"end"\t"$1
        }}' {input} > {output.bed}
        """


rule integration_intersect_snps_precursors:
    input:
        snp_vcf=rules.snp_extract_pass.output.vcf,     # from snp_postprocessing.smk
        precursor_bed=rules.mirna_precursor_bed.output.bed,
    output:
        f"{INT_OUT}/miRNA_precursor_SNPs.txt"
    conda:
        "../envs/integration/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.snp_vcf} -b {input.precursor_bed} -wa -wb > {output}"


rule integration_build_mutant_precursor:
    # NOTE: precursor id / relative SNP positions / base changes are
    # candidate-locus-specific (yours: NC_056004.1_21206). Edit
    # workflow/scripts/integration/build_mutant_precursor.py per candidate,
    # or parameterize it to read directly from the intersect output above.
    input:
        rules.integration_intersect_snps_precursors.output
    output:
        wt=f"{INT_OUT}/WT_precursor.fa",
        mut=f"{INT_OUT}/MUT_precursor.fa",
    conda:
        "../envs/integration/python.yaml"
    script:
        "../scripts/integration/build_mutant_precursor.py"


rule integration_rnafold:
    input:
        wt=rules.integration_build_mutant_precursor.output.wt,
        mut=rules.integration_build_mutant_precursor.output.mut,
    output:
        wt_fold=f"{INT_OUT}/WT_precursor.fold",
        mut_fold=f"{INT_OUT}/MUT_precursor.fold",
    conda:
        "../envs/integration/viennarna.yaml"
    shell:
        """
        RNAfold -p --noPS < {input.wt}  > {output.wt_fold}
        RNAfold -p --noPS < {input.mut} > {output.mut_fold}
        """


rule integration_seed_region_check:
    input:
        wt=rules.integration_build_mutant_precursor.output.wt
    output:
        f"{INT_OUT}/seed_region_analysis.txt"
    conda:
        "../envs/integration/python.yaml"
    script:
        "../scripts/integration/seed_region_check.py"


rule integration_build_mature_fastas:
    input:
        wt=rules.integration_build_mutant_precursor.output.wt,
        mut=rules.integration_build_mutant_precursor.output.mut,
    output:
        wt_mature=f"{INT_OUT}/WT_mature.fa",
        mut_mature=f"{INT_OUT}/MUT_mature.fa",
    conda:
        "../envs/integration/python.yaml"
    script:
        "../scripts/integration/extract_mature_from_precursor.py"


# --- [MANUAL CHECKPOINT] ---------------------------------------------------
# Submit WT_mature.fa and MUT_mature.fa to psRNATarget against your
# transcriptome, and save the results as:
#   results/integration/WT_targets.txt
#   results/integration/MUT_targets.txt
# before running the rules below. This cannot be automated (external web tool).
# -----------------------------------------------------------------------------


rule integration_compare_target_repertoire:
    input:
        wt=f"{INT_OUT}/WT_targets.txt",     # produced manually via psRNATarget
        mut=f"{INT_OUT}/MUT_targets.txt",   # produced manually via psRNATarget
    output:
        common=f"{INT_OUT}/common_genes.txt",
        lost=f"{INT_OUT}/lost_genes.txt",
        gained=f"{INT_OUT}/gained_genes.txt",
    shell:
        """
        grep -oE 'LOC[0-9]+' {input.wt}  | sort -u > {INT_OUT}/WT_genes.txt
        grep -oE 'LOC[0-9]+' {input.mut} | sort -u > {INT_OUT}/MUT_genes.txt
        comm -12 {INT_OUT}/WT_genes.txt {INT_OUT}/MUT_genes.txt > {output.common}
        comm -23 {INT_OUT}/WT_genes.txt {INT_OUT}/MUT_genes.txt > {output.lost}
        comm -13 {INT_OUT}/WT_genes.txt {INT_OUT}/MUT_genes.txt > {output.gained}
        """


rule integration_target_proteins:
    input:
        wt_genes=f"{INT_OUT}/WT_genes.txt",
        mut_genes=f"{INT_OUT}/MUT_genes.txt",
        gene_protein_map=config["snp"]["gene_protein_map"],
        proteome=config["snp"]["proteome_fasta"],
    output:
        wt_faa=f"{INT_OUT}/WT_proteins.faa",
        mut_faa=f"{INT_OUT}/MUT_proteins.faa",
    conda:
        "../envs/integration/seqtk.yaml"
    shell:
        """
        grep -Ff {input.wt_genes}  {input.gene_protein_map} | cut -f2 | sort -u > {INT_OUT}/WT_proteins.ids
        grep -Ff {input.mut_genes} {input.gene_protein_map} | cut -f2 | sort -u > {INT_OUT}/MUT_proteins.ids
        seqtk subseq {input.proteome} {INT_OUT}/WT_proteins.ids  > {output.wt_faa}
        seqtk subseq {input.proteome} {INT_OUT}/MUT_proteins.ids > {output.mut_faa}
        """


rule integration_target_interproscan:
    input:
        wt=rules.integration_target_proteins.output.wt_faa,
        mut=rules.integration_target_proteins.output.mut_faa,
    output:
        wt_tsv=f"{INT_OUT}/WT_interpro.tsv",
        mut_tsv=f"{INT_OUT}/MUT_interpro.tsv",
    conda:
        "../envs/integration/interproscan.yaml"
    shell:
        """
        interproscan.sh -i {input.wt}  -f tsv -goterms -pa -o {output.wt_tsv}
        interproscan.sh -i {input.mut} -f tsv -goterms -pa -o {output.mut_tsv}
        """


rule integration_go_kegg_comparison:
    input:
        wt=rules.integration_target_interproscan.output.wt_tsv,
        mut=rules.integration_target_interproscan.output.mut_tsv,
    output:
        f"{INT_OUT}/go_kegg_comparison.tsv"
    conda:
        "../envs/integration/r_clusterprofiler.yaml"
    script:
        "../scripts/integration/go_kegg_comparison.R"
