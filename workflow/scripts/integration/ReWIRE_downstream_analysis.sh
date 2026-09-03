#!/usr/bin/env bash
###############################################################################
# ReWIRE — Module III Downstream Analysis (merged script)
#   Part A: Genome-wide SNP characterization  
#   Part B: SNP-in-precursor rewiring case study 
#
# Run from a directory containing merged_PASS.vcf.gz, all_samples_merged.csv,
# and the reference genome/annotation files referenced below.
#
# Steps marked [MANUAL] use an external web tool (psRNATarget) or a GUI
# (IGV) and cannot be scripted — they are left as clearly marked checkpoints.
###############################################################################

set -euo pipefail

# ---------------------------------------------------------------------------
# 0. CONFIG — edit these paths for your run
# ---------------------------------------------------------------------------
WORKDIR="."
SNP_VCF="merged_PASS.vcf.gz"                 # Module I output
MIRNA_CSV="all_samples_merged.csv"           # Module II output
GENOME_FASTA="GCF_018446385.1_Zo_v1.1_genomic.fna"
GENOME_GFF="GCF_018446385.1_Zo_v1.1_genomic.gff"
PROTEOME_FASTA="GCF_018446385.1_Zo_v1.1_protein.faa"
SNPEFF_GENOME="ginger"                        # name you registered in snpEff.config
OUTDIR="${WORKDIR}/results/integration"
mkdir -p "${OUTDIR}"
cd "${WORKDIR}"

###############################################################################
# PART A — GENOME-WIDE SNP CHARACTERIZATION (Objectives 12-19)
###############################################################################

# --- A1. Whole-genome variant statistics (Objective 12) --------------------
bcftools stats "${SNP_VCF}" > "${OUTDIR}/variant_stats.txt"
plot-vcfstats -p "${OUTDIR}/stats_plots" "${OUTDIR}/variant_stats.txt" || true

bcftools view -v snps "${SNP_VCF}" -Oz -o "${OUTDIR}/SNP_only.vcf.gz"
bcftools index "${OUTDIR}/SNP_only.vcf.gz"

# --- A2. Chromosome-wise SNP distribution (Objective 13) -------------------
bcftools query -f '%CHROM\n' "${OUTDIR}/SNP_only.vcf.gz" \
  | sort | uniq -c > "${OUTDIR}/snp_density_only.txt"

grep "NC_" "${OUTDIR}/snp_density_only.txt" \
  | awk '$1>1000' > "${OUTDIR}/chromosome_snp_final.txt"

awk '{sum+=$1} END {print "Total SNPs on major chromosomes:", sum}' \
  "${OUTDIR}/chromosome_snp_final.txt"

# --- A3. Transition/Transversion ratio (Objective 14) -----------------------
bcftools stats "${OUTDIR}/SNP_only.vcf.gz" | grep "^TSTV" \
  > "${OUTDIR}/tstv_stats.txt"

# --- A4. SnpEff functional annotation (Objective 15) ------------------------
# One-time database build (skip if already built):
#   mkdir -p ~/snpeff/data/${SNPEFF_GENOME}
#   ln -s "$(pwd)/${GENOME_FASTA}" ~/snpeff/data/${SNPEFF_GENOME}/sequences.fa
#   ln -s "$(pwd)/${GENOME_GFF}"   ~/snpeff/data/${SNPEFF_GENOME}/genes.gff
#   cd ~/snpeff && snpEff build -gff3 -noCheckCds -noCheckProtein -v ${SNPEFF_GENOME}

snpEff -c ~/snpeff/snpEff.config -noStats -v "${SNPEFF_GENOME}" \
  "${OUTDIR}/SNP_only.vcf.gz" > "${OUTDIR}/SNP_annotated.vcf"

{
  echo "===== SNP EFFECT SUMMARY ====="
  for eff in missense_variant synonymous_variant stop_gained stop_lost \
             upstream_gene_variant downstream_gene_variant intergenic_region \
             intron_variant 5_prime_UTR_variant 3_prime_UTR_variant; do
    echo -n "${eff}: "
    grep -o "${eff}" "${OUTDIR}/SNP_annotated.vcf" | wc -l
  done
  echo
  echo "===== IMPACT CATEGORIES ====="
  for impact in HIGH MODERATE LOW MODIFIER; do
    echo -n "${impact} impact: "
    grep -o "|${impact}|" "${OUTDIR}/SNP_annotated.vcf" | wc -l
  done
} > "${OUTDIR}/SNP_annotation_statistics.txt"

# --- A5. Extract HIGH-impact SNPs, transcripts, proteins (Objectives 16, 18) -
grep "|HIGH|" "${OUTDIR}/SNP_annotated.vcf" > "${OUTDIR}/HIGH_impact_variants.vcf"

awk -F'|' '$3=="HIGH"{print $7}' "${OUTDIR}/SNP_annotated.vcf" \
  | sort -u > "${OUTDIR}/HIGH_only_transcripts.txt"

# transcript -> protein map built from the GFF (adjust the extraction to your
# GFF's actual attribute format if it differs)
awk -F'\t' '$3=="mRNA"' "${GENOME_GFF}" \
  | grep -oE 'Genbank:[A-Z]{2}_[0-9]+\.[0-9]+,Genbank:[A-Z]{2}_[0-9]+\.[0-9]+' \
  > "${OUTDIR}/transcript_protein_map.txt" || true

grep -Ff "${OUTDIR}/HIGH_only_transcripts.txt" "${OUTDIR}/transcript_protein_map.txt" \
  | grep -oE '[A-Z]{2}_[0-9]+\.[0-9]+' | sort -u \
  > "${OUTDIR}/HIGH_impact_proteins.ids"

# extract those protein sequences from the full proteome (needs seqtk or samtools faidx)
seqtk subseq "${PROTEOME_FASTA}" "${OUTDIR}/HIGH_impact_proteins.ids" \
  > "${OUTDIR}/HIGH_impact_proteins.faa"

# --- A6. [MANUAL] IGV validation (Objective 17) -----------------------------
# Load merged_PASS.vcf.gz + the reference genome/GFF into IGV and visually
# confirm read support (depth, allele balance) for the HIGH-impact SNPs
# listed in HIGH_impact_variants.vcf.

# --- A7. InterProScan on HIGH-impact proteins (Objective 19) ---------------
interproscan.sh -i "${OUTDIR}/HIGH_impact_proteins.faa" \
  -f tsv -goterms -pa -o "${OUTDIR}/HIGH_impact_proteins_interpro.tsv"

###############################################################################
# PART B — SNP-IN-PRECURSOR REWIRING CASE STUDY (Objectives 8-11, 21-23)
###############################################################################

# --- B1. Precursor FASTA + BED from miRNA discovery output (Objective 8) ---
awk -F',' 'NR>1 {print ">"$1"\n"$16}' "${MIRNA_CSV}" \
  > "${OUTDIR}/precursor_miRNAs.fa"

awk -F',' '
NR>1 {
  split($17,a,":")
  chr=a[1]
  split(a[2],b,"\\.\\.")
  start=b[1]-1
  end=b[2]
  print chr"\t"start"\t"end"\t"$1
}' "${MIRNA_CSV}" > "${OUTDIR}/precursor_miRNA.bed"

# --- B2. Intersect SNPs with precursor coordinates (Objective 9) -----------
bedtools intersect -a "${SNP_VCF}" -b "${OUTDIR}/precursor_miRNA.bed" \
  -wa -wb > "${OUTDIR}/miRNA_precursor_SNPs.txt"

echo "Precursors with overlapping SNPs / total SNPs found:"
cut -f4 "${OUTDIR}/miRNA_precursor_SNPs.txt" | sort -u | wc -l
wc -l < "${OUTDIR}/miRNA_precursor_SNPs.txt"

# --- B3. Compute strand-aware relative SNP positions (Objective 5/9) -------
# EDIT these four values for the candidate precursor identified above.
python3 - <<'PYEOF'
precursor_id   = "NC_056004.1_21206"
precursor_seq  = "CGCACTTTGGCTCGCCTCATGGCCGGCCATCGTCGTGA"   # from all_samples_merged.csv
precursor_start = 74212877          # genomic start of precursor (from coords)
strand = "-"
snp_positions = [74212877, 74212884, 74212888, 74212893]

L = len(precursor_seq)
print(f"Precursor {precursor_id} length = {L}")

corrected = {}
for s in snp_positions:
    if strand == "-":
        rel = L - (s - precursor_start)
    else:
        rel = (s - precursor_start) + 1
    corrected[s] = rel
    print(f"{s} -> relative position {rel}")

import json
with open("results/integration/corrected_snp_positions.json", "w") as fh:
    json.dump({"precursor_id": precursor_id, "precursor_seq": precursor_seq,
               "corrected_positions": corrected}, fh, indent=2)
PYEOF

# --- B4. Build WT and MUT precursor sequences (Objectives 6-7) -------------
# EDIT `changes` to the corrected relative positions -> ALT base for your locus.
python3 - <<'PYEOF'
wt = "CGCACTTTGGCTCGCCTCATGGCCGGCCATCGTCGTGA"
changes = {38: "T", 31: "A", 27: "T", 22: "T"}   # relative_pos: alt_base

seq = list(wt)
for pos, base in changes.items():
    seq[pos - 1] = base
mut = "".join(seq)

with open("results/integration/WT_precursor.fa", "w") as f:
    f.write(">WT_precursor\n" + wt.replace("T", "U") + "\n")   # RNA form for RNAfold
with open("results/integration/MUT_precursor.fa", "w") as f:
    f.write(">MUT_precursor\n" + mut.replace("T", "U") + "\n")

print("WT :", wt)
print("MUT:", mut)
PYEOF

# --- B5. RNA secondary structure prediction, WT vs MUT (Objective 8) -------
RNAfold -p --noPS < "${OUTDIR}/WT_precursor.fa"  > "${OUTDIR}/WT_precursor.fold"
RNAfold -p --noPS < "${OUTDIR}/MUT_precursor.fa" > "${OUTDIR}/MUT_precursor.fold"
# .fold files contain MFE structure, MFE (kcal/mol), ensemble diversity,
# and frequency of the MFE structure in the ensemble — compare these values
# directly between the two files.

# --- B6. Seed-region overlap analysis (Objective 9) -------------------------
python3 - <<'PYEOF'
precursor = "CGCACTTTGGCTCGCCTCATGGCCGGCCATCGTCGTGA"
mature    = "GGCCGGCCATCGTCGTGA"

start = precursor.find(mature)
mature_start, mature_end = start + 1, start + len(mature)
seed_start, seed_end = mature_start + 1, mature_start + 7   # positions 2-8 of mature

print(f"Mature miRNA: {mature_start}-{mature_end}  |  Seed: {seed_start}-{seed_end}")

corrected_positions = {38: 74212877, 31: 74212884, 27: 74212888, 22: 74212893}
for rel_pos, genomic in corrected_positions.items():
    in_mature = mature_start <= rel_pos <= mature_end
    in_seed = seed_start <= rel_pos <= seed_end
    print(f"SNP {genomic} (rel {rel_pos}): in mature={in_mature}, in seed={in_seed}")
PYEOF

# --- B7. [MANUAL] Build WT/MUT mature miRNA FASTAs and re-run psRNATarget --
cat > "${OUTDIR}/WT_mature.fa"  <<'EOF'
>WT_miRNA
GGCCGGCCAUCGUCGUGA
EOF
# Build MUT_mature.fa the same way using the mutant mature sequence, then
# submit WT_mature.fa and MUT_mature.fa separately to psRNATarget
# (https://www.zhaolab.org/psRNATarget/) against your transcriptome FASTA.
# Save results as WT_targets.txt and MUT_targets.txt in ${OUTDIR}.

# --- B8. Compare WT vs MUT target repertoires (Objective 10) ---------------
grep -oE 'LOC[0-9]+' "${OUTDIR}/WT_targets.txt"  | sort -u > "${OUTDIR}/WT_genes.txt"
grep -oE 'LOC[0-9]+' "${OUTDIR}/MUT_targets.txt" | sort -u > "${OUTDIR}/MUT_genes.txt"

comm -12 "${OUTDIR}/WT_genes.txt" "${OUTDIR}/MUT_genes.txt" > "${OUTDIR}/common_genes.txt"
comm -23 "${OUTDIR}/WT_genes.txt" "${OUTDIR}/MUT_genes.txt" > "${OUTDIR}/lost_genes.txt"
comm -13 "${OUTDIR}/WT_genes.txt" "${OUTDIR}/MUT_genes.txt" > "${OUTDIR}/gained_genes.txt"

echo "Common / lost / gained gene counts:"
wc -l "${OUTDIR}"/{common,lost,gained}_genes.txt

# --- B9. Gene -> protein mapping and functional annotation (Objective 21) --
# Build gene_protein_map.txt from the GFF (LOC id -> protein id), then:
grep -Ff "${OUTDIR}/WT_genes.txt"  "${OUTDIR}/gene_protein_map.txt" > "${OUTDIR}/WT_gene_protein.txt"
grep -Ff "${OUTDIR}/MUT_genes.txt" "${OUTDIR}/gene_protein_map.txt" > "${OUTDIR}/MUT_gene_protein.txt"

cut -f2 "${OUTDIR}/WT_gene_protein.txt"  | sort -u > "${OUTDIR}/WT_proteins.ids"
cut -f2 "${OUTDIR}/MUT_gene_protein.txt" | sort -u > "${OUTDIR}/MUT_proteins.ids"

seqtk subseq "${PROTEOME_FASTA}" "${OUTDIR}/WT_proteins.ids"  > "${OUTDIR}/WT_proteins.faa"
seqtk subseq "${PROTEOME_FASTA}" "${OUTDIR}/MUT_proteins.ids" > "${OUTDIR}/MUT_proteins.faa"

interproscan.sh -i "${OUTDIR}/WT_proteins.faa"  -f tsv -goterms -pa -o "${OUTDIR}/WT_interpro.tsv"
interproscan.sh -i "${OUTDIR}/MUT_proteins.faa" -f tsv -goterms -pa -o "${OUTDIR}/MUT_interpro.tsv"

# --- B10. Comparative GO / KEGG enrichment, WT vs MUT (Objectives 22-23) ---
# Run in R (clusterProfiler), comparing WT_interpro.tsv-derived GO terms
# against MUT_interpro.tsv-derived GO terms as two gene lists.
Rscript - <<'REOF'
library(clusterProfiler)
# Adjust these read/parse steps to match your actual InterProScan TSV columns
wt_go  <- read.delim("results/integration/WT_interpro.tsv",  header = FALSE)
mut_go <- read.delim("results/integration/MUT_interpro.tsv", header = FALSE)

gene_lists <- list(WT = unique(wt_go$V1), MUT = unique(mut_go$V1))

# Requires an org.db or a custom TERM2GENE mapping for Zingiber officinale
# ck <- compareCluster(geneCluster = gene_lists, fun = "enrichGO", ...)
# write.csv(as.data.frame(ck), "results/integration/go_kegg_comparison.tsv")
REOF

echo "Downstream analysis complete. Outputs in ${OUTDIR}"
