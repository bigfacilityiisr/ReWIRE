# ReWIRE

## Regulatory Workflow for Identifying Rewiring Effects

**An automated Snakemake-based computational workflow for identifying
SNP-induced alterations in plant microRNA-mediated regulatory networks**

------------------------------------------------------------------------

## Overview

**ReWIRE** is an integrated and reproducible computational workflow
designed to connect **genome-wide SNP discovery** with **miRNA
structural and regulatory analysis**.

The workflow combines:

-   Whole Genome Sequencing (**WGS**) for high-confidence SNP discovery
-   Small RNA sequencing (**sRNA-seq**) for miRNA discovery
-   SNP--miRNA precursor genomic intersection
-   Wild-type (**WT**) and mutant (**MUT**) precursor sequence
    generation
-   RNA secondary-structure comparison
-   WT/MUT miRNA target prediction
-   Functional characterization of affected target genes
-   Regulatory network rewiring analysis

The workflow is implemented using **Snakemake**, enabling modular
execution, reproducibility, scalability, and workflow visualization.

------------------------------------------------------------------------

## Pipeline Architecture

```{=html}
<p align="center">
```
`<img src="docs/ReWIRE_workflow.png" alt="ReWIRE integrated workflow" width="100%">`{=html}
```{=html}
</p>
```
**Figure.** Overview of the ReWIRE workflow. Module A performs WGS-based
SNP discovery and annotation; Module B performs small RNA processing and
miRNA discovery; Module C integrates SNPs with miRNA precursor loci and
compares WT/MUT structures; Module D evaluates changes in miRNA--target
interactions; and Module E characterizes the functional and regulatory
consequences.

> **Note:** The workflow figure is intentionally retained as the
> original project figure, including its icons and visual language.
> Place the supplied figure in the repository as
> `docs/ReWIRE_workflow.png`.

------------------------------------------------------------------------

# 1. Workflow Objectives

  ----------------------------------------------------------------------------------
  Module                  Objective               Major outputs
  ----------------------- ----------------------- ----------------------------------
  **A. WGS**              Genome-wide SNP         High-confidence SNP/INDEL dataset
                          discovery               

  **B. sRNA-seq**         miRNA discovery and     Known and candidate miRNA catalog
                          annotation              

  **C. SNP--miRNA         Identify SNPs within    SNP-associated miRNA precursors
  integration**           miRNA precursor loci    

  **D. Regulatory         Compare WT and MUT      Lost/gained target interactions
  impact**                miRNA targeting         

  **E. Functional         Characterize affected   GO/KEGG/protein-function/network
  analysis**              regulatory genes        results
  ----------------------------------------------------------------------------------

------------------------------------------------------------------------

# 2. Module A --- Whole Genome Sequencing (WGS)

### SNP Discovery and Variant Calling

The WGS module implements a Snakemake-based variant discovery workflow
following a GATK-oriented processing strategy.

### Processing workflow

``` text
Raw WGS reads
      │
      ├── SRA download / FASTQ generation
      │
      ├── FastQC + MultiQC
      │
      ├── Reference genome preparation
      │
      ├── BWA-MEM alignment
      │
      ├── Sorting / indexing
      │
      ├── Duplicate marking
      │
      ├── Read-group assignment
      │
      ├── GATK HaplotypeCaller
      │
      ├── Iterative BQSR
      │
      ├── Joint genotyping
      │
      ├── Variant filtering
      │
      └── High-confidence variants
```

### Principal tools

  Step                    Tool
  ----------------------- -----------------------------------
  Quality control         FastQC, MultiQC
  Alignment               BWA-MEM
  BAM processing          SAMtools, Picard
  Variant calling         GATK HaplotypeCaller
  Recalibration           GATK BaseRecalibrator / ApplyBQSR
  Variant processing      GATK, BCFtools
  Functional annotation   SnpEff / VEP

### Variant filtering

The workflow applies standard variant-quality filters including:

-   Quality by Depth (**QD**)
-   Fisher Strand bias (**FS**)
-   Mapping Quality (**MQ**)
-   Strand Odds Ratio (**SOR**)
-   Mapping Quality Rank Sum
-   Read Position Rank Sum

### WGS output

The final filtered dataset contained:

**6,720,824 high-confidence variants**

with SNP-level statistics subsequently generated using BCFtools.

------------------------------------------------------------------------

# 3. Module B --- Small RNA Sequencing and miRNA Discovery

The sRNA-seq module processes raw small RNA sequencing data and
identifies candidate miRNAs using **miRDeep2**.

### Processing workflow

``` text
Raw sRNA-seq reads
      │
      ├── FastQC
      │
      ├── Adapter trimming
      │
      ├── Size filtering
      │
      ├── FASTQ → FASTA conversion
      │
      ├── rRNA/tRNA/non-miRNA filtering
      │
      ├── Genome preparation + indexing
      │
      ├── Read collapsing
      │
      ├── miRDeep2 prediction
      │
      ├── Candidate filtering
      │
      ├── Redundancy removal
      │
      └── miRNA catalog
```

### Principal tools

-   FastQC
-   MultiQC
-   Cutadapt
-   Seqtk / Seqkit
-   Bowtie
-   miRDeep2
-   miRBase
-   BLASTn-short

### miRNA discovery

The workflow initially identified **91 candidate miRNAs**. After
redundancy reduction, **26 unique miRNA candidates** were retained for
downstream analyses.

Predicted mature miRNAs were compared against **miRBase** to assess
sequence homology and identify conserved plant miRNA families.

------------------------------------------------------------------------

# 4. Module C --- SNP--miRNA Integration

The central integration step connects genome-wide variation with miRNA
loci.

### Workflow

``` text
High-confidence SNPs
        │
        ├──────────────┐
        │              │
        ▼              ▼
SNP coordinates    miRNA precursor coordinates
        │              │
        └──────┬───────┘
               ▼
       BEDTools intersection
               │
               ▼
     SNP-associated precursors
               │
               ▼
       WT precursor sequence
               │
               ▼
       MUT precursor sequence
```

### Key tools

-   BEDTools
-   AWK
-   Python
-   FASTA/BED/VCF processing utilities

### Identified SNP-associated precursor

One candidate miRNA precursor was identified with overlapping SNPs:

**Candidate:** `NC_056004.1_21206`

**Coordinates:** `NC_056004.1:74212877..74212915:-`

Four SNPs were identified within the precursor locus.

The precursor was located on the **negative strand**, requiring
strand-aware coordinate conversion before downstream sequence and
seed-region analysis.

------------------------------------------------------------------------

# 5. Module D --- WT/MUT miRNA Structural Analysis

SNP alleles were incorporated into the precursor sequence to generate a
simulated mutant precursor.

### Sequence comparison

``` text
WT precursor
     │
     ├── SNP allele incorporation
     ▼
MUT precursor
     │
     └── RNA secondary-structure comparison
```

### RNAfold analysis

The ViennaRNA **RNAfold** program was used to compare WT and MUT
precursor structures.

  Feature                                      WT                  MUT
  ------------------------- --------------------- --------------------
  MFE                         **−16.90 kcal/mol**   **−8.10 kcal/mol**
  Ensemble diversity                     **1.24**             **7.27**
  MFE structure frequency               **0.780**            **0.148**

The mutant precursor showed reduced thermodynamic stability together
with increased structural heterogeneity.

### Interpretation

These results indicate that SNP accumulation can substantially alter
precursor secondary structure and may affect miRNA precursor processing
and maturation. Experimental validation would be required to establish
the biological effect.

------------------------------------------------------------------------

# 6. Mature miRNA and Seed-Region Analysis

The mature miRNA was mapped within the precursor sequence to determine
whether the identified SNPs affected mature or seed regions.

**Mature miRNA:** positions **21--38**

**Seed region:** positions **22--28**

Four SNPs mapped to precursor-relative positions:

``` text
Genomic SNP     Relative position
74212877   →    38
74212884   →    31
74212888   →    27
74212893   →    22
```

Thus:

-   All four SNPs occur within the mature miRNA region.
-   Two SNPs occur within the seed region.

This provides a mechanistic basis for evaluating changes in target
recognition.

------------------------------------------------------------------------

# 7. Module D --- Differential miRNA--Target Interactions

WT and MUT mature miRNA sequences were independently analyzed using
**psRNATarget**.

### Target prediction workflow

``` text
WT mature miRNA ──► psRNATarget ──► WT target repertoire
                                           │
                                           │
                                      comparison
                                           │
                                           ▼
MUT mature miRNA ──► psRNATarget ──► MUT target repertoire
```

### Comparative target results

  Parameter                                          WT       MUT
  ------------------------------------------- --------- ---------
  Unique target genes                           **185**   **154**
  High-confidence targets (Expectation ≤ 3)      **48**    **11**
  Common target genes                             **0**     **0**

The analysis identified **185 WT-specific target genes** and **154
MUT-specific target genes**, with no shared target genes in the compared
sets.

### Regulatory interpretation

The change in mature miRNA sequence resulted in a markedly different
predicted target repertoire, consistent with **SNP-induced target
rewiring**.

The predicted interaction profile also showed a reduction in
high-confidence interactions in the mutant sequence, supporting the
hypothesis that sequence variation within the mature/seed region can
alter miRNA--target recognition.

------------------------------------------------------------------------

# 8. Module E --- Functional Impact Analysis

The WT- and MUT-associated target sets were subjected to functional
characterization.

### Analyses

-   Gene Ontology (**GO**) enrichment
-   KEGG pathway enrichment
-   Protein/domain annotation using **InterProScan**
-   Functional categorization
-   Regulatory network construction
-   Comparative regulatory interpretation

### Functional comparison

The WT target repertoire included genes associated with functions such
as:

-   Transport
-   Protein folding
-   Organelle-associated regulation
-   Cell-wall-associated processes
-   Developmental regulation

The MUT target repertoire included genes associated with:

-   Stress adaptation
-   Signal transduction
-   Phosphate metabolism
-   Translational regulation
-   Other regulatory processes

These differences support a model in which SNP-induced changes in miRNA
sequence can redirect post-transcriptional regulatory interactions.

------------------------------------------------------------------------

# 9. End-to-End ReWIRE Workflow

``` text
                    ReWIRE
                       │
        ┌──────────────┴──────────────┐
        │                             │
        ▼                             ▼
   WGS MODULE                    sRNA-seq MODULE
        │                             │
   SNP discovery                 miRNA discovery
        │                             │
   SNP filtering                 miRNA catalog
        │                             │
        └──────────────┬──────────────┘
                       ▼
              SNP–miRNA integration
                       │
                       ▼
             WT/MUT precursor design
                       │
                       ▼
              RNAfold structure analysis
                       │
                       ▼
             WT/MUT mature miRNAs
                       │
                       ▼
            psRNATarget comparison
                       │
                       ▼
              Target rewiring analysis
                       │
                       ▼
            GO / KEGG / InterProScan
                       │
                       ▼
             Regulatory network impact
```

------------------------------------------------------------------------

# 10. Key Outputs

ReWIRE generates the following major output classes:

  Output                      Description
  --------------------------- ------------------------------------------
  `*.vcf.gz`                  Filtered genomic variants
  `*.bam`                     Processed alignment files
  `miRNA` prediction tables   Candidate miRNA predictions
  `miRNA FASTA`               Mature and precursor sequences
  `BED` files                 miRNA precursor genomic coordinates
  SNP--miRNA intersections    SNPs overlapping precursor loci
  WT/MUT FASTA                Simulated wild-type and mutant sequences
  RNAfold results             Secondary structures and MFE values
  WT/MUT target tables        Predicted miRNA--target interactions
  Functional enrichment       GO/KEGG results
  Network files               Comparative regulatory networks

------------------------------------------------------------------------

# 11. Repository Structure

ReWIRE is organized as one top-level Snakemake project that ties three
sub-projects together. `Pipeline1_SNP_calling/` (Module A) and
`Pipeline2_miRNA_discovery/` (Module B) are each complete, independently
runnable Snakemake workflows in their own right; `workflow/` holds the
post-processing and integration rules (Modules C–E) that connect them.

``` text
ReWIRE_merged/
│
├── README.md
├── CHANGES.md                     # engineering log of fixes applied when the
│                                   #   three original modules were unified
├── Snakefile                      # top-level Snakefile — runs Modules A–E together
├── config/
│   └── config.yaml                # single unified configuration (see §13)
│
├── Pipeline1_SNP_calling/         # Module A — SNP calling (standalone-runnable)
│   ├── Snakefile
│   ├── config/config.yaml
│   ├── envs/                      # gatk.yaml, mapping.yaml, qc.yaml, sra.yaml
│   └── rules/                     # download, fasterq, fastqc_se/pe, multiqc,
│                                   #   trimming, reference, alignment, dedup,
│                                   #   gvcf_r1, bqsr(1-3), gvcf_r2, final_genotype
│
├── Pipeline2_miRNA_discovery/     # Module B — miRNA discovery (standalone-runnable)
│   ├── Snakefile
│   ├── config/config.yaml
│   ├── envs/                      # download, preprocessing, mapping, mirdeepP2, blast
│   ├── scripts/                   # csv_to_fasta.py, filter_mirdeep2.py, merge_mirna.py
│   └── rules/                     # download, fasterq, qc, trimming, fasta,
│                                   #   filter_ncRNA, clean_genome, collapse,
│                                   #   mirdeepP2, filter_mirdeep2, merge_mirna,
│                                   #   mirna_to_fasta, mirbase_db, blast_mirna
│
├── workflow/                      # Module A post-processing + Modules C/D/E (integration)
│   ├── envs/integration/          # bcftools, bedtools, gatk, snpeff, seqtk,
│   │                               #   python, viennarna, interproscan, r_clusterprofiler
│   ├── rules/
│   │   ├── snp_postprocessing.smk # hard-filter → PASS → snpEff → HIGH-impact
│   │   └── integration.smk        # SNP × precursor intersection → WT/MUT
│   │                               #   precursor → RNAfold → seed check →
│   │                               #   target comparison → GO/KEGG
│   └── scripts/integration/       # build_mutant_precursor.py, seed_region_check.py,
│                                   #   extract_mature_from_precursor.py, go_kegg_comparison.R
│
├── data/
│   └── reference/                 # genome FASTA/GFF, proteome, Rfam, miRBase,
│                                   #   transcriptome — see below
│
├── results/                       # all pipeline outputs (created on first run):
│   ├── variants/                  #   Module A raw joint-genotyped VCF
│   ├── snp/                       #   Module A filtered/annotated SNPs
│   ├── mirna/                     #   Module B miRNA predictions
│   └── integration/               #   Modules C–E WT/MUT structure and targets
│
└── docs/
    └── ReWIRE_workflow.png
```

Large sequencing datasets and reference genomes are **not** committed to
GitHub. `config/config.yaml` points at local paths under `data/` for
these files — populate them locally (or via a data-retrieval script) and
document their accession numbers instead of committing the files
themselves.

------------------------------------------------------------------------

# 12. Reproducibility and Step-by-Step Execution

ReWIRE is designed around **Snakemake** to provide:

-   Modular workflow execution
-   Automatic dependency tracking
-   Reproducible analysis
-   Parallel execution
-   Intermediate-file management
-   Workflow visualization
-   Environment management

Every rule declares its own `conda:` environment file, so no tool needs
to be installed by hand — `--use-conda` builds and activates each one
automatically the first time it's needed.

### Step 1 — Clone the repository and install Snakemake

``` bash
git clone <your-repo-url> ReWIRE
cd ReWIRE
conda create -n rewire -c bioconda -c conda-forge snakemake
conda activate rewire
```

### Step 2 — Provide reference data

Place these files where `config/config.yaml` expects them (or edit the
paths to wherever you keep them):

  Reference file             Used by
  --------------------------- -------------------------------------------------
  Genome FASTA                Module A alignment/calling; Module B genome cleaning
  Genome annotation (GFF3)    HIGH-impact transcript/protein + gene→protein map
  Proteome FASTA               HIGH-impact and target-protein extraction
  Rfam FASTA                   Module B non-miRNA filtering
  miRBase FASTA                Module B homology search (BLAST)
  Transcriptome FASTA          manual psRNATarget submission (Step 7 below)
  snpEff genome database       registered locally in `snpEff.config`, for SNP annotation

### Step 3 — Set sample accessions

In `config/config.yaml`, set:

-   `se_samples` / `pe_samples` — WGS single-end / paired-end SRA
    accessions (Module A)
-   `samples` — small-RNA-seq SRA accessions (Module B)

### Step 4 — Dry-run the full DAG

Always validate before executing:

``` bash
snakemake -n
```

This resolves every rule and file dependency with nothing actually run.
Confirm the job table ends in `all` with no `MissingInputException` or
`AmbiguousRuleException` beyond files you haven't supplied yet.

Workflow visualization can still be generated the same way as before:

``` bash
snakemake --rulegraph | dot -Tpng > rulegraph.png
```

### Step 5 — Run Modules A and B

These two modules are independent of each other and can run in the same
invocation:

``` bash
snakemake --use-conda --cores 18 -p
```

Snakemake will build everything through Module B's non-redundant miRNA
catalogue and Module A's annotated/HIGH-impact SNP set, then stop
naturally at the Module C–D manual checkpoint described next (the
psRNATarget files it needs don't exist yet).

### Step 6 — Identify a SNP-in-precursor candidate (manual)

Once Modules A and B have both finished, Module C runs automatically
through:

``` text
results/integration/miRNA_precursor_SNPs.txt
```

Inspect this file and choose the SNP/precursor pair to study as your
WT vs. MUT case (see §4 for the example candidate used in this
analysis). Record the choice in `config.yaml`:

``` yaml
integration:
  candidate_locus:
    precursor_id: "NC_056004.1_21206"
    strand: "-"
    changes:                       # relative_position_in_precursor: alt_base
      38: "T"
      31: "A"
      27: "T"
      22: "T"
    mature_sequence: "GGCCGGCCATCGTCGTGA"
```

### Step 7 — Build WT/MUT precursor, fold, and check the seed region

``` bash
snakemake --use-conda --cores 4 -p \
    results/integration/WT_precursor.fold \
    results/integration/MUT_precursor.fold \
    results/integration/seed_region_analysis.txt \
    results/integration/WT_mature.fa \
    results/integration/MUT_mature.fa
```

### Step 8 — Manual checkpoint: psRNATarget

Submit `WT_mature.fa` and `MUT_mature.fa` to
[psRNATarget](https://www.zhaolab.org/psRNATarget/) against your
transcriptome FASTA, then save the results as:

``` text
results/integration/WT_targets.txt
results/integration/MUT_targets.txt
```

psRNATarget is a hosted external web tool with no stable API for
automated submission, so this step stays manual by design.

### Step 9 — Finish Modules D and E

``` bash
snakemake --use-conda --cores 8 -p results/integration/go_kegg_comparison.tsv
```

This produces the shared/lost/gained WT-vs-MUT target gene lists, their
InterProScan domain annotations, and the comparative GO/KEGG summary.

### Running a single module standalone

`Pipeline1_SNP_calling/` and `Pipeline2_miRNA_discovery/` are also
complete, independent Snakemake projects with their own
`config/config.yaml` — useful for isolated testing or reuse outside the
full ReWIRE workflow:

``` bash
cd Pipeline1_SNP_calling
snakemake --use-conda --cores 16 -n
snakemake --use-conda --cores 16 -p

cd ../Pipeline2_miRNA_discovery
snakemake --use-conda --cores 16 -n
snakemake --use-conda --cores 16 -p
```

------------------------------------------------------------------------

# 13. Configuration

All project-specific parameters live in one file, **`config/config.yaml`**
at the repository root — nothing is hard-coded into the workflow rules.
It uses YAML anchors (`&name` / `*name`) so a value defined once (e.g. a
sample list) is reused wherever it's needed, without the two copies
drifting apart.

  Section                        Consumed by                                   Key fields
  ------------------------------ ---------------------------------------------- --------------------------------------------------------------
  Flat top-level keys            `Pipeline1_SNP_calling/rules/*.smk`            `reference`, `data_dir`, `se_samples`, `pe_samples`
  `snp:` / `snpeff:`             `workflow/rules/snp_postprocessing.smk`        `output_dir`, `annotation_gff`, `proteome_fasta`, `gene_protein_map`, `genome_name`, `config_path`
  Flat top-level keys            `Pipeline2_miRNA_discovery/rules/*.smk`        `project_dir`, `samples`, `adapter`, `raw_genome_fasta`, `target_chr_id`, `genome_fasta`, `rfam_fasta`, `mirbase`
  `integration:`                 `workflow/rules/integration.smk`               `output_dir`, `transcriptome_fasta`, `candidate_locus`

Excerpt (see the repository's `config/config.yaml` for the complete,
annotated file):

``` yaml
# ---------------- Module A: SNP calling (flat keys) ------------------------
reference: &reference "data/reference/genome.fna"
data_dir: &data_dir "data"
se_samples: &se_samples
  - SRR12340661
pe_samples: &pe_samples
  - SRR36323019
  - SRR12329118

# ---------------- Module A post-processing (nested) ------------------------
snp:
  output_dir: "results/snp"
  se_samples: *se_samples
  pe_samples: *pe_samples
  reference: *reference
  annotation_gff: "data/reference/genome.gff"
  proteome_fasta: "data/reference/proteome.faa"
  gene_protein_map: "results/integration/gene_protein_map.txt"

snpeff:
  genome_name: "your_snpeff_genome_id"       # registered in your snpEff.config
  config_path: "/path/to/snpEff/snpEff.config"

# ---------------- Module B: miRNA discovery (flat keys) --------------------
project_dir: &project_dir "results/mirna"
samples: &mirna_samples
  - SRR12080165
  - SRR12080166
adapter: "AGATCGGAAGAGC"
raw_genome_fasta: "data/reference/genome_raw.fna"
target_chr_id: "NC_055986.1"
genome_fasta: "data/reference/GCF_full_clean.fa"
rfam_fasta: "data/reference/rfam.fa"
mirbase: "data/reference/mirbase.fa"

# ---------------- Modules C–E: integration (nested) -------------------------
integration:
  output_dir: "results/integration"
  transcriptome_fasta: "data/reference/transcriptome.fasta"
  candidate_locus:                 # filled in after inspecting
    precursor_id: "NC_056004.1_21206"   # miRNA_precursor_SNPs.txt — see §12, Step 6
    strand: "-"
    changes:
      38: "T"
      31: "A"
    mature_sequence: "GGCCGGCCATCGTCGTGA"
```

`Pipeline1_SNP_calling/config/config.yaml` and
`Pipeline2_miRNA_discovery/config/config.yaml` hold matching,
simplified copies of the flat keys for standalone runs of each module —
edit those instead of the root config if you're only running one module
in isolation (§12, "Running a single module standalone").

Before running the workflow, update every path above and the sample
accession lists to match your local dataset.

------------------------------------------------------------------------

# 14. Software Requirements

The workflow uses or interfaces with:

-   Snakemake
-   Conda
-   FastQC
-   MultiQC
-   BWA
-   SAMtools
-   Picard
-   GATK
-   BCFtools
-   SnpEff / VEP
-   Cutadapt
-   Seqtk / Seqkit
-   Bowtie
-   miRDeep2
-   BLAST+
-   BEDTools
-   ViennaRNA / RNAfold
-   psRNATarget
-   InterProScan
-   R-based functional enrichment tools, where applicable

Exact versions should be recorded in the repository environment files
for reproducibility.

------------------------------------------------------------------------

# 15. Scientific Significance

ReWIRE provides a computational framework for moving from **genomic
variation** to potential **regulatory consequences**:

``` text
SNP
 ↓
miRNA precursor
 ↓
precursor structure
 ↓
mature miRNA sequence
 ↓
seed-region recognition
 ↓
target repertoire
 ↓
functional pathways
 ↓
regulatory network rewiring
```

The workflow therefore enables systematic investigation of how genomic
variants may influence miRNA-mediated post-transcriptional regulation.

------------------------------------------------------------------------

# 16. Key Findings from the Demonstration Analysis

The integrated analysis demonstrated:

-   **6,720,824** high-confidence genomic variants were retained after
    filtering.
-   **91** miRNA candidates were initially predicted.
-   **26** non-redundant miRNA candidates were retained.
-   SNP--precursor intersection identified a candidate miRNA precursor
    containing **four SNPs**.
-   All four SNPs mapped within the mature miRNA sequence, including
    **two within the seed region**.
-   WT precursor MFE was **−16.90 kcal/mol**, compared with **−8.10
    kcal/mol** for the mutant precursor.
-   Ensemble diversity increased from **1.24 to 7.27** after mutation.
-   MFE structure frequency decreased from **0.780 to 0.148**.
-   WT and MUT target sets showed extensive differences, with **185
    WT-specific** and **154 MUT-specific** target genes in the analyzed
    comparison.
-   High-confidence predicted targets decreased from **48 to 11**.

Together, these observations demonstrate the potential of ReWIRE to
identify candidate **SNP-induced structural and regulatory rewiring
events** in plant miRNA networks.

------------------------------------------------------------------------

# 17. Limitations and Interpretation

ReWIRE is a **computational prediction framework**. Predicted changes in
RNA structure, miRNA processing, target recognition, and regulatory
networks represent hypotheses that require experimental validation.

Recommended validation approaches include:

-   miRNA expression analysis
-   Stem-loop RT-qPCR
-   Small RNA northern blotting
-   5′ RLM-RACE / cleavage validation
-   Reporter assays for WT and MUT target interactions
-   Genotype--phenotype association analysis

------------------------------------------------------------------------

# 18. Citation

If you use ReWIRE in your research, please cite the associated pipeline
paper:

> **ReWIRE: Regulatory Workflow for Identifying Rewiring Effects --- An
> automated Snakemake-based computational workflow for identifying
> SNP-induced alterations in plant microRNA-mediated regulatory
> networks.**

A formal citation will be added after publication.

## ReWIRE at a Glance

**Input**

`WGS FASTQ + sRNA-seq FASTQ + reference genome + annotation`

**Processing**

`SNP discovery → miRNA discovery → SNP–miRNA integration → WT/MUT structure → target prediction → functional analysis`

**Output**

`High-confidence SNPs + candidate miRNAs + structural changes + altered target repertoire + functional/regulatory impact`

------------------------------------------------------------------------

**ReWIRE**\
*From genomic variation to miRNA regulatory rewiring.*
