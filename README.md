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

A recommended GitHub organization is:

``` text
ReWIRE/
│
├── README.md
├── LICENSE
├── CITATION.cff
│
├── workflow/
│   ├── Snakefile
│   ├── config/
│   ├── rules/
│   ├── scripts/
│   └── envs/
│
├── data/
│   └── README.md
│
├── reference/
│   └── README.md
│
├── results/
│   ├── variants/
│   ├── miRNA/
│   ├── structure/
│   ├── targets/
│   └── enrichment/
│
├── docs/
│   └── ReWIRE_workflow.png
│
└── LICENSE
```

Large sequencing datasets and reference genomes should **not** be
committed directly to GitHub. Instead, document their accession numbers
and download instructions.

------------------------------------------------------------------------

# 12. Reproducibility

ReWIRE is designed around **Snakemake** to provide:

-   Modular workflow execution
-   Automatic dependency tracking
-   Reproducible analysis
-   Parallel execution
-   Intermediate-file management
-   Workflow visualization
-   Environment management

A typical workflow execution is:

``` bash
snakemake --use-conda --cores 18 -p
```

Workflow visualization can be generated using:

``` bash
snakemake --rulegraph | dot -Tpng > rulegraph.png
```

------------------------------------------------------------------------

# 13. Configuration

Project-specific parameters should be maintained in a configuration file
rather than hard-coded into workflow rules.

Example:

``` yaml
project_dir: "/path/to/ReWIRE"

reference:
  genome: "/path/to/reference/genome.fa"
  annotation: "/path/to/reference/annotation.gff"

miRNA:
  adapter: "AACTGTAGGCACCAT"

resources:
  threads: 18
```

Before running the workflow, update paths and sample information
according to the local dataset.

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

------------------------------------------------------------------------

# 19. License

Add the appropriate open-source license for the repository, for example
**MIT**, **GPL-3.0**, or another license appropriate to the project.

------------------------------------------------------------------------

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

