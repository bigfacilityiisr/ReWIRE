# Module II — miRNA Discovery Pipeline (sRNA-seq)

A complete, modular, and reproducible Snakemake workflow for miRNA identification from raw small RNA-seq (sRNA-seq) SRA data, using miRDeep2.

## Tools used

- SRA Toolkit (`prefetch`, `fasterq-dump`)
- FastQC / MultiQC
- cutadapt (adapter trimming)
- seqtk (format conversion)
- Bowtie + seqkit (non-miRNA / rRNA-tRNA filtering against Rfam)
- miRDeep2 (`mapper.pl`, core prediction, plant mode)
- BLAST (`makeblastdb`, `blastn-short`) for miRBase homology search
- psRNATarget (miRNA target prediction, external web tool)

## Workflow

```
SRA download → FASTQ conversion → FastQC/MultiQC (raw)
 → adapter trimming (cutadapt) → FastQC/MultiQC (trimmed)
 → FASTQ→FASTA → non-miRNA filtering (Bowtie vs Rfam)
 → genome indexing → read collapsing (mapper.pl)
 → miRNA prediction (miRDeep2, plant mode)
 → filtering (score/read-count thresholds) → merge & de-duplicate across samples
 → miRBase homology search (BLASTn-short) → conserved vs novel classification
 → target prediction (psRNATarget) → precursor BED file generation
```

## Project structure

```
pipeline2_mirna_discovery/
├── Snakefile
├── config/
│   └── config.yaml
├── workflow/
│   ├── rules/
│   ├── scripts/
│   └── envs/
└── test_data/
```

## Configuration

`config/config.yaml` specifies the sample accessions, adapter sequence, Rfam reference FASTA, and genome FASTA, e.g.:

```yaml
project_dir: "<path>"
samples:
  - SRR_ACCESSION_1
  - SRR_ACCESSION_2
adapter: "<adapter sequence>"
rfam_fasta: "<path to Rfam reference>"
genome_fasta: "<path to reference genome>"
```

## Running the pipeline

```bash
snakemake --use-conda --cores <N> -p
```

Visualize the rule graph:
```bash
snakemake --rulegraph | dot -Tpng > rulegraph.png
```

## Key steps

1. **Download & FASTQ conversion** — `prefetch` + `fasterq-dump`.
2. **QC (raw)** — FastQC + MultiQC.
3. **Adapter trimming** — `cutadapt`, restricted to the miRNA size range (`-m 17 -M 35 -q 20`).
4. **QC (trimmed)** — confirms quality improvement post-trimming.
5. **FASTQ → FASTA** — `seqtk`, required for miRDeep2.
6. **Non-miRNA filtering** — reads mapped to Rfam (rRNA/tRNA/etc.) with Bowtie; unmapped reads retained as the miRNA-enriched set.
7. **Genome preparation & indexing** — reference cleaned and indexed with `bowtie-build`.
8. **Read collapsing** — `mapper.pl` (removes duplicates, junk, low-complexity reads; minimum length 18 nt).
9. **miRNA prediction** — miRDeep2 core algorithm, plant mode (`-t plant`), using collapsed reads, genome, and alignment file.
10. **Filtering** — high-confidence candidates retained by score (≥5), mature read count (≥10), and star read count (≥1).
11. **Merging** — per-sample predictions combined and de-duplicated to a non-redundant miRNA set (`merge_mirna.py`).
12. **Homology search** — mature miRNAs BLASTed (`blastn-short`) against a local miRBase database to classify conserved vs. putative novel candidates.
13. **Target prediction** — non-redundant miRNAs submitted to psRNATarget against the reference transcriptome.
14. **Precursor BED generation** — precursor genomic coordinates extracted from the merged prediction table into BED format, for intersection with SNPs in Module III.

## Output

- `all_samples_merged.csv` — non-redundant predicted miRNAs with precursor/mature sequences and genomic coordinates
- `precursor_miRNA.bed` — precursor coordinates (input to Module III)
- `mirna_vs_miRBase.txt` — homology search results
- `merged_psRNATarget.txt` — predicted targets

## Notes

Full sequencing datasets are not included in this repository. A minimal `test_data/` subset is provided to verify the pipeline runs end-to-end. The Rfam and miRBase reference databases must be downloaded separately (links in `config/config.yaml` comments).
