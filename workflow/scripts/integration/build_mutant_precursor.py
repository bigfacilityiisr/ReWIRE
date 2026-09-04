"""Build WT and MUT precursor FASTAs for the SNP-in-precursor rewiring case
study (workflow/rules/integration.smk: rule integration_build_mutant_precursor).

Which precursor to use and which relative positions change to which bases is
a manual judgement call, made after inspecting {input.intersect}
(miRNA_precursor_SNPs.txt) to pick the SNP-in-precursor candidate you want to
study. Set that choice in config["integration"]["candidate_locus"]:

    candidate_locus:
      precursor_id: "NC_056004.1_21206"
      strand: "-"
      changes:            # relative_position_in_precursor: alt_base
        38: "T"
        31: "A"
        27: "T"
        22: "T"

This script then looks up the WT sequence for precursor_id in the merged
precursor FASTA and applies the substitutions to build MUT.
"""
import sys

candidate = snakemake.config["integration"]["candidate_locus"]
precursor_id = candidate["precursor_id"]
changes = {int(k): v for k, v in candidate["changes"].items()}

# Pull the WT precursor sequence out of the merged precursor FASTA
wt_seq = None
with open(snakemake.input.precursor_fasta) as fh:
    header = None
    seq_lines = []
    for line in fh:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if header == precursor_id:
                wt_seq = "".join(seq_lines)
                break
            header = line[1:].split()[0]
            seq_lines = []
        else:
            seq_lines.append(line)
    else:
        if wt_seq is None and header == precursor_id:
            wt_seq = "".join(seq_lines)

if wt_seq is None:
    sys.exit(
        f"Precursor '{precursor_id}' not found in {snakemake.input.precursor_fasta}. "
        "Update config['integration']['candidate_locus']['precursor_id'] to match "
        "an id from that FASTA (cross-check against "
        f"{snakemake.input.intersect})."
    )

mut_seq = list(wt_seq)
for rel_pos, alt_base in changes.items():
    if not (1 <= rel_pos <= len(mut_seq)):
        sys.exit(
            f"Relative position {rel_pos} is out of range for precursor "
            f"{precursor_id} (length {len(mut_seq)})."
        )
    mut_seq[rel_pos - 1] = alt_base
mut_seq = "".join(mut_seq)

# RNA-form FASTAs (T -> U) for RNAfold
with open(snakemake.output.wt, "w") as f:
    f.write(f">WT_{precursor_id}\n{wt_seq.replace('T', 'U')}\n")
with open(snakemake.output.mut, "w") as f:
    f.write(f">MUT_{precursor_id}\n{mut_seq.replace('T', 'U')}\n")

print(f"WT : {wt_seq}")
print(f"MUT: {mut_seq}")
