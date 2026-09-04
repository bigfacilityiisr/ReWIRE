"""Derive WT and MUT mature-miRNA FASTAs from the WT/MUT precursor FASTAs
(workflow/rules/integration.smk: rule integration_build_mature_fastas).

Locates config["integration"]["candidate_locus"]["mature_sequence"] inside
the WT precursor, then extracts the same coordinate window from the MUT
precursor. This only works when WT and MUT precursors are the same length
(true for single-base substitutions, which is what
integration_build_mutant_precursor produces) — for indels, build the MUT
mature FASTA manually instead.
"""
candidate = snakemake.config["integration"]["candidate_locus"]
mature = candidate["mature_sequence"].upper().replace("U", "T")


def read_fasta_seq(path):
    with open(path) as fh:
        return "".join(line.strip() for line in fh if not line.startswith(">")).replace("U", "T")


wt_precursor = read_fasta_seq(snakemake.input.wt)
mut_precursor = read_fasta_seq(snakemake.input.mut)

start = wt_precursor.find(mature)
if start == -1:
    raise SystemExit(
        f"Mature sequence not found in WT precursor for {candidate['precursor_id']} — "
        "check config['integration']['candidate_locus']['mature_sequence']."
    )
end = start + len(mature)

if len(wt_precursor) != len(mut_precursor):
    raise SystemExit(
        "WT and MUT precursors differ in length; the fixed-offset extraction "
        "used here only works for single-base substitutions. Build "
        "MUT_mature.fa manually for this candidate instead."
    )

wt_mature = wt_precursor[start:end]
mut_mature = mut_precursor[start:end]

with open(snakemake.output.wt_mature, "w") as f:
    f.write(f">WT_mature_{candidate['precursor_id']}\n{wt_mature.replace('T', 'U')}\n")
with open(snakemake.output.mut_mature, "w") as f:
    f.write(f">MUT_mature_{candidate['precursor_id']}\n{mut_mature.replace('T', 'U')}\n")

print(f"WT mature : {wt_mature}")
print(f"MUT mature: {mut_mature}")
