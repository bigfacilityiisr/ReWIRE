"""Check whether the candidate locus's SNP positions fall inside the mature
miRNA region and/or the seed region (mature positions 2-8) of the WT
precursor (workflow/rules/integration.smk: rule integration_seed_region_check).

Uses config["integration"]["candidate_locus"]["mature_sequence"] to locate
the mature miRNA within the WT precursor sequence.
"""
candidate = snakemake.config["integration"]["candidate_locus"]
mature = candidate["mature_sequence"].upper().replace("U", "T")
changes = {int(k): v for k, v in candidate["changes"].items()}

with open(snakemake.input.wt) as fh:
    precursor_rna = "".join(line.strip() for line in fh if not line.startswith(">"))
precursor = precursor_rna.replace("U", "T")

start = precursor.find(mature)

with open(snakemake.output[0], "w") as out:
    out.write(f"Precursor: {candidate['precursor_id']}\n")
    if start == -1:
        out.write(
            f"Mature sequence '{mature}' not found in this precursor — "
            "check config['integration']['candidate_locus']['mature_sequence'].\n"
        )
    else:
        mature_start, mature_end = start + 1, start + len(mature)
        seed_start, seed_end = mature_start + 1, mature_start + 7  # mature positions 2-8
        out.write(f"Mature miRNA: {mature_start}-{mature_end}  |  Seed: {seed_start}-{seed_end}\n\n")
        for rel_pos, alt_base in sorted(changes.items()):
            in_mature = mature_start <= rel_pos <= mature_end
            in_seed = seed_start <= rel_pos <= seed_end
            out.write(
                f"SNP at relative position {rel_pos} (-> {alt_base}): "
                f"in_mature={in_mature}, in_seed={in_seed}\n"
            )

print(open(snakemake.output[0]).read())
