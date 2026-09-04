# Comparative GO/KEGG-style comparison between the WT and MUT target-gene
# protein sets (workflow/rules/integration.smk: rule integration_go_kegg_comparison).
#
# NOTE: full GO/KEGG enrichment via clusterProfiler::compareCluster requires
# an organism annotation database (org.db) or a custom TERM2GENE mapping,
# which isn't available for every non-model organism. Adjust the section
# below once you have one for your species. Until then, this writes a
# side-by-side comparison of the annotated protein IDs from each InterProScan
# TSV so the rule still produces a usable output.

wt_go  <- read.delim(snakemake@input[["wt"]],  header = FALSE)
mut_go <- read.delim(snakemake@input[["mut"]], header = FALSE)

gene_lists <- list(WT = unique(wt_go$V1), MUT = unique(mut_go$V1))

# --- Organism-specific enrichment (fill in and uncomment once you have a
#     TERM2GENE mapping or org.db for your species) --------------------------
# library(clusterProfiler)
# ck <- compareCluster(geneCluster = gene_lists, fun = "enrichGO", ...)
# write.csv(as.data.frame(ck), snakemake@output[[1]], row.names = FALSE)
# -----------------------------------------------------------------------------

max_len <- max(length(gene_lists$WT), length(gene_lists$MUT))
length(gene_lists$WT)  <- max_len
length(gene_lists$MUT) <- max_len

write.csv(
  data.frame(WT = gene_lists$WT, MUT = gene_lists$MUT),
  snakemake@output[[1]], row.names = FALSE, na = ""
)
