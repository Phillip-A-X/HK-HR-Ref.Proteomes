# HPA-GTEx non-HK gene loss QC
#
# Purpose:
#   Check why the number of non-HK genes decreases after mapping the
#   full Ensembl 115 HRT Atlas background to HPA-GTEx and applying
#   the nTPM >= 1 expression threshold.
#
# Inputs:
#   data/processed/human_all_genes_with_hrt_atlas_hk_status_rescued.tsv
#   data/processed/hpa_gtex/human_hpa_gtex_all_mapped_gene_tissue_pairs.tsv
#   data/processed/hpa_gtex/HumanProteinCodingGenes_bytissue_onchromosomes.tsv
#
# Outputs:
#   results/hpa_gtex_hk/tables/nonhk_gene_loss_after_hpa_gtex_mapping_and_ntpm_filter.tsv
#   results/hpa_gtex_hk/tables/nonhk_not_mapped_to_hpa_gtex.tsv
#   results/hpa_gtex_hk/tables/nonhk_mapped_but_never_ntpm_ge_1.tsv
#   results/hpa_gtex_hk/tables/hk_nonhk_gene_counts_across_hpa_gtex_steps.tsv

# ------------------------------------------------------------
# Setup
# ------------------------------------------------------------

table_output_dir <- "results/hpa_gtex_hk/tables"
dir.create(table_output_dir, recursive = TRUE, showWarnings = FALSE)

background_file <- "data/processed/human_all_genes_with_hrt_atlas_hk_status_rescued.tsv"

mapped_all_file <- "data/processed/hpa_gtex/human_hpa_gtex_all_mapped_gene_tissue_pairs.tsv"

expressed_file <- "data/processed/hpa_gtex/HumanProteinCodingGenes_bytissue_onchromosomes.tsv"

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

background <- read.delim(
  background_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

mapped_all <- read.delim(
  mapped_all_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

expressed <- read.delim(
  expressed_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ------------------------------------------------------------
# Check required columns
# ------------------------------------------------------------

required_background_cols <- c("Gene stable ID", "hk_status")
required_mapped_cols <- c("gene_id", "hk_status")
required_expressed_cols <- c("gene_id", "hk_status")

missing_background <- setdiff(required_background_cols, colnames(background))
missing_mapped <- setdiff(required_mapped_cols, colnames(mapped_all))
missing_expressed <- setdiff(required_expressed_cols, colnames(expressed))

if (length(missing_background) > 0) {
  stop(paste("Missing columns in background:", paste(missing_background, collapse = ", ")))
}

if (length(missing_mapped) > 0) {
  stop(paste("Missing columns in mapped_all:", paste(missing_mapped, collapse = ", ")))
}

if (length(missing_expressed) > 0) {
  stop(paste("Missing columns in expressed:", paste(missing_expressed, collapse = ", ")))
}

# ------------------------------------------------------------
# Build unique gene sets
# ------------------------------------------------------------

background_genes <- unique(background[, c("Gene stable ID", "hk_status")])
colnames(background_genes)[1] <- "gene_id"

mapped_genes <- unique(mapped_all[, c("gene_id", "hk_status")])

expressed_genes <- unique(expressed[, c("gene_id", "hk_status")])

# ------------------------------------------------------------
# Count HK / non-HK genes across steps
# ------------------------------------------------------------

count_by_status <- function(df, step_name) {
  out <- as.data.frame(table(df$hk_status), stringsAsFactors = FALSE)
  colnames(out) <- c("hk_status", "n_genes")
  out$step <- step_name
  out <- out[, c("step", "hk_status", "n_genes")]
  return(out)
}

step_counts <- rbind(
  count_by_status(background_genes, "full_ensembl115_hrt_background"),
  count_by_status(mapped_genes, "successfully_mapped_to_hpa_gtex"),
  count_by_status(expressed_genes, "expressed_in_at_least_one_tissue_ntpm_ge_1")
)

write.table(
  step_counts,
  file = file.path(table_output_dir, "hk_nonhk_gene_counts_across_hpa_gtex_steps.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------
# non-HK loss categories
# ------------------------------------------------------------

background_nonhk <- background_genes$gene_id[background_genes$hk_status == "non-HK"]
mapped_nonhk <- mapped_genes$gene_id[mapped_genes$hk_status == "non-HK"]
expressed_nonhk <- expressed_genes$gene_id[expressed_genes$hk_status == "non-HK"]

nonhk_not_mapped_to_hpa <- setdiff(background_nonhk, mapped_nonhk)

nonhk_mapped_but_never_expressed <- setdiff(mapped_nonhk, expressed_nonhk)

nonhk_loss_summary <- data.frame(
  category = c(
    "non-HK in full Ensembl115-HRT background",
    "non-HK successfully mapped to HPA-GTEx",
    "non-HK expressed in at least one tissue nTPM >= 1",
    "non-HK not mapped to HPA-GTEx",
    "non-HK mapped but never nTPM >= 1"
  ),
  n_genes = c(
    length(background_nonhk),
    length(mapped_nonhk),
    length(expressed_nonhk),
    length(nonhk_not_mapped_to_hpa),
    length(nonhk_mapped_but_never_expressed)
  )
)

write.table(
  nonhk_loss_summary,
  file = file.path(table_output_dir, "nonhk_gene_loss_after_hpa_gtex_mapping_and_ntpm_filter.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  data.frame(gene_id = nonhk_not_mapped_to_hpa),
  file = file.path(table_output_dir, "nonhk_not_mapped_to_hpa_gtex.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  data.frame(gene_id = nonhk_mapped_but_never_expressed),
  file = file.path(table_output_dir, "nonhk_mapped_but_never_ntpm_ge_1.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------
# Same check for HK genes, just for completeness
# ------------------------------------------------------------

background_hk <- background_genes$gene_id[background_genes$hk_status == "HK"]
mapped_hk <- mapped_genes$gene_id[mapped_genes$hk_status == "HK"]
expressed_hk <- expressed_genes$gene_id[expressed_genes$hk_status == "HK"]

hk_not_mapped_to_hpa <- setdiff(background_hk, mapped_hk)
hk_mapped_but_never_expressed <- setdiff(mapped_hk, expressed_hk)

hk_loss_summary <- data.frame(
  category = c(
    "HK in full Ensembl115-HRT background",
    "HK successfully mapped to HPA-GTEx",
    "HK expressed in at least one tissue nTPM >= 1",
    "HK not mapped to HPA-GTEx",
    "HK mapped but never nTPM >= 1"
  ),
  n_genes = c(
    length(background_hk),
    length(mapped_hk),
    length(expressed_hk),
    length(hk_not_mapped_to_hpa),
    length(hk_mapped_but_never_expressed)
  )
)

write.table(
  hk_loss_summary,
  file = file.path(table_output_dir, "hk_gene_loss_after_hpa_gtex_mapping_and_ntpm_filter.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------
# Console output
# ------------------------------------------------------------

cat("\nDone.\n")

cat("\nHK / non-HK gene counts across steps:\n")
print(step_counts)

cat("\nnon-HK loss summary:\n")
print(nonhk_loss_summary)

cat("\nHK loss summary:\n")
print(hk_loss_summary)

cat("\nSanity check:\n")
cat("non-HK missing from expressed subset should equal full non-HK minus expressed non-HK:\n")
cat(length(background_nonhk) - length(expressed_nonhk), "\n")

cat("\nnon-HK not mapped + non-HK mapped but never expressed:\n")
cat(length(nonhk_not_mapped_to_hpa) + length(nonhk_mapped_but_never_expressed), "\n")

cat("\nTables written to:\n")
cat(table_output_dir, "\n")