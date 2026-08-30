# ============================================================
# HPA-GTEx:
# Number of expressed genes by tissue across nTPM cutoffs
# Tissues ordered by ascending mean gene length at nTPM >= 1
# ============================================================

library(ggplot2)
library(dplyr)
library(scales)

input_file <- "data/processed/hpa_gtex/HumanProteinCodingGenes_bytissue_onchromosomes.tsv"

table_output_dir <- "results/hpa_gtex/tables"
figure_output_dir <- "results/hpa_gtex/figures"

dir.create(table_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

df <- read.delim(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c("gene_id", "Tissue", "nTPM", "gene_length_bp")

missing_columns <- setdiff(required_columns, colnames(df))

if (length(missing_columns) > 0) {
  stop(
    paste(
      "Missing required columns:",
      paste(missing_columns, collapse = ", ")
    )
  )
}

df$nTPM <- as.numeric(df$nTPM)
df$gene_length_bp <- as.numeric(df$gene_length_bp)

df <- df %>%
  filter(
    !is.na(gene_id),
    !is.na(Tissue),
    !is.na(nTPM),
    !is.na(gene_length_bp),
    gene_length_bp > 0
  ) %>%
  distinct(Tissue, gene_id, .keep_all = TRUE)

# ------------------------------------------------------------
# Cutoffs
# ------------------------------------------------------------

cutoffs <- c(1, 5, 10, 20)

# ------------------------------------------------------------
# Summary function
# ------------------------------------------------------------

summarise_cutoff <- function(data, cutoff_value) {
  data %>%
    filter(nTPM >= cutoff_value) %>%
    group_by(Tissue) %>%
    summarise(
      n_expressed_genes = n(),
      mean_gene_length_bp = mean(gene_length_bp, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(cutoff = cutoff_value)
}

# ------------------------------------------------------------
# Run for all cutoffs
# ------------------------------------------------------------

all_tissue_results <- bind_rows(
  lapply(cutoffs, function(x) summarise_cutoff(df, x))
)

# ------------------------------------------------------------
# Define tissue order by ascending mean gene length at nTPM >= 1
# ------------------------------------------------------------

tissue_order_table <- all_tissue_results %>%
  filter(cutoff == 1) %>%
  arrange(mean_gene_length_bp) %>%
  mutate(order_rank = row_number())

tissue_order <- tissue_order_table$Tissue

all_tissue_results$Tissue <- factor(
  all_tissue_results$Tissue,
  levels = tissue_order
)

all_tissue_results$cutoff_label <- factor(
  paste0("nTPM >= ", all_tissue_results$cutoff),
  levels = paste0("nTPM >= ", cutoffs)
)

# ------------------------------------------------------------
# Save ordered table
# ------------------------------------------------------------

write.table(
  tissue_order_table,
  file = file.path(
    table_output_dir,
    "hpa_gtex_tissue_order_by_mean_gene_length_ntpm1.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  all_tissue_results,
  file = file.path(
    table_output_dir,
    "hpa_gtex_tissue_gene_counts_cutoff_ordered_by_ntpm1.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

p_gene_counts_ordered <- ggplot(
  all_tissue_results,
  aes(
    x = Tissue,
    y = n_expressed_genes,
    group = cutoff_label
  )
) +
  geom_line(
    aes(linetype = cutoff_label),
    linewidth = 0.7
  ) +
  geom_point(
    aes(shape = cutoff_label),
    size = 2
  ) +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Number of expressed genes by GTEx tissue across nTPM cutoffs",
    subtitle = "Tissues ordered by ascending mean gene length at nTPM >= 1",
    x = "GTEx tissue",
    y = "Number of expressed genes",
    linetype = "Expression cutoff",
    shape = "Expression cutoff"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 60, hjust = 1)
  )

ggsave(
  filename = file.path(
    figure_output_dir,
    "hpa_gtex_cutoff_sensitivity_gene_counts_ordered_by_ntpm1.png"
  ),
  plot = p_gene_counts_ordered,
  width = 14,
  height = 7,
  dpi = 300
)

cat("\nDone.\n")

cat("\nFigure written to:\n")
cat(
  file.path(
    figure_output_dir,
    "hpa_gtex_cutoff_sensitivity_gene_counts_ordered_by_ntpm1.png"
  ),
  "\n"
)

cat("\nTables written to:\n")
cat(
  file.path(
    table_output_dir,
    "hpa_gtex_tissue_order_by_mean_gene_length_ntpm1.tsv"
  ),
  "\n"
)
cat(
  file.path(
    table_output_dir,
    "hpa_gtex_tissue_gene_counts_cutoff_ordered_by_ntpm1.tsv"
  ),
  "\n"
)