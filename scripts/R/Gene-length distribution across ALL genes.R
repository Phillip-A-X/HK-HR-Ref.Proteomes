# ============================================================
# Human Ensembl 115 genomic background
# Gene-length distribution across ALL protein-coding genes
#
# Background:
#   Ensembl 115
#   protein-coding
#   chromosomes 1-22, X, Y
#
# Input:
#   data/processed/human_all_genes_with_hrt_atlas_hk_status_rescued.tsv
#
# Outputs:
#   results/figures/human_gene_background_gene_length_density_log10.png
#   results/figures/human_gene_background_gene_length_histogram_log10.png
# ============================================================


# ------------------------------------------------------------
# Setup
# ------------------------------------------------------------

library(ggplot2)
library(dplyr)

setwd("//wsl$/Ubuntu/home/phillip/projects/HK-HR-Ref.Proteomes")

input_file <-
  "data/processed/human_all_genes_with_hrt_atlas_hk_status_rescued.tsv"

figure_output_dir <- "results/figures"

dir.create(
  figure_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

df <- read.delim(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# Keep one row per gene and valid gene lengths
# ------------------------------------------------------------

df <- df %>%
  filter(
    !is.na(`Gene stable ID`),
    !is.na(gene_length_bp),
    gene_length_bp > 0
  ) %>%
  
  distinct(
    `Gene stable ID`,
    .keep_all = TRUE
  ) %>%
  
  mutate(
    log10_gene_length_bp = log10(gene_length_bp)
  )


# ------------------------------------------------------------
# QC
# ------------------------------------------------------------

cat(
  "\nNumber of unique genes:",
  n_distinct(df$`Gene stable ID`),
  "\n"
)

cat(
  "Minimum gene length:",
  min(df$gene_length_bp),
  "bp\n"
)

cat(
  "Median gene length:",
  median(df$gene_length_bp),
  "bp\n"
)

cat(
  "Mean gene length:",
  mean(df$gene_length_bp),
  "bp\n"
)

cat(
  "Maximum gene length:",
  max(df$gene_length_bp),
  "bp\n"
)


# ============================================================
# 1. Density plot
# ============================================================

p_density <- ggplot(
  df,
  aes(x = log10_gene_length_bp)
) +
  
  geom_density(
    fill = "grey50",
    color = "grey20",
    alpha = 0.8,
    linewidth = 0.5
  ) +
  
  labs(
    title = "Gene-length distribution of the human genomic background",
    subtitle = "Ensembl 115 protein-coding genes on chromosomes 1-22, X and Y",
    x = "log10 gene length (bp)",
    y = "Density"
  ) +
  
  theme_bw() +
  
  theme(
    plot.title = element_text(face = "bold")
  )


ggsave(
  filename = file.path(
    figure_output_dir,
    "human_gene_background_gene_length_density_log10.png"
  ),
  plot = p_density,
  width = 8,
  height = 6,
  dpi = 300
)


# ============================================================
# 2. Histogram / actual number of genes
# ============================================================

p_histogram <- ggplot(
  df,
  aes(x = log10_gene_length_bp)
) +
  
  geom_histogram(
    binwidth = 0.05,
    boundary = 0,
    fill = "grey55",
    color = "grey20",
    linewidth = 0.2
  ) +
  
  labs(
    title = "Gene-length distribution of the human genomic background",
    subtitle = "Ensembl 115 protein-coding genes on chromosomes 1-22, X and Y",
    x = "log10 gene length (bp)",
    y = "Number of genes"
  ) +
  
  theme_bw() +
  
  theme(
    plot.title = element_text(face = "bold")
  )


ggsave(
  filename = file.path(
    figure_output_dir,
    "human_gene_background_gene_length_histogram_log10.png"
  ),
  plot = p_histogram,
  width = 8,
  height = 6,
  dpi = 300
)


# ------------------------------------------------------------
# Finished
# ------------------------------------------------------------

cat(
  "\nFigures written to:\n",
  file.path(
    figure_output_dir,
    "human_gene_background_gene_length_density_log10.png"
  ),
  "\n",
  file.path(
    figure_output_dir,
    "human_gene_background_gene_length_histogram_log10.png"
  ),
  "\n"
)