library(dplyr)
library(ggplot2)

# ------------------------------------------------------------
# Human mitochondrial protein-coding gene length analysis
#
# Input:
#   Human nuclear + mitochondrial protein-coding gene table
#
# MT genes:
#   MT = TRUE
#   hk_status = NA
#
# Outputs:
#   Summary table
#   MT genes ordered by gene length
#   Ordered gene-length plot
#   Gene-length distribution with theoretical lognormal curve
# ------------------------------------------------------------


# Input
input_file <-
  "data/processed/human_all_genes_with_hrt_atlas_hk_status_rescued_plus_mt.tsv"

# Output directories
table_output_dir <- "results/tables"
figure_output_dir <- "results/figures"

dir.create(
  table_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figure_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# Read data
# ------------------------------------------------------------

human_genes <- read.delim(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# Extract MT genes
# ------------------------------------------------------------

mt_genes <- human_genes %>%
  filter(MT == TRUE) %>%
  mutate(
    gene_id_short = paste0(
      "...",
      substr(
        `Gene stable ID`,
        nchar(`Gene stable ID`) - 5,
        nchar(`Gene stable ID`)
      )
    )
  )


# QC
cat("\n--- MT gene QC ---\n")
cat("Number of MT genes:", nrow(mt_genes), "\n")
cat(
  "All MT hk_status NA:",
  all(is.na(mt_genes$hk_status)),
  "\n"
)


# ------------------------------------------------------------
# Summary statistics
# ------------------------------------------------------------

mt_summary <- mt_genes %>%
  summarise(
    n_genes = n(),
    mean_gene_length_bp = mean(gene_length_bp),
    median_gene_length_bp = median(gene_length_bp),
    min_gene_length_bp = min(gene_length_bp),
    max_gene_length_bp = max(gene_length_bp)
  )


cat("\n--- MT gene length summary ---\n")
print(mt_summary)


write.table(
  mt_summary,
  file = file.path(
    table_output_dir,
    "human_mt_gene_length_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# Order MT genes by gene length
# ------------------------------------------------------------

mt_genes_ordered <- mt_genes %>%
  arrange(gene_length_bp)


write.table(
  mt_genes_ordered,
  file = file.path(
    table_output_dir,
    "human_mt_genes_ordered_by_gene_length.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)


cat("\n--- MT genes ordered by gene length ---\n")

print(
  mt_genes_ordered %>%
    select(
      `Gene stable ID`,
      gene_id_short,
      gene_length_bp
    )
)


# ------------------------------------------------------------
# Figure 1
# MT genes ordered by increasing gene length
# ------------------------------------------------------------

mt_genes_ordered <- mt_genes_ordered %>%
  mutate(
    gene_id_short = factor(
      gene_id_short,
      levels = gene_id_short
    )
  )


p_ordered <- ggplot(
  mt_genes_ordered,
  aes(
    x = gene_id_short,
    y = gene_length_bp
  )
) +
  geom_col(
    width = 0.7
  ) +
  labs(
    x = "Gene ID",
    y = "Gene length (bp)",
    title = "Human mitochondrial protein-coding gene lengths"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )


ggsave(
  filename = file.path(
    figure_output_dir,
    "human_mt_gene_lengths_ordered.png"
  ),
  plot = p_ordered,
  width = 8,
  height = 5,
  dpi = 300
)


# ------------------------------------------------------------
# Figure 2
# MT gene-length distribution
#
# Histogram = observed number of genes
# Dashed line = theoretical lognormal distribution
# fitted from observed log gene lengths
# ------------------------------------------------------------

binwidth <- 200

log_mean <- mean(
  log(mt_genes$gene_length_bp)
)

log_sd <- sd(
  log(mt_genes$gene_length_bp)
)


x_values <- seq(
  min(mt_genes$gene_length_bp),
  max(mt_genes$gene_length_bp),
  length.out = 500
)


lognormal_curve <- data.frame(
  gene_length_bp = x_values,
  expected_count =
    dlnorm(
      x_values,
      meanlog = log_mean,
      sdlog = log_sd
    ) *
    nrow(mt_genes) *
    binwidth
)


p_distribution <- ggplot(
  mt_genes,
  aes(x = gene_length_bp)
) +
  geom_histogram(
    binwidth = binwidth,
    boundary = 0,
    color = "black",
    fill = "grey70"
  ) +
  geom_line(
    data = lognormal_curve,
    aes(
      x = gene_length_bp,
      y = expected_count
    ),
    inherit.aes = FALSE,
    linewidth = 1,
    linetype = "dashed"
  ) +
  labs(
    x = "Gene length (bp)",
    y = "Number of genes",
    title = "Human mitochondrial protein-coding gene-length distribution"
  ) +
  theme_classic(base_size = 12)


ggsave(
  filename = file.path(
    figure_output_dir,
    "human_mt_gene_length_distribution_lognormal.png"
  ),
  plot = p_distribution,
  width = 7,
  height = 5,
  dpi = 300
)


# ------------------------------------------------------------
# Done
# ------------------------------------------------------------

cat("\nSaved tables:\n")

cat(
  file.path(
    table_output_dir,
    "human_mt_gene_length_summary.tsv"
  ),
  "\n"
)

cat(
  file.path(
    table_output_dir,
    "human_mt_genes_ordered_by_gene_length.tsv"
  ),
  "\n"
)


cat("\nSaved figures:\n")

cat(
  file.path(
    figure_output_dir,
    "human_mt_gene_lengths_ordered.png"
  ),
  "\n"
)

cat(
  file.path(
    figure_output_dir,
    "human_mt_gene_length_distribution_lognormal.png"
  ),
  "\n"
)