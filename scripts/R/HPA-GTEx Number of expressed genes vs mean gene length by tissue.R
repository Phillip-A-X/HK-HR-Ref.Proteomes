# ------------------------------------------------------------
# HPA-GTEx:
# Number of expressed genes vs mean gene length by tissue
#
# x = number of expressed genes per tissue
# y = mean gene length per tissue
#
# expressed = nTPM >= 1
# ------------------------------------------------------------

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

df$nTPM <- as.numeric(df$nTPM)
df$gene_length_bp <- as.numeric(df$gene_length_bp)

df <- df %>%
  filter(
    !is.na(gene_id),
    !is.na(Tissue),
    !is.na(nTPM),
    !is.na(gene_length_bp),
    gene_length_bp > 0,
    nTPM >= 1
  ) %>%
  distinct(Tissue, gene_id, .keep_all = TRUE)


# ------------------------------------------------------------
# Tissue summary
# ------------------------------------------------------------

tissue_gene_count <- df %>%
  group_by(Tissue) %>%
  summarise(
    n_expressed_genes = n(),
    mean_gene_length_bp = mean(gene_length_bp, na.rm = TRUE),
    .groups = "drop"
  )


# ------------------------------------------------------------
# Correlations
# ------------------------------------------------------------

pearson_test <- cor.test(
  tissue_gene_count$n_expressed_genes,
  tissue_gene_count$mean_gene_length_bp,
  method = "pearson"
)

spearman_test <- cor.test(
  tissue_gene_count$n_expressed_genes,
  tissue_gene_count$mean_gene_length_bp,
  method = "spearman",
  exact = FALSE
)

cat("\nNumber expressed genes vs mean gene length:\n")

cat(
  "Pearson r =",
  round(unname(pearson_test$estimate), 3),
  "\np =",
  format.pval(pearson_test$p.value, digits = 3),
  "\n"
)

cat(
  "Spearman rho =",
  round(unname(spearman_test$estimate), 3),
  "\np =",
  format.pval(spearman_test$p.value, digits = 3),
  "\n"
)


# ------------------------------------------------------------
# Save table
# ------------------------------------------------------------

write.table(
  tissue_gene_count,
  file = file.path(
    table_output_dir,
    "hpa_gtex_tissue_gene_count_vs_mean_gene_length.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

plot_label <- paste0(
  "Pearson r = ",
  round(unname(pearson_test$estimate), 3),
  "\np = ",
  format.pval(pearson_test$p.value, digits = 3)
)

p_gene_count <- ggplot(
  tissue_gene_count,
  aes(
    x = n_expressed_genes,
    y = mean_gene_length_bp
  )
) +
  geom_point(size = 3) +
  #geom_smooth(
   # method = "lm",
  #  se = TRUE
#  ) +
  geom_text(
    aes(label = Tissue),
    vjust = -0.7,
    size = 3,
    check_overlap = TRUE
  ) +
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = plot_label,
    hjust = 1.1,
    vjust = -0.5,
    size = 4
  ) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Number of expressed genes vs mean gene length across GTEx tissues",
    subtitle = "Human protein-coding genes; expressed = nTPM >= 1",
    x = "Number of expressed genes",
    y = "Mean gene length (bp)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold")
  )


ggsave(
  filename = file.path(
    figure_output_dir,
    "hpa_gtex_number_expressed_genes_vs_mean_gene_length.png"
  ),
  plot = p_gene_count,
  width = 9,
  height = 7,
  dpi = 300
)

cat(
  "\nFigure written to:\n",
  file.path(
    figure_output_dir,
    "hpa_gtex_number_expressed_genes_vs_mean_gene_length.png"
  ),
  "\n"
)