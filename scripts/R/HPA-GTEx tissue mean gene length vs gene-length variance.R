# HPA-GTEx tissue mean gene length vs gene-length variance
#
# Input:
#   data/processed/hpa_gtex/HumanProteinCodingGenes_bytissue_onchromosomes.tsv
#
# One point per GTEx tissue:
#   x = mean gene length of expressed genes
#   y = variance of gene length of expressed genes
#
# Outputs:
#   results/hpa_gtex_hk/tables/hpa_gtex_tissue_gene_length_mean_variance.tsv
#   results/hpa_gtex_hk/tables/hpa_gtex_tissue_gene_length_mean_variance_correlations.tsv
#   results/hpa_gtex_hk/figures/hpa_gtex_tissue_mean_vs_variance_raw.png
#   results/hpa_gtex_hk/figures/hpa_gtex_tissue_mean_vs_variance_log10.png


# ------------------------------------------------------------
# Setup
# ------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(scales)

input_file <- "data/processed/hpa_gtex/HumanProteinCodingGenes_bytissue_onchromosomes.tsv"

table_output_dir <- "results/hpa_gtex_hk/tables"
figure_output_dir <- "results/hpa_gtex_hk/figures"

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

required_columns <- c(
  "gene_id",
  "Tissue",
  "nTPM",
  "gene_length_bp"
)

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
    gene_length_bp > 0,
    nTPM >= 1
  )


# ------------------------------------------------------------
# Safety check:
# one gene only once per tissue
# ------------------------------------------------------------

df_tissue_unique <- df %>%
  distinct(Tissue, gene_id, .keep_all = TRUE)


# ------------------------------------------------------------
# Tissue-level summary
# ------------------------------------------------------------

tissue_summary <- df_tissue_unique %>%
  group_by(Tissue) %>%
  summarise(
    n_genes = n(),
    mean_gene_length_bp = mean(gene_length_bp, na.rm = TRUE),
    variance_gene_length_bp2 = var(gene_length_bp, na.rm = TRUE),
    sd_gene_length_bp = sd(gene_length_bp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    log10_mean_gene_length_bp = log10(mean_gene_length_bp),
    log10_variance_gene_length_bp2 = log10(variance_gene_length_bp2)
  )


# ------------------------------------------------------------
# Save tissue summary
# ------------------------------------------------------------

write.table(
  tissue_summary,
  file = file.path(
    table_output_dir,
    "hpa_gtex_tissue_gene_length_mean_variance.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# Correlation tests
# ------------------------------------------------------------

raw_pearson <- cor.test(
  tissue_summary$mean_gene_length_bp,
  tissue_summary$variance_gene_length_bp2,
  method = "pearson"
)

raw_spearman <- cor.test(
  tissue_summary$mean_gene_length_bp,
  tissue_summary$variance_gene_length_bp2,
  method = "spearman",
  exact = FALSE
)

log_pearson <- cor.test(
  tissue_summary$log10_mean_gene_length_bp,
  tissue_summary$log10_variance_gene_length_bp2,
  method = "pearson"
)

log_spearman <- cor.test(
  tissue_summary$log10_mean_gene_length_bp,
  tissue_summary$log10_variance_gene_length_bp2,
  method = "spearman",
  exact = FALSE
)


correlation_results <- data.frame(
  scale = c(
    "raw",
    "raw",
    "log10-log10",
    "log10-log10"
  ),
  method = c(
    "Pearson",
    "Spearman",
    "Pearson",
    "Spearman"
  ),
  correlation = c(
    unname(raw_pearson$estimate),
    unname(raw_spearman$estimate),
    unname(log_pearson$estimate),
    unname(log_spearman$estimate)
  ),
  p_value = c(
    raw_pearson$p.value,
    raw_spearman$p.value,
    log_pearson$p.value,
    log_spearman$p.value
  )
)


write.table(
  correlation_results,
  file = file.path(
    table_output_dir,
    "hpa_gtex_tissue_gene_length_mean_variance_correlations.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# Labels for plots
# ------------------------------------------------------------

raw_label <- paste0(
  "Pearson r = ",
  round(unname(raw_pearson$estimate), 3),
  "\np = ",
  format.pval(raw_pearson$p.value, digits = 3)
)

log_label <- paste0(
  "Pearson r = ",
  round(unname(log_pearson$estimate), 3),
  "\np = ",
  format.pval(log_pearson$p.value, digits = 3)
)


# ------------------------------------------------------------
# Plot 1: raw scale
# ------------------------------------------------------------

p_raw <- ggplot(
  tissue_summary,
  aes(
    x = mean_gene_length_bp,
    y = variance_gene_length_bp2
  )
) +
  geom_point(size = 3) +
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
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
    label = raw_label,
    hjust = 1.1,
    vjust = -0.5,
    size = 4
  ) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = scientific) +
  labs(
    title = "Mean gene length vs gene-length variance across GTEx tissues",
    subtitle = "Each point represents one tissue; expressed = nTPM >= 1",
    x = "Mean gene length (bp)",
    y = expression("Gene-length variance (bp"^2*")")
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = file.path(
    figure_output_dir,
    "hpa_gtex_tissue_mean_vs_variance_raw.png"
  ),
  plot = p_raw,
  width = 9,
  height = 7,
  dpi = 300
)


# ------------------------------------------------------------
# Plot 2: log10-log10 scale
# IMPORTANT PLOT
# ------------------------------------------------------------

p_log <- ggplot(
  tissue_summary,
  aes(
    x = log10_mean_gene_length_bp,
    y = log10_variance_gene_length_bp2
  )
) +
  geom_point(size = 3) +
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
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
    label = log_label,
    hjust = 1.1,
    vjust = -0.5,
    size = 4
  ) +
  labs(
    title = "Mean gene length vs gene-length variance across GTEx tissues",
    subtitle = "log10-log10 scale; each point represents one tissue; expressed = nTPM >= 1",
    x = "log10(mean gene length bp)",
    y = expression("log10(gene-length variance bp"^2*")")
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = file.path(
    figure_output_dir,
    "hpa_gtex_tissue_mean_vs_variance_log10.png"
  ),
  plot = p_log,
  width = 9,
  height = 7,
  dpi = 300
)


# ------------------------------------------------------------
# Console output
# ------------------------------------------------------------

cat("\nDone.\n")

cat("\nTissue summary:\n")
print(tissue_summary)

cat("\nCorrelation results:\n")
print(correlation_results)

cat("\nPrimary result: log10-log10 Pearson correlation\n")
print(log_pearson)

cat("\nFigures written to:\n")
cat(figure_output_dir, "\n")

# ============================================================
# Additional check:
# Mean-variance scaling and relative gene-length variability
# ============================================================


# ------------------------------------------------------------
# 1. Log-log regression slope
# ------------------------------------------------------------

log_model <- lm(
  log10_variance_gene_length_bp2 ~ log10_mean_gene_length_bp,
  data = tissue_summary
)

log_slope <- coef(log_model)[["log10_mean_gene_length_bp"]]

log_slope_p <- summary(log_model)$coefficients[
  "log10_mean_gene_length_bp",
  "Pr(>|t|)"
]

log_r_squared <- summary(log_model)$r.squared


cat("\n----------------------------------------\n")
cat("LOG-LOG REGRESSION\n")
cat("----------------------------------------\n")

cat(
  "Slope =",
  round(log_slope, 3),
  "\nR-squared =",
  round(log_r_squared, 3),
  "\np =",
  format.pval(log_slope_p, digits = 3),
  "\n"
)


# ------------------------------------------------------------
# 2. Calculate coefficient of variation
# CV = SD / mean
# ------------------------------------------------------------

tissue_summary <- tissue_summary %>%
  mutate(
    cv_gene_length = sd_gene_length_bp / mean_gene_length_bp
  )


# ------------------------------------------------------------
# 3. Correlation:
# mean gene length vs coefficient of variation
# ------------------------------------------------------------

cv_pearson <- cor.test(
  tissue_summary$mean_gene_length_bp,
  tissue_summary$cv_gene_length,
  method = "pearson"
)

cv_spearman <- cor.test(
  tissue_summary$mean_gene_length_bp,
  tissue_summary$cv_gene_length,
  method = "spearman",
  exact = FALSE
)


cat("\n----------------------------------------\n")
cat("MEAN GENE LENGTH VS COEFFICIENT OF VARIATION\n")
cat("----------------------------------------\n")

cat(
  "Pearson r =",
  round(unname(cv_pearson$estimate), 3),
  "\np =",
  format.pval(cv_pearson$p.value, digits = 3),
  "\n"
)

cat(
  "Spearman rho =",
  round(unname(cv_spearman$estimate), 3),
  "\np =",
  format.pval(cv_spearman$p.value, digits = 3),
  "\n"
)


# ------------------------------------------------------------
# 4. Plot mean gene length vs CV
# ------------------------------------------------------------

cv_label <- paste0(
  "Pearson r = ",
  round(unname(cv_pearson$estimate), 3),
  "\np = ",
  format.pval(cv_pearson$p.value, digits = 3)
)


p_cv <- ggplot(
  tissue_summary,
  aes(
    x = mean_gene_length_bp,
    y = cv_gene_length
  )
) +
  geom_point(size = 3) +
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
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
    label = cv_label,
    hjust = 1.1,
    vjust = -0.5,
    size = 4
  ) +
  scale_x_continuous(labels = scales::comma) +
  labs(
    title = "Mean gene length vs relative gene-length variability across GTEx tissues",
    subtitle = "Coefficient of variation = SD / mean; expressed = nTPM >= 1",
    x = "Mean gene length (bp)",
    y = "Coefficient of variation of gene length"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold")
  )


ggsave(
  filename = file.path(
    figure_output_dir,
    "hpa_gtex_tissue_mean_vs_gene_length_cv.png"
  ),
  plot = p_cv,
  width = 9,
  height = 7,
  dpi = 300
)


# ------------------------------------------------------------
# 5. Save updated tissue summary including CV
# ------------------------------------------------------------

write.table(
  tissue_summary,
  file = file.path(
    table_output_dir,
    "hpa_gtex_tissue_gene_length_mean_variance_cv.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


cat("\nAdditional analysis finished.\n")
cat(
  "CV figure:",
  file.path(
    figure_output_dir,
    "hpa_gtex_tissue_mean_vs_gene_length_cv.png"
  ),
  "\n"
)