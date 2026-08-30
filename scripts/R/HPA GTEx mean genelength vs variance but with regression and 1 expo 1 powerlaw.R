# ============================================================
# HPA-GTEx mean gene length vs gene-length variance
#
# 1) log10-log10 regression:
#    log10(variance) = intercept + slope * log10(mean)
#
# 2) exponential regression:
#    ln(variance) = intercept + slope * mean
#
#    back-transformed:
#    variance = a * exp(slope * mean)
#
# expressed = nTPM >= 1
# ============================================================


# ------------------------------------------------------------
# Setup
# ------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(scales)

input_file <- "data/processed/hpa_gtex/HumanProteinCodingGenes_bytissue_onchromosomes.tsv"

figure_output_dir <- "results/hpa_gtex/figures"
table_output_dir <- "results/hpa_gtex/tables"

dir.create(
  figure_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  table_output_dir,
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

df$nTPM <- as.numeric(df$nTPM)
df$gene_length_bp <- as.numeric(df$gene_length_bp)


# ------------------------------------------------------------
# Filter expressed genes
# ------------------------------------------------------------

df <- df %>%
  filter(
    !is.na(gene_id),
    !is.na(Tissue),
    !is.na(nTPM),
    !is.na(gene_length_bp),
    gene_length_bp > 0,
    nTPM >= 1
  ) %>%
  distinct(
    Tissue,
    gene_id,
    .keep_all = TRUE
  )


# ------------------------------------------------------------
# Tissue-level summary
# ------------------------------------------------------------

tissue_summary <- df %>%
  group_by(Tissue) %>%
  summarise(
    n_genes = n(),
    mean_gene_length_bp = mean(gene_length_bp),
    variance_gene_length_bp2 = var(gene_length_bp),
    .groups = "drop"
  ) %>%
  mutate(
    log10_mean_gene_length_bp =
      log10(mean_gene_length_bp),
    
    log10_variance_gene_length_bp2 =
      log10(variance_gene_length_bp2),
    
    ln_variance_gene_length_bp2 =
      log(variance_gene_length_bp2)
  )


# ============================================================
# PART 1
# log10-log10 regression
# ============================================================

model_loglog <- lm(
  log10_variance_gene_length_bp2 ~
    log10_mean_gene_length_bp,
  data = tissue_summary
)

loglog_summary <- summary(model_loglog)

loglog_intercept <- coef(model_loglog)[1]
loglog_slope <- coef(model_loglog)[2]
loglog_r2 <- loglog_summary$r.squared


# ------------------------------------------------------------
# Label for plot
# ------------------------------------------------------------

loglog_label <- paste0(
  "log10(variance) = ",
  round(loglog_intercept, 3),
  " + ",
  round(loglog_slope, 3),
  " × log10(mean)",
  "\nSlope = ",
  round(loglog_slope, 3),
  "\nY-intercept = ",
  round(loglog_intercept, 3),
  "\nR² = ",
  round(loglog_r2, 3)
)


# ------------------------------------------------------------
# Plot log10-log10
# ------------------------------------------------------------

p_loglog <- ggplot(
  tissue_summary,
  aes(
    x = log10_mean_gene_length_bp,
    y = log10_variance_gene_length_bp2
  )
) +
  
  geom_point(
    size = 3
  ) +
  
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 1
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
    label = loglog_label,
    hjust = 1.05,
    vjust = -0.3,
    size = 3.8
  ) +
  
  labs(
    title = "Mean gene length vs gene-length variance across GTEx tissues",
    subtitle = "log10-log10 regression; expressed = nTPM >= 1",
    x = "log10 mean gene length (bp)",
    y = expression(
      "log10 gene-length variance (bp"^2*")"
    )
  ) +
  
  theme_bw() +
  
  theme(
    plot.title = element_text(face = "bold")
  )


ggsave(
  filename = file.path(
    figure_output_dir,
    "hpa_gtex_tissue_mean_vs_variance_log10_regression_fit.png"
  ),
  plot = p_loglog,
  width = 9,
  height = 7,
  dpi = 300
)


# ============================================================
# PART 2
# Exponential regression
#
# ln(variance) = intercept + slope * mean
#
# therefore:
#
# variance = exp(intercept) *
#            exp(slope * mean)
# ============================================================

model_exponential <- lm(
  ln_variance_gene_length_bp2 ~
    mean_gene_length_bp,
  data = tissue_summary
)

exp_summary <- summary(model_exponential)

exp_intercept <- coef(model_exponential)[1]
exp_slope <- coef(model_exponential)[2]
exp_r2 <- exp_summary$r.squared

exp_a <- exp(exp_intercept)


# ------------------------------------------------------------
# Generate predicted exponential curve
# ------------------------------------------------------------

exp_curve <- data.frame(
  mean_gene_length_bp = seq(
    min(tissue_summary$mean_gene_length_bp),
    max(tissue_summary$mean_gene_length_bp),
    length.out = 500
  )
)

exp_curve$predicted_variance_bp2 <-
  exp_a *
  exp(
    exp_slope *
      exp_curve$mean_gene_length_bp
  )


# ------------------------------------------------------------
# Label
# ------------------------------------------------------------

exp_label <- paste0(
  "ln(variance) = ",
  round(exp_intercept, 3),
  " + ",
  format(
    exp_slope,
    scientific = TRUE,
    digits = 3
  ),
  " × mean",
  "\nR² = ",
  round(exp_r2, 3),
  "\n",
  "variance = ",
  format(
    exp_a,
    scientific = TRUE,
    digits = 3
  ),
  " × exp(",
  format(
    exp_slope,
    scientific = TRUE,
    digits = 3
  ),
  " × mean)"
)


# ------------------------------------------------------------
# Raw-scale plot with predicted exponential
# ------------------------------------------------------------

p_exponential <- ggplot() +
  
  geom_point(
    data = tissue_summary,
    aes(
      x = mean_gene_length_bp,
      y = variance_gene_length_bp2
    ),
    size = 3
  ) +
  
  geom_line(
    data = exp_curve,
    aes(
      x = mean_gene_length_bp,
      y = predicted_variance_bp2
    ),
    linewidth = 1
  ) +
  
  geom_text(
    data = tissue_summary,
    aes(
      x = mean_gene_length_bp,
      y = variance_gene_length_bp2,
      label = Tissue
    ),
    vjust = -0.7,
    size = 3,
    check_overlap = TRUE
  ) +
  
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = exp_label,
    hjust = 1.05,
    vjust = -0.3,
    size = 3.6
  ) +
  
  scale_x_continuous(
    labels = comma
  ) +
  
  scale_y_continuous(
    labels = scientific
  ) +
  
  labs(
    title = "Mean gene length vs gene-length variance across GTEx tissues",
    subtitle = "Raw scale with predicted exponential fit; expressed = nTPM >= 1",
    x = "Mean gene length (bp)",
    y = expression(
      "Gene-length variance (bp"^2*")"
    )
  ) +
  
  theme_bw() +
  
  theme(
    plot.title = element_text(face = "bold")
  )


ggsave(
  filename = file.path(
    figure_output_dir,
    "hpa_gtex_tissue_mean_vs_variance_exponential_fit.png"
  ),
  plot = p_exponential,
  width = 9,
  height = 7,
  dpi = 300
)


# ============================================================
# Save regression statistics
# ============================================================

regression_summary <- data.frame(
  
  model = c(
    "log10-log10",
    "exponential"
  ),
  
  intercept = c(
    loglog_intercept,
    exp_intercept
  ),
  
  slope = c(
    loglog_slope,
    exp_slope
  ),
  
  R2 = c(
    loglog_r2,
    exp_r2
  )
)


write.table(
  regression_summary,
  file = file.path(
    table_output_dir,
    "hpa_gtex_mean_variance_regression_fits.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# Console output
# ============================================================

cat("\n========================================\n")
cat("LOG10-LOG10 REGRESSION\n")
cat("========================================\n")

cat(
  "\nSlope:",
  loglog_slope
)

cat(
  "\nY-intercept:",
  loglog_intercept
)

cat(
  "\nR²:",
  loglog_r2,
  "\n"
)

cat(
  "\nEquation:\nlog10(variance) =",
  loglog_intercept,
  "+",
  loglog_slope,
  "* log10(mean)\n"
)


cat("\n========================================\n")
cat("EXPONENTIAL REGRESSION\n")
cat("========================================\n")

cat(
  "\nSlope:",
  exp_slope
)

cat(
  "\nY-intercept of ln regression:",
  exp_intercept
)

cat(
  "\nR²:",
  exp_r2,
  "\n"
)

cat(
  "\nRaw-scale equation:\nvariance =",
  exp_a,
  "* exp(",
  exp_slope,
  "* mean)\n"
)


cat("\n========================================\n")
cat("DONE\n")
cat("========================================\n")