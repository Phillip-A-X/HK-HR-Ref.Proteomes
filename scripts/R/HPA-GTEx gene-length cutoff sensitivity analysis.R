# ============================================================
# HPA-GTEx gene-length cutoff sensitivity analysis
#
# Cutoffs tested:
#   nTPM >= 1
#   nTPM >= 5
#   nTPM >= 10
#   nTPM >= 20
#
# Input:
#   data/processed/hpa_gtex/
#   HumanProteinCodingGenes_bytissue_onchromosomes.tsv
#
# Outputs:
#   results/hpa_gtex/tables/
#   results/hpa_gtex/figures/
# ============================================================


# ------------------------------------------------------------
# Setup
# ------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(scales)

input_file <- paste0(
  "data/processed/hpa_gtex/",
  "HumanProteinCodingGenes_bytissue_onchromosomes.tsv"
)

table_output_dir <- "results/hpa_gtex/tables"
figure_output_dir <- "results/hpa_gtex/figures"

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

missing_columns <- setdiff(
  required_columns,
  colnames(df)
)

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
  distinct(
    Tissue,
    gene_id,
    .keep_all = TRUE
  )


# ------------------------------------------------------------
# Cutoffs
# ------------------------------------------------------------

cutoffs <- c(1, 5, 10, 20)


# ------------------------------------------------------------
# Functions
# ------------------------------------------------------------

analyse_cutoff <- function(data, cutoff_value) {
  
  df_cutoff <- data %>%
    filter(
      nTPM >= cutoff_value
    )
  
  tissue_summary <- df_cutoff %>%
    group_by(Tissue) %>%
    summarise(
      n_expressed_genes = n(),
      mean_gene_length_bp =
        mean(gene_length_bp, na.rm = TRUE),
      variance_gene_length_bp2 =
        var(gene_length_bp, na.rm = TRUE),
      sd_gene_length_bp =
        sd(gene_length_bp, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      cutoff = cutoff_value,
      cv_gene_length =
        sd_gene_length_bp / mean_gene_length_bp,
      log10_mean_gene_length_bp =
        log10(mean_gene_length_bp),
      log10_variance_gene_length_bp2 =
        log10(variance_gene_length_bp2)
    )
  
  return(tissue_summary)
}


# ------------------------------------------------------------
# Run tissue-level analysis for all cutoffs
# ------------------------------------------------------------

all_tissue_results <- bind_rows(
  lapply(
    cutoffs,
    function(x) analyse_cutoff(df, x)
  )
)


# ------------------------------------------------------------
# Save complete tissue-level table
# ------------------------------------------------------------

write.table(
  all_tissue_results,
  file = file.path(
    table_output_dir,
    "hpa_gtex_cutoff_sensitivity_tissue_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# Correlation / regression analysis
# ============================================================

cutoff_statistics <- data.frame()


for (cutoff_i in cutoffs) {
  
  d <- all_tissue_results %>%
    filter(cutoff == cutoff_i)
  
  # ----------------------------------------------------------
  # Mean vs raw variance
  # ----------------------------------------------------------
  
  raw_pearson <- cor.test(
    d$mean_gene_length_bp,
    d$variance_gene_length_bp2,
    method = "pearson"
  )
  
  raw_spearman <- cor.test(
    d$mean_gene_length_bp,
    d$variance_gene_length_bp2,
    method = "spearman",
    exact = FALSE
  )
  
  
  # ----------------------------------------------------------
  # log10(mean) vs log10(variance)
  # ----------------------------------------------------------
  
  log_pearson <- cor.test(
    d$log10_mean_gene_length_bp,
    d$log10_variance_gene_length_bp2,
    method = "pearson"
  )
  
  log_spearman <- cor.test(
    d$log10_mean_gene_length_bp,
    d$log10_variance_gene_length_bp2,
    method = "spearman",
    exact = FALSE
  )
  
  
  # ----------------------------------------------------------
  # log-log linear model
  # ----------------------------------------------------------
  
  log_model <- lm(
    log10_variance_gene_length_bp2 ~
      log10_mean_gene_length_bp,
    data = d
  )
  
  log_model_summary <- summary(log_model)
  
  log_slope <- coef(log_model)[[
    "log10_mean_gene_length_bp"
  ]]
  
  log_r_squared <-
    log_model_summary$r.squared
  
  log_slope_p <-
    log_model_summary$coefficients[
      "log10_mean_gene_length_bp",
      "Pr(>|t|)"
    ]
  
  
  # ----------------------------------------------------------
  # Mean vs CV
  # ----------------------------------------------------------
  
  cv_pearson <- cor.test(
    d$mean_gene_length_bp,
    d$cv_gene_length,
    method = "pearson"
  )
  
  cv_spearman <- cor.test(
    d$mean_gene_length_bp,
    d$cv_gene_length,
    method = "spearman",
    exact = FALSE
  )
  
  
  # ----------------------------------------------------------
  # Gene count vs mean gene length
  # ----------------------------------------------------------
  
  count_pearson <- cor.test(
    d$n_expressed_genes,
    d$mean_gene_length_bp,
    method = "pearson"
  )
  
  count_spearman <- cor.test(
    d$n_expressed_genes,
    d$mean_gene_length_bp,
    method = "spearman",
    exact = FALSE
  )
  
  
  # ----------------------------------------------------------
  # Add results
  # ----------------------------------------------------------
  
  cutoff_statistics <- bind_rows(
    cutoff_statistics,
    data.frame(
      cutoff = cutoff_i,
      
      raw_mean_variance_pearson_r =
        unname(raw_pearson$estimate),
      raw_mean_variance_pearson_p =
        raw_pearson$p.value,
      
      raw_mean_variance_spearman_rho =
        unname(raw_spearman$estimate),
      raw_mean_variance_spearman_p =
        raw_spearman$p.value,
      
      log_mean_variance_pearson_r =
        unname(log_pearson$estimate),
      log_mean_variance_pearson_p =
        log_pearson$p.value,
      
      log_mean_variance_spearman_rho =
        unname(log_spearman$estimate),
      log_mean_variance_spearman_p =
        log_spearman$p.value,
      
      log_log_slope =
        log_slope,
      
      log_log_R2 =
        log_r_squared,
      
      log_log_slope_p =
        log_slope_p,
      
      mean_cv_pearson_r =
        unname(cv_pearson$estimate),
      mean_cv_pearson_p =
        cv_pearson$p.value,
      
      mean_cv_spearman_rho =
        unname(cv_spearman$estimate),
      mean_cv_spearman_p =
        cv_spearman$p.value,
      
      count_mean_pearson_r =
        unname(count_pearson$estimate),
      count_mean_pearson_p =
        count_pearson$p.value,
      
      count_mean_spearman_rho =
        unname(count_spearman$estimate),
      count_mean_spearman_p =
        count_spearman$p.value
    )
  )
}


# ------------------------------------------------------------
# Save statistics
# ------------------------------------------------------------

write.table(
  cutoff_statistics,
  file = file.path(
    table_output_dir,
    "hpa_gtex_cutoff_sensitivity_statistics.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# Additional summary:
# gene retention by cutoff
# ============================================================

gene_retention <- all_tissue_results %>%
  select(
    Tissue,
    cutoff,
    n_expressed_genes
  ) %>%
  group_by(Tissue) %>%
  mutate(
    genes_at_cutoff1 =
      n_expressed_genes[cutoff == 1],
    fraction_vs_cutoff1 =
      n_expressed_genes / genes_at_cutoff1,
    percent_vs_cutoff1 =
      100 * fraction_vs_cutoff1
  ) %>%
  ungroup()


write.table(
  gene_retention,
  file = file.path(
    table_output_dir,
    "hpa_gtex_cutoff_sensitivity_gene_retention.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# PLOTS
# ============================================================

all_tissue_results$cutoff_label <- factor(
  paste0("nTPM >= ", all_tissue_results$cutoff),
  levels = paste0("nTPM >= ", cutoffs)
)


# ------------------------------------------------------------
# Plot 1:
# Number of expressed genes by cutoff
# ------------------------------------------------------------

p_gene_count_cutoff <- ggplot(
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
  scale_y_continuous(
    labels = comma
  ) +
  labs(
    title =
      "Number of expressed genes by GTEx tissue and nTPM cutoff",
    subtitle =
      "Human protein-coding genes",
    x = "GTEx tissue",
    y = "Number of expressed genes",
    linetype = "Expression cutoff",
    shape = "Expression cutoff"
  ) +
  theme_bw() +
  theme(
    plot.title =
      element_text(face = "bold"),
    axis.text.x =
      element_text(
        angle = 60,
        hjust = 1
      )
  )

ggsave(
  filename = file.path(
    figure_output_dir,
    "hpa_gtex_cutoff_sensitivity_gene_counts.png"
  ),
  plot = p_gene_count_cutoff,
  width = 14,
  height = 7,
  dpi = 300
)


# ------------------------------------------------------------
# Plot 2:
# Mean gene length by cutoff
# ------------------------------------------------------------

p_mean_cutoff <- ggplot(
  all_tissue_results,
  aes(
    x = Tissue,
    y = mean_gene_length_bp,
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
  scale_y_continuous(
    labels = comma
  ) +
  labs(
    title =
      "Mean gene length by GTEx tissue and nTPM cutoff",
    subtitle =
      "Human protein-coding genes",
    x = "GTEx tissue",
    y = "Mean gene length (bp)",
    linetype = "Expression cutoff",
    shape = "Expression cutoff"
  ) +
  theme_bw() +
  theme(
    plot.title =
      element_text(face = "bold"),
    axis.text.x =
      element_text(
        angle = 60,
        hjust = 1
      )
  )

ggsave(
  filename = file.path(
    figure_output_dir,
    "hpa_gtex_cutoff_sensitivity_mean_gene_length.png"
  ),
  plot = p_mean_cutoff,
  width = 14,
  height = 7,
  dpi = 300
)


# ------------------------------------------------------------
# Plot 3:
# log mean vs log variance by cutoff
# ------------------------------------------------------------

p_log_variance_cutoff <- ggplot(
  all_tissue_results,
  aes(
    x = log10_mean_gene_length_bp,
    y = log10_variance_gene_length_bp2
  )
) +
  geom_point(size = 2) +
  geom_smooth(
    method = "lm",
    se = FALSE
  ) +
  facet_wrap(
    ~ cutoff_label,
    ncol = 2
  ) +
  geom_text(
    aes(label = Tissue),
    size = 2.5,
    vjust = -0.6,
    check_overlap = TRUE
  ) +
  labs(
    title =
      "Mean gene length vs gene-length variance across nTPM cutoffs",
    subtitle =
      "log10-log10 scale; each point represents one GTEx tissue",
    x =
      "log10(mean gene length bp)",
    y =
      expression(
        "log10(gene-length variance bp"^2*")"
      )
  ) +
  theme_bw() +
  theme(
    plot.title =
      element_text(face = "bold")
  )

ggsave(
  filename = file.path(
    figure_output_dir,
    "hpa_gtex_cutoff_sensitivity_log_mean_vs_variance.png"
  ),
  plot = p_log_variance_cutoff,
  width = 12,
  height = 9,
  dpi = 300
)


# ------------------------------------------------------------
# Plot 4:
# Mean gene length vs CV by cutoff
# ------------------------------------------------------------

p_cv_cutoff <- ggplot(
  all_tissue_results,
  aes(
    x = mean_gene_length_bp,
    y = cv_gene_length
  )
) +
  geom_point(size = 2) +
  geom_smooth(
    method = "lm",
    se = FALSE
  ) +
  facet_wrap(
    ~ cutoff_label,
    ncol = 2
  ) +
  geom_text(
    aes(label = Tissue),
    size = 2.5,
    vjust = -0.6,
    check_overlap = TRUE
  ) +
  labs(
    title =
      "Mean gene length vs relative gene-length variability across nTPM cutoffs",
    subtitle =
      "CV = SD / mean; each point represents one GTEx tissue",
    x =
      "Mean gene length (bp)",
    y =
      "Coefficient of variation of gene length"
  ) +
  theme_bw() +
  theme(
    plot.title =
      element_text(face = "bold")
  )

ggsave(
  filename = file.path(
    figure_output_dir,
    "hpa_gtex_cutoff_sensitivity_mean_vs_cv.png"
  ),
  plot = p_cv_cutoff,
  width = 12,
  height = 9,
  dpi = 300
)


# ------------------------------------------------------------
# Plot 5:
# Number expressed genes vs mean gene length
# ------------------------------------------------------------

p_count_mean_cutoff <- ggplot(
  all_tissue_results,
  aes(
    x = n_expressed_genes,
    y = mean_gene_length_bp
  )
) +
  geom_point(size = 2) +
  geom_smooth(
    method = "lm",
    se = FALSE
  ) +
  facet_wrap(
    ~ cutoff_label,
    ncol = 2
  ) +
  geom_text(
    aes(label = Tissue),
    size = 2.5,
    vjust = -0.6,
    check_overlap = TRUE
  ) +
  scale_x_continuous(
    labels = comma
  ) +
  scale_y_continuous(
    labels = comma
  ) +
  labs(
    title =
      "Number of expressed genes vs mean gene length across nTPM cutoffs",
    subtitle =
      "Each point represents one GTEx tissue",
    x =
      "Number of expressed genes",
    y =
      "Mean gene length (bp)"
  ) +
  theme_bw() +
  theme(
    plot.title =
      element_text(face = "bold")
  )

ggsave(
  filename = file.path(
    figure_output_dir,
    "hpa_gtex_cutoff_sensitivity_gene_count_vs_mean.png"
  ),
  plot = p_count_mean_cutoff,
  width = 12,
  height = 9,
  dpi = 300
)


# ============================================================
# Console output
# ============================================================

cat("\n========================================\n")
cat("CUTOFF SENSITIVITY ANALYSIS FINISHED\n")
cat("========================================\n")

cat("\nCutoffs tested:\n")
print(cutoffs)

cat("\nMain statistics:\n")
print(cutoff_statistics)

cat("\nGene retention preview:\n")
print(head(gene_retention, 20))

cat("\nTables written to:\n")
cat(table_output_dir, "\n")

cat("\nFigures written to:\n")
cat(figure_output_dir, "\n")