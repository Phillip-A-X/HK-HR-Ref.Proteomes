# ============================================================
# HPA-GTEx:
# Separate mean gene length vs gene-length variance plots
# for different nTPM cutoffs
#
# One plot per cutoff:
#   x = log10(mean gene length)
#   y = log10(gene-length variance)
#
# Cutoffs:
#   nTPM >= 1, 5, 10, 20
#
# Output:
#   results/hpa_gtex/figures/
# ============================================================


# ------------------------------------------------------------
# Setup
# ------------------------------------------------------------

library(ggplot2)
library(dplyr)

# If needed:
# setwd("//wsl$/Ubuntu/home/phillip/projects/HK-HR-Ref.Proteomes")

input_file <- "data/processed/hpa_gtex/HumanProteinCodingGenes_bytissue_onchromosomes.tsv"

figure_output_dir <- "results/hpa_gtex/figures"
table_output_dir  <- "results/hpa_gtex/tables"

dir.create(figure_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_output_dir, recursive = TRUE, showWarnings = FALSE)


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
# Storage for regression summary
# ------------------------------------------------------------

regression_summary <- data.frame()


# ------------------------------------------------------------
# Loop over cutoffs
# ------------------------------------------------------------

for (cutoff_i in cutoffs) {
  
  df_cutoff <- df %>%
    filter(nTPM >= cutoff_i)
  
  tissue_summary <- df_cutoff %>%
    group_by(Tissue) %>%
    summarise(
      n_expressed_genes = n(),
      mean_gene_length_bp = mean(gene_length_bp, na.rm = TRUE),
      variance_gene_length_bp2 = var(gene_length_bp, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      log10_mean_gene_length_bp = log10(mean_gene_length_bp),
      log10_variance_gene_length_bp2 = log10(variance_gene_length_bp2)
    )
  
  # Regression
  model <- lm(
    log10_variance_gene_length_bp2 ~ log10_mean_gene_length_bp,
    data = tissue_summary
  )
  
  model_summary <- summary(model)
  
  slope <- coef(model)[["log10_mean_gene_length_bp"]]
  intercept <- coef(model)[["(Intercept)"]]
  r_squared <- model_summary$r.squared
  p_value <- model_summary$coefficients[
    "log10_mean_gene_length_bp",
    "Pr(>|t|)"
  ]
  
  regression_summary <- bind_rows(
    regression_summary,
    data.frame(
      cutoff = cutoff_i,
      intercept = intercept,
      slope = slope,
      R2 = r_squared,
      p_value = p_value
    )
  )
  
  plot_label <- paste0(
    "Slope = ", round(slope, 3),
    "\nY-intercept = ", round(intercept, 3),
    "\nR² = ", round(r_squared, 3),
    "\np = ", format.pval(p_value, digits = 3)
  )
  
  p <- ggplot(
    tissue_summary,
    aes(
      x = log10_mean_gene_length_bp,
      y = log10_variance_gene_length_bp2
    )
  ) +
    geom_point(size = 3) +
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
      label = plot_label,
      hjust = 1.05,
      vjust = -0.3,
      size = 3.8
    ) +
    labs(
      title = "Mean gene length vs gene-length variance across GTEx tissues",
      subtitle = paste0("log10-log10 scale; expressed = nTPM >= ", cutoff_i),
      x = "log10(mean gene length bp)",
      y = expression("log10(gene-length variance bp"^2*")")
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold")
    )
  
  output_file <- file.path(
    figure_output_dir,
    paste0("hpa_gtex_mean_vs_variance_log10_ntpm", cutoff_i, ".png")
  )
  
  ggsave(
    filename = output_file,
    plot = p,
    width = 9,
    height = 7,
    dpi = 300
  )
  
  cat("Written:", output_file, "\n")
}


# ------------------------------------------------------------
# Save regression summary table
# ------------------------------------------------------------

write.table(
  regression_summary,
  file = file.path(
    table_output_dir,
    "hpa_gtex_mean_vs_variance_log10_regression_by_cutoff.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# Done
# ------------------------------------------------------------

cat("\nDone.\n")

cat("\nRegression summary written to:\n")
cat(
  file.path(
    table_output_dir,
    "hpa_gtex_mean_vs_variance_log10_regression_by_cutoff.tsv"
  ),
  "\n"
)