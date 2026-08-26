# HPA-GTEx mean gene length vs gene-length variance
# Raw-scale plot with power-law fit
#
# Input:
#   data/processed/hpa_gtex/HumanProteinCodingGenes_bytissue_onchromosomes.tsv
#
# Expression definition:
#   expressed in tissue = nTPM >= 1
#
# Output:
#   results/hpa_gtex/figures/
#   hpa_gtex_tissue_mean_vs_variance_raw_powerlaw.png


# ------------------------------------------------------------
# Setup
# ------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(scales)

# Only needed if working directory is not already the project root:
# setwd("//wsl$/Ubuntu/home/phillip/projects/HK-HR-Ref.Proteomes")

input_file <- "data/processed/hpa_gtex/HumanProteinCodingGenes_bytissue_onchromosomes.tsv"

figure_output_dir <- "results/hpa_gtex/figures"

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
# Tissue-level mean and variance
# ------------------------------------------------------------

tissue_summary <- df %>%
  group_by(Tissue) %>%
  summarise(
    n_genes = n(),
    mean_gene_length_bp = mean(
      gene_length_bp,
      na.rm = TRUE
    ),
    variance_gene_length_bp2 = var(
      gene_length_bp,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    log10_mean_gene_length_bp =
      log10(mean_gene_length_bp),
    
    log10_variance_gene_length_bp2 =
      log10(variance_gene_length_bp2)
  )


# ------------------------------------------------------------
# Fit log-log model
#
# log10(variance) =
# intercept + slope * log10(mean)
# ------------------------------------------------------------

log_model <- lm(
  log10_variance_gene_length_bp2 ~
    log10_mean_gene_length_bp,
  data = tissue_summary
)

model_summary <- summary(log_model)

log_intercept <- coef(log_model)[["(Intercept)"]]

log_slope <- coef(log_model)[[
  "log10_mean_gene_length_bp"
]]

r_squared <- model_summary$r.squared

p_value <- model_summary$coefficients[
  "log10_mean_gene_length_bp",
  "Pr(>|t|)"
]


# ------------------------------------------------------------
# Back-transform to raw scale
#
# variance = a * mean^b
# ------------------------------------------------------------

power_a <- 10^log_intercept
power_b <- log_slope


# ------------------------------------------------------------
# Generate fitted power-law curve
# ------------------------------------------------------------

power_curve <- data.frame(
  mean_gene_length_bp = seq(
    min(tissue_summary$mean_gene_length_bp),
    max(tissue_summary$mean_gene_length_bp),
    length.out = 500
  )
)

power_curve$predicted_variance_bp2 <-
  power_a *
  power_curve$mean_gene_length_bp^power_b


# ------------------------------------------------------------
# Plot label
# ------------------------------------------------------------

power_label <- paste0(
  "Power-law exponent = ",
  round(power_b, 3),
  "\nR² = ",
  round(r_squared, 3),
  "\np = ",
  format.pval(p_value, digits = 3)
)


# ------------------------------------------------------------
# Plot:
# raw mean gene length vs raw variance
# ------------------------------------------------------------

p_power <- ggplot() +
  
  geom_point(
    data = tissue_summary,
    aes(
      x = mean_gene_length_bp,
      y = variance_gene_length_bp2
    ),
    size = 3
  ) +
  
  geom_line(
    data = power_curve,
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
    label = power_label,
    hjust = 1.1,
    vjust = -0.5,
    size = 4
  ) +
  
  scale_x_continuous(
    labels = comma
  ) +
  
  scale_y_continuous(
    labels = scientific
  ) +
  
  labs(
    title = "Mean gene length vs gene-length variance across GTEx tissues",
    subtitle = "Raw scale with power-law fit; expressed = nTPM >= 1",
    x = "Mean gene length (bp)",
    y = expression(
      "Gene-length variance (bp"^2*")"
    )
  ) +
  
  theme_bw() +
  
  theme(
    plot.title = element_text(face = "bold")
  )


# ------------------------------------------------------------
# Save figure
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    figure_output_dir,
    "hpa_gtex_tissue_mean_vs_variance_raw_powerlaw.png"
  ),
  plot = p_power,
  width = 9,
  height = 7,
  dpi = 300
)


# ------------------------------------------------------------
# Console output
# ------------------------------------------------------------

cat("\nDone.\n")

cat(
  "\nPower-law exponent =",
  round(power_b, 3),
  "\n"
)

cat(
  "R-squared =",
  round(r_squared, 3),
  "\n"
)

cat(
  "p =",
  format.pval(p_value, digits = 3),
  "\n"
)

cat(
  "\nPower-law equation:\nvariance =",
  format(power_a, scientific = TRUE),
  "* mean^",
  round(power_b, 3),
  "\n"
)

cat(
  "\nFigure written to:\n",
  file.path(
    figure_output_dir,
    "hpa_gtex_tissue_mean_vs_variance_raw_powerlaw.png"
  ),
  "\n"
)