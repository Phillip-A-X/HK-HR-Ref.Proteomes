# ============================================================
# HPA-GTEx:
# Compare linear, power-law and exponential models
# for mean gene length vs gene-length variance
#
# One point = one GTEx tissue
# expressed = nTPM >= 1
#
# Models:
#   Linear:       variance = a + b * mean
#   Power law:    variance = a * mean^b
#   Exponential:  variance = a * exp(b * mean)
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

# If necessary:
# setwd("//wsl$/Ubuntu/home/phillip/projects/HK-HR-Ref.Proteomes")

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
  distinct(
    Tissue,
    gene_id,
    .keep_all = TRUE
  )


# ------------------------------------------------------------
# Tissue summary
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
    # Scaling only for numerical stability
    mean_gene_length_kb = mean_gene_length_bp / 1000,
    variance_1e10 = variance_gene_length_bp2 / 1e10
  )


# ============================================================
# MODEL 1: Linear
# ============================================================

model_linear <- lm(
  variance_1e10 ~ mean_gene_length_kb,
  data = tissue_summary
)


# ============================================================
# MODEL 2: Power law
#
# y = a * x^b
#
# Starting values obtained from log-log model
# ============================================================

power_start <- lm(
  log(variance_1e10) ~ log(mean_gene_length_kb),
  data = tissue_summary
)

start_a_power <- exp(coef(power_start)[1])
start_b_power <- coef(power_start)[2]

model_power <- nls(
  variance_1e10 ~ a * mean_gene_length_kb^b,
  data = tissue_summary,
  start = list(
    a = start_a_power,
    b = start_b_power
  )
)


# ============================================================
# MODEL 3: Exponential
#
# y = a * exp(b*x)
#
# Starting values obtained from log(y) ~ x
# ============================================================

exp_start <- lm(
  log(variance_1e10) ~ mean_gene_length_kb,
  data = tissue_summary
)

start_a_exp <- exp(coef(exp_start)[1])
start_b_exp <- coef(exp_start)[2]

model_exponential <- nls(
  variance_1e10 ~ a * exp(b * mean_gene_length_kb),
  data = tissue_summary,
  start = list(
    a = start_a_exp,
    b = start_b_exp
  )
)


# ============================================================
# Model comparison
# ============================================================

observed <- tissue_summary$variance_1e10

pred_linear <- predict(model_linear)
pred_power <- predict(model_power)
pred_exponential <- predict(model_exponential)

sst <- sum(
  (observed - mean(observed))^2
)


# ------------------------------------------------------------
# Function for raw-scale R2 and RMSE
# ------------------------------------------------------------

model_metrics <- function(observed, predicted) {
  
  rss <- sum(
    (observed - predicted)^2
  )
  
  r_squared <- 1 - rss / sst
  
  rmse <- sqrt(
    mean(
      (observed - predicted)^2
    )
  )
  
  return(
    c(
      RSS = rss,
      R2 = r_squared,
      RMSE = rmse
    )
  )
}


metrics_linear <- model_metrics(
  observed,
  pred_linear
)

metrics_power <- model_metrics(
  observed,
  pred_power
)

metrics_exponential <- model_metrics(
  observed,
  pred_exponential
)


# ------------------------------------------------------------
# Comparison table
# ------------------------------------------------------------

model_comparison <- data.frame(
  model = c(
    "Linear",
    "Power law",
    "Exponential"
  ),
  
  R2 = c(
    metrics_linear["R2"],
    metrics_power["R2"],
    metrics_exponential["R2"]
  ),
  
  RMSE = c(
    metrics_linear["RMSE"],
    metrics_power["RMSE"],
    metrics_exponential["RMSE"]
  ),
  
  AIC = c(
    AIC(model_linear),
    AIC(model_power),
    AIC(model_exponential)
  )
)

model_comparison <- model_comparison %>%
  arrange(AIC)


# ------------------------------------------------------------
# Save model comparison
# ------------------------------------------------------------

write.table(
  model_comparison,
  file = file.path(
    table_output_dir,
    "hpa_gtex_mean_variance_model_comparison.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# Model equations
# ============================================================

linear_coef <- coef(model_linear)
power_coef <- coef(model_power)
exp_coef <- coef(model_exponential)

cat("\n========================================\n")
cat("MODEL EQUATIONS\n")
cat("========================================\n")

cat(
  "\nLinear:\n",
  "variance =",
  round(linear_coef[1], 4),
  "+",
  round(linear_coef[2], 4),
  "* mean\n"
)

cat(
  "\nPower law:\n",
  "variance =",
  signif(power_coef["a"], 4),
  "* mean^",
  round(power_coef["b"], 4),
  "\n"
)

cat(
  "\nExponential:\n",
  "variance =",
  signif(exp_coef["a"], 4),
  "* exp(",
  round(exp_coef["b"], 6),
  "* mean)\n"
)


# ============================================================
# Generate fitted curves
# ============================================================

curve_df <- data.frame(
  mean_gene_length_kb = seq(
    min(tissue_summary$mean_gene_length_kb),
    max(tissue_summary$mean_gene_length_kb),
    length.out = 500
  )
)

curve_df$Linear <- predict(
  model_linear,
  newdata = curve_df
)

curve_df$Power_law <- predict(
  model_power,
  newdata = curve_df
)

curve_df$Exponential <- predict(
  model_exponential,
  newdata = curve_df
)


# Convert to long format without additional packages
curve_long <- bind_rows(
  data.frame(
    mean_gene_length_kb = curve_df$mean_gene_length_kb,
    predicted_variance = curve_df$Linear,
    Model = "Linear"
  ),
  data.frame(
    mean_gene_length_kb = curve_df$mean_gene_length_kb,
    predicted_variance = curve_df$Power_law,
    Model = "Power law"
  ),
  data.frame(
    mean_gene_length_kb = curve_df$mean_gene_length_kb,
    predicted_variance = curve_df$Exponential,
    Model = "Exponential"
  )
)


# ============================================================
# Plot all three models
# ============================================================

p_models <- ggplot() +
  
  geom_point(
    data = tissue_summary,
    aes(
      x = mean_gene_length_kb,
      y = variance_1e10
    ),
    size = 3
  ) +
  
  geom_text(
    data = tissue_summary,
    aes(
      x = mean_gene_length_kb,
      y = variance_1e10,
      label = Tissue
    ),
    vjust = -0.7,
    size = 3,
    check_overlap = TRUE
  ) +
  
  geom_line(
    data = curve_long,
    aes(
      x = mean_gene_length_kb,
      y = predicted_variance,
      linetype = Model
    ),
    linewidth = 1
  ) +
  
  labs(
    title = "Mean gene length vs gene-length variance across GTEx tissues",
    subtitle = "Comparison of linear, power-law and exponential fits; expressed = nTPM >= 1",
    x = "Mean gene length (kb)",
    y = expression(
      "Gene-length variance ("*10^10*" bp"^2*")"
    ),
    linetype = "Model"
  ) +
  
  theme_bw() +
  
  theme(
    plot.title = element_text(face = "bold")
  )


ggsave(
  filename = file.path(
    figure_output_dir,
    "hpa_gtex_mean_variance_linear_power_exponential_comparison.png"
  ),
  plot = p_models,
  width = 9,
  height = 7,
  dpi = 300
)


# ============================================================
# Console output
# ============================================================

cat("\n========================================\n")
cat("MODEL COMPARISON\n")
cat("========================================\n\n")

print(model_comparison)

cat(
  "\nPower-law exponent:",
  round(power_coef["b"], 4),
  "\n"
)

cat(
  "\nBest model by lowest AIC:",
  model_comparison$model[1],
  "\n"
)

cat(
  "\nFigure written to:\n",
  file.path(
    figure_output_dir,
    "hpa_gtex_mean_variance_linear_power_exponential_comparison.png"
  ),
  "\n"
)