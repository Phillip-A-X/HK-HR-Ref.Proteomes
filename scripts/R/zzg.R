# ------------------------------------------------------------
# Density plots of log10 gene length by GTEx tissue
# ------------------------------------------------------------

library(ggplot2)
library(dplyr)

input_file <- "data/processed/hpa_gtex/HumanProteinCodingGenes_bytissue_onchromosomes.tsv"
figure_output_dir <- "results/hpa_gtex/figures"

dir.create(figure_output_dir, recursive = TRUE, showWarnings = FALSE)

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
    gene_length_bp > 0,
    nTPM >= 1
  ) %>%
  distinct(Tissue, gene_id, .keep_all = TRUE) %>%
  mutate(
    log10_gene_length_bp = log10(gene_length_bp)
  )

# optional:
# Tissue-Reihenfolge nach mean gene length
tissue_order <- df %>%
  group_by(Tissue) %>%
  summarise(
    mean_gene_length_bp = mean(gene_length_bp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mean_gene_length_bp) %>%
  pull(Tissue)

df$Tissue <- factor(df$Tissue, levels = tissue_order)

p_density_facets <- ggplot(
  df,
  aes(
    x = log10_gene_length_bp,
    y = after_stat(density)
  )
) +
  geom_density(
    fill = "grey40",
    color = "grey20",
    alpha = 0.85,
    linewidth = 0.3
  ) +
  facet_wrap(
    ~ Tissue,
    scales = "fixed",
    ncol = 6
  ) +
  labs(
    title = "Density plots of log10 gene length by GTEx tissue",
    subtitle = "Human protein-coding genes; expressed = nTPM >= 1",
    x = "log10(gene length bp)",
    y = "Density"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8)
  )

ggsave(
  filename = file.path(
    figure_output_dir,
    "hpa_gtex_log10_gene_length_density_by_tissue.png"
  ),
  plot = p_density_facets,
  width = 16,
  height = 12,
  dpi = 300
)

cat(
  "Figure written to:\n",
  file.path(
    figure_output_dir,
    "hpa_gtex_log10_gene_length_density_by_tissue.png"
  ),
  "\n"
)

# ------------------------------------------------------------
# Add theoretical normal distribution to each tissue density
# ------------------------------------------------------------

# Calculate tissue-specific mean and SD of log10 gene length
normal_parameters <- df %>%
  group_by(Tissue) %>%
  summarise(
    mean_log10_gene_length = mean(log10_gene_length_bp, na.rm = TRUE),
    sd_log10_gene_length = sd(log10_gene_length_bp, na.rm = TRUE),
    .groups = "drop"
  )


# Generate theoretical normal density for each tissue
theoretical_density <- do.call(
  rbind,
  lapply(seq_len(nrow(normal_parameters)), function(i) {
    
    tissue_i <- as.character(normal_parameters$Tissue[i])
    
    mu <- normal_parameters$mean_log10_gene_length[i]
    sigma <- normal_parameters$sd_log10_gene_length[i]
    
    observed_values <- df$log10_gene_length_bp[
      as.character(df$Tissue) == tissue_i
    ]
    
    x_values <- seq(
      from = min(observed_values, na.rm = TRUE),
      to = max(observed_values, na.rm = TRUE),
      length.out = 500
    )
    
    data.frame(
      Tissue = tissue_i,
      log10_gene_length_bp = x_values,
      theoretical_density = dnorm(
        x_values,
        mean = mu,
        sd = sigma
      )
    )
  })
)


# Keep same tissue order as original density plot
theoretical_density$Tissue <- factor(
  theoretical_density$Tissue,
  levels = tissue_order
)


# ------------------------------------------------------------
# Observed density + theoretical normal distribution
# ------------------------------------------------------------

p_density_theoretical <- ggplot() +
  
  # Observed gene-length density
  geom_density(
    data = df,
    aes(
      x = log10_gene_length_bp,
      y = after_stat(density)
    ),
    fill = "grey40",
    color = "grey20",
    alpha = 0.65,
    linewidth = 0.3
  ) +
  
  # Theoretical normal distribution
  geom_line(
    data = theoretical_density,
    aes(
      x = log10_gene_length_bp,
      y = theoretical_density
    ),
    linetype = "dashed",
    linewidth = 0.9
  ) +
  
  facet_wrap(
    ~ Tissue,
    scales = "fixed",
    ncol = 6
  ) +
  
  labs(
    title = "Observed vs theoretical log10 gene-length distributions by GTEx tissue",
    subtitle = "Dashed line = theoretical normal distribution based on tissue-specific mean and SD; expressed = nTPM >= 1",
    x = "log10(gene length bp)",
    y = "Density"
  ) +
  
  theme_bw() +
  
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8)
  )


# ------------------------------------------------------------
# Save new figure
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    figure_output_dir,
    "hpa_gtex_log10_gene_length_density_by_tissue_with_theoretical_normal.png"
  ),
  plot = p_density_theoretical,
  width = 16,
  height = 12,
  dpi = 300
)


cat(
  "Observed + theoretical density figure written to:\n",
  file.path(
    figure_output_dir,
    "hpa_gtex_log10_gene_length_density_by_tissue_with_theoretical_normal.png"
  ),
  "\n"
)