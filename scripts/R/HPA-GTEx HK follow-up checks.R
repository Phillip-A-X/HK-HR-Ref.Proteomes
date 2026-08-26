# HPA-GTEx HK follow-up checks
#
# Purpose:
#   1. Effect-size table for HK vs broad non-HK
#   2. HK detection per tissue
#   3. Fixed-cutoff sensitivity check for non-HK breadth groups
#
# Input:
#   data/processed/hpa_gtex/HumanProteinCodingGenes_bytissue_onchromosomes.tsv
#
# Definition:
#   expressed in tissue = nTPM >= 1
#
# Outputs:
#   results/hpa_gtex_hk/tables/
#   results/hpa_gtex_hk/figures/

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
  "gene_length_bp",
  "log10_gene_length_bp",
  "hk_status"
)

missing_columns <- setdiff(required_columns, colnames(df))

if (length(missing_columns) > 0) {
  stop(paste("Missing required columns:", paste(missing_columns, collapse = ", ")))
}

df$nTPM <- as.numeric(df$nTPM)
df$gene_length_bp <- as.numeric(df$gene_length_bp)
df$log10_gene_length_bp <- as.numeric(df$log10_gene_length_bp)

# ------------------------------------------------------------
# Gene-level summary
# ------------------------------------------------------------

gene_summary <- df %>%
  group_by(gene_id, hk_status) %>%
  summarise(
    n_tissues_expressed = n_distinct(Tissue),
    mean_nTPM_across_expressed_tissues = mean(nTPM, na.rm = TRUE),
    median_nTPM_across_expressed_tissues = median(nTPM, na.rm = TRUE),
    variance_log10_nTPM_plus1 = var(log10(nTPM + 1), na.rm = TRUE),
    sd_log10_nTPM_plus1 = sd(log10(nTPM + 1), na.rm = TRUE),
    gene_length_bp = first(gene_length_bp),
    log10_gene_length_bp = first(log10_gene_length_bp),
    .groups = "drop"
  )

gene_summary$variance_log10_nTPM_plus1[is.na(gene_summary$variance_log10_nTPM_plus1)] <- 0
gene_summary$sd_log10_nTPM_plus1[is.na(gene_summary$sd_log10_nTPM_plus1)] <- 0

# ------------------------------------------------------------
# 1. Effect sizes: HK vs broad non-HK
# ------------------------------------------------------------

# Quartile-based broad non-HK definition:
# broad non-HK = upper quartile of non-HK expression breadth.
# In the current dataset Q3 was 36, but it is recalculated here.

nonhk_breadth <- gene_summary %>%
  filter(hk_status == "non-HK") %>%
  pull(n_tissues_expressed)

q3_nonhk <- as.numeric(quantile(nonhk_breadth, 0.75, na.rm = TRUE, type = 7))

effect_data <- gene_summary %>%
  mutate(
    effect_group = case_when(
      hk_status == "HK" ~ "HK",
      hk_status == "non-HK" & n_tissues_expressed >= q3_nonhk ~ "broad non-HK",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(effect_group))

effect_summary <- effect_data %>%
  group_by(effect_group) %>%
  summarise(
    n_genes = n(),
    mean_gene_length_bp = mean(gene_length_bp, na.rm = TRUE),
    median_gene_length_bp = median(gene_length_bp, na.rm = TRUE),
    mean_log10_gene_length_bp = mean(log10_gene_length_bp, na.rm = TRUE),
    median_log10_gene_length_bp = median(log10_gene_length_bp, na.rm = TRUE),
    mean_expression_breadth = mean(n_tissues_expressed, na.rm = TRUE),
    median_expression_breadth = median(n_tissues_expressed, na.rm = TRUE),
    mean_sd_log10_nTPM_plus1 = mean(sd_log10_nTPM_plus1, na.rm = TRUE),
    median_sd_log10_nTPM_plus1 = median(sd_log10_nTPM_plus1, na.rm = TRUE),
    mean_nTPM_across_expressed_tissues = mean(mean_nTPM_across_expressed_tissues, na.rm = TRUE),
    median_nTPM_across_expressed_tissues = median(median_nTPM_across_expressed_tissues, na.rm = TRUE),
    .groups = "drop"
  )

hk_row <- effect_summary %>% filter(effect_group == "HK")
broad_row <- effect_summary %>% filter(effect_group == "broad non-HK")

effect_size_table <- data.frame(
  comparison = "HK vs broad non-HK",
  variable = c(
    "mean_gene_length_bp",
    "median_gene_length_bp",
    "mean_log10_gene_length_bp",
    "median_log10_gene_length_bp",
    "mean_expression_breadth",
    "median_expression_breadth",
    "mean_sd_log10_nTPM_plus1",
    "median_sd_log10_nTPM_plus1",
    "mean_nTPM_across_expressed_tissues",
    "median_nTPM_across_expressed_tissues"
  ),
  HK = c(
    hk_row$mean_gene_length_bp,
    hk_row$median_gene_length_bp,
    hk_row$mean_log10_gene_length_bp,
    hk_row$median_log10_gene_length_bp,
    hk_row$mean_expression_breadth,
    hk_row$median_expression_breadth,
    hk_row$mean_sd_log10_nTPM_plus1,
    hk_row$median_sd_log10_nTPM_plus1,
    hk_row$mean_nTPM_across_expressed_tissues,
    hk_row$median_nTPM_across_expressed_tissues
  ),
  broad_nonHK = c(
    broad_row$mean_gene_length_bp,
    broad_row$median_gene_length_bp,
    broad_row$mean_log10_gene_length_bp,
    broad_row$median_log10_gene_length_bp,
    broad_row$mean_expression_breadth,
    broad_row$median_expression_breadth,
    broad_row$mean_sd_log10_nTPM_plus1,
    broad_row$median_sd_log10_nTPM_plus1,
    broad_row$mean_nTPM_across_expressed_tissues,
    broad_row$median_nTPM_across_expressed_tissues
  )
)

effect_size_table$difference_HK_minus_broad_nonHK <- effect_size_table$HK - effect_size_table$broad_nonHK
effect_size_table$ratio_HK_div_broad_nonHK <- effect_size_table$HK / effect_size_table$broad_nonHK

write.table(
  effect_summary,
  file = file.path(table_output_dir, "hpa_gtex_hk_vs_broad_nonhk_effect_group_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  effect_size_table,
  file = file.path(table_output_dir, "hpa_gtex_hk_vs_broad_nonhk_effect_sizes.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------
# 2. HK detection per tissue
# ------------------------------------------------------------

total_hk_genes <- gene_summary %>%
  filter(hk_status == "HK") %>%
  summarise(n = n_distinct(gene_id)) %>%
  pull(n)

hk_detection_per_tissue <- df %>%
  filter(hk_status == "HK") %>%
  group_by(Tissue) %>%
  summarise(
    hk_genes_detected = n_distinct(gene_id),
    .groups = "drop"
  ) %>%
  mutate(
    total_hk_genes = total_hk_genes,
    percent_hk_genes_detected = 100 * hk_genes_detected / total_hk_genes
  ) %>%
  arrange(desc(percent_hk_genes_detected))

write.table(
  hk_detection_per_tissue,
  file = file.path(table_output_dir, "hpa_gtex_hk_detection_per_tissue.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

p_hk_detection <- ggplot(
  hk_detection_per_tissue,
  aes(
    x = reorder(Tissue, percent_hk_genes_detected),
    y = percent_hk_genes_detected
  )
) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Detection of HRT Atlas HK genes across GTEx tissues",
    subtitle = "Detected = nTPM >= 1 in the respective tissue",
    x = "GTEx tissue",
    y = "Detected HRT-HK genes (%)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = file.path(figure_output_dir, "hpa_gtex_hk_detection_per_tissue.png"),
  plot = p_hk_detection,
  width = 8,
  height = 9,
  dpi = 300
)

# ------------------------------------------------------------
# 3. Fixed-cutoff sensitivity check
# ------------------------------------------------------------

# Fixed cutoffs:
# restricted non-HK: expressed in <= 5 tissues
# intermediate non-HK: expressed in 6-29 tissues
# broad non-HK: expressed in >= 30 tissues
# HK: HRT Atlas HK genes

fixed_summary <- gene_summary %>%
  mutate(
    fixed_breadth_group = case_when(
      hk_status == "HK" ~ "HK",
      hk_status == "non-HK" & n_tissues_expressed <= 5 ~ "restricted non-HK fixed",
      hk_status == "non-HK" & n_tissues_expressed >= 30 ~ "broad non-HK fixed",
      hk_status == "non-HK" ~ "intermediate non-HK fixed",
      TRUE ~ "unknown"
    )
  )

fixed_summary$fixed_breadth_group <- factor(
  fixed_summary$fixed_breadth_group,
  levels = c(
    "restricted non-HK fixed",
    "intermediate non-HK fixed",
    "broad non-HK fixed",
    "HK"
  )
)

fixed_group_summary <- fixed_summary %>%
  group_by(fixed_breadth_group) %>%
  summarise(
    n_genes = n(),
    mean_expression_breadth = mean(n_tissues_expressed, na.rm = TRUE),
    median_expression_breadth = median(n_tissues_expressed, na.rm = TRUE),
    min_expression_breadth = min(n_tissues_expressed, na.rm = TRUE),
    max_expression_breadth = max(n_tissues_expressed, na.rm = TRUE),
    mean_gene_length_bp = mean(gene_length_bp, na.rm = TRUE),
    median_gene_length_bp = median(gene_length_bp, na.rm = TRUE),
    sd_gene_length_bp = sd(gene_length_bp, na.rm = TRUE),
    mean_log10_gene_length_bp = mean(log10_gene_length_bp, na.rm = TRUE),
    median_log10_gene_length_bp = median(log10_gene_length_bp, na.rm = TRUE),
    sd_log10_gene_length_bp = sd(log10_gene_length_bp, na.rm = TRUE),
    mean_sd_log10_nTPM_plus1 = mean(sd_log10_nTPM_plus1, na.rm = TRUE),
    median_sd_log10_nTPM_plus1 = median(sd_log10_nTPM_plus1, na.rm = TRUE),
    mean_nTPM_across_expressed_tissues = mean(mean_nTPM_across_expressed_tissues, na.rm = TRUE),
    median_nTPM_across_expressed_tissues = median(median_nTPM_across_expressed_tissues, na.rm = TRUE),
    .groups = "drop"
  )

write.table(
  fixed_summary,
  file = file.path(table_output_dir, "hpa_gtex_gene_summary_with_fixed_nonhk_breadth_groups.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  fixed_group_summary,
  file = file.path(table_output_dir, "hpa_gtex_fixed_nonhk_breadth_group_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

p_fixed_length <- ggplot(
  fixed_summary,
  aes(
    x = fixed_breadth_group,
    y = log10_gene_length_bp,
    fill = fixed_breadth_group
  )
) +
  geom_boxplot(outlier.alpha = 0.2) +
  labs(
    title = "Gene length by fixed HK/non-HK breadth group",
    subtitle = "restricted non-HK <= 5 tissues; broad non-HK >= 30 tissues",
    x = "Gene group",
    y = "log10(gene length bp)",
    fill = "Gene group"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = file.path(figure_output_dir, "hpa_gtex_fixed_nonhk_breadth_groups_gene_length_log10_boxplot.png"),
  plot = p_fixed_length,
  width = 9,
  height = 6,
  dpi = 300
)

p_fixed_variability <- ggplot(
  fixed_summary,
  aes(
    x = fixed_breadth_group,
    y = sd_log10_nTPM_plus1,
    fill = fixed_breadth_group
  )
) +
  geom_boxplot(outlier.alpha = 0.2) +
  labs(
    title = "Expression variability by fixed HK/non-HK breadth group",
    subtitle = "Variability = SD of log10(nTPM + 1) across expressed tissues",
    x = "Gene group",
    y = "SD log10(nTPM + 1)",
    fill = "Gene group"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = file.path(figure_output_dir, "hpa_gtex_fixed_nonhk_breadth_groups_expression_variability_boxplot.png"),
  plot = p_fixed_variability,
  width = 9,
  height = 6,
  dpi = 300
)

p_fixed_density <- ggplot(
  fixed_summary,
  aes(
    x = log10_gene_length_bp,
    color = fixed_breadth_group,
    linetype = fixed_breadth_group
  )
) +
  geom_density(linewidth = 1.1) +
  labs(
    title = "Gene length distributions by fixed HK/non-HK breadth group",
    subtitle = "Density of log10(gene length bp)",
    x = "log10(gene length bp)",
    y = "Density",
    color = "Gene group",
    linetype = "Gene group"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(
  filename = file.path(figure_output_dir, "hpa_gtex_fixed_nonhk_breadth_groups_gene_length_density.png"),
  plot = p_fixed_density,
  width = 9,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# Console output
# ------------------------------------------------------------

cat("\nDone.\n")

cat("\nEffect-size summary: HK vs broad non-HK\n")
print(effect_size_table)

cat("\nHK detection per tissue preview:\n")
print(head(hk_detection_per_tissue, 10))

cat("\nFixed-cutoff group summary:\n")
print(fixed_group_summary)

cat("\nTables written to:\n")
cat(table_output_dir, "\n")

cat("\nFigures written to:\n")
cat(figure_output_dir, "\n")