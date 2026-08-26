# HPA-GTEx non-HK expression breadth quartile analysis
#
# Purpose:
#   Split non-HK genes into restricted / intermediate / broad groups
#   based on quartiles of expression breadth.
#
# Input:
#   data/processed/hpa_gtex/HumanProteinCodingGenes_bytissue_onchromosomes.tsv
#
# Definition:
#   expressed in tissue = nTPM >= 1
#   expression breadth = number of GTEx tissues with nTPM >= 1
#
# Groups:
#   HK:
#       HRT Atlas HK genes
#   restricted non-HK:
#       non-HK genes with expression breadth <= Q1 of non-HK breadth distribution
#   intermediate non-HK:
#       non-HK genes between Q1 and Q3
#   broad non-HK:
#       non-HK genes with expression breadth >= Q3 of non-HK breadth distribution
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

df <- df %>%
  filter(
    !is.na(gene_id),
    !is.na(Tissue),
    !is.na(nTPM),
    !is.na(gene_length_bp),
    !is.na(log10_gene_length_bp),
    !is.na(hk_status)
  )

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
# Calculate non-HK expression-breadth quartiles
# ------------------------------------------------------------

nonhk_breadth <- gene_summary %>%
  filter(hk_status == "non-HK") %>%
  pull(n_tissues_expressed)

q1_nonhk <- as.numeric(quantile(nonhk_breadth, 0.25, na.rm = TRUE, type = 7))
q3_nonhk <- as.numeric(quantile(nonhk_breadth, 0.75, na.rm = TRUE, type = 7))

quartile_info <- data.frame(
  metric = c(
    "nonHK_expression_breadth_Q1",
    "nonHK_expression_breadth_Q3",
    "number_of_GTEx_tissues"
  ),
  value = c(
    q1_nonhk,
    q3_nonhk,
    length(unique(df$Tissue))
  )
)

write.table(
  quartile_info,
  file = file.path(table_output_dir, "nonhk_expression_breadth_quartile_cutoffs.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------
# Assign analysis group
# ------------------------------------------------------------

gene_summary <- gene_summary %>%
  mutate(
    breadth_group = case_when(
      hk_status == "HK" ~ "HK",
      hk_status == "non-HK" & n_tissues_expressed <= q1_nonhk ~ "restricted non-HK",
      hk_status == "non-HK" & n_tissues_expressed >= q3_nonhk ~ "broad non-HK",
      hk_status == "non-HK" ~ "intermediate non-HK",
      TRUE ~ "unknown"
    )
  )

gene_summary$breadth_group <- factor(
  gene_summary$breadth_group,
  levels = c(
    "restricted non-HK",
    "intermediate non-HK",
    "broad non-HK",
    "HK"
  )
)

write.table(
  gene_summary,
  file = file.path(table_output_dir, "hpa_gtex_gene_summary_with_nonhk_breadth_groups.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------
# Group summary
# ------------------------------------------------------------

group_summary <- gene_summary %>%
  group_by(breadth_group) %>%
  summarise(
    n_genes = n(),
    mean_expression_breadth = mean(n_tissues_expressed, na.rm = TRUE),
    median_expression_breadth = median(n_tissues_expressed, na.rm = TRUE),
    min_expression_breadth = min(n_tissues_expressed, na.rm = TRUE),
    max_expression_breadth = max(n_tissues_expressed, na.rm = TRUE),
    mean_gene_length_bp = mean(gene_length_bp, na.rm = TRUE),
    median_gene_length_bp = median(gene_length_bp, na.rm = TRUE),
    variance_gene_length_bp = var(gene_length_bp, na.rm = TRUE),
    sd_gene_length_bp = sd(gene_length_bp, na.rm = TRUE),
    mean_log10_gene_length_bp = mean(log10_gene_length_bp, na.rm = TRUE),
    median_log10_gene_length_bp = median(log10_gene_length_bp, na.rm = TRUE),
    variance_log10_gene_length_bp = var(log10_gene_length_bp, na.rm = TRUE),
    sd_log10_gene_length_bp = sd(log10_gene_length_bp, na.rm = TRUE),
    mean_sd_log10_nTPM_plus1 = mean(sd_log10_nTPM_plus1, na.rm = TRUE),
    median_sd_log10_nTPM_plus1 = median(sd_log10_nTPM_plus1, na.rm = TRUE),
    mean_variance_log10_nTPM_plus1 = mean(variance_log10_nTPM_plus1, na.rm = TRUE),
    median_variance_log10_nTPM_plus1 = median(variance_log10_nTPM_plus1, na.rm = TRUE),
    mean_nTPM_across_expressed_tissues = mean(mean_nTPM_across_expressed_tissues, na.rm = TRUE),
    median_nTPM_across_expressed_tissues = median(median_nTPM_across_expressed_tissues, na.rm = TRUE),
    .groups = "drop"
  )

write.table(
  group_summary,
  file = file.path(table_output_dir, "hpa_gtex_nonhk_breadth_group_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------
# Statistical tests
# ------------------------------------------------------------

test_results <- data.frame(
  comparison = character(),
  variable = character(),
  test = character(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

add_test <- function(comparison, variable, test, p_value) {
  data.frame(
    comparison = comparison,
    variable = variable,
    test = test,
    p_value = p_value,
    stringsAsFactors = FALSE
  )
}

run_pairwise_tests <- function(data, group_a, group_b, variable_name) {
  sub <- data %>%
    filter(breadth_group %in% c(group_a, group_b)) %>%
    droplevels()
  
  sub$breadth_group <- factor(sub$breadth_group, levels = c(group_a, group_b))
  
  w <- wilcox.test(as.formula(paste(variable_name, "~ breadth_group")), data = sub)
  t <- t.test(as.formula(paste(variable_name, "~ breadth_group")), data = sub)
  
  bind_rows(
    add_test(paste(group_a, "vs", group_b), variable_name, "Wilcoxon rank-sum", w$p.value),
    add_test(paste(group_a, "vs", group_b), variable_name, "Welch t-test", t$p.value)
  )
}

pairwise_comparisons <- list(
  c("HK", "broad non-HK"),
  c("HK", "restricted non-HK"),
  c("broad non-HK", "restricted non-HK"),
  c("broad non-HK", "intermediate non-HK"),
  c("restricted non-HK", "intermediate non-HK")
)

variables_to_test <- c(
  "log10_gene_length_bp",
  "gene_length_bp",
  "n_tissues_expressed",
  "sd_log10_nTPM_plus1",
  "mean_nTPM_across_expressed_tissues"
)

for (comp in pairwise_comparisons) {
  for (var in variables_to_test) {
    test_results <- bind_rows(
      test_results,
      run_pairwise_tests(gene_summary, comp[1], comp[2], var)
    )
  }
}

write.table(
  test_results,
  file = file.path(table_output_dir, "hpa_gtex_nonhk_breadth_group_pairwise_tests.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------
# Plot 1: expression breadth by group
# ------------------------------------------------------------

p_breadth <- ggplot(
  gene_summary,
  aes(x = breadth_group, y = n_tissues_expressed, fill = breadth_group)
) +
  geom_boxplot(outlier.alpha = 0.2) +
  labs(
    title = "Expression breadth by HK/non-HK breadth group",
    subtitle = "Expression breadth = number of GTEx tissues with nTPM >= 1",
    x = "Gene group",
    y = "Number of tissues expressed",
    fill = "Gene group"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = file.path(figure_output_dir, "hpa_gtex_nonhk_breadth_groups_expression_breadth_boxplot.png"),
  plot = p_breadth,
  width = 9,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# Plot 2: log10 gene length by group
# ------------------------------------------------------------

p_length <- ggplot(
  gene_summary,
  aes(x = breadth_group, y = log10_gene_length_bp, fill = breadth_group)
) +
  geom_boxplot(outlier.alpha = 0.2) +
  labs(
    title = "Gene length by HK/non-HK breadth group",
    subtitle = "non-HK groups are based on expression-breadth quartiles",
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
  filename = file.path(figure_output_dir, "hpa_gtex_nonhk_breadth_groups_gene_length_log10_boxplot.png"),
  plot = p_length,
  width = 9,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# Plot 3: expression variability by group
# ------------------------------------------------------------

p_var <- ggplot(
  gene_summary,
  aes(x = breadth_group, y = sd_log10_nTPM_plus1, fill = breadth_group)
) +
  geom_boxplot(outlier.alpha = 0.2) +
  labs(
    title = "Expression variability by HK/non-HK breadth group",
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
  filename = file.path(figure_output_dir, "hpa_gtex_nonhk_breadth_groups_expression_variability_boxplot.png"),
  plot = p_var,
  width = 9,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# Plot 4: density of gene length by group
# ------------------------------------------------------------

p_density <- ggplot(
  gene_summary,
  aes(
    x = log10_gene_length_bp,
    color = breadth_group,
    linetype = breadth_group
  )
) +
  geom_density(linewidth = 1.1) +
  labs(
    title = "Gene length distributions by HK/non-HK breadth group",
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
  filename = file.path(figure_output_dir, "hpa_gtex_nonhk_breadth_groups_gene_length_density.png"),
  plot = p_density,
  width = 9,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# Plot 5: expression breadth vs gene length by group
# ------------------------------------------------------------

p_scatter <- ggplot(
  gene_summary,
  aes(
    x = n_tissues_expressed,
    y = log10_gene_length_bp,
    color = breadth_group
  )
) +
  geom_point(alpha = 0.35, size = 0.8) +
  geom_smooth(method = "loess", se = TRUE) +
  labs(
    title = "Expression breadth and gene length by HK/non-HK breadth group",
    subtitle = "non-HK groups are based on expression-breadth quartiles",
    x = "Number of tissues expressed",
    y = "log10(gene length bp)",
    color = "Gene group"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(
  filename = file.path(figure_output_dir, "hpa_gtex_nonhk_breadth_groups_breadth_vs_gene_length.png"),
  plot = p_scatter,
  width = 9,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# Plot 6: mean nTPM by group
# ------------------------------------------------------------

p_mean_ntpm <- ggplot(
  gene_summary,
  aes(x = breadth_group, y = mean_nTPM_across_expressed_tissues, fill = breadth_group)
) +
  geom_boxplot(outlier.alpha = 0.2) +
  scale_y_log10(labels = comma) +
  labs(
    title = "Mean nTPM across expressed tissues by gene group",
    subtitle = "Y-axis log10-scaled",
    x = "Gene group",
    y = "Mean nTPM across expressed tissues",
    fill = "Gene group"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = file.path(figure_output_dir, "hpa_gtex_nonhk_breadth_groups_mean_ntpm_boxplot.png"),
  plot = p_mean_ntpm,
  width = 9,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# Console output
# ------------------------------------------------------------

cat("\nDone.\n")

cat("\nnon-HK expression-breadth quartile cutoffs:\n")
print(quartile_info)

cat("\nGroup counts:\n")
print(table(gene_summary$breadth_group))

cat("\nGroup summary:\n")
print(group_summary)

cat("\nPairwise test preview:\n")
print(head(test_results, 20))

cat("\nTables written to:\n")
cat(table_output_dir, "\n")

cat("\nFigures written to:\n")
cat(figure_output_dir, "\n")