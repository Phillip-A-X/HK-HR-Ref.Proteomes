# HPA GTEx HRT Atlas HK vs non-HK expression breadth and gene length analysis
#
# Input:
#   data/processed/hpa_gtex/HumanProteinCodingGenes_bytissue_onchromosomes.tsv
#
# This table contains HPA-GTEx gene-tissue pairs with nTPM >= 1,
# mapped to the Ensembl 115 human protein-coding main-chromosome background.
#
# Aim:
#   Compare HRT Atlas HK genes and non-HK genes regarding:
#       expression breadth
#       expression variability
#       gene length
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

df <- read.delim(input_file, stringsAsFactors = FALSE, check.names = FALSE)

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
# Gene-level summary:
# one row per gene
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

# If a gene is expressed in only one tissue, variance/sd are NA.
# Set these to 0 because there is no across-tissue variation among expressed tissues.
gene_summary$variance_log10_nTPM_plus1[is.na(gene_summary$variance_log10_nTPM_plus1)] <- 0
gene_summary$sd_log10_nTPM_plus1[is.na(gene_summary$sd_log10_nTPM_plus1)] <- 0

write.table(
  gene_summary,
  file = file.path(table_output_dir, "hpa_gtex_hrt_hk_gene_expression_breadth_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------
# HK vs non-HK summary statistics
# ------------------------------------------------------------

hk_summary <- gene_summary %>%
  group_by(hk_status) %>%
  summarise(
    n_genes = n(),
    mean_expression_breadth = mean(n_tissues_expressed, na.rm = TRUE),
    median_expression_breadth = median(n_tissues_expressed, na.rm = TRUE),
    variance_expression_breadth = var(n_tissues_expressed, na.rm = TRUE),
    sd_expression_breadth = sd(n_tissues_expressed, na.rm = TRUE),
    mean_gene_length_bp = mean(gene_length_bp, na.rm = TRUE),
    median_gene_length_bp = median(gene_length_bp, na.rm = TRUE),
    variance_gene_length_bp = var(gene_length_bp, na.rm = TRUE),
    sd_gene_length_bp = sd(gene_length_bp, na.rm = TRUE),
    mean_log10_gene_length_bp = mean(log10_gene_length_bp, na.rm = TRUE),
    median_log10_gene_length_bp = median(log10_gene_length_bp, na.rm = TRUE),
    mean_sd_log10_nTPM_plus1 = mean(sd_log10_nTPM_plus1, na.rm = TRUE),
    median_sd_log10_nTPM_plus1 = median(sd_log10_nTPM_plus1, na.rm = TRUE),
    mean_variance_log10_nTPM_plus1 = mean(variance_log10_nTPM_plus1, na.rm = TRUE),
    median_variance_log10_nTPM_plus1 = median(variance_log10_nTPM_plus1, na.rm = TRUE),
    .groups = "drop"
  )

write.table(
  hk_summary,
  file = file.path(table_output_dir, "hpa_gtex_hrt_hk_vs_nonhk_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------
# Statistical tests
# ------------------------------------------------------------

hk_test_results <- data.frame(
  comparison = character(),
  test = character(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

add_test <- function(comparison_name, test_name, p_value) {
  data.frame(
    comparison = comparison_name,
    test = test_name,
    p_value = p_value,
    stringsAsFactors = FALSE
  )
}

# Wilcoxon tests
test_breadth <- wilcox.test(n_tissues_expressed ~ hk_status, data = gene_summary)
test_length <- wilcox.test(gene_length_bp ~ hk_status, data = gene_summary)
test_log_length <- wilcox.test(log10_gene_length_bp ~ hk_status, data = gene_summary)
test_expr_var <- wilcox.test(sd_log10_nTPM_plus1 ~ hk_status, data = gene_summary)

# Welch t-tests
ttest_breadth <- t.test(n_tissues_expressed ~ hk_status, data = gene_summary)
ttest_length <- t.test(gene_length_bp ~ hk_status, data = gene_summary)
ttest_log_length <- t.test(log10_gene_length_bp ~ hk_status, data = gene_summary)
ttest_expr_var <- t.test(sd_log10_nTPM_plus1 ~ hk_status, data = gene_summary)

hk_test_results <- bind_rows(
  add_test("expression_breadth", "Wilcoxon rank-sum", test_breadth$p.value),
  add_test("gene_length_bp", "Wilcoxon rank-sum", test_length$p.value),
  add_test("log10_gene_length_bp", "Wilcoxon rank-sum", test_log_length$p.value),
  add_test("sd_log10_nTPM_plus1", "Wilcoxon rank-sum", test_expr_var$p.value),
  add_test("expression_breadth", "Welch t-test", ttest_breadth$p.value),
  add_test("gene_length_bp", "Welch t-test", ttest_length$p.value),
  add_test("log10_gene_length_bp", "Welch t-test", ttest_log_length$p.value),
  add_test("sd_log10_nTPM_plus1", "Welch t-test", ttest_expr_var$p.value)
)

write.table(
  hk_test_results,
  file = file.path(table_output_dir, "hpa_gtex_hrt_hk_vs_nonhk_tests.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------
# Plot 1: expression breadth HK vs non-HK
# ------------------------------------------------------------

p_breadth_box <- ggplot(gene_summary, aes(x = hk_status, y = n_tissues_expressed)) +
  geom_boxplot(outlier.alpha = 0.2) +
  labs(
    title = "Expression breadth of HRT Atlas HK vs non-HK genes",
    subtitle = "Expression breadth = number of GTEx tissues with nTPM >= 1",
    x = "HRT Atlas status",
    y = "Number of tissues expressed"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = file.path(figure_output_dir, "hpa_gtex_hrt_hk_vs_nonhk_expression_breadth_boxplot.png"),
  plot = p_breadth_box,
  width = 7,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# Plot 2: gene length HK vs non-HK
# ------------------------------------------------------------

p_length_box <- ggplot(gene_summary, aes(x = hk_status, y = log10_gene_length_bp)) +
  geom_boxplot(outlier.alpha = 0.2) +
  labs(
    title = "Gene length of HRT Atlas HK vs non-HK genes",
    subtitle = "Human protein-coding genes mapped to HPA-GTEx; y-axis log10-transformed",
    x = "HRT Atlas status",
    y = "log10(gene length bp)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = file.path(figure_output_dir, "hpa_gtex_hrt_hk_vs_nonhk_gene_length_log10_boxplot.png"),
  plot = p_length_box,
  width = 7,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# Plot 3: expression variability HK vs non-HK
# ------------------------------------------------------------

p_var_box <- ggplot(gene_summary, aes(x = hk_status, y = sd_log10_nTPM_plus1)) +
  geom_boxplot(outlier.alpha = 0.2) +
  labs(
    title = "Expression variability of HRT Atlas HK vs non-HK genes",
    subtitle = "Variability = SD of log10(nTPM + 1) across expressed tissues",
    x = "HRT Atlas status",
    y = "SD log10(nTPM + 1)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = file.path(figure_output_dir, "hpa_gtex_hrt_hk_vs_nonhk_expression_variability_boxplot.png"),
  plot = p_var_box,
  width = 7,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# Plot 4: expression breadth vs gene length, colored by HK status
# ------------------------------------------------------------

p_breadth_length <- ggplot(
  gene_summary,
  aes(x = n_tissues_expressed, y = log10_gene_length_bp, color = hk_status)
) +
  geom_point(alpha = 0.35, size = 0.8) +
  geom_smooth(method = "loess", se = TRUE) +
  labs(
    title = "Expression breadth and gene length by HRT Atlas status",
    subtitle = "Expression breadth = number of GTEx tissues with nTPM >= 1",
    x = "Number of tissues expressed",
    y = "log10(gene length bp)",
    color = "HRT Atlas status"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = file.path(figure_output_dir, "hpa_gtex_expression_breadth_vs_gene_length_by_hrt_status.png"),
  plot = p_breadth_length,
  width = 9,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# Plot 5: density of gene length by HK status
# ------------------------------------------------------------

p_density_length <- ggplot(
  gene_summary,
  aes(
    x = log10_gene_length_bp,
    color = hk_status,
    linetype = hk_status
  )
) +
  geom_density(linewidth = 1.1) +
  labs(
    title = "Gene length distributions of HRT Atlas HK vs non-HK genes",
    subtitle = "Density of log10(gene length bp)",
    x = "log10(gene length bp)",
    y = "Density",
    color = "HRT Atlas status",
    linetype = "HRT Atlas status"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(
  filename = file.path(figure_output_dir, "hpa_gtex_hrt_hk_vs_nonhk_gene_length_density.png"),
  plot = p_density_length,
  width = 8,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# Console output
# ------------------------------------------------------------

cat("\nDone.\n")

cat("\nTables written to:\n")
cat(table_output_dir, "\n")

cat("\nFigures written to:\n")
cat(figure_output_dir, "\n")

cat("\nHK status counts:\n")
print(table(gene_summary$hk_status))

cat("\nHK vs non-HK summary:\n")
print(hk_summary)

cat("\nTest results:\n")
print(hk_test_results)