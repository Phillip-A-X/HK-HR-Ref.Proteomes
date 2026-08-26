# HPA-GTEx protein length analysis by HK/non-HK expression breadth group
#
# Input:
#   data/processed/hpa_gtex_protein/human_reference_proteome_longest_protein_per_gene_with_hpa_gtex_breadth.tsv
#
# This table contains:
#   one representative protein per ENSG
#   representative = longest mapped UniProt reference-proteome protein sequence per ENSG
#   HPA-GTEx expression breadth groups from RNA data
#
# Output:
#   results/hpa_gtex_protein/tables/
#   results/hpa_gtex_protein/figures/

# ------------------------------------------------------------
# Setup
# ------------------------------------------------------------

library(ggplot2)
library(dplyr)

input_file <- "data/processed/hpa_gtex_protein/human_reference_proteome_longest_protein_per_gene_with_hpa_gtex_breadth.tsv"

table_output_dir <- "results/hpa_gtex_protein/tables"
figure_output_dir <- "results/hpa_gtex_protein/figures"

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
  "uniprot_accession",
  "protein_length_aa",
  "hk_status_gene_level",
  "n_tissues_expressed",
  "breadth_group"
)

missing_columns <- setdiff(required_columns, colnames(df))

if (length(missing_columns) > 0) {
  stop(paste("Missing required columns:", paste(missing_columns, collapse = ", ")))
}

df$protein_length_aa <- as.numeric(df$protein_length_aa)
df$n_tissues_expressed <- as.numeric(df$n_tissues_expressed)
df$log10_protein_length_aa <- log10(df$protein_length_aa)

df <- df %>%
  filter(
    !is.na(gene_id),
    !is.na(uniprot_accession),
    !is.na(protein_length_aa),
    !is.na(log10_protein_length_aa),
    !is.na(breadth_group)
  )

df$breadth_group <- factor(
  df$breadth_group,
  levels = c(
    "restricted non-HK",
    "intermediate non-HK",
    "broad non-HK",
    "HK"
  )
)

# ------------------------------------------------------------
# QC table
# ------------------------------------------------------------

qc <- data.frame(
  metric = c(
    "input_rows",
    "unique_genes",
    "unique_uniprot_accessions",
    "min_protein_length_aa",
    "max_protein_length_aa"
  ),
  value = c(
    nrow(df),
    length(unique(df$gene_id)),
    length(unique(df$uniprot_accession)),
    min(df$protein_length_aa, na.rm = TRUE),
    max(df$protein_length_aa, na.rm = TRUE)
  )
)

write.table(
  qc,
  file = file.path(table_output_dir, "hpa_gtex_protein_length_breadth_group_qc.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------
# Group summary
# ------------------------------------------------------------

group_summary <- df %>%
  group_by(breadth_group) %>%
  summarise(
    n_genes = n(),
    mean_expression_breadth = mean(n_tissues_expressed, na.rm = TRUE),
    median_expression_breadth = median(n_tissues_expressed, na.rm = TRUE),
    min_expression_breadth = min(n_tissues_expressed, na.rm = TRUE),
    max_expression_breadth = max(n_tissues_expressed, na.rm = TRUE),
    mean_protein_length_aa = mean(protein_length_aa, na.rm = TRUE),
    median_protein_length_aa = median(protein_length_aa, na.rm = TRUE),
    variance_protein_length_aa = var(protein_length_aa, na.rm = TRUE),
    sd_protein_length_aa = sd(protein_length_aa, na.rm = TRUE),
    q1_protein_length_aa = quantile(protein_length_aa, 0.25, na.rm = TRUE),
    q3_protein_length_aa = quantile(protein_length_aa, 0.75, na.rm = TRUE),
    iqr_protein_length_aa = IQR(protein_length_aa, na.rm = TRUE),
    min_protein_length_aa = min(protein_length_aa, na.rm = TRUE),
    max_protein_length_aa = max(protein_length_aa, na.rm = TRUE),
    mean_log10_protein_length_aa = mean(log10_protein_length_aa, na.rm = TRUE),
    median_log10_protein_length_aa = median(log10_protein_length_aa, na.rm = TRUE),
    sd_log10_protein_length_aa = sd(log10_protein_length_aa, na.rm = TRUE),
    .groups = "drop"
  )

write.table(
  group_summary,
  file = file.path(table_output_dir, "hpa_gtex_protein_length_breadth_group_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------
# Pairwise statistical tests
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
  "protein_length_aa",
  "log10_protein_length_aa",
  "n_tissues_expressed"
)

for (comp in pairwise_comparisons) {
  for (var in variables_to_test) {
    test_results <- bind_rows(
      test_results,
      run_pairwise_tests(df, comp[1], comp[2], var)
    )
  }
}

write.table(
  test_results,
  file = file.path(table_output_dir, "hpa_gtex_protein_length_breadth_group_pairwise_tests.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------
# Effect sizes: HK vs broad non-HK
# ------------------------------------------------------------

effect_summary <- df %>%
  filter(breadth_group %in% c("HK", "broad non-HK")) %>%
  group_by(breadth_group) %>%
  summarise(
    n_genes = n(),
    mean_protein_length_aa = mean(protein_length_aa, na.rm = TRUE),
    median_protein_length_aa = median(protein_length_aa, na.rm = TRUE),
    mean_log10_protein_length_aa = mean(log10_protein_length_aa, na.rm = TRUE),
    median_log10_protein_length_aa = median(log10_protein_length_aa, na.rm = TRUE),
    mean_expression_breadth = mean(n_tissues_expressed, na.rm = TRUE),
    median_expression_breadth = median(n_tissues_expressed, na.rm = TRUE),
    .groups = "drop"
  )

hk_row <- effect_summary %>% filter(breadth_group == "HK")
broad_row <- effect_summary %>% filter(breadth_group == "broad non-HK")

effect_size_table <- data.frame(
  comparison = "HK vs broad non-HK",
  variable = c(
    "mean_protein_length_aa",
    "median_protein_length_aa",
    "mean_log10_protein_length_aa",
    "median_log10_protein_length_aa",
    "mean_expression_breadth",
    "median_expression_breadth"
  ),
  HK = c(
    hk_row$mean_protein_length_aa,
    hk_row$median_protein_length_aa,
    hk_row$mean_log10_protein_length_aa,
    hk_row$median_log10_protein_length_aa,
    hk_row$mean_expression_breadth,
    hk_row$median_expression_breadth
  ),
  broad_nonHK = c(
    broad_row$mean_protein_length_aa,
    broad_row$median_protein_length_aa,
    broad_row$mean_log10_protein_length_aa,
    broad_row$median_log10_protein_length_aa,
    broad_row$mean_expression_breadth,
    broad_row$median_expression_breadth
  )
)

effect_size_table$difference_HK_minus_broad_nonHK <- effect_size_table$HK - effect_size_table$broad_nonHK
effect_size_table$ratio_HK_div_broad_nonHK <- effect_size_table$HK / effect_size_table$broad_nonHK

write.table(
  effect_size_table,
  file = file.path(table_output_dir, "hpa_gtex_protein_length_hk_vs_broad_nonhk_effect_sizes.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------
# Plot 1: protein length boxplot
# ------------------------------------------------------------

p_box <- ggplot(
  df,
  aes(x = breadth_group, y = log10_protein_length_aa, fill = breadth_group)
) +
  geom_boxplot(outlier.alpha = 0.2) +
  labs(
    title = "Protein length by HK/non-HK expression breadth group",
    subtitle = "One longest mapped UniProt reference-proteome protein retained per ENSG",
    x = "Gene group",
    y = "log10(protein length aa)",
    fill = "Gene group"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = file.path(figure_output_dir, "hpa_gtex_protein_length_breadth_group_log10_boxplot.png"),
  plot = p_box,
  width = 9,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# Plot 2: protein length density
# ------------------------------------------------------------

p_density <- ggplot(
  df,
  aes(
    x = log10_protein_length_aa,
    color = breadth_group,
    linetype = breadth_group
  )
) +
  geom_density(linewidth = 1.1) +
  labs(
    title = "Protein length distributions by HK/non-HK expression breadth group",
    subtitle = "Density of log10(protein length aa)",
    x = "log10(protein length aa)",
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
  filename = file.path(figure_output_dir, "hpa_gtex_protein_length_breadth_group_density.png"),
  plot = p_density,
  width = 9,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# Plot 3: expression breadth vs protein length
# ------------------------------------------------------------

p_scatter <- ggplot(
  df,
  aes(
    x = n_tissues_expressed,
    y = log10_protein_length_aa,
    color = breadth_group
  )
) +
  geom_point(alpha = 0.35, size = 0.8) +
  #geom_smooth(method = "loess", se = TRUE) +
  labs(
    title = "Expression breadth and protein length",
    subtitle = "Protein length from longest mapped UniProt reference-proteome sequence per ENSG",
    x = "Number of tissues expressed",
    y = "log10(protein length aa)",
    color = "Gene group"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(
  filename = file.path(figure_output_dir, "hpa_gtex_expression_breadth_vs_protein_length.png"),
  plot = p_scatter,
  width = 9,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# Console output
# ------------------------------------------------------------

cat("\nDone.\n")

cat("\nQC:\n")
print(qc)

cat("\nGroup summary:\n")
print(group_summary)

cat("\nEffect sizes HK vs broad non-HK:\n")
print(effect_size_table)

cat("\nPairwise tests preview:\n")
print(head(test_results, 20))

cat("\nTables written to:\n")
cat(table_output_dir, "\n")

cat("\nFigures written to:\n")
cat(figure_output_dir, "\n")