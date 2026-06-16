# E. coli HK/essential-proxy vs non-HK gene length analysis
# Source: Gerdes et al. 2003 essential-gene proxy
# Input: protein-coding E. coli genes with assigned HK/non-HK status

library(ggplot2)
library(dplyr)
library(scales)

# Optional: set working directory manually if needed
# setwd("~/projects/HK-HR-Ref.Proteomes")

input_file <- "data/processed/gene_lengths_with_hk_status/escherichia_coli_k12_mg1655_gene_lengths_with_hk_status.tsv"

table_output_dir <- "results/ecoli/tables"
figure_output_dir <- "results/ecoli/figures"

dir.create(table_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_output_dir, recursive = TRUE, showWarnings = FALSE)

df <- read.delim(input_file, stringsAsFactors = FALSE, check.names = FALSE)

required_columns <- c("gene_id", "gene_name", "gene_length_bp", "hk_status")

missing_columns <- setdiff(required_columns, colnames(df))

if (length(missing_columns) > 0) {
  stop(paste("Missing required columns:", paste(missing_columns, collapse = ", ")))
}

df <- df %>%
  filter(!is.na(gene_length_bp)) %>%
  filter(gene_length_bp > 0) %>%
  filter(hk_status %in% c("HK", "non-HK"))

df$hk_status <- factor(df$hk_status, levels = c("HK", "non-HK"))

hk_lengths <- df$gene_length_bp[df$hk_status == "HK"]
nonhk_lengths <- df$gene_length_bp[df$hk_status == "non-HK"]

hk_log_lengths <- log10(hk_lengths)
nonhk_log_lengths <- log10(nonhk_lengths)

# Summary statistics
summary_table <- df %>%
  group_by(hk_status) %>%
  summarise(
    n = n(),
    mean_gene_length_bp = mean(gene_length_bp),
    median_gene_length_bp = median(gene_length_bp),
    sd_gene_length_bp = sd(gene_length_bp),
    q1_gene_length_bp = quantile(gene_length_bp, 0.25),
    q3_gene_length_bp = quantile(gene_length_bp, 0.75),
    iqr_gene_length_bp = IQR(gene_length_bp),
    min_gene_length_bp = min(gene_length_bp),
    max_gene_length_bp = max(gene_length_bp),
    mean_log10_gene_length = mean(log10(gene_length_bp)),
    median_log10_gene_length = median(log10(gene_length_bp)),
    .groups = "drop"
  )

print(summary_table)

write.table(
  summary_table,
  file = file.path(table_output_dir, "ecoli_gerdes_gene_length_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Wilcoxon rank-sum test on raw gene lengths
wilcox_result <- wilcox.test(
  hk_lengths,
  nonhk_lengths,
  exact = FALSE
)

wilcox_output <- data.frame(
  test = "Wilcoxon rank-sum test",
  comparison = "HK/essential-proxy vs non-HK gene_length_bp",
  statistic = unname(wilcox_result$statistic),
  p_value = wilcox_result$p.value
)

print(wilcox_output)

write.table(
  wilcox_output,
  file = file.path(table_output_dir, "ecoli_gerdes_gene_length_wilcox_test.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Welch t-test on raw gene lengths
welch_raw_result <- t.test(
  hk_lengths,
  nonhk_lengths
)

welch_raw_output <- data.frame(
  test = "Welch t-test",
  comparison = "HK/essential-proxy vs non-HK raw gene_length_bp",
  mean_HK_bp = mean(hk_lengths),
  mean_nonHK_bp = mean(nonhk_lengths),
  mean_difference_bp = mean(hk_lengths) - mean(nonhk_lengths),
  statistic = unname(welch_raw_result$statistic),
  df = unname(welch_raw_result$parameter),
  p_value = welch_raw_result$p.value,
  conf_int_low = welch_raw_result$conf.int[1],
  conf_int_high = welch_raw_result$conf.int[2]
)

print(welch_raw_output)

write.table(
  welch_raw_output,
  file = file.path(table_output_dir, "ecoli_gerdes_gene_length_welch_t_test_raw.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Welch t-test on log10-transformed gene lengths
welch_log_result <- t.test(
  hk_log_lengths,
  nonhk_log_lengths
)

welch_log_output <- data.frame(
  test = "Welch t-test",
  comparison = "HK/essential-proxy vs non-HK log10(gene_length_bp)",
  mean_log10_HK = mean(hk_log_lengths),
  mean_log10_nonHK = mean(nonhk_log_lengths),
  mean_log10_difference = mean(hk_log_lengths) - mean(nonhk_log_lengths),
  statistic = unname(welch_log_result$statistic),
  df = unname(welch_log_result$parameter),
  p_value = welch_log_result$p.value,
  conf_int_low = welch_log_result$conf.int[1],
  conf_int_high = welch_log_result$conf.int[2]
)

print(welch_log_output)

write.table(
  welch_log_output,
  file = file.path(table_output_dir, "ecoli_gerdes_gene_length_welch_t_test_log10.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Boxplot with log10-scaled y-axis
png(
  filename = file.path(figure_output_dir, "ecoli_gerdes_gene_length_boxplot_log10_y.png"),
  width = 1200,
  height = 1000,
  res = 200
)

print(
  ggplot(df, aes(x = hk_status, y = gene_length_bp)) +
    geom_boxplot(outlier.alpha = 0.25) +
    scale_y_log10(labels = label_comma()) +
    labs(
      title = "E. coli Gerdes essential-proxy vs non-HK gene lengths",
      x = "Gene status",
      y = "Gene length (bp, log10 scale)"
    ) +
    theme_classic()
)

dev.off()

# Mean gene length barplot
mean_barplot_table <- summary_table %>%
  select(hk_status, mean_gene_length_bp)

png(
  filename = file.path(figure_output_dir, "ecoli_gerdes_mean_gene_length_barplot.png"),
  width = 1200,
  height = 1000,
  res = 200
)

print(
  ggplot(mean_barplot_table, aes(x = hk_status, y = mean_gene_length_bp)) +
    geom_col(width = 0.65) +
    geom_text(
      aes(label = round(mean_gene_length_bp, 0)),
      vjust = -0.4,
      size = 4
    ) +
    scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, 0.12))) +
    labs(
      title = "Mean E. coli gene length",
      x = "Gene status",
      y = "Mean gene length (bp)"
    ) +
    theme_classic()
)

dev.off()

# Density plot on log10 gene length
png(
  filename = file.path(figure_output_dir, "ecoli_gerdes_gene_length_distribution_density_log10.png"),
  width = 1400,
  height = 1000,
  res = 200
)

print(
  ggplot(df, aes(x = log10(gene_length_bp), fill = hk_status)) +
    geom_density(alpha = 0.45) +
    labs(
      title = "E. coli gene length distribution",
      x = "log10 gene length (bp)",
      y = "Density",
      fill = "Gene status"
    ) +
    theme_classic()
)

dev.off()

# Frequency polygon on log10 gene length
png(
  filename = file.path(figure_output_dir, "ecoli_gerdes_gene_length_distribution_freqpoly_log10.png"),
  width = 1400,
  height = 1000,
  res = 200
)

print(
  ggplot(df, aes(x = log10(gene_length_bp), color = hk_status)) +
    geom_freqpoly(bins = 50, linewidth = 1.1) +
    labs(
      title = "E. coli gene length distribution",
      x = "log10 gene length (bp)",
      y = "Number of genes",
      color = "Gene status"
    ) +
    theme_classic()
)

dev.off()

# Faceted histogram on log10 gene length
png(
  filename = file.path(figure_output_dir, "ecoli_gerdes_gene_length_distribution_histogram_log10_faceted.png"),
  width = 1400,
  height = 1100,
  res = 200
)

print(
  ggplot(df, aes(x = log10(gene_length_bp))) +
    geom_histogram(bins = 50) +
    facet_wrap(~ hk_status, ncol = 1, scales = "free_y") +
    labs(
      title = "E. coli gene length distribution by gene status",
      x = "log10 gene length (bp)",
      y = "Number of genes"
    ) +
    theme_classic()
)

dev.off()

cat("\nAnalysis complete.\n")
cat("Tables written to:", table_output_dir, "\n")
cat("Figures written to:", figure_output_dir, "\n")
