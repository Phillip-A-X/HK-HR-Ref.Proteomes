# scripts/R/humanHK_eisenberg_genelengthcompare.R


# Human Eisenberg HK vs non-HK gene length analysis

# Input:
# C:\Users\phill\OneDrive\Desktop\BSc\R\human_all_genes_with_eisenberg_hk_status.tsv

# Expected columns:
# Gene stable ID
# Gene start (bp)
# Gene end (bp)
# Chromosome/scaffold name
# gene_length_bp
# hk_status

# Output:
# eisenberg/results/tables/human_eisenberg_gene_length_summary.tsv
# eisenberg/results/tables/human_eisenberg_gene_length_wilcox_test.tsv
# eisenberg/results/figures/human_eisenberg_gene_length_boxplot.png
# eisenberg/results/figures/human_eisenberg_mean_gene_length_barplot.png


# Setup

input_file <- file.choose() #"data/processed/human_all_genes_with_eisenberg_hk_status.tsv"

table_output_dir <- "eisenberg/results/tables"
figure_output_dir <- "eisenberg/results/figures"

dir.create(table_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_output_dir, recursive = TRUE, showWarnings = FALSE)


# Libraries

library(ggplot2)
library(scales)


# Load data

df <- read.delim(input_file, stringsAsFactors = FALSE, check.names = FALSE)


# Basic checks

required_columns <- c("Gene stable ID", "gene_length_bp", "hk_status")
missing_columns <- setdiff(required_columns, colnames(df))

if (length(missing_columns) > 0) {
  stop(paste("Missing required columns:", paste(missing_columns, collapse = ", ")))
}

df$hk_status <- factor(df$hk_status, levels = c("HK", "non-HK"))


# Descriptive statistics

human_eisenberg_genelengthcompstat_table <- aggregate(
  gene_length_bp ~ hk_status,
  data = df,
  FUN = function(x) {
    c(
      n = length(x),
      mean = mean(x),
      median = median(x),
      sd = sd(x),
      min = min(x),
      q1 = quantile(x, 0.25),
      q3 = quantile(x, 0.75),
      iqr = IQR(x),
      max = max(x)
    )
  }
)

# Convert aggregate matrix-column into normal data frame
human_eisenberg_genelengthcompstat_table <- do.call(data.frame, human_eisenberg_genelengthcompstat_table)

# Clean column names
colnames(human_eisenberg_genelengthcompstat_table) <- c(
  "hk_status",
  "n",
  "mean_gene_length_bp",
  "median_gene_length_bp",
  "sd_gene_length_bp",
  "min_gene_length_bp",
  "q1_gene_length_bp",
  "q3_gene_length_bp",
  "iqr_gene_length_bp",
  "max_gene_length_bp"
)

# Write summary table
write.table(
  human_eisenberg_genelengthcompstat_table,
  file = file.path(table_output_dir, "human_eisenberg_gene_length_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

print(human_eisenberg_genelengthcompstat_table)


# Boxplot

# Gene lengths are often right-skewed.
# Therefore, the y-axis is log10-scaled for readability.

png(
  filename = file.path(figure_output_dir, "human_eisenberg_gene_length_boxplot.png"),
  width = 1200,
  height = 1000,
  res = 200
)

print(
  ggplot(df, aes(x = hk_status, y = gene_length_bp)) +
    geom_boxplot(outlier.alpha = 0.25) +
    scale_y_log10(labels = label_comma()) +
    labs(
      title = "Human Eisenberg HK vs non-HK gene lengths",
      x = "Gene status",
      y = "Gene length (bp, log10 scale)"
    ) +
    theme_classic()
)

dev.off()


# Mean gene length barplot

png(
  filename = file.path(figure_output_dir, "human_eisenberg_mean_gene_length_barplot.png"),
  width = 1200,
  height = 1000,
  res = 200
)

print(
  ggplot(
    human_eisenberg_genelengthcompstat_table,
    aes(x = hk_status, y = mean_gene_length_bp)
  ) +
    geom_col(width = 0.6) +
    geom_text(
      aes(label = comma(round(mean_gene_length_bp, 0))),
      vjust = -0.5,
      size = 4
    ) +
    scale_y_continuous(labels = label_comma()) +
    labs(
      title = "Human Eisenberg HK vs non-HK mean gene length",
      x = "Gene status",
      y = "Mean gene length (bp)"
    ) +
    theme_classic()
)

dev.off()


# Statistical test for significance

hk_lengths <- df$gene_length_bp[df$hk_status == "HK"]
nonhk_lengths <- df$gene_length_bp[df$hk_status == "non-HK"]

wilcox_result <- wilcox.test(hk_lengths, nonhk_lengths)

print(wilcox_result)

wilcox_output <- data.frame(
  test = "Wilcoxon rank-sum test",
  statistic = unname(wilcox_result$statistic),
  p_value = wilcox_result$p.value
)

write.table(
  wilcox_output,
  file = file.path(table_output_dir, "human_eisenberg_gene_length_wilcox_test.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Welch t-test for comparing mean gene lengths

t_test_result <- t.test(hk_lengths, nonhk_lengths)

print(t_test_result)

t_test_output <- data.frame(
  test = "Welch t-test",
  comparison = "HK vs non-HK gene_length_bp",
  mean_HK_bp = mean(hk_lengths),
  mean_nonHK_bp = mean(nonhk_lengths),
  mean_difference_bp = mean(hk_lengths) - mean(nonhk_lengths),
  statistic = unname(t_test_result$statistic),
  df = unname(t_test_result$parameter),
  p_value = t_test_result$p.value,
  conf_int_low = t_test_result$conf.int[1],
  conf_int_high = t_test_result$conf.int[2]
)

write.table(
  t_test_output,
  file = file.path(table_output_dir, "human_eisenberg_gene_length_welch_t_test.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
