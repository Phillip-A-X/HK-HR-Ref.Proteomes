# scripts/R/compare_humanHK_vsnonHK_hrt_protein_lengths.R


# Human HRT Atlas HK vs non-HK protein length analysis

# Input:
# data/processed/human_reference_proteome_with_hrt_atlas_hk_status.tsv

# Expected columns:
# uniprot_accession
# entry_name
# description
# protein_length_aa
# hk_status
# hk_source

# Output:
# results/tables/human_hrt_protein_length_summary.tsv
# results/tables/human_hrt_protein_length_wilcox_test.tsv
# results/tables/human_hrt_protein_length_welch_t_test.tsv
# results/figures/human_hrt_protein_length_boxplot.png
# results/figures/human_hrt_mean_protein_length_barplot.png


# Setup

input_file <- file.choose() #"data/processed/human_reference_proteome_with_hrt_atlas_hk_status.tsv"

table_output_dir <- "results/tables"
figure_output_dir <- "results/figures"

dir.create(table_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_output_dir, recursive = TRUE, showWarnings = FALSE)


# Libraries

library(ggplot2)
library(scales)


# Load data

df <- read.delim(input_file, stringsAsFactors = FALSE, check.names = FALSE)


# Basic checks

required_columns <- c("protein_length_aa", "hk_status")
missing_columns <- setdiff(required_columns, colnames(df))

if (length(missing_columns) > 0) {
  stop(paste("Missing required columns:", paste(missing_columns, collapse = ", ")))
}

df <- df[!is.na(df$protein_length_aa) & !is.na(df$hk_status), ]

df$hk_status <- factor(df$hk_status, levels = c("HK", "non-HK"))


# Descriptive statistics

human_hrt_proteinlengthcompstat_table <- aggregate(
  protein_length_aa ~ hk_status,
  data = df,
  FUN = function(x) {
    c(
      n = length(x),
      mean = mean(x),
      median = median(x),
      sd = sd(x),
      min = min(x),
      q1 = as.numeric(quantile(x, 0.25)),
      q3 = as.numeric(quantile(x, 0.75)),
      iqr = IQR(x),
      max = max(x)
    )
  }
)

# Convert aggregate matrix-column into normal data frame
human_hrt_proteinlengthcompstat_table <- do.call(data.frame, human_hrt_proteinlengthcompstat_table)

# Clean column names
colnames(human_hrt_proteinlengthcompstat_table) <- c(
  "hk_status",
  "n",
  "mean_protein_length_aa",
  "median_protein_length_aa",
  "sd_protein_length_aa",
  "min_protein_length_aa",
  "q1_protein_length_aa",
  "q3_protein_length_aa",
  "iqr_protein_length_aa",
  "max_protein_length_aa"
)

# Write summary table
write.table(
  human_hrt_proteinlengthcompstat_table,
  file = file.path(table_output_dir, "human_hrt_protein_length_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

print(human_hrt_proteinlengthcompstat_table)


# Boxplot

# Protein lengths are also right-skewed.
# Therefore, the y-axis is log10-scaled for readability.

png(
  filename = file.path(figure_output_dir, "human_hrt_protein_length_boxplot.png"),
  width = 1200,
  height = 1000,
  res = 200
)

print(
  ggplot(df, aes(x = hk_status, y = protein_length_aa)) +
    geom_boxplot(outlier.alpha = 0.25) +
    scale_y_log10(labels = label_comma()) +
    labs(
      title = "Human HRT Atlas HK vs non-HK protein lengths",
      x = "Protein status",
      y = "Protein length (aa, log10 scale)"
    ) +
    theme_classic()
)

dev.off()


# Mean protein length barplot

png(
  filename = file.path(figure_output_dir, "human_hrt_mean_protein_length_barplot.png"),
  width = 1200,
  height = 1000,
  res = 200
)

print(
  ggplot(
    human_hrt_proteinlengthcompstat_table,
    aes(x = hk_status, y = mean_protein_length_aa)
  ) +
    geom_col(width = 0.6) +
    geom_text(
      aes(label = comma(round(mean_protein_length_aa, 0))),
      vjust = -0.5,
      size = 4
    ) +
    scale_y_continuous(labels = label_comma()) +
    labs(
      title = "Human HRT Atlas HK vs non-HK mean protein length",
      x = "Protein status",
      y = "Mean protein length (aa)"
    ) +
    theme_classic()
)

dev.off()


# Wilcoxon rank-sum test

hk_lengths <- df$protein_length_aa[df$hk_status == "HK"]
nonhk_lengths <- df$protein_length_aa[df$hk_status == "non-HK"]

wilcox_result <- wilcox.test(hk_lengths, nonhk_lengths)

print(wilcox_result)

wilcox_output <- data.frame(
  test = "Wilcoxon rank-sum test",
  comparison = "HK vs non-HK protein_length_aa",
  statistic = unname(wilcox_result$statistic),
  p_value = wilcox_result$p.value
)

write.table(
  wilcox_output,
  file = file.path(table_output_dir, "human_hrt_protein_length_wilcox_test.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# Welch t-test for comparing mean protein lengths

t_test_result <- t.test(hk_lengths, nonhk_lengths)

print(t_test_result)

t_test_output <- data.frame(
  test = "Welch t-test",
  comparison = "HK vs non-HK protein_length_aa",
  mean_HK_aa = mean(hk_lengths),
  mean_nonHK_aa = mean(nonhk_lengths),
  mean_difference_aa = mean(hk_lengths) - mean(nonhk_lengths),
  statistic = unname(t_test_result$statistic),
  df = unname(t_test_result$parameter),
  p_value = t_test_result$p.value,
  conf_int_low = t_test_result$conf.int[1],
  conf_int_high = t_test_result$conf.int[2]
)

write.table(
  t_test_output,
  file = file.path(table_output_dir, "human_hrt_protein_length_welch_t_test.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
