# ============================================================
# HPA-GTEx log10 gene-length histograms by tissue
# with theoretical normal distribution
#
# Input:
#   data/processed/hpa_gtex/
#   HumanProteinCodingGenes_bytissue_onchromosomes.tsv
#
# expressed = nTPM >= 1
#
# Output:
#   results/hpa_gtex/figures/
#   hpa_gtex_log10_gene_length_histogram_by_tissue_with_theoretical_normal.png
# ============================================================


# ------------------------------------------------------------
# Setup
# ------------------------------------------------------------

library(ggplot2)
library(dplyr)

input_file <-
  "data/processed/hpa_gtex/HumanProteinCodingGenes_bytissue_onchromosomes.tsv"

figure_output_dir <-
  "results/hpa_gtex/figures"

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

df$nTPM <- as.numeric(df$nTPM)

df$gene_length_bp <-
  as.numeric(df$gene_length_bp)


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
  ) %>%
  
  mutate(
    log10_gene_length_bp =
      log10(gene_length_bp)
  )


# ------------------------------------------------------------
# Histogram settings
#
# Smaller bin width gives more resolution and makes
# small secondary peaks easier to detect.
# ------------------------------------------------------------

binwidth <- 0.05


# ------------------------------------------------------------
# Order tissues by mean gene length
# ------------------------------------------------------------

tissue_order <- df %>%
  group_by(Tissue) %>%
  summarise(
    mean_gene_length_bp =
      mean(gene_length_bp),
    .groups = "drop"
  ) %>%
  arrange(mean_gene_length_bp) %>%
  pull(Tissue)

df$Tissue <- factor(
  df$Tissue,
  levels = tissue_order
)


# ------------------------------------------------------------
# Calculate tissue-specific parameters
# for theoretical normal distribution
# ------------------------------------------------------------

normal_parameters <- df %>%
  group_by(Tissue) %>%
  summarise(
    
    n_genes = n(),
    
    mean_log10 =
      mean(log10_gene_length_bp),
    
    sd_log10 =
      sd(log10_gene_length_bp),
    
    min_log10 =
      min(log10_gene_length_bp),
    
    max_log10 =
      max(log10_gene_length_bp),
    
    .groups = "drop"
  )


# ------------------------------------------------------------
# Generate theoretical normal curves
#
# Histogram uses COUNTS rather than density.
#
# Therefore theoretical density has to be converted to:
#
# expected count =
# density * number of genes * bin width
# ------------------------------------------------------------

normal_curve <- bind_rows(
  
  lapply(
    seq_len(nrow(normal_parameters)),
    function(i) {
      
      tissue_name <-
        normal_parameters$Tissue[i]
      
      n <-
        normal_parameters$n_genes[i]
      
      mu <-
        normal_parameters$mean_log10[i]
      
      sigma <-
        normal_parameters$sd_log10[i]
      
      xmin <-
        normal_parameters$min_log10[i]
      
      xmax <-
        normal_parameters$max_log10[i]
      
      
      x_values <- seq(
        xmin,
        xmax,
        length.out = 500
      )
      
      
      data.frame(
        
        Tissue = tissue_name,
        
        log10_gene_length_bp =
          x_values,
        
        theoretical_count =
          dnorm(
            x_values,
            mean = mu,
            sd = sigma
          ) *
          n *
          binwidth
      )
    }
  )
)


normal_curve$Tissue <- factor(
  normal_curve$Tissue,
  levels = tissue_order
)


# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

p_hist_theoretical <- ggplot(
  df,
  aes(
    x = log10_gene_length_bp
  )
) +
  
  geom_histogram(
    binwidth = binwidth,
    boundary = 0,
    fill = "grey55",
    color = "grey20",
    linewidth = 0.2
  ) +
  
  geom_line(
    data = normal_curve,
    aes(
      x = log10_gene_length_bp,
      y = theoretical_count
    ),
    inherit.aes = FALSE,
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  
  # Reference for ~1000 bp:
  geom_vline(
    xintercept = 3,
    linetype = "dotted",
    linewidth = 0.4
  ) +
  
  facet_wrap(
    ~ Tissue,
    scales = "fixed",
    ncol = 6
  ) +
  
  labs(
    title =
      "Log10 gene-length distributions across GTEx tissues",
    
    subtitle =
      "Histogram with theoretical normal distribution; dotted line = 1,000 bp",
    
    x =
      "log10 gene length (bp)",
    
    y =
      "Number of genes"
  ) +
  
  theme_bw() +
  
  theme(
    
    plot.title =
      element_text(
        face = "bold"
      ),
    
    strip.text =
      element_text(
        size = 8
      ),
    
    axis.text =
      element_text(
        size = 7
      ),
    
    axis.title =
      element_text(
        size = 10
      )
  )


# ------------------------------------------------------------
# Save figure
# ------------------------------------------------------------

output_file <- file.path(
  figure_output_dir,
  "hpa_gtex_log10_gene_length_histogram_by_tissue_with_theoretical_normal.png"
)

ggsave(
  filename = output_file,
  plot = p_hist_theoretical,
  width = 14,
  height = 10,
  dpi = 300
)


# ------------------------------------------------------------
# Console output
# ------------------------------------------------------------

cat(
  "\nDone.\n",
  "Histogram bin width = ",
  binwidth,
  "\n",
  "1000 bp corresponds to log10 = 3\n",
  "Figure written to:\n",
  output_file,
  "\n",
  sep = ""
)