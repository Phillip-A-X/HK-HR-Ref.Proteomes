# ============================================================
# HPA-GTEx unique expressed genes per nTPM cutoff
# ============================================================

library(dplyr)

input_file <- "data/processed/hpa_gtex/HumanProteinCodingGenes_bytissue_onchromosomes.tsv"

output_file <- "results/hpa_gtex/tables/hpa_gtex_cutoff_unique_gene_counts.tsv"

df <- read.delim(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

df$nTPM <- as.numeric(df$nTPM)

cutoffs <- c(1, 5, 10, 20)

cutoff_unique_gene_counts <- bind_rows(
  lapply(cutoffs, function(cutoff_value) {
    
    filtered <- df %>%
      filter(
        !is.na(gene_id),
        !is.na(Tissue),
        !is.na(nTPM),
        nTPM >= cutoff_value
      )
    
    data.frame(
      nTPM_cutoff = cutoff_value,
      unique_expressed_genes = n_distinct(filtered$gene_id),
      expressed_gene_tissue_pairs = nrow(
        distinct(filtered, Tissue, gene_id)
      ),
      tissues_retained = n_distinct(filtered$Tissue)
    )
  })
)

write.table(
  cutoff_unique_gene_counts,
  output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

print(cutoff_unique_gene_counts)

cat("\nWritten to:\n", output_file, "\n")