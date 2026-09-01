library(dplyr)

# ------------------------------------------------------------
# Add mitochondrial protein-coding genes to human gene table
#
# Existing nuclear background:
#   20,084 protein-coding genes on chr 1-22, X, Y
#
# MT genes:
#   extracted from original Ensembl 115 BioMart export
#   hk_status = NA
#   MT = TRUE
#
# Existing nuclear genes:
#   existing hk_status retained
#   MT = FALSE
#
# Original table remains unchanged.
# ------------------------------------------------------------


# Existing processed nuclear gene table
nuclear_file <- "data/processed/human_all_genes_with_hrt_atlas_hk_status_rescued.tsv"

# Select original BioMart export:
# biomart_export_HumanAllGenes_220526.tsv
biomart_file <- "//wsl$/Ubuntu/home/phillip/projects/HK-HR-Ref.Proteomes/data/raw/biomart/biomart_export_HumanAllGenes_220526.tsv"


# Read files
nuclear <- read.delim(
  nuclear_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

biomart <- read.delim(
  biomart_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# Extract mitochondrial protein-coding genes
# ------------------------------------------------------------

mt_genes <- biomart %>%
  filter(`Chromosome/scaffold name` == "MT") %>%
  mutate(
    gene_length_bp =
      `Gene end (bp)` - `Gene start (bp)` + 1,
    hk_status = NA_character_,
    MT = TRUE
  )


# Mark all existing nuclear genes as non-mitochondrial
nuclear <- nuclear %>%
  mutate(
    MT = FALSE
  )


# ------------------------------------------------------------
# Make MT table compatible with existing processed table
# ------------------------------------------------------------

# Add any columns present in nuclear table but absent from MT table
missing_cols <- setdiff(
  names(nuclear),
  names(mt_genes)
)

for (col in missing_cols) {
  mt_genes[[col]] <- NA
}

# Keep exactly the same columns and order
mt_genes <- mt_genes[
  ,
  names(nuclear),
  drop = FALSE
]


# ------------------------------------------------------------
# Combine nuclear + mitochondrial genes
# ------------------------------------------------------------

human_plus_mt <- bind_rows(
  nuclear,
  mt_genes
)


# ------------------------------------------------------------
# QC
# ------------------------------------------------------------

cat("\n--- QC ---\n")

cat(
  "Original nuclear genes:",
  nrow(nuclear),
  "\n"
)

cat(
  "MT genes added:",
  sum(human_plus_mt$MT),
  "\n"
)

cat(
  "Total genes:",
  nrow(human_plus_mt),
  "\n"
)

cat(
  "MT genes with hk_status NA:",
  sum(
    human_plus_mt$MT &
      is.na(human_plus_mt$hk_status)
  ),
  "\n"
)

cat(
  "Duplicated Gene stable IDs:",
  sum(
    duplicated(
      human_plus_mt$`Gene stable ID`
    )
  ),
  "\n"
)


# Show MT genes
print(
  human_plus_mt %>%
    filter(MT == TRUE) %>%
    select(
      `Gene stable ID`,
      `Gene start (bp)`,
      `Gene end (bp)`,
      `Chromosome/scaffold name`,
      gene_length_bp,
      hk_status,
      MT
    )
)


# ------------------------------------------------------------
# Save NEW table
# ------------------------------------------------------------

output_file <-
  "data/processed/human_all_genes_with_hrt_atlas_hk_status_rescued_plus_mt.tsv"

write.table(
  human_plus_mt,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

cat(
  "\nSaved:",
  output_file,
  "\n"
)