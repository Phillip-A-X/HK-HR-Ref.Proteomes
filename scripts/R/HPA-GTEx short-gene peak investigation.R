# ============================================================
# HPA-GTEx short-gene peak investigation
#
# Peak region:
#   log10(gene length) = 2.85 - 3.15
#   approximately 708 - 1413 bp
#
# Aim:
#   Identify genes forming the short-gene peak and investigate
#   their HPA-GTEx tissue expression.
#
# Outputs:
#   results/hpa_gtex/tables/
# ============================================================


# ------------------------------------------------------------
# Setup
# ------------------------------------------------------------

library(dplyr)

setwd("//wsl$/Ubuntu/home/phillip/projects/HK-HR-Ref.Proteomes")

background_file <-
  "data/processed/human_all_genes_with_hrt_atlas_hk_status_rescued.tsv"

hpa_file <-
  "data/processed/hpa_gtex/human_hpa_gtex_all_mapped_gene_tissue_pairs.tsv"

output_dir <-
  "results/hpa_gtex/tables"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# Peak definition
# ------------------------------------------------------------

peak_min <- 2.85
peak_max <- 3.15


# ============================================================
# 1. Load genomic background
# ============================================================

background <- read.delim(
  background_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_background_columns <- c(
  "Gene stable ID",
  "gene_length_bp"
)

missing_background_columns <- setdiff(
  required_background_columns,
  colnames(background)
)

if (length(missing_background_columns) > 0) {
  stop(
    paste(
      "Missing background columns:",
      paste(missing_background_columns, collapse = ", ")
    )
  )
}


background <- background %>%
  filter(
    !is.na(`Gene stable ID`),
    !is.na(gene_length_bp),
    gene_length_bp > 0
  ) %>%
  distinct(
    `Gene stable ID`,
    .keep_all = TRUE
  ) %>%
  mutate(
    gene_id = `Gene stable ID`,
    log10_gene_length_bp = log10(gene_length_bp)
  )


# ------------------------------------------------------------
# Select genes forming the short-gene peak
# ------------------------------------------------------------

peak_genes <- background %>%
  filter(
    log10_gene_length_bp >= peak_min,
    log10_gene_length_bp <= peak_max
  ) %>%
  select(
    gene_id,
    gene_length_bp,
    log10_gene_length_bp
  ) %>%
  arrange(gene_length_bp)


cat(
  "\nGenes in peak region:",
  nrow(peak_genes),
  "\n"
)

cat(
  "Peak region:",
  peak_min,
  "-",
  peak_max,
  "log10(bp)\n"
)

cat(
  "Approximate bp range:",
  round(10^peak_min),
  "-",
  round(10^peak_max),
  "bp\n"
)


# ============================================================
# 2. Load ALL mapped HPA-GTEx gene-tissue pairs
# ============================================================

hpa <- read.delim(
  hpa_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# Detect relevant columns automatically
# ------------------------------------------------------------

find_column <- function(df, candidates) {
  
  found <- candidates[
    candidates %in% colnames(df)
  ]
  
  if (length(found) == 0) {
    return(NA_character_)
  }
  
  found[1]
}


gene_id_col <- find_column(
  hpa,
  c(
    "gene_id",
    "Gene",
    "Gene stable ID",
    "Gene.stable.ID"
  )
)

tissue_col <- find_column(
  hpa,
  c(
    "Tissue",
    "tissue"
  )
)

ntpm_col <- find_column(
  hpa,
  c(
    "nTPM",
    "ntpm"
  )
)

gene_name_col <- find_column(
  hpa,
  c(
    "gene_name",
    "Gene name",
    "Gene.name"
  )
)

description_col <- find_column(
  hpa,
  c(
    "description",
    "Description",
    "Gene description",
    "Gene.description"
  )
)


if (is.na(gene_id_col)) {
  stop("Could not identify gene ID column in HPA table.")
}

if (is.na(tissue_col)) {
  stop("Could not identify Tissue column in HPA table.")
}

if (is.na(ntpm_col)) {
  stop("Could not identify nTPM column in HPA table.")
}


# ------------------------------------------------------------
# Standardize column names
# ------------------------------------------------------------

hpa_standard <- data.frame(
  gene_id = as.character(hpa[[gene_id_col]]),
  Tissue = as.character(hpa[[tissue_col]]),
  nTPM = as.numeric(hpa[[ntpm_col]]),
  stringsAsFactors = FALSE
)

# Remove possible Ensembl version suffix:
# ENSG00000123456.15 -> ENSG00000123456

hpa_standard$gene_id <-
  sub(
    "\\..*$",
    "",
    hpa_standard$gene_id
  )


# ------------------------------------------------------------
# Add gene name if available
# ------------------------------------------------------------

if (!is.na(gene_name_col)) {
  
  hpa_standard$gene_name <-
    as.character(
      hpa[[gene_name_col]]
    )
  
} else {
  
  hpa_standard$gene_name <- NA_character_
  
  warning(
    "No gene-name column found in HPA table."
  )
}


# ------------------------------------------------------------
# Add description if available
# ------------------------------------------------------------

if (!is.na(description_col)) {
  
  hpa_standard$description <-
    as.character(
      hpa[[description_col]]
    )
  
} else {
  
  hpa_standard$description <- NA_character_
  
  warning(
    "No description column found. Description will be NA."
  )
}


# ============================================================
# 3. Restrict HPA data to peak genes
# ============================================================

peak_hpa_long <- peak_genes %>%
  left_join(
    hpa_standard,
    by = "gene_id"
  ) %>%
  mutate(
    expressed_ntpm_ge_1 =
      !is.na(nTPM) & nTPM >= 1
  )


# ------------------------------------------------------------
# Save detailed gene x tissue table
# ------------------------------------------------------------

write.table(
  peak_hpa_long,
  file = file.path(
    output_dir,
    "hpa_gtex_short_gene_peak_expression_by_tissue.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 4. Gene-level expression summary
# ============================================================


# ------------------------------------------------------------
# Gene name / description lookup
# ------------------------------------------------------------

annotation_summary <- hpa_standard %>%
  filter(
    gene_id %in% peak_genes$gene_id
  ) %>%
  group_by(gene_id) %>%
  summarise(
    
    gene_name = {
      x <- gene_name[
        !is.na(gene_name) &
          gene_name != ""
      ]
      
      if (length(x) == 0) {
        NA_character_
      } else {
        x[1]
      }
    },
    
    description = {
      x <- description[
        !is.na(description) &
          description != ""
      ]
      
      if (length(x) == 0) {
        NA_character_
      } else {
        x[1]
      }
    },
    
    .groups = "drop"
  )


# ------------------------------------------------------------
# Expression summary
# ------------------------------------------------------------

expression_summary <- hpa_standard %>%
  filter(
    gene_id %in% peak_genes$gene_id
  ) %>%
  group_by(gene_id) %>%
  summarise(
    
    HPA_present = TRUE,
    
    max_nTPM_total =
      if (
        all(is.na(nTPM))
      ) {
        NA_real_
      } else {
        max(
          nTPM,
          na.rm = TRUE
        )
      },
    
    number_expressed_tissues =
      n_distinct(
        Tissue[
          !is.na(nTPM) &
            nTPM >= 1
        ]
      ),
    
    tissue_with_max_expression = {
      
      if (
        all(is.na(nTPM)) ||
        max(nTPM, na.rm = TRUE) < 1
      ) {
        
        NA_character_
        
      } else {
        
        max_value <-
          max(
            nTPM,
            na.rm = TRUE
          )
        
        max_tissues <-
          unique(
            Tissue[
              !is.na(nTPM) &
                nTPM == max_value
            ]
          )
        
        paste(
          max_tissues,
          collapse = "; "
        )
      }
    },
    
    .groups = "drop"
  )

# ============================================================
# 5. Combine peak genes with expression information
# ============================================================

peak_gene_summary <- peak_genes %>%
  
  left_join(
    annotation_summary,
    by = "gene_id"
  ) %>%
  
  left_join(
    expression_summary,
    by = "gene_id"
  ) %>%
  
  mutate(
    
    HPA_present =
      ifelse(
        is.na(HPA_present),
        FALSE,
        HPA_present
      ),
    
    number_expressed_tissues =
      ifelse(
        is.na(number_expressed_tissues),
        0,
        number_expressed_tissues
      )
  ) %>%
  
  select(
    gene_id,
    gene_name,
    gene_length_bp,
    log10_gene_length_bp,
    description,
    HPA_present,
    max_nTPM_total,
    number_expressed_tissues,
    tissue_with_max_expression
  ) %>%
  
  arrange(
    gene_length_bp
  )


# ------------------------------------------------------------
# Save main summary table
# ------------------------------------------------------------

write.table(
  peak_gene_summary,
  file = file.path(
    output_dir,
    "hpa_gtex_short_gene_peak_expression_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 6. QC summary
# ============================================================

qc_summary <- data.frame(
  
  peak_log10_min =
    peak_min,
  
  peak_log10_max =
    peak_max,
  
  peak_bp_min =
    round(10^peak_min),
  
  peak_bp_max =
    round(10^peak_max),
  
  n_peak_genes =
    nrow(peak_gene_summary),
  
  n_HPA_present =
    sum(
      peak_gene_summary$HPA_present
    ),
  
  n_HPA_absent =
    sum(
      !peak_gene_summary$HPA_present
    ),
  
  n_expressed_in_at_least_one_tissue =
    sum(
      peak_gene_summary$number_expressed_tissues > 0
    ),
  
  n_never_ntpm_ge_1 =
    sum(
      peak_gene_summary$HPA_present &
        peak_gene_summary$number_expressed_tissues == 0
    )
)


write.table(
  qc_summary,
  file = file.path(
    output_dir,
    "hpa_gtex_short_gene_peak_expression_qc.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# Console output
# ============================================================

cat("\n========================================\n")
cat("SHORT-GENE PEAK ANALYSIS\n")
cat("========================================\n\n")

print(qc_summary)

cat(
  "\nMain summary:\n",
  file.path(
    output_dir,
    "hpa_gtex_short_gene_peak_expression_summary.tsv"
  ),
  "\n"
)

cat(
  "\nGene x tissue table:\n",
  file.path(
    output_dir,
    "hpa_gtex_short_gene_peak_expression_by_tissue.tsv"
  ),
  "\n"
)

cat(
  "\nQC table:\n",
  file.path(
    output_dir,
    "hpa_gtex_short_gene_peak_expression_qc.tsv"
  ),
  "\n"
)

cat("\nFirst peak genes:\n")

print(
  head(
    peak_gene_summary,
    20
  )
)