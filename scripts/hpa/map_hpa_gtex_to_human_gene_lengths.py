from pathlib import Path
import pandas as pd
import numpy as np

# ------------------------------------------------------------
# HPA GTEx tissue-expression mapping to human Ensembl 115
# protein-coding main-chromosome gene-length background
#
# Input:
#   data/raw/hpa/rna_tissue_gtex.tsv
#   data/processed/human_all_genes_with_hrt_atlas_hk_status_rescued.tsv
#
# Main output:
#   data/processed/hpa_gtex/HumanProteinCodingGenes_bytissue_onchromosomes.tsv
#
# Definition:
#   expressed in tissue = nTPM >= 1
# ------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[2]

HPA_FILE = PROJECT_ROOT / "data" / "raw" / "hpa" / "rna_tissue_gtex.tsv"

GENE_LENGTH_FILE = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "human_all_genes_with_hrt_atlas_hk_status_rescued.tsv"
)

OUTPUT_DIR = PROJECT_ROOT / "data" / "processed" / "hpa_gtex"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Main final output requested
FINAL_OUTPUT = OUTPUT_DIR / "HumanProteinCodingGenes_bytissue_onchromosomes.tsv"

# Additional QC / summary outputs
MAPPED_ALL_OUTPUT = OUTPUT_DIR / "human_hpa_gtex_all_mapped_gene_tissue_pairs.tsv"
TISSUE_SUMMARY_OUTPUT = OUTPUT_DIR / "human_hpa_gtex_tissue_gene_length_summary.tsv"
QC_OUTPUT = OUTPUT_DIR / "human_hpa_gtex_mapping_qc.tsv"
UNMATCHED_HPA_OUTPUT = OUTPUT_DIR / "human_hpa_gtex_unmatched_hpa_genes.tsv"
ENSEMBL_WITHOUT_HPA_OUTPUT = OUTPUT_DIR / "human_ensembl115_genes_without_hpa_gtex.tsv"

EXPRESSION_THRESHOLD = 1.0


def remove_ensembl_version(x):
    """
    Removes Ensembl version suffix if present.

    Example:
        ENSG00000123456.15 -> ENSG00000123456
        ENSG00000123456    -> ENSG00000123456
    """
    if pd.isna(x):
        return x
    return str(x).split(".")[0]


# ------------------------------------------------------------
# Load HPA GTEx tissue-expression data
# ------------------------------------------------------------

print(f"Reading HPA GTEx file:\n{HPA_FILE}")

if not HPA_FILE.exists():
    raise FileNotFoundError(f"HPA input file not found: {HPA_FILE}")

hpa = pd.read_csv(HPA_FILE, sep="\t", dtype=str)

print("\nHPA columns:")
print(list(hpa.columns))

required_hpa_columns = ["Gene", "Tissue", "nTPM"]
missing_hpa_columns = [col for col in required_hpa_columns if col not in hpa.columns]

if missing_hpa_columns:
    raise ValueError(f"Missing required HPA columns: {missing_hpa_columns}")

hpa["gene_id"] = hpa["Gene"].apply(remove_ensembl_version)
hpa["nTPM"] = pd.to_numeric(hpa["nTPM"], errors="coerce")

# Optional: keep original HPA gene symbol if available
if "Gene name" not in hpa.columns:
    hpa["Gene name"] = pd.NA


# ------------------------------------------------------------
# Load Ensembl 115 human protein-coding main-chromosome
# gene-length background
# ------------------------------------------------------------

print(f"\nReading human gene-length background:\n{GENE_LENGTH_FILE}")

if not GENE_LENGTH_FILE.exists():
    raise FileNotFoundError(f"Gene-length file not found: {GENE_LENGTH_FILE}")

genes = pd.read_csv(GENE_LENGTH_FILE, sep="\t", dtype=str)

print("\nGene-length columns:")
print(list(genes.columns))

required_gene_columns = ["Gene stable ID", "gene_length_bp"]
missing_gene_columns = [col for col in required_gene_columns if col not in genes.columns]

if missing_gene_columns:
    raise ValueError(f"Missing required gene-length columns: {missing_gene_columns}")

genes["gene_id"] = genes["Gene stable ID"].apply(remove_ensembl_version)
genes["gene_length_bp"] = pd.to_numeric(genes["gene_length_bp"], errors="coerce")

# Keep one row per Ensembl gene
genes_unique = genes.drop_duplicates(subset=["gene_id"]).copy()


# ------------------------------------------------------------
# Map HPA GTEx genes to Ensembl 115 protein-coding background
# ------------------------------------------------------------

mapped = hpa.merge(
    genes_unique,
    on="gene_id",
    how="inner",
    suffixes=("_hpa", "_ensembl115")
)

mapped["log10_gene_length_bp"] = np.log10(mapped["gene_length_bp"])

# Main biological filter:
# expressed in tissue = nTPM >= 1
expressed = mapped[mapped["nTPM"] >= EXPRESSION_THRESHOLD].copy()


# ------------------------------------------------------------
# Select and order columns for final table
# ------------------------------------------------------------

preferred_columns = [
    "gene_id",
    "Gene",
    "Gene name",
    "Tissue",
    "nTPM",
    "gene_length_bp",
    "log10_gene_length_bp",
    "Gene stable ID",
    "Gene start (bp)",
    "Gene end (bp)",
    "Chromosome/scaffold name",
    "hk_status",
]

final_columns = [col for col in preferred_columns if col in expressed.columns]
remaining_columns = [col for col in expressed.columns if col not in final_columns]

expressed_final = expressed[final_columns + remaining_columns].copy()


# ------------------------------------------------------------
# QC
# ------------------------------------------------------------

hpa_genes = set(hpa["gene_id"].dropna().unique())
ensembl_genes = set(genes_unique["gene_id"].dropna().unique())
mapped_genes = set(mapped["gene_id"].dropna().unique())
expressed_genes = set(expressed["gene_id"].dropna().unique())

unmatched_hpa_genes = sorted(hpa_genes - ensembl_genes)
ensembl_without_hpa = sorted(ensembl_genes - hpa_genes)

qc = pd.DataFrame(
    [
        ["hpa_rows_gene_tissue_pairs", len(hpa)],
        ["hpa_unique_genes", len(hpa_genes)],
        ["ensembl115_background_genes", len(ensembl_genes)],
        ["mapped_unique_genes", len(mapped_genes)],
        ["expressed_unique_genes_ntpm_ge_1", len(expressed_genes)],
        ["unmatched_hpa_genes", len(unmatched_hpa_genes)],
        ["ensembl115_genes_without_hpa_gtex", len(ensembl_without_hpa)],
        [
            "mapping_rate_hpa_to_ensembl115_percent",
            round(len(mapped_genes) / len(hpa_genes) * 100, 3),
        ],
        ["expression_threshold_nTPM", EXPRESSION_THRESHOLD],
        ["mapped_rows_gene_tissue_pairs", len(mapped)],
        ["expressed_rows_gene_tissue_pairs_ntpm_ge_1", len(expressed)],
        ["number_of_tissues_hpa", hpa["Tissue"].nunique()],
        ["number_of_tissues_after_mapping", mapped["Tissue"].nunique()],
        ["number_of_tissues_after_ntpm_filter", expressed["Tissue"].nunique()],
    ],
    columns=["metric", "value"],
)


# ------------------------------------------------------------
# Tissue-level gene-length summary
# ------------------------------------------------------------

tissue_summary = (
    expressed
    .groupby("Tissue", as_index=False)
    .agg(
        n_expressed_genes=("gene_id", "nunique"),
        mean_gene_length_bp=("gene_length_bp", "mean"),
        median_gene_length_bp=("gene_length_bp", "median"),
        variance_gene_length_bp=("gene_length_bp", "var"),
        sd_gene_length_bp=("gene_length_bp", "std"),
        min_gene_length_bp=("gene_length_bp", "min"),
        q1_gene_length_bp=("gene_length_bp", lambda x: x.quantile(0.25)),
        q3_gene_length_bp=("gene_length_bp", lambda x: x.quantile(0.75)),
        max_gene_length_bp=("gene_length_bp", "max"),
        mean_log10_gene_length_bp=("log10_gene_length_bp", "mean"),
        median_log10_gene_length_bp=("log10_gene_length_bp", "median"),
        variance_log10_gene_length_bp=("log10_gene_length_bp", "var"),
        sd_log10_gene_length_bp=("log10_gene_length_bp", "std"),
        mean_nTPM=("nTPM", "mean"),
        median_nTPM=("nTPM", "median"),
    )
    .sort_values("Tissue")
)


# ------------------------------------------------------------
# Write outputs
# ------------------------------------------------------------

mapped.to_csv(MAPPED_ALL_OUTPUT, sep="\t", index=False)
expressed_final.to_csv(FINAL_OUTPUT, sep="\t", index=False)
tissue_summary.to_csv(TISSUE_SUMMARY_OUTPUT, sep="\t", index=False)
qc.to_csv(QC_OUTPUT, sep="\t", index=False)

pd.DataFrame({"gene_id": unmatched_hpa_genes}).to_csv(
    UNMATCHED_HPA_OUTPUT,
    sep="\t",
    index=False,
)

pd.DataFrame({"gene_id": ensembl_without_hpa}).to_csv(
    ENSEMBL_WITHOUT_HPA_OUTPUT,
    sep="\t",
    index=False,
)


# ------------------------------------------------------------
# Console output
# ------------------------------------------------------------

print("\nDone.")

print("\nMain final output:")
print(FINAL_OUTPUT)

print("\nAdditional outputs:")
print(MAPPED_ALL_OUTPUT)
print(TISSUE_SUMMARY_OUTPUT)
print(QC_OUTPUT)
print(UNMATCHED_HPA_OUTPUT)
print(ENSEMBL_WITHOUT_HPA_OUTPUT)

print("\nQC:")
print(qc.to_string(index=False))

print("\nTissue summary preview:")
print(tissue_summary.head(10).to_string(index=False))
