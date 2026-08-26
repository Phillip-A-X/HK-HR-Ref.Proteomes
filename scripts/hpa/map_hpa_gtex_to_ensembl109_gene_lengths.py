from pathlib import Path
import pandas as pd
import numpy as np

# ------------------------------------------------------------
# Robustness analysis:
# Map HPA-GTEx tissue expression data to an Ensembl 109
# human protein-coding main-chromosome gene-length background.
#
# Purpose:
#   Check whether using Ensembl 109 improves HPA-GTEx mapping
#   compared with the main Ensembl 115 analysis.
#
# Input:
#   data/raw/hpa/rna_tissue_gtex.tsv
#   data/raw/biomart/biomart_export_HumanAllGenes_Ensembl109.tsv
#
# Main output:
#   data/processed/hpa_gtex_ensembl109/HumanProteinCodingGenes_bytissue_onchromosomes_ensembl109.tsv
# ------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[2]

HPA_FILE = PROJECT_ROOT / "data" / "raw" / "hpa" / "rna_tissue_gtex.tsv"

ENSEMBL109_FILE = (
    PROJECT_ROOT
    / "data"
    / "raw"
    / "biomart"
    / "biomart_export_HumanAllGenes_Ensembl109.tsv"
)

OUTPUT_DIR = PROJECT_ROOT / "data" / "processed" / "hpa_gtex_ensembl109"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

GENE_LENGTH_OUTPUT = OUTPUT_DIR / "human_ensembl109_protein_coding_main_chromosome_gene_lengths.tsv"
FINAL_OUTPUT = OUTPUT_DIR / "HumanProteinCodingGenes_bytissue_onchromosomes_ensembl109.tsv"
MAPPED_ALL_OUTPUT = OUTPUT_DIR / "human_hpa_gtex_all_mapped_gene_tissue_pairs_ensembl109.tsv"
TISSUE_SUMMARY_OUTPUT = OUTPUT_DIR / "human_hpa_gtex_tissue_gene_length_summary_ensembl109.tsv"
QC_OUTPUT = OUTPUT_DIR / "human_hpa_gtex_mapping_qc_ensembl109.tsv"
UNMATCHED_HPA_OUTPUT = OUTPUT_DIR / "human_hpa_gtex_unmatched_hpa_genes_ensembl109.tsv"
ENSEMBL_WITHOUT_HPA_OUTPUT = OUTPUT_DIR / "human_ensembl109_genes_without_hpa_gtex.tsv"

EXPRESSION_THRESHOLD = 1.0
MAIN_CHROMOSOMES = {str(i) for i in range(1, 23)} | {"X", "Y"}


def remove_ensembl_version(x):
    if pd.isna(x):
        return x
    return str(x).split(".")[0]


# ------------------------------------------------------------
# Load HPA-GTEx
# ------------------------------------------------------------

print(f"Reading HPA GTEx file:\n{HPA_FILE}")

if not HPA_FILE.exists():
    raise FileNotFoundError(f"HPA file not found: {HPA_FILE}")

hpa = pd.read_csv(HPA_FILE, sep="\t", dtype=str)

required_hpa_columns = ["Gene", "Tissue", "nTPM"]
missing_hpa = [col for col in required_hpa_columns if col not in hpa.columns]

if missing_hpa:
    raise ValueError(f"Missing HPA columns: {missing_hpa}")

hpa["gene_id"] = hpa["Gene"].apply(remove_ensembl_version)
hpa["nTPM"] = pd.to_numeric(hpa["nTPM"], errors="coerce")

if "Gene name" not in hpa.columns:
    hpa["Gene name"] = pd.NA


# ------------------------------------------------------------
# Load Ensembl 109 BioMart export and create gene lengths
# ------------------------------------------------------------

print(f"\nReading Ensembl 109 BioMart file:\n{ENSEMBL109_FILE}")

if not ENSEMBL109_FILE.exists():
    raise FileNotFoundError(f"Ensembl 109 file not found: {ENSEMBL109_FILE}")

genes = pd.read_csv(ENSEMBL109_FILE, sep="\t", dtype=str)

required_gene_columns = [
    "Gene stable ID",
    "Gene start (bp)",
    "Gene end (bp)",
    "Chromosome/scaffold name",
]
missing_genes = [col for col in required_gene_columns if col not in genes.columns]

if missing_genes:
    raise ValueError(f"Missing Ensembl 109 columns: {missing_genes}")

genes["gene_id"] = genes["Gene stable ID"].apply(remove_ensembl_version)
genes["Gene start (bp)"] = pd.to_numeric(genes["Gene start (bp)"], errors="coerce")
genes["Gene end (bp)"] = pd.to_numeric(genes["Gene end (bp)"], errors="coerce")

genes["gene_length_bp"] = genes["Gene end (bp)"] - genes["Gene start (bp)"] + 1

genes_main = genes[genes["Chromosome/scaffold name"].isin(MAIN_CHROMOSOMES)].copy()
genes_main = genes_main.drop_duplicates(subset=["gene_id"]).copy()

genes_main.to_csv(GENE_LENGTH_OUTPUT, sep="\t", index=False)


# ------------------------------------------------------------
# Map HPA-GTEx to Ensembl 109 background
# ------------------------------------------------------------

mapped = hpa.merge(
    genes_main,
    on="gene_id",
    how="inner",
    suffixes=("_hpa", "_ensembl109"),
)

mapped["log10_gene_length_bp"] = np.log10(mapped["gene_length_bp"])

expressed = mapped[mapped["nTPM"] >= EXPRESSION_THRESHOLD].copy()


# ------------------------------------------------------------
# Final table
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
]

final_columns = [col for col in preferred_columns if col in expressed.columns]
remaining_columns = [col for col in expressed.columns if col not in final_columns]

expressed_final = expressed[final_columns + remaining_columns].copy()


# ------------------------------------------------------------
# QC
# ------------------------------------------------------------

hpa_genes = set(hpa["gene_id"].dropna().unique())
ensembl109_genes = set(genes_main["gene_id"].dropna().unique())
mapped_genes = set(mapped["gene_id"].dropna().unique())
expressed_genes = set(expressed["gene_id"].dropna().unique())

unmatched_hpa_genes = sorted(hpa_genes - ensembl109_genes)
ensembl_without_hpa = sorted(ensembl109_genes - hpa_genes)

qc = pd.DataFrame(
    [
        ["hpa_rows_gene_tissue_pairs", len(hpa)],
        ["hpa_unique_genes", len(hpa_genes)],
        ["ensembl109_raw_rows_with_header_removed", len(genes)],
        ["ensembl109_main_chromosome_background_genes", len(ensembl109_genes)],
        ["mapped_unique_genes", len(mapped_genes)],
        ["expressed_unique_genes_ntpm_ge_1", len(expressed_genes)],
        ["unmatched_hpa_genes", len(unmatched_hpa_genes)],
        ["ensembl109_genes_without_hpa_gtex", len(ensembl_without_hpa)],
        [
            "mapping_rate_hpa_to_ensembl109_percent",
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
# Tissue summary
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

print("\nQC:")
print(qc.to_string(index=False))

print("\nTissue summary preview:")
print(tissue_summary.head(10).to_string(index=False))