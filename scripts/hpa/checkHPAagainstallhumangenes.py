from pathlib import Path
import pandas as pd

# ------------------------------------------------------------
# QC:
# Check whether HPA-GTEx genes that did not map to the final
# Ensembl 115 main-chromosome protein-coding background
# are present in the unfiltered BioMart export.
#
# Inputs:
#   data/processed/hpa_gtex/human_hpa_gtex_unmatched_hpa_genes.tsv
#   data/raw/biomart/biomart_export_HumanAllGenes_220526.tsv
#
# Outputs:
#   data/processed/hpa_gtex/hpa_unmatched_vs_unfiltered_biomart.tsv
#   data/processed/hpa_gtex/hpa_unmatched_vs_unfiltered_biomart_summary.tsv
# ------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[2]

HPA_UNMATCHED_FILE = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "hpa_gtex"
    / "human_hpa_gtex_unmatched_hpa_genes.tsv"
)

UNFILTERED_BIOMART_FILE = (
    PROJECT_ROOT
    / "data"
    / "raw"
    / "biomart"
    / "biomart_export_HumanAllGenes_220526.tsv"
)

OUTPUT_DIR = PROJECT_ROOT / "data" / "processed" / "hpa_gtex"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

DETAIL_OUTPUT = OUTPUT_DIR / "hpa_unmatched_vs_unfiltered_biomart.tsv"
SUMMARY_OUTPUT = OUTPUT_DIR / "hpa_unmatched_vs_unfiltered_biomart_summary.tsv"

MAIN_CHROMOSOMES = {str(i) for i in range(1, 23)} | {"X", "Y"}


def remove_ensembl_version(x):
    if pd.isna(x):
        return x
    return str(x).split(".")[0]


# ------------------------------------------------------------
# Load HPA unmatched genes
# ------------------------------------------------------------

if not HPA_UNMATCHED_FILE.exists():
    raise FileNotFoundError(f"Missing file: {HPA_UNMATCHED_FILE}")

hpa_unmatched = pd.read_csv(HPA_UNMATCHED_FILE, sep="\t", dtype=str)

if "gene_id" not in hpa_unmatched.columns:
    raise ValueError("Expected column 'gene_id' in HPA unmatched file.")

hpa_unmatched["gene_id"] = hpa_unmatched["gene_id"].apply(remove_ensembl_version)

hpa_unmatched_ids = set(hpa_unmatched["gene_id"].dropna().unique())


# ------------------------------------------------------------
# Load unfiltered BioMart export
# ------------------------------------------------------------

if not UNFILTERED_BIOMART_FILE.exists():
    raise FileNotFoundError(f"Missing file: {UNFILTERED_BIOMART_FILE}")

biomart = pd.read_csv(UNFILTERED_BIOMART_FILE, sep="\t", dtype=str)

required_cols = [
    "Gene stable ID",
    "Gene start (bp)",
    "Gene end (bp)",
    "Chromosome/scaffold name",
]

missing = [col for col in required_cols if col not in biomart.columns]
if missing:
    raise ValueError(f"Missing required BioMart columns: {missing}")

biomart["gene_id"] = biomart["Gene stable ID"].apply(remove_ensembl_version)

# Keep one row per gene_id if duplicates exist.
# If the same gene appears multiple times on different scaffolds, keep all for detail first.
biomart_hits = biomart[biomart["gene_id"].isin(hpa_unmatched_ids)].copy()

biomart_hits["is_main_chromosome"] = biomart_hits["Chromosome/scaffold name"].isin(MAIN_CHROMOSOMES)

# Collapse to one row per HPA-unmatched gene for summary/detail
collapsed = (
    biomart_hits
    .groupby("gene_id", as_index=False)
    .agg(
        present_in_unfiltered_biomart=("gene_id", "size"),
        chromosomes_or_scaffolds=("Chromosome/scaffold name", lambda x: ";".join(sorted(set(x.dropna())))),
        any_main_chromosome=("is_main_chromosome", "any"),
        all_non_main_chromosome=("is_main_chromosome", lambda x: not any(x)),
    )
)

collapsed["present_in_unfiltered_biomart"] = collapsed["present_in_unfiltered_biomart"] > 0

# Add HPA-unmatched genes that were not found at all in unfiltered BioMart
all_hpa_unmatched = pd.DataFrame({"gene_id": sorted(hpa_unmatched_ids)})

detail = all_hpa_unmatched.merge(
    collapsed,
    on="gene_id",
    how="left",
)

detail["present_in_unfiltered_biomart"] = detail["present_in_unfiltered_biomart"].fillna(False)
detail["chromosomes_or_scaffolds"] = detail["chromosomes_or_scaffolds"].fillna("not_present_in_unfiltered_biomart")
detail["any_main_chromosome"] = detail["any_main_chromosome"].fillna(False)
detail["all_non_main_chromosome"] = detail["all_non_main_chromosome"].fillna(False)

def classify(row):
    if not row["present_in_unfiltered_biomart"]:
        return "absent_from_unfiltered_biomart"
    if row["all_non_main_chromosome"]:
        return "present_in_unfiltered_biomart_but_non_main_chromosome"
    if row["any_main_chromosome"]:
        return "present_in_unfiltered_biomart_on_main_chromosome"
    return "unclear"

detail["qc_class"] = detail.apply(classify, axis=1)

summary = (
    detail
    .groupby("qc_class", as_index=False)
    .agg(n_genes=("gene_id", "nunique"))
    .sort_values("qc_class")
)

# Add total row
total_row = pd.DataFrame(
    [{"qc_class": "total_hpa_unmatched_genes", "n_genes": len(hpa_unmatched_ids)}]
)

summary = pd.concat([summary, total_row], ignore_index=True)

# ------------------------------------------------------------
# Write outputs
# ------------------------------------------------------------

detail.to_csv(DETAIL_OUTPUT, sep="\t", index=False)
summary.to_csv(SUMMARY_OUTPUT, sep="\t", index=False)

print("\nDone.")

print("\nInput counts:")
print(f"HPA unmatched unique genes: {len(hpa_unmatched_ids)}")
print(f"Unfiltered BioMart rows: {len(biomart)}")
print(f"Unfiltered BioMart unique genes: {biomart['gene_id'].nunique()}")

print("\nSummary:")
print(summary.to_string(index=False))

print("\nOutputs:")
print(DETAIL_OUTPUT)
print(SUMMARY_OUTPUT)