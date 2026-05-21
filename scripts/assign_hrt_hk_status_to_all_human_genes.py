from pathlib import Path
import sys

import pandas as pd


# Allows importing from lib/ when running this script from the project root
sys.path.append(str(Path(__file__).resolve().parents[1]))

from lib.constants import BIOMART_DIR, PROCESSED_DATA_DIR


ALL_GENES_FILE = BIOMART_DIR / "biomart_export_HumanAllGenes_220526.tsv"

HRT_HK_GENES_FILE = (
    BIOMART_DIR / "biomart_export_HRTatlashkgenelist_ENSGlength180526.tsv"
)

OUTPUT_FILE = (
    PROCESSED_DATA_DIR / "human_all_genes_with_hrt_atlas_hk_status.tsv"
)


def calculate_gene_length(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add gene_length_bp column using inclusive genomic coordinates.
    """
    df = df.copy()

    df["Gene start (bp)"] = pd.to_numeric(df["Gene start (bp)"], errors="raise")
    df["Gene end (bp)"] = pd.to_numeric(df["Gene end (bp)"], errors="raise")

    df["gene_length_bp"] = df["Gene end (bp)"] - df["Gene start (bp)"] + 1

    return df


def main() -> None:
    PROCESSED_DATA_DIR.mkdir(parents=True, exist_ok=True)

    if not ALL_GENES_FILE.exists():
        raise FileNotFoundError(f"All genes file not found: {ALL_GENES_FILE}")

    if not HRT_HK_GENES_FILE.exists():
        raise FileNotFoundError(f"HRT Atlas HK genes file not found: {HRT_HK_GENES_FILE}")

    all_genes = pd.read_csv(ALL_GENES_FILE, sep="\t")
    hrt_hk_genes = pd.read_csv(HRT_HK_GENES_FILE, sep="\t")

    required_columns = [
        "Gene stable ID",
        "Gene start (bp)",
        "Gene end (bp)",
        "Chromosome/scaffold name",
    ]

    for table_name, df in [
        ("all_genes", all_genes),
        ("hrt_hk_genes", hrt_hk_genes),
    ]:
        missing_columns = [col for col in required_columns if col not in df.columns]

        if missing_columns:
            raise ValueError(
                f"Missing columns in {table_name}: {missing_columns}\n"
                f"Available columns: {list(df.columns)}"
            )

    # Clean IDs
    all_genes["Gene stable ID"] = all_genes["Gene stable ID"].astype(str).str.strip()
    hrt_hk_genes["Gene stable ID"] = hrt_hk_genes["Gene stable ID"].astype(str).str.strip()

    # Calculate gene length for all genes
    all_genes = calculate_gene_length(all_genes)

    # Remove duplicate gene entries if BioMart returned repeated ENSG rows
    all_genes = all_genes.drop_duplicates(subset=["Gene stable ID"]).copy()
    hrt_hk_genes = hrt_hk_genes.drop_duplicates(subset=["Gene stable ID"]).copy()

    # Optional: keep only main chromosomes.
    # Current setting excludes MT and non-standard scaffolds/patches.
    main_chromosomes = [str(i) for i in range(1, 23)] + ["X", "Y"]

    all_genes = all_genes[
        all_genes["Chromosome/scaffold name"].astype(str).isin(main_chromosomes)
    ].copy()

    hrt_hk_gene_ids = set(hrt_hk_genes["Gene stable ID"])

    all_genes["hk_status"] = all_genes["Gene stable ID"].isin(hrt_hk_gene_ids).map(
        {
            True: "HK",
            False: "non-HK",
        }
    )

    all_genes.to_csv(OUTPUT_FILE, sep="\t", index=False)

    print("HRT Atlas HK status assignment to all human genes")
    print(f"All genes input rows: {len(pd.read_csv(ALL_GENES_FILE, sep='\t')):,}")
    print(f"All genes after unique ENSG + main chromosome filter: {len(all_genes):,}")
    print(f"HRT Atlas HK genes input rows: {len(hrt_hk_genes):,}")
    print()
    print("HK status counts:")
    print(all_genes["hk_status"].value_counts())
    print()
    print(f"Output written to: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()