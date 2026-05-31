from pathlib import Path
import sys

import pandas as pd


# Allows importing from lib/ when running this script from the project root
sys.path.append(str(Path(__file__).resolve().parents[1]))

from lib.constants import BIOMART_DIR, PROCESSED_DATA_DIR


INPUT_FILE = (
    BIOMART_DIR
    / "biomart_export_humanHK_eisenberg.tsv"
)

OUTPUT_FILE = (
    PROCESSED_DATA_DIR
    / "eisenberg_human_hk_gene_lengths_main_chromosomes.tsv"
)


def main() -> None:
    PROCESSED_DATA_DIR.mkdir(parents=True, exist_ok=True)

    if not INPUT_FILE.exists():
        raise FileNotFoundError(f"Input file not found: {INPUT_FILE}")

    df = pd.read_csv(INPUT_FILE, sep="\t")

    required_columns = [
        "HGNC symbol",
        "Gene stable ID",
        "Gene start (bp)",
        "Gene end (bp)",
        "Chromosome/scaffold name",
    ]

    missing_columns = [col for col in required_columns if col not in df.columns]

    if missing_columns:
        raise ValueError(
            f"Missing required columns: {missing_columns}\n"
            f"Available columns: {list(df.columns)}"
        )

    df["Gene start (bp)"] = pd.to_numeric(df["Gene start (bp)"], errors="raise")
    df["Gene end (bp)"] = pd.to_numeric(df["Gene end (bp)"], errors="raise")

    df["gene_length_bp"] = df["Gene end (bp)"] - df["Gene start (bp)"] + 1

    # Keep only main nuclear chromosomes, exclude MT, patches, scaffolds, haplotypes etc.
    main_chromosomes = [str(i) for i in range(1, 23)] + ["X", "Y"]

    df_main = df[
        df["Chromosome/scaffold name"].astype(str).isin(main_chromosomes)
    ].copy()

    # One row per Ensembl gene ID
    df_main = df_main.drop_duplicates(subset=["Gene stable ID"]).copy()

    df_main.to_csv(OUTPUT_FILE, sep="\t", index=False)

    print("Eisenberg HK gene length processing")
    print(f"Input rows: {len(df):,}")
    print(f"Rows on main chromosomes: {len(df_main):,}")
    print(f"Unique ENSG genes on main chromosomes: {df_main['Gene stable ID'].nunique():,}")
    print(f"Output written to: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()


    print(df["Gene stable ID"].duplicated().sum())