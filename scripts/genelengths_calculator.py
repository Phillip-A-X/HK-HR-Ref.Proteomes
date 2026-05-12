from pathlib import Path
import sys

import pandas as pd


# libimport
sys.path.append(str(Path(__file__).resolve().parents[1]))

from lib.constants import BIOMART_EISENBERG_GENE_LENGTH_FILE, PROCESSED_DATA_DIR


OUTPUT_FILE = PROCESSED_DATA_DIR / "eisenberg_gene_lengths_with_length_bp.tsv"


def main() -> None:
    PROCESSED_DATA_DIR.mkdir(parents=True, exist_ok=True)

    if not BIOMART_EISENBERG_GENE_LENGTH_FILE.exists():
        raise FileNotFoundError(
            f"Input file not found: {BIOMART_EISENBERG_GENE_LENGTH_FILE}"
        )

    df = pd.read_csv(BIOMART_EISENBERG_GENE_LENGTH_FILE, sep="\t")

    required_columns = [
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
#accounts for counting of difference; eg if positions are counted one bp goes missing so one is added after the subtraction!
    df["gene_length_bp"] = df["Gene end (bp)"] - df["Gene start (bp)"] + 1

    df.to_csv(OUTPUT_FILE, sep="\t", index=False)

    print(f"Input rows: {len(df):,}")
    print(f"Output written to: {OUTPUT_FILE}")
    print(df.head())


if __name__ == "__main__":
    main()
