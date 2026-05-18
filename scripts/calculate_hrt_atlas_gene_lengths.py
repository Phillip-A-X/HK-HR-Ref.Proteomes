from pathlib import Path
import sys

import pandas as pd


# Allows importing from lib/ when running this script from the project root
sys.path.append(str(Path(__file__).resolve().parents[1]))

from lib.constants import (
    BIOMART_HRT_ATLAS_HUMAN_HK_TRANSCRIPTS_TO_ENSG_FILE,
    HRT_ATLAS_HUMAN_HK_GENE_LENGTHS_FILE,
    HRT_ATLAS_HUMAN_HK_TRANSCRIPT_TO_GENE_MAPPING_QC_FILE,
    HRT_ATLAS_HUMAN_HK_DUPLICATED_GENES_FILE,
    PROCESSED_DATA_DIR,
)


def main() -> None:
    PROCESSED_DATA_DIR.mkdir(parents=True, exist_ok=True)

    input_file = BIOMART_HRT_ATLAS_HUMAN_HK_TRANSCRIPTS_TO_ENSG_FILE

    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")

    df = pd.read_csv(input_file, sep="\t")

    required_columns = [
        "Transcript stable ID",
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

    transcripts_per_gene = (
        df.groupby("Gene stable ID")["Transcript stable ID"]
        .nunique()
        .reset_index(name="n_hk_transcripts")
        .sort_values("n_hk_transcripts", ascending=False)
    )

    duplicated_genes = transcripts_per_gene[
        transcripts_per_gene["n_hk_transcripts"] > 1
    ].copy()

    gene_lengths = (
        df[
            [
                "Gene stable ID",
                "Gene start (bp)",
                "Gene end (bp)",
                "Chromosome/scaffold name",
                "gene_length_bp",
            ]
        ]
        .drop_duplicates(subset=["Gene stable ID"])
        .sort_values("Gene stable ID")
        .copy()
    )

    gene_lengths.to_csv(
        HRT_ATLAS_HUMAN_HK_GENE_LENGTHS_FILE,
        sep="\t",
        index=False,
    )

    transcripts_per_gene.to_csv(
        HRT_ATLAS_HUMAN_HK_TRANSCRIPT_TO_GENE_MAPPING_QC_FILE,
        sep="\t",
        index=False,
    )

    duplicated_genes.to_csv(
        HRT_ATLAS_HUMAN_HK_DUPLICATED_GENES_FILE,
        sep="\t",
        index=False,
    )

    print("HRT Atlas human HK transcript-to-gene mapping")
    print(f"Input rows: {len(df):,}")
    print(f"Unique transcripts: {df['Transcript stable ID'].nunique():,}")
    print(f"Unique genes: {df['Gene stable ID'].nunique():,}")
    print(f"Genes with >1 HK transcript: {len(duplicated_genes):,}")
    print()
    print(f"Gene length table written to: {HRT_ATLAS_HUMAN_HK_GENE_LENGTHS_FILE}")
    print(f"Transcript-to-gene QC written to: {HRT_ATLAS_HUMAN_HK_TRANSCRIPT_TO_GENE_MAPPING_QC_FILE}")
    print(f"Duplicated genes table written to: {HRT_ATLAS_HUMAN_HK_DUPLICATED_GENES_FILE}")


if __name__ == "__main__":
    main()