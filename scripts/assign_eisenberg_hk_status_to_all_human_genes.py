from pathlib import Path
import sys

import pandas as pd


# Allows importing from lib/ when running this script from the project root
sys.path.append(str(Path(__file__).resolve().parents[1]))

from lib.constants import PROCESSED_DATA_DIR


# Human protein-coding main-chromosome background table
# The existing hk_status column from HRT Atlas is ignored here.
ALL_HUMAN_GENES_FILE = (
    PROCESSED_DATA_DIR / "human_all_genes_with_hrt_atlas_hk_status_rescued.tsv"
)

# Eisenberg HK genes mapped to ENSG and filtered to main chromosomes
EISENBERG_HK_GENES_FILE = (
    PROCESSED_DATA_DIR / "eisenberg_human_hk_gene_lengths_main_chromosomes.tsv"
)

OUTPUT_FILE = (
    PROCESSED_DATA_DIR / "human_all_genes_with_eisenberg_hk_status.tsv"
)


def main() -> None:
    PROCESSED_DATA_DIR.mkdir(parents=True, exist_ok=True)

    if not ALL_HUMAN_GENES_FILE.exists():
        raise FileNotFoundError(f"All human genes file not found: {ALL_HUMAN_GENES_FILE}")

    if not EISENBERG_HK_GENES_FILE.exists():
        raise FileNotFoundError(f"Eisenberg HK genes file not found: {EISENBERG_HK_GENES_FILE}")

    all_genes = pd.read_csv(ALL_HUMAN_GENES_FILE, sep="\t")
    eisenberg_hk = pd.read_csv(EISENBERG_HK_GENES_FILE, sep="\t")

    required_all_gene_columns = [
        "Gene stable ID",
        "Gene start (bp)",
        "Gene end (bp)",
        "Chromosome/scaffold name",
        "gene_length_bp",
    ]

    required_eisenberg_columns = [
        "Gene stable ID",
    ]

    missing_all = [col for col in required_all_gene_columns if col not in all_genes.columns]
    missing_eisenberg = [col for col in required_eisenberg_columns if col not in eisenberg_hk.columns]

    if missing_all:
        raise ValueError(
            f"Missing columns in all human genes table: {missing_all}\n"
            f"Available columns: {list(all_genes.columns)}"
        )

    if missing_eisenberg:
        raise ValueError(
            f"Missing columns in Eisenberg HK table: {missing_eisenberg}\n"
            f"Available columns: {list(eisenberg_hk.columns)}"
        )

    all_genes["Gene stable ID"] = all_genes["Gene stable ID"].astype(str).str.strip()
    eisenberg_hk["Gene stable ID"] = eisenberg_hk["Gene stable ID"].astype(str).str.strip()

    all_genes = all_genes.drop_duplicates(subset=["Gene stable ID"]).copy()
    eisenberg_hk = eisenberg_hk.drop_duplicates(subset=["Gene stable ID"]).copy()

    eisenberg_hk_ids = set(eisenberg_hk["Gene stable ID"])

    # Remove HRT-specific columns if present, because this table is Eisenberg-based
    columns_to_drop = [
        "hk_status",
        "hrt_rescue_applied",
    ]

    all_genes = all_genes.drop(
        columns=[col for col in columns_to_drop if col in all_genes.columns]
    )

    all_genes["hk_status"] = all_genes["Gene stable ID"].isin(eisenberg_hk_ids).map(
        {
            True: "HK",
            False: "non-HK",
        }
    )

    all_genes.to_csv(OUTPUT_FILE, sep="\t", index=False)

    print("Eisenberg HK status assignment to all human genes")
    print(f"All human background genes: {len(all_genes):,}")
    print(f"Eisenberg HK genes input: {len(eisenberg_hk):,}")
    print()
    print("HK status counts:")
    print(all_genes["hk_status"].value_counts())
    print()
    print(f"Output written to: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()