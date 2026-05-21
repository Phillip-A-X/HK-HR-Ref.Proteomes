from pathlib import Path
import sys

import pandas as pd


# Allows importing from lib/ when running this script from the project root
sys.path.append(str(Path(__file__).resolve().parents[1]))

from lib.constants import BIOMART_DIR, PROCESSED_DATA_DIR


RESCUE_MAPPING_FILE = (
    BIOMART_DIR
    / "biomart_export_HRTAtlas_35unmappedENST_HGNCsymbol_rescueMapping_toENSG_geneCoords_200526.tsv"
)

FINAL_HRT_HK_STATUS_FILE = (
    PROCESSED_DATA_DIR
    / "human_all_genes_with_hrt_atlas_hk_status.tsv"
)

OUTPUT_FILE = (
    PROCESSED_DATA_DIR
    / "hrt_atlas_rescue_mapping_check_against_final_table.tsv"
)


def main() -> None:
    PROCESSED_DATA_DIR.mkdir(parents=True, exist_ok=True)

    if not RESCUE_MAPPING_FILE.exists():
        raise FileNotFoundError(f"Rescue mapping file not found: {RESCUE_MAPPING_FILE}")

    if not FINAL_HRT_HK_STATUS_FILE.exists():
        raise FileNotFoundError(f"Final HK status table not found: {FINAL_HRT_HK_STATUS_FILE}")

    rescue = pd.read_csv(RESCUE_MAPPING_FILE, sep="\t")
    final = pd.read_csv(FINAL_HRT_HK_STATUS_FILE, sep="\t")

    required_rescue_columns = [
        "Gene stable ID",
        "Gene start (bp)",
        "Gene end (bp)",
        "Chromosome/scaffold name",
    ]

    required_final_columns = [
        "Gene stable ID",
        "Gene start (bp)",
        "Gene end (bp)",
        "Chromosome/scaffold name",
        "gene_length_bp",
        "hk_status",
    ]

    missing_rescue = [col for col in required_rescue_columns if col not in rescue.columns]
    missing_final = [col for col in required_final_columns if col not in final.columns]

    if missing_rescue:
        raise ValueError(
            f"Missing columns in rescue mapping file: {missing_rescue}\n"
            f"Available columns: {list(rescue.columns)}"
        )

    if missing_final:
        raise ValueError(
            f"Missing columns in final HK status table: {missing_final}\n"
            f"Available columns: {list(final.columns)}"
        )

    rescue["Gene stable ID"] = rescue["Gene stable ID"].astype(str).str.strip()
    final["Gene stable ID"] = final["Gene stable ID"].astype(str).str.strip()

    # Remove exact duplicate rescue ENSG rows if present
    rescue_unique = rescue.drop_duplicates(subset=["Gene stable ID"]).copy()

    final_subset = final[
        [
            "Gene stable ID",
            "Gene start (bp)",
            "Gene end (bp)",
            "Chromosome/scaffold name",
            "gene_length_bp",
            "hk_status",
        ]
    ].copy()

    merged = rescue_unique.merge(
        final_subset,
        on="Gene stable ID",
        how="left",
        suffixes=("_rescue", "_final"),
    )

    def classify(row):
        if pd.isna(row["hk_status"]):
            return "absent_from_final_table_do_not_add_automatically"

        if row["hk_status"] == "HK":
            return "already_present_as_HK_no_change_needed"

        if row["hk_status"] == "non-HK":
            return "present_as_nonHK_possible_rescue_candidate_check_manually"

        return "present_with_unknown_status_check_manually"

    merged["rescue_decision"] = merged.apply(classify, axis=1)

    merged.to_csv(OUTPUT_FILE, sep="\t", index=False)

    print("HRT Atlas rescue mapping check")
    print()
    print(f"Rescue mapping input rows: {len(rescue):,}")
    print(f"Unique rescue ENSG IDs: {rescue_unique['Gene stable ID'].nunique():,}")
    print(f"Final table rows: {len(final):,}")
    print()
    print("Decision counts:")
    print(merged["rescue_decision"].value_counts())
    print()
    print(f"Output written to: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()