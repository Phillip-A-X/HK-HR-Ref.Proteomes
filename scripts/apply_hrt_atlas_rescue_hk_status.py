from pathlib import Path
import sys

import pandas as pd


# Allows importing from lib/ when running this script from the project root
sys.path.append(str(Path(__file__).resolve().parents[1]))

from lib.constants import PROCESSED_DATA_DIR


FINAL_GENE_TABLE_FILE = (
    PROCESSED_DATA_DIR / "human_all_genes_with_hrt_atlas_hk_status.tsv"
)

RESCUE_CHECK_FILE = (
    PROCESSED_DATA_DIR / "hrt_atlas_rescue_refseq_ccds_match_check_gene_level.tsv"
)

OUTPUT_FILE = (
    PROCESSED_DATA_DIR / "human_all_genes_with_hrt_atlas_hk_status_rescued.tsv"
)

RESCUE_APPLIED_LOG_FILE = (
    PROCESSED_DATA_DIR / "hrt_atlas_rescue_applied_log.tsv"
)


ELIGIBLE_RESCUE_DECISION = (
    "eligible for rescue - unambiguous match to one original HRT Atlas entry"
)


def main() -> None:
    PROCESSED_DATA_DIR.mkdir(parents=True, exist_ok=True)

    if not FINAL_GENE_TABLE_FILE.exists():
        raise FileNotFoundError(f"Final gene table not found: {FINAL_GENE_TABLE_FILE}")

    if not RESCUE_CHECK_FILE.exists():
        raise FileNotFoundError(f"Rescue check file not found: {RESCUE_CHECK_FILE}")

    final = pd.read_csv(FINAL_GENE_TABLE_FILE, sep="\t")
    rescue = pd.read_csv(RESCUE_CHECK_FILE, sep="\t")

    required_final_columns = ["Gene stable ID", "hk_status"]
    required_rescue_columns = [
        "rescue_ensembl_gene_id",
        "gene_level_rescue_decision",
    ]

    missing_final = [col for col in required_final_columns if col not in final.columns]
    missing_rescue = [col for col in required_rescue_columns if col not in rescue.columns]

    if missing_final:
        raise ValueError(
            f"Missing columns in final table: {missing_final}\n"
            f"Available columns: {list(final.columns)}"
        )

    if missing_rescue:
        raise ValueError(
            f"Missing columns in rescue check table: {missing_rescue}\n"
            f"Available columns: {list(rescue.columns)}"
        )

    final["Gene stable ID"] = final["Gene stable ID"].astype(str).str.strip()
    rescue["rescue_ensembl_gene_id"] = rescue["rescue_ensembl_gene_id"].astype(str).str.strip()

    eligible_rescue = rescue[
        rescue["gene_level_rescue_decision"] == ELIGIBLE_RESCUE_DECISION
    ].copy()

    eligible_rescue_ids = set(eligible_rescue["rescue_ensembl_gene_id"])

    print(f"Eligible rescue genes: {len(eligible_rescue_ids):,}")

    rescue_present_mask = final["Gene stable ID"].isin(eligible_rescue_ids)
    rescue_present = final[rescue_present_mask].copy()

    rescue_not_in_final = sorted(eligible_rescue_ids - set(final["Gene stable ID"]))

    print(f"Eligible rescue genes present in final table: {len(rescue_present):,}")
    print(f"Eligible rescue genes missing from final table: {len(rescue_not_in_final):,}")

    print("Status before rescue:")
    print(rescue_present["hk_status"].value_counts())
    print()

    # Only change genes that are currently non-HK
    change_mask = rescue_present_mask & (final["hk_status"] == "non-HK")
    n_changed = int(change_mask.sum())

    final.loc[change_mask, "hk_status"] = "HK"

    # Optional trace column
    final["hrt_rescue_applied"] = False
    final.loc[change_mask, "hrt_rescue_applied"] = True

    rescue_log = final.loc[
        change_mask,
        [
            "Gene stable ID",
            "Gene start (bp)",
            "Gene end (bp)",
            "Chromosome/scaffold name",
            "gene_length_bp",
            "hk_status",
            "hrt_rescue_applied",
        ],
    ].copy()

    final.to_csv(OUTPUT_FILE, sep="\t", index=False)
    rescue_log.to_csv(RESCUE_APPLIED_LOG_FILE, sep="\t", index=False)

    print(f"Changed from non-HK to HK: {n_changed:,}")
    print()
    print("Final HK status counts after rescue:")
    print(final["hk_status"].value_counts())
    print()
    print(f"Updated final table written to: {OUTPUT_FILE}")
    print(f"Rescue applied log written to: {RESCUE_APPLIED_LOG_FILE}")

    if rescue_not_in_final:
        print()
        print("Warning: Some eligible rescue genes were not present in the final table:")
        for gene_id in rescue_not_in_final:
            print(gene_id)


if __name__ == "__main__":
    main()