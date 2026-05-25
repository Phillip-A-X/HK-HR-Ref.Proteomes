from pathlib import Path
import sys

import pandas as pd


# Allows importing from lib/ when running this script from the project root
sys.path.append(str(Path(__file__).resolve().parents[1]))

from lib.constants import RAW_DATA_DIR, BIOMART_DIR, PROCESSED_DATA_DIR


HRT_ATLAS_ORIGINAL_FILE = RAW_DATA_DIR / "Housekeeping_GenesHuman_HRTatlas.tsv"

RESCUE_ANNOTATION_FILE = (
    BIOMART_DIR
    / "biomart_export_humanHK_rescuemappedENSG_toENST_CCD_RefSeq250526.tsv"
)

OUTPUT_FILE = (
    PROCESSED_DATA_DIR
    / "hrt_atlas_rescue_refseq_ccds_match_check_gene_level.tsv"
)

SUMMARY_FILE = (
    PROCESSED_DATA_DIR
    / "hrt_atlas_rescue_refseq_ccds_match_summary_gene_level.tsv"
)


def clean_id(value) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def main() -> None:
    PROCESSED_DATA_DIR.mkdir(parents=True, exist_ok=True)

    if not HRT_ATLAS_ORIGINAL_FILE.exists():
        raise FileNotFoundError(f"HRT Atlas original file not found: {HRT_ATLAS_ORIGINAL_FILE}")

    if not RESCUE_ANNOTATION_FILE.exists():
        raise FileNotFoundError(f"Rescue annotation file not found: {RESCUE_ANNOTATION_FILE}")

    # Original HRT Atlas file has no header:
    # column 0 = ENST
    # column 1 = gene symbol
    # column 2 = RefSeq mRNA
    # column 3 = CCDS ID
    hrt = pd.read_csv(
        HRT_ATLAS_ORIGINAL_FILE,
        sep="\t",
        header=None,
        names=["hrt_enst", "hrt_gene_symbol", "hrt_refseq_mrna", "hrt_ccds_id"],
    )

    rescue = pd.read_csv(RESCUE_ANNOTATION_FILE, sep="\t")

    required_rescue_columns = [
        "Gene stable ID",
        "Transcript stable ID",
        "CCDS ID",
        "RefSeq mRNA ID",
    ]

    missing_columns = [col for col in required_rescue_columns if col not in rescue.columns]

    if missing_columns:
        raise ValueError(
            f"Missing required columns in rescue annotation file: {missing_columns}\n"
            f"Available columns: {list(rescue.columns)}"
        )

    # Clean columns
    for col in ["hrt_enst", "hrt_gene_symbol", "hrt_refseq_mrna", "hrt_ccds_id"]:
        hrt[col] = hrt[col].apply(clean_id)

    for col in required_rescue_columns:
        rescue[col] = rescue[col].apply(clean_id)

    output_rows = []

    # Decision unit = Gene stable ID / ENSG
    for ensg, group in rescue.groupby("Gene stable ID"):
        group = group.copy()

        rescue_ensts = sorted(set(group["Transcript stable ID"]) - {""})
        rescue_refseqs = sorted(set(group["RefSeq mRNA ID"]) - {""})
        rescue_ccds = sorted(set(group["CCDS ID"]) - {""})

        # Match HRT entries by RefSeq and CCDS
        refseq_matches = hrt[
            (hrt["hrt_refseq_mrna"] != "")
            & (hrt["hrt_refseq_mrna"].isin(rescue_refseqs))
        ].copy()

        ccds_matches = hrt[
            (hrt["hrt_ccds_id"] != "")
            & (hrt["hrt_ccds_id"].isin(rescue_ccds))
        ].copy()

        matched_indices = set(refseq_matches.index).union(set(ccds_matches.index))
        matched_hrt = hrt.loc[sorted(matched_indices)].copy()

        matched_hrt_ensts = sorted(set(matched_hrt["hrt_enst"]) - {""})
        matched_hrt_symbols = sorted(set(matched_hrt["hrt_gene_symbol"]) - {""})
        matched_hrt_refseqs = sorted(set(matched_hrt["hrt_refseq_mrna"]) - {""})
        matched_hrt_ccds = sorted(set(matched_hrt["hrt_ccds_id"]) - {""})

        n_matched_hrt_entries = len(matched_hrt)

        if n_matched_hrt_entries == 0:
            gene_level_rescue_decision = "exclude - no RefSeq or CCDS match to original HRT Atlas entry"
        elif n_matched_hrt_entries == 1:
            gene_level_rescue_decision = "eligible for rescue - unambiguous match to one original HRT Atlas entry"
        else:
            gene_level_rescue_decision = "exclude or manual review - ambiguous match to multiple original HRT Atlas entries"

        output_rows.append(
            {
                "rescue_ensembl_gene_id": ensg,
                "n_biomart_rows_for_gene": len(group),
                "n_rescue_transcripts_for_gene": len(rescue_ensts),
                "rescue_transcript_ids": ";".join(rescue_ensts),
                "rescue_refseq_mrna_ids": ";".join(rescue_refseqs),
                "rescue_ccds_ids": ";".join(rescue_ccds),
                "n_matched_original_hrt_entries": n_matched_hrt_entries,
                "matched_hrt_enst_ids": ";".join(matched_hrt_ensts),
                "matched_hrt_gene_symbols": ";".join(matched_hrt_symbols),
                "matched_hrt_refseq_mrna_ids": ";".join(matched_hrt_refseqs),
                "matched_hrt_ccds_ids": ";".join(matched_hrt_ccds),
                "gene_level_rescue_decision": gene_level_rescue_decision,
            }
        )

    result = pd.DataFrame(output_rows)

    summary = (
        result.groupby("gene_level_rescue_decision")
        .size()
        .reset_index(name="n_genes")
        .sort_values("n_genes", ascending=False)
    )

    result.to_csv(OUTPUT_FILE, sep="\t", index=False)
    summary.to_csv(SUMMARY_FILE, sep="\t", index=False)

    print("HRT Atlas rescue RefSeq/CCDS match check on gene level")
    print()
    print(f"Original HRT Atlas rows: {len(hrt):,}")
    print(f"Rescue BioMart rows: {len(rescue):,}")
    print(f"Unique rescue ENSG genes: {result['rescue_ensembl_gene_id'].nunique():,}")
    print()
    print("Gene-level decision summary:")
    print(summary)
    print()
    print(f"Detailed gene-level output written to: {OUTPUT_FILE}")
    print(f"Summary output written to: {SUMMARY_FILE}")


if __name__ == "__main__":
    main()