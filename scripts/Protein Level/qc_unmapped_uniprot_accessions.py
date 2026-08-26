from pathlib import Path
import gzip
import pandas as pd


PROJECT_DIR = Path(".")
RAW_UNIPROT_DIR = PROJECT_DIR / "data" / "raw" / "uniprot"
PROCESSED_DATA_DIR = PROJECT_DIR / "data" / "processed"

REFERENCE_FASTA = RAW_UNIPROT_DIR / "UP000005640_9606.fasta.gz"

PROTEIN_TABLE = (
    PROCESSED_DATA_DIR / "human_reference_proteome_with_hrt_atlas_hk_status.tsv"
)

UNIPROT_TO_ENSG_MAPPING = (
    RAW_UNIPROT_DIR / "human_reference_proteome_uniprot_to_ensg_idmapping_2026_08_13.tsv"
)

OUTPUT_DIR = PROCESSED_DATA_DIR / "hpa_gtex_protein"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_UNMAPPED_TABLE = (
    OUTPUT_DIR / "uniprot_accessions_not_mapped_to_ensg_qc.tsv"
)

OUTPUT_SUMMARY = (
    OUTPUT_DIR / "uniprot_accessions_not_mapped_to_ensg_summary.tsv"
)

OUTPUT_MAPPED_VS_UNMAPPED_LENGTH_SUMMARY = (
    OUTPUT_DIR / "uniprot_mapped_vs_unmapped_length_summary.tsv"
)


def parse_fasta_header_types(fasta_file: Path) -> pd.DataFrame:
    """
    Extract UniProt accession and database prefix from FASTA headers.

    Example:
        >sp|O60762|DPM1_HUMAN ...
        >tr|A0A...|...
    """
    records = []

    with gzip.open(fasta_file, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                header = line.strip()[1:]
                id_part = header.split(" ", 1)[0]
                fields = id_part.split("|")

                if len(fields) >= 3:
                    db_prefix = fields[0]
                    accession = fields[1]
                    entry_name = fields[2]

                    records.append(
                        {
                            "uniprot_accession": accession,
                            "uniprot_db_prefix": db_prefix,
                            "entry_name_from_fasta": entry_name,
                        }
                    )

    return pd.DataFrame(records)


def summarize_lengths(df: pd.DataFrame, group_col: str) -> pd.DataFrame:
    return (
        df.groupby(group_col)
        .agg(
            n=("uniprot_accession", "nunique"),
            mean_protein_length_aa=("protein_length_aa", "mean"),
            median_protein_length_aa=("protein_length_aa", "median"),
            sd_protein_length_aa=("protein_length_aa", "std"),
            min_protein_length_aa=("protein_length_aa", "min"),
            max_protein_length_aa=("protein_length_aa", "max"),
        )
        .reset_index()
    )


def main() -> None:
    protein = pd.read_csv(PROTEIN_TABLE, sep="\t")
    mapping = pd.read_csv(UNIPROT_TO_ENSG_MAPPING, sep="\t")
    fasta_types = parse_fasta_header_types(REFERENCE_FASTA)

    required_protein_cols = {
        "uniprot_accession",
        "entry_name",
        "description",
        "protein_length_aa",
        "hk_status",
    }

    missing = required_protein_cols - set(protein.columns)
    if missing:
        raise ValueError(f"Missing protein table columns: {missing}")

    protein["uniprot_accession"] = protein["uniprot_accession"].astype(str).str.strip()
    mapping["From"] = mapping["From"].astype(str).str.strip()

    mapped_accessions = set(mapping["From"].dropna().astype(str).str.strip())
    all_accessions = set(protein["uniprot_accession"])

    unmapped_accessions = all_accessions - mapped_accessions

    protein["ensg_mapping_status"] = protein["uniprot_accession"].apply(
        lambda x: "mapped_to_ensg" if x in mapped_accessions else "not_mapped_to_ensg"
    )

    protein = protein.merge(
        fasta_types,
        on="uniprot_accession",
        how="left",
    )

    protein["review_status_from_fasta"] = protein["uniprot_db_prefix"].map(
        {
            "sp": "Swiss-Prot reviewed",
            "tr": "TrEMBL unreviewed",
        }
    )

    unmapped = protein[protein["ensg_mapping_status"] == "not_mapped_to_ensg"].copy()

    unmapped.to_csv(
        OUTPUT_UNMAPPED_TABLE,
        sep="\t",
        index=False,
    )

    # Summary metrics
    summary_rows = []

    def add_metric(metric, value):
        summary_rows.append({"metric": metric, "value": value})

    add_metric("reference_proteome_total_accessions", len(all_accessions))
    add_metric("mapped_to_ensg_accessions", len(mapped_accessions & all_accessions))
    add_metric("not_mapped_to_ensg_accessions", len(unmapped_accessions))
    add_metric(
        "not_mapped_to_ensg_percent",
        100 * len(unmapped_accessions) / len(all_accessions),
    )

    # HK/non-HK affected
    hk_counts = (
        unmapped["hk_status"]
        .value_counts(dropna=False)
        .rename_axis("hk_status")
        .reset_index(name="n")
    )

    for _, row in hk_counts.iterrows():
        add_metric(f"unmapped_hk_status_{row['hk_status']}", row["n"])

    # sp/tr affected
    review_counts = (
        unmapped["review_status_from_fasta"]
        .value_counts(dropna=False)
        .rename_axis("review_status_from_fasta")
        .reset_index(name="n")
    )

    for _, row in review_counts.iterrows():
        add_metric(f"unmapped_review_status_{row['review_status_from_fasta']}", row["n"])

    # Length stats of unmapped
    add_metric("unmapped_mean_protein_length_aa", unmapped["protein_length_aa"].mean())
    add_metric("unmapped_median_protein_length_aa", unmapped["protein_length_aa"].median())
    add_metric("unmapped_min_protein_length_aa", unmapped["protein_length_aa"].min())
    add_metric("unmapped_max_protein_length_aa", unmapped["protein_length_aa"].max())

    summary = pd.DataFrame(summary_rows)

    summary.to_csv(
        OUTPUT_SUMMARY,
        sep="\t",
        index=False,
    )

    mapped_vs_unmapped_summary = summarize_lengths(
        protein,
        "ensg_mapping_status",
    )

    mapped_vs_unmapped_summary.to_csv(
        OUTPUT_MAPPED_VS_UNMAPPED_LENGTH_SUMMARY,
        sep="\t",
        index=False,
    )

    print("Done.")
    print()
    print("Main mapping QC:")
    print(summary.to_string(index=False))
    print()
    print("Mapped vs unmapped length summary:")
    print(mapped_vs_unmapped_summary.to_string(index=False))
    print()
    print(f"Wrote: {OUTPUT_UNMAPPED_TABLE}")
    print(f"Wrote: {OUTPUT_SUMMARY}")
    print(f"Wrote: {OUTPUT_MAPPED_VS_UNMAPPED_LENGTH_SUMMARY}")


if __name__ == "__main__":
    main()