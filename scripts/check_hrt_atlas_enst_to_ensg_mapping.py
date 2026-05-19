from pathlib import Path
import sys

import pandas as pd



sys.path.append(str(Path(__file__).resolve().parents[1]))

from lib.constants import (
    HRT_ATLAS_HUMAN_HK_GENES_FILE,
    BIOMART_HRT_ATLAS_HUMAN_HK_ENST_TO_ENSG_FILE,
    HRT_ATLAS_HUMAN_HK_GENE_LENGTHS_FILE,
    HRT_ATLAS_TRANSCRIPTS_PER_GENE_QC_FILE,
    HRT_ATLAS_HUMAN_HK_DUPLICATED_GENES_FILE,
    HRT_ATLAS_UNMAPPED_ENST_FILE,
    HRT_ATLAS_ENST_TO_ENSG_MAPPING_SUMMARY_FILE,
    PROCESSED_DATA_DIR,
)


def main() -> None:
    PROCESSED_DATA_DIR.mkdir(parents=True, exist_ok=True)

    if not HRT_ATLAS_HUMAN_HK_GENES_FILE.exists():
        raise FileNotFoundError(
            f"Original HRT Atlas file not found: {HRT_ATLAS_HUMAN_HK_GENES_FILE}"
        )

    if not BIOMART_HRT_ATLAS_HUMAN_HK_ENST_TO_ENSG_FILE.exists():
        raise FileNotFoundError(
            f"BioMart mapping file not found: {BIOMART_HRT_ATLAS_HUMAN_HK_ENST_TO_ENSG_FILE}"
        )

    # Original HRT Atlas file has no header:
    # column 0 = Ensembl transcript ID (ENST...)
    # column 1 = Gene name
    # column 2 = RefSeq ID
    # column 3 = CCDS ID
    hrt = pd.read_csv(
        HRT_ATLAS_HUMAN_HK_GENES_FILE,
        sep="\t",
        header=None,
        names=["hrt_enst", "gene_name", "refseq_id", "ccds_id"],
    )

    biomart = pd.read_csv(
        BIOMART_HRT_ATLAS_HUMAN_HK_ENST_TO_ENSG_FILE,
        sep="\t",
    )

    required_columns = [
        "Transcript stable ID",
        "Gene stable ID",
        "Gene start (bp)",
        "Gene end (bp)",
        "Chromosome/scaffold name",
    ]

    missing_columns = [col for col in required_columns if col not in biomart.columns]

    if missing_columns:
        raise ValueError(
            f"Missing required BioMart columns: {missing_columns}\n"
            f"Available columns: {list(biomart.columns)}"
        )

    # Clean IDs
    hrt["hrt_enst"] = hrt["hrt_enst"].astype(str).str.strip()
    biomart["Transcript stable ID"] = biomart["Transcript stable ID"].astype(str).str.strip()
    biomart["Gene stable ID"] = biomart["Gene stable ID"].astype(str).str.strip()

    # Calculate gene length
    biomart["Gene start (bp)"] = pd.to_numeric(
        biomart["Gene start (bp)"],
        errors="raise",
    )

    biomart["Gene end (bp)"] = pd.to_numeric(
        biomart["Gene end (bp)"],
        errors="raise",
    )

    biomart["gene_length_bp"] = (
        biomart["Gene end (bp)"] - biomart["Gene start (bp)"] + 1
    )

    # Compare original HRT ENST IDs with BioMart mapped ENST IDs
    input_enst = set(hrt["hrt_enst"].dropna())
    mapped_enst = set(biomart["Transcript stable ID"].dropna())

    unmapped_enst = sorted(input_enst - mapped_enst)

    # Count how many HRT transcripts map to each ENSG
    transcripts_per_gene = (
        biomart
        .dropna(subset=["Transcript stable ID", "Gene stable ID"])
        .groupby("Gene stable ID")["Transcript stable ID"]
        .nunique()
        .reset_index(name="n_hk_transcripts")
        .sort_values("n_hk_transcripts", ascending=False)
    )

    duplicated_genes = transcripts_per_gene[
        transcripts_per_gene["n_hk_transcripts"] > 1
    ].copy()

    # One row per unique ENSG for gene-length analysis
    gene_lengths = (
        biomart[
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

    # Summary table
    summary = pd.DataFrame(
        [
            {
                "metric": "original_hrt_rows",
                "value": len(hrt),
            },
            {
                "metric": "unique_original_hrt_enst",
                "value": len(input_enst),
            },
            {
                "metric": "biomart_rows",
                "value": len(biomart),
            },
            {
                "metric": "unique_mapped_enst_in_biomart",
                "value": len(mapped_enst),
            },
            {
                "metric": "unmapped_enst",
                "value": len(unmapped_enst),
            },
            {
                "metric": "unique_mapped_ensg_genes",
                "value": biomart["Gene stable ID"].nunique(),
            },
            {
                "metric": "genes_with_more_than_one_hk_transcript",
                "value": len(duplicated_genes),
            },
        ]
    )

    # Write outputs
    gene_lengths.to_csv(
        HRT_ATLAS_HUMAN_HK_GENE_LENGTHS_FILE,
        sep="\t",
        index=False,
    )

    transcripts_per_gene.to_csv(
        HRT_ATLAS_TRANSCRIPTS_PER_GENE_QC_FILE,
        sep="\t",
        index=False,
    )

    duplicated_genes.to_csv(
        HRT_ATLAS_HUMAN_HK_DUPLICATED_GENES_FILE,
        sep="\t",
        index=False,
    )

    pd.Series(unmapped_enst, name="unmapped_enst").to_csv(
        HRT_ATLAS_UNMAPPED_ENST_FILE,
        sep="\t",
        index=False,
    )

    summary.to_csv(
        HRT_ATLAS_ENST_TO_ENSG_MAPPING_SUMMARY_FILE,
        sep="\t",
        index=False,
    )

    print("HRT Atlas ENST to ENSG mapping QC")
    print()
    print(f"Original HRT Atlas rows: {len(hrt):,}")
    print(f"Unique original HRT ENST IDs: {len(input_enst):,}")
    print(f"BioMart rows: {len(biomart):,}")
    print(f"Unique mapped ENST IDs in BioMart: {len(mapped_enst):,}")
    print(f"Unmapped ENST IDs: {len(unmapped_enst):,}")
    print(f"Unique mapped ENSG genes: {biomart['Gene stable ID'].nunique():,}")
    print(f"Genes with >1 HRT transcript: {len(duplicated_genes):,}")
    print()
    print("Output files:")
    print(f"- Gene lengths: {HRT_ATLAS_HUMAN_HK_GENE_LENGTHS_FILE}")
    print(f"- Transcripts per gene QC: {HRT_ATLAS_TRANSCRIPTS_PER_GENE_QC_FILE}")
    print(f"- Genes with multiple transcripts: {HRT_ATLAS_HUMAN_HK_DUPLICATED_GENES_FILE}")
    print(f"- Unmapped ENSTs: {HRT_ATLAS_UNMAPPED_ENST_FILE}")
    print(f"- Summary: {HRT_ATLAS_ENST_TO_ENSG_MAPPING_SUMMARY_FILE}")

    if len(unmapped_enst) == 0:
        print()
        print(
            "Interpretation: No original HRT Atlas ENST IDs are missing from the BioMart export. "
            "The reduction from transcript entries to unique genes is therefore explained by "
            "multiple transcripts mapping to the same Ensembl gene."
        )
    else:
        print()
        print(
            "Interpretation: Some original HRT Atlas ENST IDs were not returned by BioMart. "
            "The reduction is therefore caused by both transcript-to-gene collapsing and unmapped ENST IDs."
        )


if __name__ == "__main__":
    main()