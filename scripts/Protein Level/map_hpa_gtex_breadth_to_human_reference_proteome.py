from pathlib import Path
import pandas as pd


PROJECT_DIR = Path(".")
PROCESSED_DATA_DIR = PROJECT_DIR / "data" / "processed"
RAW_UNIPROT_DIR = PROJECT_DIR / "data" / "raw" / "uniprot"

PROTEIN_TABLE = (
    PROCESSED_DATA_DIR / "human_reference_proteome_with_hrt_atlas_hk_status.tsv"
)

UNIPROT_TO_ENSG_MAPPING = (
    RAW_UNIPROT_DIR / "human_reference_proteome_uniprot_to_ensg_idmapping_2026_08_13.tsv"
)

HPA_GTEX_EXPRESSED_TABLE = (
    PROCESSED_DATA_DIR / "hpa_gtex" / "HumanProteinCodingGenes_bytissue_onchromosomes.tsv"
)

OUTPUT_DIR = PROCESSED_DATA_DIR / "hpa_gtex_protein"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_ALL_MAPPED = (
    OUTPUT_DIR / "human_reference_proteome_with_hpa_gtex_breadth_all_mapped.tsv"
)

OUTPUT_LONGEST_PER_GENE = (
    OUTPUT_DIR / "human_reference_proteome_longest_protein_per_gene_with_hpa_gtex_breadth.tsv"
)

OUTPUT_QC = (
    OUTPUT_DIR / "human_reference_proteome_hpa_gtex_breadth_mapping_qc.tsv"
)


def remove_ensembl_version(x):
    if pd.isna(x):
        return x
    return str(x).split(".")[0]


def main():
    protein = pd.read_csv(PROTEIN_TABLE, sep="\t")
    mapping = pd.read_csv(UNIPROT_TO_ENSG_MAPPING, sep="\t")
    hpa = pd.read_csv(HPA_GTEX_EXPRESSED_TABLE, sep="\t")

    required_protein_cols = {
        "uniprot_accession",
        "entry_name",
        "description",
        "protein_length_aa",
        "hk_status",
    }

    required_mapping_cols = {"From", "To"}
    required_hpa_cols = {"gene_id", "Tissue", "nTPM", "hk_status"}

    missing_protein = required_protein_cols - set(protein.columns)
    missing_mapping = required_mapping_cols - set(mapping.columns)
    missing_hpa = required_hpa_cols - set(hpa.columns)

    if missing_protein:
        raise ValueError(f"Missing protein columns: {missing_protein}")
    if missing_mapping:
        raise ValueError(f"Missing mapping columns: {missing_mapping}")
    if missing_hpa:
        raise ValueError(f"Missing HPA-GTEx columns: {missing_hpa}")

    protein["uniprot_accession"] = protein["uniprot_accession"].astype(str).str.strip()
    protein["protein_length_aa"] = pd.to_numeric(protein["protein_length_aa"], errors="coerce")

    mapping["uniprot_accession"] = mapping["From"].astype(str).str.strip()
    mapping["gene_id"] = mapping["To"].apply(remove_ensembl_version)

    # Keep only Ensembl gene IDs
    mapping = mapping[mapping["gene_id"].astype(str).str.startswith("ENSG")].copy()

    hpa["gene_id"] = hpa["gene_id"].apply(remove_ensembl_version)

    # Gene-level expression breadth from expressed HPA-GTEx table.
    # This table already contains only gene-tissue pairs with nTPM >= 1.
    gene_summary = (
        hpa.groupby(["gene_id", "hk_status"], as_index=False)
        .agg(
            n_tissues_expressed=("Tissue", "nunique"),
            mean_nTPM_across_expressed_tissues=("nTPM", "mean"),
            median_nTPM_across_expressed_tissues=("nTPM", "median"),
        )
    )

    nonhk_breadth = gene_summary.loc[
        gene_summary["hk_status"] == "non-HK",
        "n_tissues_expressed",
    ]

    q1_nonhk = nonhk_breadth.quantile(0.25)
    q3_nonhk = nonhk_breadth.quantile(0.75)

    def assign_breadth_group(row):
        if row["hk_status"] == "HK":
            return "HK"
        if row["n_tissues_expressed"] <= q1_nonhk:
            return "restricted non-HK"
        if row["n_tissues_expressed"] >= q3_nonhk:
            return "broad non-HK"
        return "intermediate non-HK"

    gene_summary["breadth_group"] = gene_summary.apply(assign_breadth_group, axis=1)

    # Merge protein lengths with ENSG mapping
    protein_mapped_to_ensg = protein.merge(
        mapping[["uniprot_accession", "gene_id"]],
        on="uniprot_accession",
        how="left",
    )

    # Merge with HPA-GTEx breadth groups
    protein_with_breadth = protein_mapped_to_ensg.merge(
        gene_summary,
        on="gene_id",
        how="inner",
        suffixes=("_protein_table", "_hpa_gtex"),
    )

    # Sanity: use hk_status from HPA-GTEx / HRT gene-level table for group definition
    # but keep protein-level status for comparison.
    protein_with_breadth = protein_with_breadth.rename(
        columns={
            "hk_status_protein_table": "hk_status_protein_level",
            "hk_status_hpa_gtex": "hk_status_gene_level",
        }
    )

    protein_with_breadth.to_csv(
        OUTPUT_ALL_MAPPED,
        sep="\t",
        index=False,
    )

    # Longest protein sequence per ENSG
    longest_per_gene = (
        protein_with_breadth
        .sort_values(
            by=["gene_id", "protein_length_aa"],
            ascending=[True, False],
        )
        .drop_duplicates(subset=["gene_id"], keep="first")
        .copy()
    )

    longest_per_gene.to_csv(
        OUTPUT_LONGEST_PER_GENE,
        sep="\t",
        index=False,
    )

    qc = pd.DataFrame(
        {
            "metric": [
                "reference_proteome_proteins",
                "unique_reference_proteome_accessions",
                "uniprot_to_ensg_mapping_rows_after_ensg_filter",
                "unique_uniprot_accessions_with_ensg_mapping",
                "unique_ensg_ids_in_mapping",
                "hpa_gtex_expressed_genes",
                "protein_entries_mapped_to_hpa_gtex_genes",
                "unique_genes_with_protein_and_hpa_gtex_breadth",
                "longest_protein_per_gene_rows",
                "nonHK_expression_breadth_Q1",
                "nonHK_expression_breadth_Q3",
            ],
            "value": [
                len(protein),
                protein["uniprot_accession"].nunique(),
                len(mapping),
                mapping["uniprot_accession"].nunique(),
                mapping["gene_id"].nunique(),
                gene_summary["gene_id"].nunique(),
                len(protein_with_breadth),
                protein_with_breadth["gene_id"].nunique(),
                len(longest_per_gene),
                q1_nonhk,
                q3_nonhk,
            ],
        }
    )

    qc.to_csv(
        OUTPUT_QC,
        sep="\t",
        index=False,
    )

    print("Done.")
    print()
    print("QC:")
    print(qc.to_string(index=False))
    print()
    print("Breadth group counts, longest protein per gene:")
    print(longest_per_gene["breadth_group"].value_counts())
    print()
    print(f"Wrote: {OUTPUT_ALL_MAPPED}")
    print(f"Wrote: {OUTPUT_LONGEST_PER_GENE}")
    print(f"Wrote: {OUTPUT_QC}")


if __name__ == "__main__":
    main()