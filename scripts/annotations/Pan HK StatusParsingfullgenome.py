#!/usr/bin/env python3

from pathlib import Path
import csv
import gzip
import re
from collections import Counter


# -------------------------------------------------------------------------
# Project paths
# -------------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[2]

GTF_FILE = (
    PROJECT_ROOT
    / "data"
    / "raw"
    / "ensembl"
    / "pan_troglodytes"
    / "release-115"
    / "Pan_troglodytes.Pan_tro_3.0.115.gtf.gz"
)

HK_LIST_FILE = (
    PROJECT_ROOT
    / "data"
    / "raw"
    / "hk_lists"
    / "pan_troglodytes"
    / "pan_hk_gene_list_joshi_2022.txt"
)

OUTPUT_DIR = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "pan_joshi_full115"
)

OUTPUT_FILE = (
    OUTPUT_DIR
    / "pan_troglodytes_Pan_tro_3.0_115_full_protein_coding_with_joshi_hk_status.tsv"
)

UNMAPPED_HK_FILE = (
    OUTPUT_DIR
    / "pan_joshi_2022_hk_ids_not_found_in_full115_protein_coding.tsv"
)

SUMMARY_FILE = (
    OUTPUT_DIR
    / "pan_joshi_full115_mapping_summary.tsv"
)


# -------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------

def parse_gtf_attributes(attribute_string):
    """
    Parse quoted key-value attributes from the ninth GTF column.
    """

    attributes = {}

    for match in re.finditer(r'(\S+) "([^"]+)"', attribute_string):
        key = match.group(1)
        value = match.group(2)
        attributes[key] = value

    return attributes


def normalise_pan_gene_id(gene_id):
    """
    Remove surrounding whitespace and optional Ensembl version suffixes.

    Example:
    ENSPTRG00000000001.2 -> ENSPTRG00000000001
    """

    gene_id = gene_id.strip()

    if "." in gene_id:
        gene_id = gene_id.split(".", 1)[0]

    return gene_id


def load_hk_gene_ids(hk_list_file):
    """
    Load Pan troglodytes Ensembl Gene IDs from the Joshi 2022 HK list.

    Empty lines, the header 'Genes' and invalid entries are ignored.
    Duplicate IDs are removed.
    """

    hk_gene_ids = set()
    invalid_entries = []

    with open(
        hk_list_file,
        "r",
        encoding="utf-8-sig"
    ) as infile:

        for line_number, line in enumerate(infile, start=1):
            raw_value = line.strip()

            if not raw_value:
                continue

            if raw_value.lower() == "genes":
                continue

            gene_id = normalise_pan_gene_id(raw_value)

            if not gene_id.startswith("ENSPTRG"):
                invalid_entries.append(
                    (line_number, raw_value)
                )
                continue

            hk_gene_ids.add(gene_id)

    if invalid_entries:
        print("\nWarning: invalid HK-list entries were ignored:")

        for line_number, value in invalid_entries[:20]:
            print(f"  line {line_number}: {value}")

        if len(invalid_entries) > 20:
            print(
                f"  ... and {len(invalid_entries) - 20} more"
            )

    return hk_gene_ids


# -------------------------------------------------------------------------
# Main processing
# -------------------------------------------------------------------------

def main():
    if not GTF_FILE.exists():
        raise FileNotFoundError(
            f"Complete Pan GTF not found: {GTF_FILE}"
        )

    if not HK_LIST_FILE.exists():
        raise FileNotFoundError(
            f"Pan Joshi HK list not found: {HK_LIST_FILE}"
        )

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    hk_gene_ids = load_hk_gene_ids(
        HK_LIST_FILE
    )

    print(f"Input GTF: {GTF_FILE}")
    print(f"Input HK list: {HK_LIST_FILE}")
    print(
        "Unique HK IDs loaded:",
        len(hk_gene_ids)
    )

    feature_counts = Counter()
    gene_biotype_counts = Counter()
    sequence_counts = Counter()
    hk_status_counts = Counter()

    all_gene_rows = 0
    protein_coding_rows = 0

    unique_all_gene_ids = set()
    unique_protein_coding_gene_ids = set()
    mapped_hk_gene_ids = set()

    missing_gene_id = 0
    missing_gene_name = 0
    missing_gene_biotype = 0

    with gzip.open(
        GTF_FILE,
        "rt",
        encoding="utf-8"
    ) as infile, open(
        OUTPUT_FILE,
        "w",
        newline="",
        encoding="utf-8"
    ) as outfile:

        writer = csv.writer(
            outfile,
            delimiter="\t"
        )

        writer.writerow([
            "species",
            "gene_id",
            "gene_name",
            "genetype",
            "chromosome",
            "start",
            "end",
            "gene_length_bp",
            "strand",
            "hk_status",
            "hk_source",
            "source_file"
        ])

        for line in infile:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 9:
                continue

            (
                seqname,
                source,
                feature,
                start,
                end,
                score,
                strand,
                frame,
                attributes_raw
            ) = fields

            feature_counts[feature] += 1

            if feature != "gene":
                continue

            all_gene_rows += 1

            attributes = parse_gtf_attributes(
                attributes_raw
            )

            gene_id = normalise_pan_gene_id(
                attributes.get("gene_id", "")
            )

            gene_name = attributes.get(
                "gene_name",
                ""
            )

            gene_biotype = attributes.get(
                "gene_biotype",
                attributes.get("gene_type", "")
            )

            if gene_id:
                unique_all_gene_ids.add(gene_id)
            else:
                missing_gene_id += 1

            if not gene_name:
                missing_gene_name += 1

            if gene_biotype:
                gene_biotype_counts[gene_biotype] += 1
            else:
                gene_biotype_counts["MISSING"] += 1
                missing_gene_biotype += 1

            sequence_counts[seqname] += 1

            if gene_biotype != "protein_coding":
                continue

            if not gene_id:
                continue

            protein_coding_rows += 1

            unique_protein_coding_gene_ids.add(
                gene_id
            )

            if gene_id in hk_gene_ids:
                hk_status = "HK"
                mapped_hk_gene_ids.add(gene_id)
            else:
                hk_status = "non-HK"

            hk_status_counts[hk_status] += 1

            start_int = int(start)
            end_int = int(end)

            gene_length_bp = (
                end_int
                - start_int
                + 1
            )

            if gene_length_bp <= 0:
                raise ValueError(
                    f"Invalid gene length for {gene_id}: "
                    f"start={start_int}, end={end_int}"
                )

            writer.writerow([
                "pan_troglodytes",
                gene_id,
                gene_name,
                gene_biotype,
                seqname,
                start_int,
                end_int,
                gene_length_bp,
                strand,
                hk_status,
                "Joshi et al. 2022",
                str(
                    GTF_FILE.relative_to(
                        PROJECT_ROOT
                    )
                )
            ])

    unmapped_hk_gene_ids = sorted(
        hk_gene_ids
        - mapped_hk_gene_ids
    )

    with open(
        UNMAPPED_HK_FILE,
        "w",
        newline="",
        encoding="utf-8"
    ) as outfile:

        writer = csv.writer(
            outfile,
            delimiter="\t"
        )

        writer.writerow([
            "gene_id",
            "reason"
        ])

        for gene_id in unmapped_hk_gene_ids:
            writer.writerow([
                gene_id,
                (
                    "Not found among complete Ensembl 115 "
                    "Pan troglodytes genes with exact biotype "
                    "protein_coding"
                )
            ])

    mapping_rate = (
        len(mapped_hk_gene_ids)
        / len(hk_gene_ids)
        * 100
        if hk_gene_ids
        else 0
    )

    summary_rows = [
        (
            "all_gene_rows_in_full_gtf",
            all_gene_rows
        ),
        (
            "unique_all_gene_ids_in_full_gtf",
            len(unique_all_gene_ids)
        ),
        (
            "protein_coding_gene_rows",
            protein_coding_rows
        ),
        (
            "unique_protein_coding_gene_ids",
            len(unique_protein_coding_gene_ids)
        ),
        (
            "unique_hk_ids_in_input_list",
            len(hk_gene_ids)
        ),
        (
            "mapped_hk_gene_ids",
            len(mapped_hk_gene_ids)
        ),
        (
            "unmapped_hk_gene_ids",
            len(unmapped_hk_gene_ids)
        ),
        (
            "hk_mapping_rate_percent",
            round(mapping_rate, 4)
        ),
        (
            "output_HK_rows",
            hk_status_counts["HK"]
        ),
        (
            "output_non_HK_rows",
            hk_status_counts["non-HK"]
        ),
        (
            "missing_gene_id_rows",
            missing_gene_id
        ),
        (
            "missing_gene_name_rows",
            missing_gene_name
        ),
        (
            "missing_gene_biotype_rows",
            missing_gene_biotype
        )
    ]

    with open(
        SUMMARY_FILE,
        "w",
        newline="",
        encoding="utf-8"
    ) as outfile:

        writer = csv.writer(
            outfile,
            delimiter="\t"
        )

        writer.writerow([
            "metric",
            "value"
        ])

        writer.writerows(
            summary_rows
        )

    print("\nGene counts")
    print(
        f"All gene rows in complete GTF: "
        f"{all_gene_rows}"
    )
    print(
        f"Unique gene IDs in complete GTF: "
        f"{len(unique_all_gene_ids)}"
    )
    print(
        f"Protein-coding gene rows: "
        f"{protein_coding_rows}"
    )
    print(
        f"Unique protein-coding gene IDs: "
        f"{len(unique_protein_coding_gene_ids)}"
    )

    print("\nJoshi 2022 HK mapping")
    print(
        f"Unique HK IDs in input list: "
        f"{len(hk_gene_ids)}"
    )
    print(
        f"Mapped HK IDs: "
        f"{len(mapped_hk_gene_ids)}"
    )
    print(
        f"Unmapped HK IDs: "
        f"{len(unmapped_hk_gene_ids)}"
    )
    print(
        f"Mapping rate: "
        f"{mapping_rate:.2f}%"
    )

    print("\nFinal classification")
    print(
        f"HK genes: "
        f"{hk_status_counts['HK']}"
    )
    print(
        f"non-HK genes: "
        f"{hk_status_counts['non-HK']}"
    )

    print("\nGene biotypes in complete GTF")
    for biotype, count in gene_biotype_counts.most_common():
        print(f"{biotype}\t{count}")

    print("\nSequences containing gene rows")
    for sequence, count in sequence_counts.most_common():
        print(f"{sequence}\t{count}")

    print("\nOutput files")
    print(OUTPUT_FILE)
    print(UNMAPPED_HK_FILE)
    print(SUMMARY_FILE)


if __name__ == "__main__":
    main()
