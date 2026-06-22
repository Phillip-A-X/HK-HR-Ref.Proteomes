#!/usr/bin/env python3

from pathlib import Path
import gzip
import re
from collections import Counter


PROJECT_ROOT = Path(__file__).resolve().parents[2]

matches = list(
    (PROJECT_ROOT / "data" / "raw").rglob(
        "Pan_troglodytes.Pan_tro_3.0.115.gtf.gz"
    )
)

if len(matches) == 0:
    raise FileNotFoundError("Full Pan troglodytes GTF not found")

if len(matches) > 1:
    raise RuntimeError(f"Multiple matching GTF files found: {matches}")

GTF_FILE = matches[0]


def parse_gtf_attributes(attribute_string):
    attributes = {}

    for match in re.finditer(r'(\S+) "([^"]+)"', attribute_string):
        key = match.group(1)
        value = match.group(2)
        attributes[key] = value

    return attributes


def main():
    if not GTF_FILE.exists():
        raise FileNotFoundError(f"GTF file not found: {GTF_FILE}")

    print(f"Using GTF file: {GTF_FILE}")
    print()

    feature_counts = Counter()
    gene_biotype_counts = Counter()
    sequence_counts = Counter()

    gene_rows = 0
    protein_coding_rows = 0

    unique_gene_ids = set()
    unique_protein_coding_gene_ids = set()

    missing_gene_id = 0
    missing_biotype = 0

    main_chromosome_gene_rows = 0
    other_sequence_gene_rows = 0

    main_chromosome_protein_coding_rows = 0
    other_sequence_protein_coding_rows = 0

    main_chromosomes = {
        "1", "2A", "2B", "3", "4", "5", "6", "7", "8", "9",
        "10", "11", "12", "13", "14", "15", "16", "17", "18",
        "19", "20", "21", "22", "X", "Y", "MT"
    }

    with gzip.open(GTF_FILE, "rt", encoding="utf-8") as infile:
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

            gene_rows += 1

            attributes = parse_gtf_attributes(attributes_raw)

            gene_id = attributes.get("gene_id", "")
            gene_biotype = attributes.get(
                "gene_biotype",
                attributes.get("gene_type", "")
            )

            if gene_id:
                unique_gene_ids.add(gene_id)
            else:
                missing_gene_id += 1

            if gene_biotype:
                gene_biotype_counts[gene_biotype] += 1
            else:
                gene_biotype_counts["MISSING"] += 1
                missing_biotype += 1

            sequence_counts[seqname] += 1

            if seqname in main_chromosomes:
                main_chromosome_gene_rows += 1
            else:
                other_sequence_gene_rows += 1

            if gene_biotype == "protein_coding":
                protein_coding_rows += 1

                if gene_id:
                    unique_protein_coding_gene_ids.add(gene_id)

                if seqname in main_chromosomes:
                    main_chromosome_protein_coding_rows += 1
                else:
                    other_sequence_protein_coding_rows += 1

    print("Gene counts")
    print(f"All gene rows: {gene_rows}")
    print(f"Unique gene IDs: {len(unique_gene_ids)}")
    print(f"Protein-coding gene rows: {protein_coding_rows}")
    print(
        "Unique protein-coding gene IDs:",
        len(unique_protein_coding_gene_ids)
    )
    print(f"Gene rows without gene_id: {missing_gene_id}")
    print(f"Gene rows without biotype: {missing_biotype}")

    print("\nChromosome versus other sequence counts")
    print(
        "All genes on main chromosomes:",
        main_chromosome_gene_rows
    )
    print(
        "All genes on other sequences:",
        other_sequence_gene_rows
    )
    print(
        "Protein-coding genes on main chromosomes:",
        main_chromosome_protein_coding_rows
    )
    print(
        "Protein-coding genes on other sequences:",
        other_sequence_protein_coding_rows
    )

    print("\nGene biotypes")
    for biotype, count in gene_biotype_counts.most_common():
        print(f"{biotype}\t{count}")

    print("\nTop sequence regions")
    for sequence, count in sequence_counts.most_common(40):
        print(f"{sequence}\t{count}")

    print("\nFeature counts")
    for feature, count in feature_counts.most_common():
        print(f"{feature}\t{count}")


if __name__ == "__main__":
    main()