#!/usr/bin/env python3

from pathlib import Path
import csv
import re


PROJECT_ROOT = Path(__file__).resolve().parents[2]

HK_RAW_FILE = PROJECT_ROOT / "data" / "raw" / "Housekeeping and Essential Gene Lists" / "ecoli" / "E Coli Essential Gene List Gerdes 2003.tsv"

GENE_LENGTH_FILE = PROJECT_ROOT / "data" / "processed" / "gene_lengths" / "escherichia_coli_k12_mg1655_gene_lengths.tsv"

NORMALIZED_HK_FILE = PROJECT_ROOT / "data" / "processed" / "hk_gene_lists_normalized" / "ecoli_k12_mg1655_essential_gene_ids.tsv"

OUTPUT_FILE = PROJECT_ROOT / "data" / "processed" / "gene_lengths_with_hk_status" / "escherichia_coli_k12_mg1655_gene_lengths_with_hk_status.tsv"


def extract_essential_gene_ids():
    essential_gene_ids = set()

    with open(HK_RAW_FILE, "r", newline="", encoding="utf-8") as infile:
        reader = csv.reader(infile, delimiter="\t")

        for row in reader:
            if len(row) < 10:
                continue

            essentiality = row[3].strip()
            blattner_id = row[9].strip()

            if essentiality == "E" and re.fullmatch(r"b\d{4}", blattner_id):
                essential_gene_ids.add(blattner_id)

    return essential_gene_ids


def write_normalized_hk_ids(essential_gene_ids):
    NORMALIZED_HK_FILE.parent.mkdir(parents=True, exist_ok=True)

    with open(NORMALIZED_HK_FILE, "w", newline="", encoding="utf-8") as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        writer.writerow(["hk_gene_id"])

        for gene_id in sorted(essential_gene_ids):
            writer.writerow([gene_id])


def assign_hk_status(essential_gene_ids):
    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)

    with open(GENE_LENGTH_FILE, "r", newline="", encoding="utf-8") as infile, open(OUTPUT_FILE, "w", newline="", encoding="utf-8") as outfile:
        reader = csv.DictReader(infile, delimiter="\t")

        output_fieldnames = reader.fieldnames + ["hk_status", "hk_source"]

        writer = csv.DictWriter(outfile, delimiter="\t", fieldnames=output_fieldnames)
        writer.writeheader()

        total_genes = 0
        hk_genes = 0
        annotation_gene_ids = set()

        for row in reader:
            total_genes += 1

            gene_id = row["gene_id"]
            annotation_gene_ids.add(gene_id)

            if gene_id in essential_gene_ids:
                row["hk_status"] = "HK"
                row["hk_source"] = "Gerdes et al. 2003 essential-gene proxy"
                hk_genes += 1
            else:
                row["hk_status"] = "non-HK"
                row["hk_source"] = "not in Gerdes et al. 2003 essential-gene proxy"

            writer.writerow(row)

    unmatched_hk_ids = essential_gene_ids - annotation_gene_ids

    print(f"Total protein-coding genes: {total_genes}")
    print(f"Essential/HK-proxy genes in source list: {len(essential_gene_ids)}")
    print(f"Matched HK genes in annotation table: {hk_genes}")
    print(f"non-HK genes: {total_genes - hk_genes}")
    print(f"Unmatched HK IDs: {len(unmatched_hk_ids)}")
    print(f"Wrote normalized HK IDs to: {NORMALIZED_HK_FILE}")
    print(f"Wrote HK-status table to: {OUTPUT_FILE}")


def main():
    essential_gene_ids = extract_essential_gene_ids()
    write_normalized_hk_ids(essential_gene_ids)
    assign_hk_status(essential_gene_ids)


if __name__ == "__main__":
    main()
