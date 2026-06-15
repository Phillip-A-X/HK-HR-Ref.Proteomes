#!/usr/bin/env python3

from pathlib import Path
import gzip
import csv
from urllib.parse import unquote


PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT_BASE = PROJECT_ROOT / "data" / "raw" / "ensembl_genomes" / "bacteria"
OUTPUT_DIR = PROJECT_ROOT / "data" / "processed" / "gene_lengths"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def parse_gff3_attributes(attribute_string):
    """
    Parse the attributes column of a GFF3 file.
    Example:
    ID=gene:b0001;Name=thrL;biotype=protein_coding;gene_id=b0001
    """
    attributes = {}

    for item in attribute_string.split(";"):
        if "=" not in item:
            continue

        key, value = item.split("=", 1)
        attributes[key] = unquote(value)

    return attributes


def extract_species_name(gff3_path):
    """
    Extract species name from path like:
    data/raw/ensembl_genomes/bacteria/escherichia_coli_k12_mg1655/release-62/file.gff3.gz
    """
    return gff3_path.parts[-3]


def process_gff3(gff3_path):
    species = extract_species_name(gff3_path)
    output_file = OUTPUT_DIR / f"{species}_gene_lengths.tsv"
    rows_written = 0

    with gzip.open(gff3_path, "rt") as infile, open(output_file, "w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t")

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
            "description",
            "source_file"
        ])

        for line in infile:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 9:
                continue

            chromosome, source, feature, start, end, score, strand, phase, attributes_raw = fields

            if feature != "gene":
                continue

            attributes = parse_gff3_attributes(attributes_raw)

            genetype = attributes.get("biotype", "")

            if genetype != "protein_coding":
                continue

            gene_id = attributes.get("gene_id", "")

            if not gene_id:
                raw_id = attributes.get("ID", "")
                gene_id = raw_id.replace("gene:", "")

            gene_name = attributes.get("Name", "")
            description = attributes.get("description", "")

            start_int = int(start)
            end_int = int(end)
            gene_length_bp = end_int - start_int + 1

            writer.writerow([
                species,
                gene_id,
                gene_name,
                genetype,
                chromosome,
                start_int,
                end_int,
                gene_length_bp,
                strand,
                description,
                str(gff3_path.relative_to(PROJECT_ROOT))
            ])

            rows_written += 1

    print(f"{species}: wrote {rows_written} protein-coding genes to {output_file}")


def main():
    gff3_files = sorted(INPUT_BASE.rglob("*.gff3.gz"))

    if not gff3_files:
        raise FileNotFoundError(f"No .gff3.gz files found under {INPUT_BASE}")

    for gff3_path in gff3_files:
        process_gff3(gff3_path)


if __name__ == "__main__":
    main()