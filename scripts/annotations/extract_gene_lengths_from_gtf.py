#!/usr/bin/env python3

from pathlib import Path
import gzip
import csv
import re


PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT_BASES = [
    PROJECT_ROOT / "data" / "raw" / "ensembl",
    PROJECT_ROOT / "data" / "raw" / "ensembl_genomes"
]
OUTPUT_DIR = PROJECT_ROOT / "data" / "processed" / "gene_lengths"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def parse_gtf_attributes(attribute_string):
    """
    Parse the attributes column of a GTF file.
    Example:
    gene_id "ENSG..."; gene_name "BRCA2"; gene_biotype "protein_coding";
    """
    attributes = {}

    for match in re.finditer(r'(\S+) "([^"]+)"', attribute_string):
        key = match.group(1)
        value = match.group(2)
        attributes[key] = value

    return attributes


def extract_species_name(gtf_path):
    """
    Extract species name from path like:
    data/raw/ensembl/gorilla_gorilla/release-115/file.gtf.gz
    """
    return gtf_path.parts[-3]


def process_gtf(gtf_path):
    species = extract_species_name(gtf_path)
    output_file = OUTPUT_DIR / f"{species}_gene_lengths.tsv"
    rows_written = 0

    with gzip.open(gtf_path, "rt") as infile, open(output_file, "w", newline="") as outfile:
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
            "source_file"
        ])

        for line in infile:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 9:
                continue

            seqname, source, feature, start, end, score, strand, frame, attributes_raw = fields

            if feature != "gene":
                continue

            attributes = parse_gtf_attributes(attributes_raw)

            gene_biotype = attributes.get("gene_biotype", attributes.get("gene_type", ""))

            if gene_biotype != "protein_coding":
                continue

            gene_id = attributes.get("gene_id", "")
            gene_name = attributes.get("gene_name", "")

            start_int = int(start)
            end_int = int(end)
            gene_length_bp = end_int - start_int + 1

            writer.writerow([
                species,
                gene_id,
                gene_name,
                gene_biotype,
                seqname,
                start_int,
                end_int,
                gene_length_bp,
                strand,
                str(gtf_path.relative_to(PROJECT_ROOT))
            ])

            rows_written += 1

    print(f"{species}: wrote {rows_written} protein-coding genes to {output_file}")


def main():
    gtf_files = []

    for input_base in INPUT_BASES:
        gtf_files.extend(input_base.rglob("*.gtf.gz"))

    gtf_files = sorted(gtf_files)

    if not gtf_files:
        raise FileNotFoundError("No .gtf.gz files found under input directories")

    for gtf_path in gtf_files:
        process_gtf(gtf_path)


if __name__ == "__main__":
    main()