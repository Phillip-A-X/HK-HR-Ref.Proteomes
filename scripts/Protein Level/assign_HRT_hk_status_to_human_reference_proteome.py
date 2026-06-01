from pathlib import Path
import gzip
import sys

import pandas as pd


# Allows importing from lib/ when running this script from the project root
sys.path.append(str(Path(__file__).resolve().parents[1]))

from lib.constants import RAW_DATA_DIR, PROCESSED_DATA_DIR


# Input 1:
# UniProt Reference Proteome FASTA, extracted from UniProt release 2026_01
REFERENCE_PROTEOME_FASTA = (
    RAW_DATA_DIR / "uniprot" / "UP000005640_9606.fasta.gz"
)

# Input 2:
# UniProt ID mapping result for final HRT Atlas HK ENSG IDs
HRT_UNIPROT_MAPPING_FILE = (
    RAW_DATA_DIR / "uniprot" / "HRTATLAS HUMAN idmapping_2026_06_01.tsv"
)

# Output:
# Final human reference proteome table with HRT Atlas HK status
OUTPUT_FILE = (
    PROCESSED_DATA_DIR / "human_reference_proteome_with_hrt_atlas_hk_status.tsv"
)

# Optional QC output:
# HRT Atlas mapped UniProt accessions that were not found in the reference proteome FASTA
MAPPED_HK_ACCESSIONS_NOT_IN_FASTA_FILE = (
    PROCESSED_DATA_DIR / "hrt_atlas_hk_uniprot_accessions_not_in_reference_proteome.tsv"
)


def parse_uniprot_fasta(fasta_file: Path) -> pd.DataFrame:
    """
    Parse UniProt FASTA file.

    Expected UniProt header example:
    >sp|O60762|DPM1_HUMAN Dolichol-phosphate mannosyltransferase subunit 1 OS=Homo sapiens ...

    Extracted:
    - uniprot_accession: O60762
    - entry_name: DPM1_HUMAN
    - protein_length_aa: length of amino acid sequence
    """

    if not fasta_file.exists():
        raise FileNotFoundError(f"Reference proteome FASTA not found: {fasta_file}")

    records = []

    current_accession = None
    current_entry_name = None
    current_description = None
    current_sequence_parts = []

    with gzip.open(fasta_file, "rt") as handle:
        for line in handle:
            line = line.strip()

            if line.startswith(">"):
                # Save previous record
                if current_accession is not None:
                    sequence = "".join(current_sequence_parts)
                    records.append(
                        {
                            "uniprot_accession": current_accession,
                            "entry_name": current_entry_name,
                            "description": current_description,
                            "protein_length_aa": len(sequence),
                        }
                    )

                # Start new record
                header = line[1:]
                header_parts = header.split(" ", 1)
                id_part = header_parts[0]
                description = header_parts[1] if len(header_parts) > 1 else ""

                # UniProt ID part: sp|O60762|DPM1_HUMAN or tr|A0A...|...
                id_fields = id_part.split("|")

                if len(id_fields) >= 3:
                    accession = id_fields[1]
                    entry_name = id_fields[2]
                else:
                    raise ValueError(f"Unexpected FASTA header format: {line}")

                current_accession = accession
                current_entry_name = entry_name
                current_description = description
                current_sequence_parts = []

            else:
                current_sequence_parts.append(line)

        # Save last record
        if current_accession is not None:
            sequence = "".join(current_sequence_parts)
            records.append(
                {
                    "uniprot_accession": current_accession,
                    "entry_name": current_entry_name,
                    "description": current_description,
                    "protein_length_aa": len(sequence),
                }
            )

    return pd.DataFrame(records)


def main() -> None:
    PROCESSED_DATA_DIR.mkdir(parents=True, exist_ok=True)

    if not HRT_UNIPROT_MAPPING_FILE.exists():
        raise FileNotFoundError(f"HRT UniProt mapping file not found: {HRT_UNIPROT_MAPPING_FILE}")

    # Load HRT Atlas HK ENSG -> UniProt mapping table
    hrt_mapping = pd.read_csv(HRT_UNIPROT_MAPPING_FILE, sep="\t")

    required_mapping_columns = ["From", "Entry"]
    missing_columns = [col for col in required_mapping_columns if col not in hrt_mapping.columns]

    if missing_columns:
        raise ValueError(
            f"Missing columns in HRT mapping file: {missing_columns}\n"
            f"Available columns: {list(hrt_mapping.columns)}"
        )

    hrt_mapping["From"] = hrt_mapping["From"].astype(str).str.strip()
    hrt_mapping["Entry"] = hrt_mapping["Entry"].astype(str).str.strip()

    # Entry = UniProt accession
    hrt_hk_uniprot_accessions = set(
        hrt_mapping["Entry"]
        .dropna()
        .astype(str)
        .str.strip()
    )

    # Remove empty or invalid placeholder values if present
    hrt_hk_uniprot_accessions = {
        acc for acc in hrt_hk_uniprot_accessions
        if acc != "" and acc.lower() != "nan"
    }

    # Parse full human reference proteome FASTA
    reference_proteome = parse_uniprot_fasta(REFERENCE_PROTEOME_FASTA)

    reference_accessions = set(reference_proteome["uniprot_accession"])

    # Assign HRT Atlas HK status on protein level
    reference_proteome["hk_status"] = reference_proteome["uniprot_accession"].isin(
        hrt_hk_uniprot_accessions
    ).map(
        {
            True: "HK",
            False: "non-HK",
        }
    )

    # Optional trace column
    reference_proteome["hk_source"] = reference_proteome["hk_status"].map(
        {
            "HK": "HRT Atlas",
            "non-HK": "",
        }
    )

    # QC: mapped HK UniProt accessions that are not in the reference proteome FASTA
    mapped_not_in_fasta = sorted(hrt_hk_uniprot_accessions - reference_accessions)

    pd.Series(
        mapped_not_in_fasta,
        name="uniprot_accession_not_in_reference_proteome",
    ).to_csv(
        MAPPED_HK_ACCESSIONS_NOT_IN_FASTA_FILE,
        sep="\t",
        index=False,
    )

    reference_proteome.to_csv(
        OUTPUT_FILE,
        sep="\t",
        index=False,
    )

    print("HRT Atlas HK status assignment to human UniProt reference proteome")
    print()
    print(f"Reference proteome FASTA: {REFERENCE_PROTEOME_FASTA}")
    print(f"HRT UniProt mapping file: {HRT_UNIPROT_MAPPING_FILE}")
    print()
    print(f"Reference proteome proteins: {len(reference_proteome):,}")
    print(f"Unique HRT mapped UniProt accessions: {len(hrt_hk_uniprot_accessions):,}")
    print(f"HRT mapped accessions found in reference proteome: {len(hrt_hk_uniprot_accessions & reference_accessions):,}")
    print(f"HRT mapped accessions not found in reference proteome: {len(mapped_not_in_fasta):,}")
    print()
    print("Protein HK status counts:")
    print(reference_proteome["hk_status"].value_counts())
    print()
    print(f"Output written to: {OUTPUT_FILE}")
    print(f"Mapped HK accessions not in FASTA written to: {MAPPED_HK_ACCESSIONS_NOT_IN_FASTA_FILE}")


if __name__ == "__main__":
    main()
