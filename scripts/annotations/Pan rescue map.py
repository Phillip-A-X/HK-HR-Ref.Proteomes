#!/usr/bin/env python3

from pathlib import Path
import csv
import gzip
import json
import re
import time
from collections import Counter, defaultdict
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen


# -------------------------------------------------------------------------
# Configuration
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
    / "pan_joshi_full115_rescued"
)

FINAL_OUTPUT_FILE = (
    OUTPUT_DIR
    / "pan_troglodytes_full115_protein_coding_with_joshi_hk_status_rescued.tsv"
)

MAPPING_AUDIT_FILE = (
    OUTPUT_DIR
    / "pan_joshi_2022_id_mapping_audit.tsv"
)

UNRESOLVED_FILE = (
    OUTPUT_DIR
    / "pan_joshi_2022_unresolved_or_ambiguous_ids.tsv"
)

SUMMARY_FILE = (
    OUTPUT_DIR
    / "pan_joshi_full115_rescued_mapping_summary.tsv"
)

API_CACHE_FILE = (
    OUTPUT_DIR
    / "ensembl_archive_api_cache.json"
)

ENSEMBL_REST_SERVER = "https://rest.ensembl.org"

REQUEST_DELAY_SECONDS = 0.12
MAX_RETRIES = 4
REQUEST_TIMEOUT_SECONDS = 30


# -------------------------------------------------------------------------
# General helper functions
# -------------------------------------------------------------------------

def normalise_gene_id(gene_id):
    """
    Remove surrounding whitespace and optional version suffix.

    Example:
    ENSPTRG00000022029.3 -> ENSPTRG00000022029
    """

    if gene_id is None:
        return ""

    gene_id = str(gene_id).strip()

    if "." in gene_id:
        gene_id = gene_id.split(".", 1)[0]

    return gene_id


def parse_gtf_attributes(attribute_string):
    """
    Parse quoted attributes from the ninth GTF column.
    """

    attributes = {}

    for match in re.finditer(r'(\S+) "([^"]+)"', attribute_string):
        key = match.group(1)
        value = match.group(2)
        attributes[key] = value

    return attributes


def load_hk_ids(hk_list_file):
    """
    Read the Joshi HK list.

    The header 'Genes', empty lines and invalid entries are ignored.
    Duplicate IDs are removed.
    """

    hk_ids = set()
    invalid_entries = []

    with open(
        hk_list_file,
        "r",
        encoding="utf-8-sig"
    ) as infile:

        for line_number, line in enumerate(infile, start=1):
            value = line.strip()

            if not value:
                continue

            if value.lower() == "genes":
                continue

            gene_id = normalise_gene_id(value)

            if not gene_id.startswith("ENSPTRG"):
                invalid_entries.append(
                    (line_number, value)
                )
                continue

            hk_ids.add(gene_id)

    if invalid_entries:
        print("\nInvalid HK-list entries ignored:")

        for line_number, value in invalid_entries[:20]:
            print(f"  line {line_number}: {value}")

    return hk_ids


# -------------------------------------------------------------------------
# Parse Ensembl 115 GTF
# -------------------------------------------------------------------------

def load_gtf_genes(gtf_file):
    """
    Parse every gene row from the complete Ensembl 115 GTF.

    Returns:
    gene_records:
        dictionary keyed by normalised gene ID

    protein_coding_ids:
        set containing exact protein_coding gene IDs
    """

    gene_records = {}
    protein_coding_ids = set()

    biotype_counts = Counter()
    sequence_counts = Counter()

    gene_rows = 0
    duplicate_gene_ids = set()

    with gzip.open(
        gtf_file,
        "rt",
        encoding="utf-8"
    ) as infile:

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

            if feature != "gene":
                continue

            gene_rows += 1

            attributes = parse_gtf_attributes(
                attributes_raw
            )

            gene_id = normalise_gene_id(
                attributes.get("gene_id", "")
            )

            if not gene_id:
                continue

            gene_name = attributes.get(
                "gene_name",
                ""
            )

            gene_biotype = attributes.get(
                "gene_biotype",
                attributes.get("gene_type", "")
            )

            start_int = int(start)
            end_int = int(end)

            gene_length_bp = end_int - start_int + 1

            if gene_length_bp <= 0:
                raise ValueError(
                    f"Invalid coordinates for {gene_id}: "
                    f"{start_int}-{end_int}"
                )

            if gene_id in gene_records:
                duplicate_gene_ids.add(gene_id)

            gene_records[gene_id] = {
                "species": "pan_troglodytes",
                "gene_id": gene_id,
                "gene_name": gene_name,
                "genetype": gene_biotype,
                "chromosome": seqname,
                "start": start_int,
                "end": end_int,
                "gene_length_bp": gene_length_bp,
                "strand": strand,
                "source_file": str(
                    gtf_file.relative_to(PROJECT_ROOT)
                )
            }

            biotype_counts[gene_biotype or "MISSING"] += 1
            sequence_counts[seqname] += 1

            if gene_biotype == "protein_coding":
                protein_coding_ids.add(gene_id)

    if duplicate_gene_ids:
        raise RuntimeError(
            f"{len(duplicate_gene_ids)} duplicated gene IDs "
            f"detected in the GTF."
        )

    return {
        "gene_records": gene_records,
        "protein_coding_ids": protein_coding_ids,
        "gene_rows": gene_rows,
        "biotype_counts": biotype_counts,
        "sequence_counts": sequence_counts
    }


# -------------------------------------------------------------------------
# Ensembl Archive API
# -------------------------------------------------------------------------

def load_api_cache(cache_file):
    """
    Load previous Archive API responses so interrupted runs can continue.
    """

    if not cache_file.exists():
        return {}

    with open(
        cache_file,
        "r",
        encoding="utf-8"
    ) as infile:
        return json.load(infile)


def save_api_cache(cache, cache_file):
    """
    Save Archive API responses after each request.
    """

    temporary_file = cache_file.with_suffix(".tmp")

    with open(
        temporary_file,
        "w",
        encoding="utf-8"
    ) as outfile:
        json.dump(
            cache,
            outfile,
            indent=2,
            sort_keys=True
        )

    temporary_file.replace(cache_file)


def query_archive_api(gene_id):
    """
    Query Ensembl GET archive/id/:id.

    Returns a dictionary containing either:
    {
        "request_status": "success",
        "data": {...}
    }

    or an error record.
    """

    encoded_id = quote(gene_id)

    url = (
        f"{ENSEMBL_REST_SERVER}"
        f"/archive/id/{encoded_id}"
        f"?content-type=application/json"
    )

    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json",
        "User-Agent": "HK-HR-RefProteomes/1.0"
    }

    for attempt in range(1, MAX_RETRIES + 1):
        request = Request(
            url,
            headers=headers,
            method="GET"
        )

        try:
            with urlopen(
                request,
                timeout=REQUEST_TIMEOUT_SECONDS
            ) as response:

                response_text = response.read().decode(
                    "utf-8"
                )

                data = json.loads(response_text)

                return {
                    "request_status": "success",
                    "http_status": response.status,
                    "data": data
                }

        except HTTPError as error:
            retry_after = error.headers.get(
                "Retry-After"
            )

            if error.code == 429:
                wait_seconds = (
                    float(retry_after)
                    if retry_after
                    else attempt * 2
                )

                print(
                    f"Rate limit for {gene_id}; "
                    f"waiting {wait_seconds} seconds."
                )

                time.sleep(wait_seconds)
                continue

            if error.code in {400, 404}:
                return {
                    "request_status": "not_found",
                    "http_status": error.code,
                    "error": str(error)
                }

            if attempt == MAX_RETRIES:
                return {
                    "request_status": "http_error",
                    "http_status": error.code,
                    "error": str(error)
                }

            time.sleep(attempt * 2)

        except (
            URLError,
            TimeoutError,
            json.JSONDecodeError
        ) as error:

            if attempt == MAX_RETRIES:
                return {
                    "request_status": "connection_error",
                    "http_status": "",
                    "error": str(error)
                }

            time.sleep(attempt * 2)

    return {
        "request_status": "unknown_error",
        "http_status": "",
        "error": "Maximum attempts exceeded"
    }


def extract_candidate_ids(api_data):
    """
    Extract possible current/replacement Ensembl gene IDs from an
    Archive API response.

    Supported fields include:
    latest
    id
    possible_replacement

    possible_replacement may contain strings or dictionaries.
    """

    candidate_ids = set()

    latest = api_data.get("latest")

    if latest:
        latest_id = normalise_gene_id(latest)

        if latest_id.startswith("ENSPTRG"):
            candidate_ids.add(latest_id)

    returned_id = api_data.get("id")

    if returned_id:
        returned_id = normalise_gene_id(returned_id)

        if returned_id.startswith("ENSPTRG"):
            candidate_ids.add(returned_id)

    possible_replacements = api_data.get(
        "possible_replacement",
        []
    )

    if isinstance(possible_replacements, dict):
        possible_replacements = [
            possible_replacements
        ]

    if isinstance(possible_replacements, str):
        possible_replacements = [
            possible_replacements
        ]

    if possible_replacements is None:
        possible_replacements = []

    for replacement in possible_replacements:
        candidate = ""

        if isinstance(replacement, str):
            candidate = replacement

        elif isinstance(replacement, dict):
            for key in (
                "id",
                "stable_id",
                "latest",
                "replacement",
                "new_id"
            ):
                if replacement.get(key):
                    candidate = replacement[key]
                    break

        candidate = normalise_gene_id(candidate)

        if candidate.startswith("ENSPTRG"):
            candidate_ids.add(candidate)

    return sorted(candidate_ids)


# -------------------------------------------------------------------------
# Resolve each Joshi HK ID
# -------------------------------------------------------------------------

def resolve_hk_ids(
    hk_ids,
    gene_records,
    protein_coding_ids,
    api_cache
):
    """
    Resolve direct and historic HK IDs.

    Only a single unambiguous protein-coding target is rescued
    automatically.
    """

    mapping_records = []
    resolved_target_to_original_ids = defaultdict(set)

    direct_count = 0
    rescued_count = 0

    hk_ids_sorted = sorted(hk_ids)

    for index, original_id in enumerate(
        hk_ids_sorted,
        start=1
    ):
        # Direct protein-coding match
        if original_id in protein_coding_ids:
            resolved_target_to_original_ids[
                original_id
            ].add(original_id)

            mapping_records.append({
                "original_hk_gene_id": original_id,
                "resolved_gene_id": original_id,
                "mapping_method": "direct",
                "mapping_status": "exact_protein_coding",
                "candidate_ids": original_id,
                "candidate_count": 1,
                "local_gene_biotype": "protein_coding",
                "api_release": "",
                "api_latest": "",
                "api_is_current": "",
                "api_request_status": "not_required",
                "notes": ""
            })

            direct_count += 1
            continue

        print(
            f"[{index}/{len(hk_ids_sorted)}] "
            f"Checking archive history for {original_id}"
        )

        if original_id not in api_cache:
            api_cache[original_id] = query_archive_api(
                original_id
            )

            save_api_cache(
                api_cache,
                API_CACHE_FILE
            )

            time.sleep(
                REQUEST_DELAY_SECONDS
            )

        cached_result = api_cache[original_id]

        request_status = cached_result.get(
            "request_status",
            "unknown"
        )

        if request_status != "success":
            mapping_records.append({
                "original_hk_gene_id": original_id,
                "resolved_gene_id": "",
                "mapping_method": "ensembl_archive",
                "mapping_status": "api_not_found_or_error",
                "candidate_ids": "",
                "candidate_count": 0,
                "local_gene_biotype": "",
                "api_release": "",
                "api_latest": "",
                "api_is_current": "",
                "api_request_status": request_status,
                "notes": cached_result.get(
                    "error",
                    ""
                )
            })

            continue

        api_data = cached_result.get(
            "data",
            {}
        )

        candidate_ids = extract_candidate_ids(
            api_data
        )

        candidates_present_locally = [
            candidate
            for candidate in candidate_ids
            if candidate in gene_records
        ]

        protein_coding_candidates = [
            candidate
            for candidate in candidates_present_locally
            if candidate in protein_coding_ids
        ]

        non_protein_coding_candidates = [
            candidate
            for candidate in candidates_present_locally
            if candidate not in protein_coding_ids
        ]

        api_release = api_data.get(
            "release",
            ""
        )

        api_latest = api_data.get(
            "latest",
            ""
        )

        api_is_current = api_data.get(
            "is_current",
            ""
        )

        if len(protein_coding_candidates) == 1:
            resolved_id = protein_coding_candidates[0]

            resolved_target_to_original_ids[
                resolved_id
            ].add(original_id)

            mapping_records.append({
                "original_hk_gene_id": original_id,
                "resolved_gene_id": resolved_id,
                "mapping_method": "ensembl_archive",
                "mapping_status": "rescued_one_to_one",
                "candidate_ids": ";".join(
                    candidate_ids
                ),
                "candidate_count": len(
                    candidate_ids
                ),
                "local_gene_biotype": "protein_coding",
                "api_release": api_release,
                "api_latest": api_latest,
                "api_is_current": api_is_current,
                "api_request_status": request_status,
                "notes": ""
            })

            rescued_count += 1

        elif len(protein_coding_candidates) > 1:
            mapping_records.append({
                "original_hk_gene_id": original_id,
                "resolved_gene_id": "",
                "mapping_method": "ensembl_archive",
                "mapping_status": "ambiguous_multiple_protein_coding",
                "candidate_ids": ";".join(
                    protein_coding_candidates
                ),
                "candidate_count": len(
                    protein_coding_candidates
                ),
                "local_gene_biotype": "",
                "api_release": api_release,
                "api_latest": api_latest,
                "api_is_current": api_is_current,
                "api_request_status": request_status,
                "notes": (
                    "Multiple local protein-coding candidates; "
                    "not rescued automatically"
                )
            })

        elif non_protein_coding_candidates:
            candidate_biotypes = sorted({
                gene_records[candidate]["genetype"]
                for candidate
                in non_protein_coding_candidates
            })

            mapping_records.append({
                "original_hk_gene_id": original_id,
                "resolved_gene_id": "",
                "mapping_method": "ensembl_archive",
                "mapping_status": "replacement_not_protein_coding",
                "candidate_ids": ";".join(
                    non_protein_coding_candidates
                ),
                "candidate_count": len(
                    non_protein_coding_candidates
                ),
                "local_gene_biotype": ";".join(
                    candidate_biotypes
                ),
                "api_release": api_release,
                "api_latest": api_latest,
                "api_is_current": api_is_current,
                "api_request_status": request_status,
                "notes": (
                    "Candidate exists in local Ensembl 115 GTF "
                    "but is not protein_coding"
                )
            })

        elif candidate_ids:
            mapping_records.append({
                "original_hk_gene_id": original_id,
                "resolved_gene_id": "",
                "mapping_method": "ensembl_archive",
                "mapping_status": "candidate_not_in_local_release115_gtf",
                "candidate_ids": ";".join(
                    candidate_ids
                ),
                "candidate_count": len(
                    candidate_ids
                ),
                "local_gene_biotype": "",
                "api_release": api_release,
                "api_latest": api_latest,
                "api_is_current": api_is_current,
                "api_request_status": request_status,
                "notes": (
                    "Archive candidate was not found in the "
                    "local Ensembl 115 GTF"
                )
            })

        else:
            mapping_records.append({
                "original_hk_gene_id": original_id,
                "resolved_gene_id": "",
                "mapping_method": "ensembl_archive",
                "mapping_status": "no_replacement_candidate",
                "candidate_ids": "",
                "candidate_count": 0,
                "local_gene_biotype": "",
                "api_release": api_release,
                "api_latest": api_latest,
                "api_is_current": api_is_current,
                "api_request_status": request_status,
                "notes": (
                    "No usable ENSPTRG replacement returned"
                )
            })

    return {
        "mapping_records": mapping_records,
        "resolved_target_to_original_ids":
            resolved_target_to_original_ids,
        "direct_count": direct_count,
        "rescued_count": rescued_count
    }


# -------------------------------------------------------------------------
# Write output files
# -------------------------------------------------------------------------

def write_mapping_audit(mapping_records):
    fieldnames = [
        "original_hk_gene_id",
        "resolved_gene_id",
        "mapping_method",
        "mapping_status",
        "candidate_ids",
        "candidate_count",
        "local_gene_biotype",
        "api_release",
        "api_latest",
        "api_is_current",
        "api_request_status",
        "notes"
    ]

    with open(
        MAPPING_AUDIT_FILE,
        "w",
        newline="",
        encoding="utf-8"
    ) as outfile:

        writer = csv.DictWriter(
            outfile,
            fieldnames=fieldnames,
            delimiter="\t"
        )

        writer.writeheader()
        writer.writerows(mapping_records)


def write_unresolved(mapping_records):
    unresolved_records = [
        record
        for record in mapping_records
        if not record["resolved_gene_id"]
    ]

    fieldnames = [
        "original_hk_gene_id",
        "mapping_status",
        "candidate_ids",
        "candidate_count",
        "local_gene_biotype",
        "api_release",
        "api_latest",
        "api_is_current",
        "api_request_status",
        "notes"
    ]

    with open(
        UNRESOLVED_FILE,
        "w",
        newline="",
        encoding="utf-8"
    ) as outfile:

        writer = csv.DictWriter(
            outfile,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore"
        )

        writer.writeheader()
        writer.writerows(unresolved_records)


def write_final_gene_table(
    gene_records,
    protein_coding_ids,
    resolved_target_to_original_ids,
    mapping_records
):
    mapping_by_target = defaultdict(list)

    for record in mapping_records:
        resolved_gene_id = record["resolved_gene_id"]

        if resolved_gene_id:
            mapping_by_target[resolved_gene_id].append(
                record
            )

    fieldnames = [
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
        "mapping_method",
        "mapping_status",
        "original_hk_gene_ids",
        "source_file"
    ]

    with open(
        FINAL_OUTPUT_FILE,
        "w",
        newline="",
        encoding="utf-8"
    ) as outfile:

        writer = csv.DictWriter(
            outfile,
            fieldnames=fieldnames,
            delimiter="\t"
        )

        writer.writeheader()

        for gene_id in sorted(
            protein_coding_ids
        ):
            gene_record = gene_records[gene_id]

            associated_mappings = mapping_by_target.get(
                gene_id,
                []
            )

            if associated_mappings:
                hk_status = "HK"
                hk_source = "Joshi et al. 2022"

                original_ids = sorted(
                    resolved_target_to_original_ids[
                        gene_id
                    ]
                )

                methods = sorted({
                    record["mapping_method"]
                    for record in associated_mappings
                })

                statuses = sorted({
                    record["mapping_status"]
                    for record in associated_mappings
                })

                mapping_method = ";".join(methods)
                mapping_status = ";".join(statuses)
                original_hk_gene_ids = ";".join(
                    original_ids
                )

            else:
                hk_status = "non-HK"
                hk_source = "Joshi et al. 2022"
                mapping_method = "not_applicable"
                mapping_status = "not_in_hk_list"
                original_hk_gene_ids = ""

            output_record = {
                **gene_record,
                "hk_status": hk_status,
                "hk_source": hk_source,
                "mapping_method": mapping_method,
                "mapping_status": mapping_status,
                "original_hk_gene_ids":
                    original_hk_gene_ids
            }

            writer.writerow(output_record)


def write_summary(
    hk_ids,
    gene_records,
    protein_coding_ids,
    mapping_records,
    resolved_target_to_original_ids
):
    status_counts = Counter(
        record["mapping_status"]
        for record in mapping_records
    )

    direct_original_ids = {
        record["original_hk_gene_id"]
        for record in mapping_records
        if record["mapping_method"] == "direct"
        and record["resolved_gene_id"]
    }

    rescued_original_ids = {
        record["original_hk_gene_id"]
        for record in mapping_records
        if record["mapping_method"] == "ensembl_archive"
        and record["resolved_gene_id"]
    }

    unresolved_original_ids = {
        record["original_hk_gene_id"]
        for record in mapping_records
        if not record["resolved_gene_id"]
    }

    final_hk_gene_ids = set(
        resolved_target_to_original_ids
    )

    direct_mapping_rate = (
        len(direct_original_ids)
        / len(hk_ids)
        * 100
    )

    rescued_mapping_rate = (
        (
            len(direct_original_ids)
            + len(rescued_original_ids)
        )
        / len(hk_ids)
        * 100
    )

    summary_rows = [
        (
            "all_gene_rows_in_full_gtf",
            len(gene_records)
        ),
        (
            "protein_coding_gene_ids_in_full_gtf",
            len(protein_coding_ids)
        ),
        (
            "unique_hk_ids_in_input_list",
            len(hk_ids)
        ),
        (
            "directly_mapped_original_hk_ids",
            len(direct_original_ids)
        ),
        (
            "archive_rescued_original_hk_ids",
            len(rescued_original_ids)
        ),
        (
            "unresolved_original_hk_ids",
            len(unresolved_original_ids)
        ),
        (
            "direct_mapping_rate_percent",
            round(direct_mapping_rate, 4)
        ),
        (
            "mapping_rate_after_rescue_percent",
            round(rescued_mapping_rate, 4)
        ),
        (
            "unique_final_hk_target_genes",
            len(final_hk_gene_ids)
        ),
        (
            "final_non_hk_genes",
            len(protein_coding_ids)
            - len(final_hk_gene_ids)
        )
    ]

    for status, count in sorted(
        status_counts.items()
    ):
        summary_rows.append(
            (
                f"mapping_status_{status}",
                count
            )
        )

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

        writer.writerows(summary_rows)

    return {
        "direct_original_ids":
            direct_original_ids,
        "rescued_original_ids":
            rescued_original_ids,
        "unresolved_original_ids":
            unresolved_original_ids,
        "final_hk_gene_ids":
            final_hk_gene_ids,
        "direct_mapping_rate":
            direct_mapping_rate,
        "rescued_mapping_rate":
            rescued_mapping_rate,
        "status_counts":
            status_counts
    }


# -------------------------------------------------------------------------
# Main
# -------------------------------------------------------------------------

def main():
    if not GTF_FILE.exists():
        raise FileNotFoundError(
            f"GTF file not found: {GTF_FILE}"
        )

    if not HK_LIST_FILE.exists():
        raise FileNotFoundError(
            f"HK list not found: {HK_LIST_FILE}"
        )

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    print("Loading Joshi 2022 HK list...")

    hk_ids = load_hk_ids(
        HK_LIST_FILE
    )

    print(
        f"Unique HK IDs loaded: {len(hk_ids)}"
    )

    print("\nParsing complete Ensembl 115 GTF...")

    gtf_result = load_gtf_genes(
        GTF_FILE
    )

    gene_records = gtf_result[
        "gene_records"
    ]

    protein_coding_ids = gtf_result[
        "protein_coding_ids"
    ]

    print(
        f"All gene IDs: {len(gene_records)}"
    )

    print(
        f"Protein-coding gene IDs: "
        f"{len(protein_coding_ids)}"
    )

    api_cache = load_api_cache(
        API_CACHE_FILE
    )

    print(
        f"Cached API responses: "
        f"{len(api_cache)}"
    )

    print("\nResolving HK IDs...")

    resolution_result = resolve_hk_ids(
        hk_ids=hk_ids,
        gene_records=gene_records,
        protein_coding_ids=protein_coding_ids,
        api_cache=api_cache
    )

    mapping_records = resolution_result[
        "mapping_records"
    ]

    resolved_target_to_original_ids = (
        resolution_result[
            "resolved_target_to_original_ids"
        ]
    )

    write_mapping_audit(
        mapping_records
    )

    write_unresolved(
        mapping_records
    )

    write_final_gene_table(
        gene_records=gene_records,
        protein_coding_ids=protein_coding_ids,
        resolved_target_to_original_ids=
            resolved_target_to_original_ids,
        mapping_records=mapping_records
    )

    summary = write_summary(
        hk_ids=hk_ids,
        gene_records=gene_records,
        protein_coding_ids=protein_coding_ids,
        mapping_records=mapping_records,
        resolved_target_to_original_ids=
            resolved_target_to_original_ids
    )

    print("\nMapping summary")

    print(
        "Input HK IDs:",
        len(hk_ids)
    )

    print(
        "Directly mapped HK IDs:",
        len(
            summary["direct_original_ids"]
        )
    )

    print(
        "Archive-rescued HK IDs:",
        len(
            summary["rescued_original_ids"]
        )
    )

    print(
        "Unresolved HK IDs:",
        len(
            summary["unresolved_original_ids"]
        )
    )

    print(
        "Direct mapping rate:",
        f"{summary['direct_mapping_rate']:.2f}%"
    )

    print(
        "Mapping rate after rescue:",
        f"{summary['rescued_mapping_rate']:.2f}%"
    )

    print(
        "Unique final HK target genes:",
        len(
            summary["final_hk_gene_ids"]
        )
    )

    print("\nMapping status counts")

    for status, count in sorted(
        summary["status_counts"].items()
    ):
        print(
            f"{status}\t{count}"
        )

    print("\nOutput files")

    print(FINAL_OUTPUT_FILE)
    print(MAPPING_AUDIT_FILE)
    print(UNRESOLVED_FILE)
    print(SUMMARY_FILE)
    print(API_CACHE_FILE)


if __name__ == "__main__":
    main()