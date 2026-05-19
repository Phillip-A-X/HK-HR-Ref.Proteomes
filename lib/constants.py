#Constant Values
#Path

from pathlib import Path


UNIPROT_RELEASE = "2026_01"

MOUNT_POINT = Path("/mnt/d/data/compressed")

REFERENCE_PROTEOME_DIR = (
    MOUNT_POINT /"ftp.uniprot.org"/"pub"/"databases"/"uniprot"/"release-2026_01"/"knowledgebase"/"reference_proteomes"
)




##Local directories

DATA_DIR = Path("data")
RAW_DATA_DIR = DATA_DIR / "raw"
PROCESSED_DATA_DIR = DATA_DIR / "processed"

BIOMART_DIR = RAW_DATA_DIR / "biomart"

RESULTS_DIR = Path("results")
TABLES_DIR = RESULTS_DIR / "tables"
FIGURES_DIR = RESULTS_DIR / "figures"




#Human

HUMAN_PROTEOME_ID = "UP000005640_9606"
HUMAN_FASTA_FILE = REFERENCE_PROTEOME_DIR / f"{HUMAN_PROTEOME_ID}.fasta.gz"


#Human Eisenberg HK Gene List
EISENBERG_HK_GENES_FILE = RAW_DATA_DIR / "eisenberg_hk_genes.tsv"

BIOMART_EISENBERG_SWISSPROT_FILE = (
    BIOMART_DIR / "biomartexport_100526_humanuniprotid.tsv"
)

BIOMART_EISENBERG_GENE_LENGTH_FILE = (
    BIOMART_DIR / "biomartexport_120526_humanhkgenelength.tsv"
)


## HRT Atlas Human HK genes

# Original HRT Atlas human HK gene/transcript list.
# No header in file.
# Columns:
# 1. Ensembl transcript ID (ENST...)
# 2. Gene name
# 3. RefSeq ID
# 4. CCDS ID
HRT_ATLAS_HUMAN_HK_GENES_FILE = (
    RAW_DATA_DIR / "Housekeeping_GenesHuman_HRTatlas.tsv"
)


# BioMart export:
# Dataset: Human genes (GRCh38.p14), Ensembl 115
# Filter: HRT Atlas ENST IDs
# Attributes:
# - Transcript stable ID
# - Gene stable ID
# - Gene start (bp)
# - Gene end (bp)
# - Chromosome/scaffold name
BIOMART_HRT_ATLAS_HUMAN_HK_ENST_TO_ENSG_FILE = (
    BIOMART_DIR / "biomart_export_HRTAtlas_HumanHKGenes_ENSTmappedtoENSG_duplicatesallowed_180526.tsv"
)


# Main processed output:
# one row per unique Ensembl gene ID, with calculated gene length
HRT_ATLAS_HUMAN_HK_GENE_LENGTHS_FILE = (
    PROCESSED_DATA_DIR / "hrt_atlas_human_hk_gene_lengths.tsv"
)


# QC output:
# number of HRT Atlas transcripts mapped to each Ensembl gene
HRT_ATLAS_TRANSCRIPTS_PER_GENE_QC_FILE = (
    PROCESSED_DATA_DIR / "hrt_atlas_transcripts_per_gene_mapping_qc.tsv"
)


# QC output:
# only genes with more than one HRT Atlas transcript
HRT_ATLAS_HUMAN_HK_DUPLICATED_GENES_FILE = (
    PROCESSED_DATA_DIR / "hrt_atlas_human_hk_genes_with_multiple_transcripts.tsv"
)


# QC output:
# ENST IDs from the original HRT Atlas file that were not returned by BioMart
HRT_ATLAS_UNMAPPED_ENST_FILE = (
    PROCESSED_DATA_DIR / "hrt_atlas_unmapped_enst_from_biomart.tsv"
)


# QC output:
# summary numbers for ENST-to-ENSG mapping
HRT_ATLAS_ENST_TO_ENSG_MAPPING_SUMMARY_FILE = (
    PROCESSED_DATA_DIR / "hrt_atlas_enst_to_ensg_mapping_summary.tsv"
)


POLYX_IDENTICAL_AA = 8
POLYX_WINDOW_SIZE = 10
