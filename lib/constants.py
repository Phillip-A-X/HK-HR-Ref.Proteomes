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

EISENBERG_HK_GENES_FILE = RAW_DATA_DIR / "eisenberg_hk_genes.tsv"

BIOMART_EISENBERG_SWISSPROT_FILE = (
    BIOMART_DIR / "biomartexport_100526_humanuniprotid.tsv"
)

BIOMART_EISENBERG_GENE_LENGTH_FILE = (
    BIOMART_DIR / "biomartexport_120526_humanhkgenelength.tsv"
)



POLYX_IDENTICAL_AA = 8
POLYX_WINDOW_SIZE = 10
