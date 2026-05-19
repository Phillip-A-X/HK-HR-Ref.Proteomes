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

#Human HRTAtlas Hk Gene List
BIOMART_HRT_ATLAS_HUMAN_HK_TRANSCRIPTS_TO_ENSG_FILE = (
    BIOMART_DIR / "biomart_export_HRTAtlas_HumanHKGenes_ENSTmappedtoENSG_duplicatesallowed_180526.tsv"
)

HRT_ATLAS_HUMAN_HK_GENE_LENGTHS_FILE = (
    PROCESSED_DATA_DIR / "hrt_atlas_human_hk_gene_lengths.tsv"
)

HRT_ATLAS_HUMAN_HK_TRANSCRIPT_TO_GENE_MAPPING_QC_FILE = (
    PROCESSED_DATA_DIR / "hrt_atlas_human_hk_transcript_to_gene_mapping_qc.tsv"
)

HRT_ATLAS_HUMAN_HK_DUPLICATED_GENES_FILE = (
    PROCESSED_DATA_DIR / "hrt_atlas_human_hk_genes_with_multiple_transcripts.tsv"
)






# HRT Atlas Human HK genes: ENST to ENSG + gene coordinates - addon with protein coding only


BIOMART_HRT_ATLAS_HUMAN_HK_ENST_TO_ENSG_PROTEIN_CODING_DUPLICATES_ALLOWED_FILE = (
    BIOMART_DIR
    / "biomart_export_HRTAtlas_HumanHKGenes_ENSTmappedtoENSG_geneCoords_proteinCoding_duplicatesAllowed_200526.tsv"
)

BIOMART_HRT_ATLAS_HUMAN_HK_ENST_TO_ENSG_PROTEIN_CODING_UNIQUE_RESULTS_FILE = (
    BIOMART_DIR
    / "biomart_export_HRTAtlas_HumanHKGenes_ENSTmappedtoENSG_geneCoords_proteinCoding_uniqueResultsOnly_200526.tsv"
)

HRT_ATLAS_HUMAN_HK_GENE_LENGTHS_PROTEIN_CODING_FILE = (
    PROCESSED_DATA_DIR / "hrt_atlas_human_hk_gene_lengths_protein_coding.tsv"
)

HRT_ATLAS_HUMAN_HK_TRANSCRIPT_TO_GENE_MAPPING_QC_PROTEIN_CODING_FILE = (
    PROCESSED_DATA_DIR / "hrt_atlas_human_hk_transcript_to_gene_mapping_qc_protein_coding.tsv"
)

HRT_ATLAS_HUMAN_HK_GENES_WITH_MULTIPLE_TRANSCRIPTS_PROTEIN_CODING_FILE = (
    PROCESSED_DATA_DIR / "hrt_atlas_human_hk_genes_with_multiple_transcripts_protein_coding.tsv"
)


POLYX_IDENTICAL_AA = 8
POLYX_WINDOW_SIZE = 10
