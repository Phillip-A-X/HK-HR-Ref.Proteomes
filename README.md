# HK-HR-Ref.Proteomes
Computational pipeline for homorepeat analysis in housekeeping proteins using UniProt reference proteomes.

Protocol (not final form)
06.05.26 14:36 [powershell] Curl of UniProt Reference Proteome Release 2026_01, 28-Jan-2026 into "C:\uniprot_reference_proteomes"

07.05.26 13:51  [powershell] wsl --install 
                [ubuntu]     sudo mkdir -p /mnt/d
                             sudo mount -t drvfs D: /mnt/d
                             ls /mnt/d
07.05.26 21:36 [powershell] robocopy "C:\uniprot_reference_proteomes" "D:\data\compressed\ftp.uniprot.org\pub\databases\uniprot\release-2026_01\knowledgebase\reference_proteomes" "Reference_Proteomes_2026_01.tar.gz" /Z
#Linux/Ubuntu path: /mnt/d/data/compressed/ftp.uniprot.org/pub/databases/uniprot/release-2026_01/knowledgebase/reference_proteomes/
07.05.26 23:12 "https://www.tau.ac.il/~elieis/HKG/HK_genes.txt" was accessed and saved as "eisenberg_hk_genes.tsv" for Eisenberg list of human hk genelist containing gene symbols and refseqmrna

Ensembl Biomart
