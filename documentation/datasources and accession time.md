#Datasources

##biomartexport_100526_humanuniprotid.tsv
- https://www.ensembl.org/biomart/martview/edfc4b2acd0521d29cc397c4cd58aaab
- 10.05.26
- Biomart Export
- Dataset: Human genes (GRCh38.p14) 3971 / 86369 Genes found
- Filters HGNC symbol(s) [e.g. A1BG]: [ID-list specified] - HGNC symbols of eisenberg_hk_genes.tsv accessed 08.05.26 Eisenberg list contains 3804 rows, probably due to certain genes not being located on numerated chromosomes but alternative scaffolds or haplotype regions. Removing these could bring the number closer together, will try [12.05.26 14:52]
- Attributes HGNC symbol UniProtKB/Swiss-Prot ID Gene stable ID Transcript stable ID Protein stable ID

##biomartexport_120526_humanhkgenelength.tsv
- https://www.ensembl.org/biomart/martview/edfc4b2acd0521d29cc397c4cd58aaab
- 12.05.26
- Biomart Export
- Dataset: Human genes (GRCh38.p14) 3971 / 86369 Genes found
- Filters HGNC symbol(s) [e.g. A1BG]: [ID-list specified] - HGNC symbols of eisenberg_hk_genes.tsv accessed 08.05.26 - Note: Eisenberg list contains 3804 rows, probably due to certain genes not being located on numerated chromosomes but alternative scaffolds or haplotype regions. Thus i added Chromosome/scaffold name into the attributes. Removing these could bring the number closer together, will try [12.05.26 14:52]
- Attributes Gene stable ID Gene start (bp) Gene end (bp) Chromosome/scaffold name


##eisenberg_hk_genes.tsv
- https://www.tau.ac.il/~elieis/HKG/HK_genes.txt
- 08.05.26
- Source: Eisenberg & Levanon housekeeping gene list  
Content: Human housekeeping genes with RefSeq mRNA accessions  
Columns:
- gene_symbol: HGNC gene symbol
- refseq_mrna: RefSeq mRNA accession, e.g. NM_015665
