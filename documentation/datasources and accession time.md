#Datasources

##biomartexport_100526_humanuniprotid.tsv
- https://www.ensembl.org/biomart/martview/edfc4b2acd0521d29cc397c4cd58aaab
- 10.05.26
- Biomart Export
- Dataset: Human genes (GRCh38.p14) 3971 / 86369 Genes found
- Filters HGNC symbol(s) [e.g. A1BG]: [ID-list specified] - HGNC symbols of eisenberg_hk_genes.tsv accessed 08.05.26 Eisenberg list contains 3804 rows, probably due to certain genes not being located on numerated chromosomes but alternative scaffolds or haplotype regions. Removing these could bring the number closer together, will try [12.05.26 14:52]
- Attributes HGNC symbol UniProtKB/Swiss-Prot ID Gene stable ID Transcript stable ID Protein stable ID

##biomartexport_120526_humanhkgenelength.tsv
- https://www.ensembl.org/biomart/martview/7a34f955b4379c0b815f0eb97152867d
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

##hrtatlas_hk_genes.tsv
- https://housekeeping.unicamp.br/?download
- 12.05.26
- Source: https://pubmed.ncbi.nlm.nih.gov/32663312/
- Content: Human housekeeping genes
- Columns: Ensembl (Transcript IDs),	Gene.name,	Refseq,	CCDS.ID

##biomart_export_HRTAtlas_HumanHKGenes_ENSTmappedtoENSG_duplicatesallowed_180526.tsv
- https://www.ensembl.org/biomart/martview/7dba1dc5838c33f434b6891cabc8bf58
- 18.05.26
- Biomart Export
- Dataset: Human genes (GRCh38.p14) Ensembl 115
- Filter: Transcript stable ID(s) [e.g. ENST00000000233]: [ID-list specified] HRT Atlas
- Attributes Transcript stable ID, Gene stable ID, Gene start (bp), Gene end (bp), Chromosome/scaffold name

##biomart_export_HRTAtlas_35unmappedENST_HGNCsymbol_rescueMapping_toENSG_geneCoords_200526.tsv
- 20.05.26
- Biomart Export
- Dataset: Human genes (GRCh38.p14) Ensembl 115
- Filter: HGNC symbol(s) [e.g. A1BG]: [ID-list specified] - unmapped HK genes from HRTAtlas
- Attributes: Gene Stable ID, Gene Start, Gene end, Chromosome/scaffold name

##biomart_export_HumanAllGenes_220526.tsv
- https://www.ensembl.org/biomart/martview/d5187c55b721f402c5ab60972d8ea13f
- 22.05.26
- Biomart Export
- Dataset: 23262 / 86369 Genes Human genes (GRCh38.p14) Ensembl 115
- Filters Gene type: protein_coding
- Attributes Gene stable ID, Gene start (bp), Gene end (bp), Chromosome/scaffold name

##biomart_export_humanHK_rescuemappedENSG_toENST_CCD_RefSeq250526.tsv
- https://www.ensembl.org/biomart/martview/99e4d9b03d9b751102db79c2179c210e
- 25.05.26
- Dataset: 20 / 86369 Genes Human genes (GRCh38.p14)
- Filters Gene stable ID(s) [e.g. ENSG00000000003]: [ID-list specified] Gene type: protein_coding
- Attributes Gene stable ID Transcript stable ID CCDS ID RefSeq mRNA ID
