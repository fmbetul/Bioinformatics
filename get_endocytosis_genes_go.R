# Get list of endocytosis list using GO Terms
library(biomaRt)

# Connect to BioMart
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Query genes related to endocytosis (example query, adjust as necessary)
genes_endocytosis_go <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'go_id', 'name_1006'),
                              filters = 'go', values = c('GO:0006897', # GO:0006897 is for endocytosis
                                                         'GO:0030139', # endocytic vesicles
                                                         'GO:0006898', # receptor mediated endocytosis
                                                         'GO:0030100' # regulation of endocytosis
                                                         ),
                              mart = ensembl)