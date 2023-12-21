# Feature Plots in one page

gene_vector <- endo82_genes # enter your gene list vector HERE
gene_vector <- unique(gene_vector)
gene_vector <- subset(gene_vector, gene_vector %in% row.names(SeuFile))

col_num <- 5 # enter desired number of FP columns in the page HERE
row_num <- ceiling(length(gene_vector)/col_num)

pdf(file = "3_output/v2.3/1_FPs_anderson_endo82_genes.pdf", width = col_num*3.5, height = row_num*3)
DefaultAssay(SeuFile) <- "RNA"
print(FeaturePlot(SeuFile, features = gene_vector, ncol = col_num))
dev.off()