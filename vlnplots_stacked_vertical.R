# Vertical Stacked VlnPlots

gene_vector <- endo82_genes # enter your gene list vector here
gene_vector <- unique(gene_vector)
gene_vector <- subset(gene_vector, gene_vector %in% row.names(SeuFile))

DefaultAssay(SeuFile) <- "RNA"
StorePlots = list()

# First loop: Process all genes except the last one
for(gene in gene_vector[-length(gene_vector)]){
  vlnA <- VlnPlot(SeuFile, features = gene, pt.size = 0, same.y.lims = F,)
  vlnA <- vlnA +  
    theme(axis.title.x=element_blank(), 
          axis.ticks.x= element_blank(), 
          axis.text.x=element_blank(), 
          axis.title.y=element_blank(), 
          axis.ticks.y = element_blank(), 
          axis.text.y = element_blank(),
          legend.position = "none", 
          plot.title = element_text(size = 12, hjust = 0)) # Align title to the left
  StorePlots[[gene]] = vlnA 
}

# Second loop: Process only the last gene
gene = gene_vector[length(gene_vector)]
vlnB <- VlnPlot(SeuFile, features = gene, pt.size = 0, same.y.lims = FALSE)
vlnB <- vlnB + 
  theme(axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        #axis.text.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        legend.position = "none", 
        plot.title = element_text(size = 12, hjust = 0))
StorePlots[[gene]] = vlnB

# Arrange VlnPlots stored in StorePlots
AllVlns <- arrangeGrob(grobs = StorePlots, 
                       heights = c(rep(unit(1, "null"), length(StorePlots) - 1), unit(2, "null")), # double last space
                       ncol = 1)

# save
ggsave(filename = "3_output/v2.3/2_Vlns_anderson_endo82_genes.pdf", 
       plot = AllVlns,
       width = length(unique(Idents(SeuFile))),
       height = 1.2*length(StorePlots)+1,
       limitsize = F)