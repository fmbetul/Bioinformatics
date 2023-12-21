# Horizontal Stacked VlnPlots

Vln_num <- 21 # enter desired number of the VlnPlots in one row HERE
# Rename your Seurat Object as SeuFile

DefaultAssay(SeuFile) <- "RNA"

# Generate a list of gene vectors, each with Vln_num genes
gene_vector <- endo82_genes # enter your gene vector HERE
gene_vector <- unique(gene_vector)
gene_vector <- subset(gene_vector, gene_vector %in% row.names(SeuFile))
gene_chunks <- split(gene_vector, ceiling(seq_along(gene_vector)/Vln_num))

# Initialize a list to store the combined plots for each chunk
combined_plots <- list()

# Generate a blank ggplot object to use as a filler
blank_plot <- ggplot() + 
  theme_void() + 
  theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# Loop through each chunk of genes
for(i in seq_along(gene_chunks)){
  gene_chunk <- gene_chunks[[i]]
  StorePlots <- list()
  
  # Loop through each gene in the chunk to create violin plots
  for(j in seq_along(gene_chunk)){
    gene <- gene_chunk[j]
    vlnPlot <- VlnPlot(SeuFile, features = gene, pt.size = 0, same.y.lims = FALSE) + 
      coord_flip()
    
    # For the first plot in the chunk, add cluster names on the y-axis
    if(j == 1){
      vlnPlot <- vlnPlot + 
        theme(axis.title.x = element_blank(), 
              axis.ticks.x = element_blank(), 
              axis.text.x = element_blank(), 
              axis.title.y = element_blank(), 
              axis.ticks.y = element_blank(), 
              # axis.text.y = element_blank(),
              legend.position = "none", 
              plot.title = element_text(size = 12, hjust = 0)) 
    } else {
      vlnPlot <- vlnPlot + 
        theme(axis.title.x = element_blank(), 
              axis.ticks.x = element_blank(), 
              axis.text.x = element_blank(), 
              axis.title.y = element_blank(), 
              axis.ticks.y = element_blank(), 
              axis.text.y = element_blank(), # Disable y-axis text for non-first plots
              legend.position = "none", 
              plot.title = element_text(size = 12, hjust = 0))
    }
    
    StorePlots[[gene]] = vlnPlot
  }
  
  # If there are fewer than Vln_num plots, fill in the remaining space with blank plots
  num_plots_to_add <- Vln_num - length(StorePlots)
  if(num_plots_to_add > 0){
    for(k in 1:num_plots_to_add){
      StorePlots[[paste0("blank", k)]] <- blank_plot
    }
  }
  
  # Arrange the violin plots horizontally for this chunk
  AllVlns <- ggarrange(plotlist = StorePlots, widths = c(2.5, rep(1, length(StorePlots) - 1)), nrow = 1)
  
  # Add the combined plot to the list
  combined_plots[[i]] <- AllVlns
}

# Now, combine all the chunks vertically
final_plot <- ggarrange(plotlist = combined_plots, ncol = 1, nrow = length(combined_plots))

# Save the final plot to a file
ggsave(filename = "3_output/v2.3/2_Vlns_horizontal_anderson_endo82_genes.pdf", 
       plot = final_plot,
       height = length(unique(Idents(SeuFile)))*length(gene_chunks), # cluster_num * gene_chunk
       width = 1*Vln_num + 1.5) # 1 for each VlnPlot plus some padding for cluster names