# FeatureBot_v2.0
# Generated: 2023/5/1
# Object: generate a tool to make scaled FeaturePlots


# Load Libraries
library(Seurat)
library(Signac)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(patchwork)



# Define the generate_feature_plots() function to generate feature plots for a list of Seurat objects and a list of genes
generate_feature_plots <- function(seurat_objects, genes) {
  
  # Find the maximum value for each gene across all Seurat objects
  max_values <- sapply(genes, function(gene) {
    max(sapply(seurat_objects, function(seurat_obj) {
      
      # Set the default assay for the Seurat object to 'RNA'
      DefaultAssay(seurat_obj) <- 'RNA'
      
      # Check if the gene is present in the RNA assay
      if (gene %in% rownames(seurat_obj@assays$RNA@counts)) {
        
        # Fetch the gene expression data
        max(Seurat::FetchData(seurat_obj, vars = gene), na.rm = TRUE)
      } else {
        NA
        }
      }))
    }, USE.NAMES = TRUE)
  
  
  
  # Loop through the Seurat objects
  for (seurat_obj_name in names(seurat_objects)) {
    
    # Get the current Seurat object
    seurat_obj <- seurat_objects[[seurat_obj_name]]
    
    # Set the default assay for the Seurat object to 'RNA'
    DefaultAssay(seurat_obj) <- 'RNA'
    
    # Create a combined feature plot for all genes in the Seurat object
    p <- FeaturePlot(seurat_obj, features = genes, combine = TRUE)
    
    # Loop through the genes and add a color scale to each feature plot
    for (i in seq_along(genes)) {
      gene <- genes[[i]]
      p[[i]] <- p[[i]] + scale_color_gradientn(colours = c("lightgrey", "navy"), limits = c(0, max_values[[gene]]))
      }
    
    # Save the feature plots as a PDF file
    if (!file.exists("output")) {
      dir.create("output")
    }
    
    pdf_filename <- paste0("output/", seurat_obj_name, "_FeaturePlots.pdf")
    pdf(file = pdf_filename, width = 15, height = 18)
    print(p)
    dev.off()
    
    }
  }


