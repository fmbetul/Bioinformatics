# Bioinformatics

Bioinformatics tools

## generate_feature_plots

-   Generates multiple FeaturePlots for the SeuratObjects and genes of interest in a scaled manner.

-   The `generate_feature_plots()` function takes Seurat objects and the list of genes you are interested in as input. It generates FeaturePlots for your specified genes, sets the upper limit according to the max value of each gene, and saves the plots as separate PDF files in the output folder, one for each Seurat object.

    -   Inputs:
        -   List of SeuratObject
        -   Vector of genes
        
    -   Output:
        -   FeaturePlots saved in PDF

-   Make sure you have your Seurat objects in your global environment.

```{r}
# Load the function:
source("https://raw.githubusercontent.com/fmbetul/Bioinformatics/main/generate_feature_plots.R")

# Set Inputs
SeuFiles <- list(File1 = SeuratObject1, File2 = SeuratObject2)  # List of Seurat objects
gene_list <- c()  # Character vector of genes you are interested in

# Generate scaled FeaturePlots using the generate_feature_plots() function
generate_feature_plots(seurat_objects = SeuFiles, genes = gene_list)
```
