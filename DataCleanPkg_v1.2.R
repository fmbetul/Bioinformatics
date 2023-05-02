# DataCleanFunc_v1.1
# 2023/04/30
# Minor developments

# Functions included in this package:
  # CheckUMAP()
  # CropCluster()
  # getLeftovers()
  # assignExtraCluster()
  # validateExtraCluster()
  # getFinalClusters()



library(Seurat)
library(Signac)
library(ggplot2)



# CheckUMAP() function
CheckUMAP = function(CheckInput, sample = SeuData){
  if(class(CheckInput) == "data.frame"){
    CheckMeta = as.data.frame(c(rep("RoI", length(row.names(CheckInput)))))
    colnames(CheckMeta) = "Pop"
    row.names(CheckMeta) = row.names(CheckInput)
  }else{
    CheckMeta = as.data.frame(c(rep("RoI", length(colnames(CheckInput)))))
    colnames(CheckMeta) = "Pop"
    row.names(CheckMeta) = colnames(CheckInput)
  }
  
  sample = AddMetaData(sample, CheckMeta, "CheckMeta")     # CheckMeta is sth like: NA NA NA
  Idents(sample) = "CheckMeta"
  DimPlot(sample, reduction="umap", cols = c("#cb181d"), na.value = "grey85")
  
  }






# CropCluster() function
CropCluster = function(SeuData, ClusterName, ClusterNo = 0, 
                       min.x = -100, max.x = 100,
                       min.y = -100, max.y = 100){
  
  
  Filename = as.character(substitute(SeuData))
  
  # get cluster names:
  Clusters <- levels(Idents(SeuData))
  
  # get the UMAP coordinates of each cell
  UMAP_Data <- as.data.frame(SeuData@reductions$umap@cell.embeddings)  
  
  
  
  
  # 1. Select a new cluster & check its properties
  
  if(ClusterNo > 0){
    ClusterName <- Clusters[ClusterNo]
    if (!is.character(ClusterName)) {
      stop("Error: This ClusterNo is not available.")
    }
  }

  if(!(ClusterName %in% Clusters)){
    stop("Error: This ClusterName is not valid.")
  }
  
  cluster.barcodes <- WhichCells(SeuData, idents = ClusterName)
  
  # Check-point
  print(paste("Current cluster:", ClusterName, sep = " "))
  print(paste("Total number of cells in this cluster:", length(cluster.barcodes), sep = " "))
  
  # get UMAP coordinates of the cells in the selected cluster
  Cluster_Coords <- subset(UMAP_Data, row.names(UMAP_Data) %in% cluster.barcodes)
  
  # look at the selected cluster on UMAP
  check.cluster <- CheckUMAP(Cluster_Coords, sample = SeuData)
  print(check.cluster)
  
  ###
  
  # set the cluster boundaries
  Cluster_Clean = subset(UMAP_Data,
                         row.names(UMAP_Data) %in% cluster.barcodes 
                         & UMAP_Data $UMAP_1 > min.x                               
                         & UMAP_Data $UMAP_1 < max.x
                         & UMAP_Data $UMAP_2 > min.y
                         & UMAP_Data $UMAP_2 < max.y)
  
  
  print(paste("Cells left in", ClusterName ,"cluster:", dim(Cluster_Clean)[1], sep = " "))
  print(paste("Current number of cells in current Leftovers must be", length(cluster.barcodes) - dim(Cluster_Clean)[1], sep = " "))
  
  # look at the new cluster on UMAP
  check.cropped.cluster <- CheckUMAP(Cluster_Clean, sample = SeuData)
  print(check.cropped.cluster)
  
  # Save the _Clean cluster
  set.name <- paste0(ClusterName, "_Clean")     # eg."POMC_Clean"
  assign(set.name, Cluster_Clean, envir = globalenv())
                         
                         
  ###
  
  # 3. Current_Leftovers (double-check)
  
  Current_Leftovers <- subset(Cluster_Coords, !(row.names(Cluster_Coords) %in% row.names(Cluster_Clean)))
  print(paste("Current number of cells in current Leftovers is", dim(Current_Leftovers)[1], sep = " "))
  
  # look at the Leftovers on UMAP
  check.leftovers <- CheckUMAP(Current_Leftovers, sample = SeuData)
  print(check.leftovers)
  
  # look at all clusters overall
  print(DimPlot(SeuData, reduction = "umap", label  =TRUE))
  
  }







# getLeftovers() function
getLeftovers <- function(SeuData){
  
  # get the UMAP coordinates of each cell
  UMAP_Data <- as.data.frame(SeuData@reductions$umap@cell.embeddings)  
  
  # get the vector containing names of the _Clean objects
  Clean_Cluster_Names <- ls(envir = globalenv())[grepl("_Clean$", ls(envir = globalenv()))]
  
  # get the cell barcodes in the _Clean objects
  assigned_barcodes <- character(0)
  for(n in Clean_Cluster_Names){
    assigned_barcodes <- c(assigned_barcodes, row.names(get(n, envir = globalenv())))
  }
  
  print(paste("Number of cells assigned:", length(assigned_barcodes), sep = " "))
  print(paste("Number of cells in Leftovers must be:", dim(UMAP_Data)[1] - length(assigned_barcodes), sep = " "))
  
  # get the cell barcodes of the Leftovers
  Leftovers <- subset(UMAP_Data, !(row.names(UMAP_Data) %in% assigned_barcodes))
  
  # Save the Leftovers object in the global environment
  assign("Leftovers", Leftovers, envir = globalenv())
  
  print(paste("Number of cells in Leftovers:", dim(Leftovers)[1], sep = " "))
  print(CheckUMAP(sample = SeuData, CheckInput = Leftovers))
  
  }





# assignExtraCluster() function
assignExtraCluster <- function(SeuData, ClusterName, min.x = -50, max.x = 50, min.y = -50, max.y = 50) {
  
  # Get the UMAP coordinates of each cell
  UMAP_Data <- as.data.frame(SeuData@reductions$umap@cell.embeddings)
  
  # Check Leftovers
  Leftovers <- get("Leftovers", envir = globalenv())
  print(paste("Number of cells in current Leftovers:", dim(Leftovers)[1], sep = " "))
  print(CheckUMAP(sample = SeuData, CheckInput = Leftovers))
  
  # Define the borders of the _Extra cluster
  ExtraCluster <- subset(UMAP_Data, row.names(UMAP_Data) %in% row.names(Leftovers) 
                         & UMAP_Data$UMAP_1 > min.x
                         & UMAP_Data$UMAP_1 < max.x
                         & UMAP_Data$UMAP_2 > min.y
                         & UMAP_Data$UMAP_2 < max.y)
  
  # Visualize the _Extra cluster
  print(paste("Number of cells assigned to", ClusterName, "_Extra cluster:", dim(ExtraCluster)[1], sep = " "))
  print(CheckUMAP(sample = SeuData, CheckInput = ExtraCluster))
  
  # Save the _Extra cluster
  set.name <- paste0(ClusterName, "_Extra")
  assign(set.name, ExtraCluster, envir = globalenv())
  
  # Visualize and save NewLeftovers in the global environment
  NewLeftovers <- subset(Leftovers, !(row.names(Leftovers) %in% row.names(ExtraCluster)))
  assign("NewLeftovers", NewLeftovers, envir = globalenv()) # <-- This line is changed
  print(paste("Number of cells in NewLeftovers must be:", dim(Leftovers)[1] - dim(ExtraCluster)[1], sep = " "))
  print(paste("Number of cells in NewLeftovers is:", dim(NewLeftovers)[1], sep = " "))
  print(CheckUMAP(sample = SeuData, CheckInput = NewLeftovers))
  
  # look at all clusters overall
  print(DimPlot(SeuData, reduction = "umap", label  =TRUE))
  
  
  }







# validateExtraCluster() function
validateExtraCluster <- function(SeuData){
  
  # Update Leftovers
  assign("Leftovers", NewLeftovers, envir = globalenv())
  
  print(paste("Number of cells in updated Leftovers must be:", dim(NewLeftovers)[1], sep = " "))
  print(paste("Number of cells in updated Leftovers:", dim(Leftovers)[1], sep = " "))
  print(CheckUMAP(sample = SeuData, CheckInput = Leftovers))
  
  }




# getFinalClusters() function
getFinalClusters <- function(SeuData) {
  
  # Get the list of object names ending with _Clean or _Extra
  object_names <- ls(envir = globalenv())[grepl("_Clean$|_Extra$", ls(envir = globalenv()))]
  
  # Loop through the object names
  for (name in object_names) {
    # Extract the cluster name by removing the _Clean or _Extra suffix
    cluster_name <- gsub("(_Clean$)|(_Extra$)", "", name)
    
    # Get the cell barcodes from the data frame
    cell_barcodes <- row.names(get(name, envir = globalenv()))
    
    # Update the Idents of the clusters in the Seurat object
    Idents(object = SeuData, cells = cell_barcodes) <- as.character(cluster_name)
  }
  
  # Update the celltype column in the Seurat object
  SeuData$celltype <- Idents(SeuData)
  
  plot <- DimPlot(SeuData, reduction = "umap", label = TRUE)
  print(plot)
  
  # Return the updated Seurat object
  return(SeuData)
  
  }











