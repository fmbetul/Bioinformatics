# Check if required packages are installed, install if not
required_packages <- c("shiny", "Seurat", "ggplot2", "scales")
installed_packages <- installed.packages()[, "Package"]
missing_packages <- required_packages[!(required_packages %in% installed_packages)]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, dependencies = TRUE)
}

# Load libraries
library(shiny)
library(Seurat)
library(ggplot2)
library(scales)

# Define UI
ui <- fluidPage(
  titlePanel("Spatial Data Visualization with Seurat"),
  selectInput("sampleSelect", "Choose a Sample:",
              choices = list.files(path = "2_input/", pattern = "*.rds", full.names = FALSE)),
  textInput("geneInput", "Enter Gene Name:", value = "DBH"),
  fluidRow(
    column(6, plotOutput("spatialFeaturePlot")),
    column(6, plotOutput("featurePlot"))
  ),
  uiOutput("clusterCheckboxes"),
  fluidRow(
    column(6, plotOutput("spatialDimPlot")),
    column(6, plotOutput("umapDimPlot"))
  ),
  tags$script(HTML('
    $(document).on("change", ".cluster-checkbox", function () {
      var cluster = $(this).val();
      if ($(this).is(":checked")) {
        Shiny.onInputChange("addCluster", cluster);
      } else {
        Shiny.onInputChange("removeCluster", cluster);
      }
    });
  '))
)

# Define server logic
server <- function(input, output, session) {
  
  # Reactive expression to read the selected Seurat object
  selectedSeuratObject <- reactive({
    req(input$sampleSelect)
    readRDS(file.path("2_input/", input$sampleSelect))
  })
  
  # Update cluster checkboxes
  output$clusterCheckboxes <- renderUI({
    seuratObject <- selectedSeuratObject()
    checkboxGroupInput("clusterSelect", "Select Clusters to Highlight:",
                       choices = levels(Idents(seuratObject)),
                       selected = levels(Idents(seuratObject))[1],
                       inline = TRUE, width = "100%")
  })
  
  # Reactive expression for cluster colors based on selection
  clusterColors <- reactive({
    seuratObject <- selectedSeuratObject()
    ident_levels <- levels(Idents(seuratObject))
    ident_colors <- hue_pal()(length(ident_levels))
    named_ident_colors <- setNames(ident_colors, ident_levels)
    clusters_to_highlight <- input$clusterSelect
    
    # Initialize all clusters to 'lightgrey'
    highlight_colors <- setNames(rep('lightgrey', length(ident_levels)), ident_levels)
    
    # Assign highlight colors only to selected clusters
    if (length(clusters_to_highlight) > 0) {
      for (cluster in clusters_to_highlight) {
        if (cluster %in% names(highlight_colors)) {
          highlight_colors[cluster] <- named_ident_colors[cluster]
        }
      }
    }
    highlight_colors
  })
  
  # Render the SpatialFeaturePlot when a new gene is entered
  output$spatialFeaturePlot <- renderPlot({
    req(input$geneInput) # Make sure the gene name input is filled
    seuratObject <- selectedSeuratObject() # Get the selected Seurat object
    gene <- input$geneInput # Get the entered gene name
    
    # Verify if the gene exists in the dataset
    if (gene %in% rownames(seuratObject@assays$SCT@data)) {
      SpatialFeaturePlot(seuratObject, features = gene, interactive = FALSE,
                         pt.size.factor = 1, alpha = c(0.3, 0.8))
    } else {
      # Provide a message if the gene is not found
      ggplot() + ggtitle(paste("Gene", gene, "not found in dataset")) + theme_void()
    }
  })
  
  # Render the FeaturePlot when a new gene is entered
  output$featurePlot <- renderPlot({
    req(input$geneInput) # Make sure the gene name input is filled
    seuratObject <- selectedSeuratObject() # Get the selected Seurat object
    gene <- input$geneInput # Get the entered gene name
    
    # Verify if the gene exists in the dataset
    if (gene %in% rownames(seuratObject@assays$SCT@data)) {
      FeaturePlot(seuratObject, features = gene)
    } else {
      # Provide a message if the gene is not found
      ggplot() + ggtitle(paste("Gene", gene, "not found in dataset")) + theme_void()
    }
  })
  
  # Render the spatial dimensionality plot highlighting selected clusters
  output$spatialDimPlot <- renderPlot({
    seuratObject <- selectedSeuratObject()
    SpatialDimPlot(seuratObject, group.by = 'ident', 
                   pt.size.factor = 1, alpha = 0.7, 
                   cols = clusterColors()) +
      theme_minimal() +
      ggtitle("Spatial Dimensionality Plot")
  })
  
  # Render the UMAP dimensionality plot highlighting selected clusters
  output$umapDimPlot <- renderPlot({
    seuratObject <- selectedSeuratObject()
    DimPlot(seuratObject, group.by = 'ident', reduction = 'umap') +
      scale_color_manual(values = clusterColors()) +
      theme_minimal() +
      ggtitle("UMAP Plot")
  })
  
  # Update cluster selection based on checkbox selection
  observeEvent(input$addCluster, {
    selected_clusters <- input$clusterSelect
    selected_clusters <- c(selected_clusters, input$addCluster)
    updateCheckboxGroupInput(session, "clusterSelect",
                             choices = levels(Idents(selectedSeuratObject())),
                             selected = selected_clusters)
  })
  
  observeEvent(input$removeCluster, {
    selected_clusters <- input$clusterSelect
    selected_clusters <- selected_clusters[selected_clusters != input$removeCluster]
    updateCheckboxGroupInput(session, "clusterSelect",
                             choices = levels(Idents(selectedSeuratObject())),
                             selected = selected_clusters)
  })
}

# Run the app
shinyApp(ui = ui, server = server)
