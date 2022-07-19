#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
getwd()

# Data
MyPalette <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(n = 100)

seurat_lognorm <- readRDS("../../../output/processed_data/SCD-VI-i003/SCD-VI-i003_lognorm_seurat_predictions.rds")
DefaultAssay(seurat_lognorm) <- "prediction.score.id"

seurat_sct <- readRDS("../../../output/processed_data/SCD-VI-i003/SCD-VI-i003_sctransform_seurat_predictions.rds")
DefaultAssay(seurat_sct) <- "prediction.score.id"

res <- readRDS("../../../output/processed_data/SCD-VI-i003/spotlight_lognorm-hvg_rawdata.rds")
spotlight <- readRDS("../../../output/processed_data/SCD-VI-i003/SCD-VI-i003_lognorm_pca_clustered.rds")
mat <- res$mat
spotlight[["predictionsspotlight"]] <- CreateAssayObject(t(mat))
DefaultAssay(spotlight) <- "predictionsspotlight"

cell2loc <- Load10X_Spatial("../../../data/SCD-VI-i001-004/spaceranger_output/SCD-VI-i003/outs",
                          filename = "filtered_feature_bc_matrix.h5",
                          slice = "slice1",
                          filter.matrix = TRUE)
cell2loc_pred <- read.csv("../../../output/cell2location/zebrafish_heart/SCD-VI-i003/cell2location_map/cell2location_prediction.csv",
                          row.names=1)
colnames(cell2loc_pred) <- c("B-cells","Bl.ves.EC (apnln)","Bl.ves.EC (lyve1)","Bl.ves.EC (plvapb)",
                             "Cardiomyocytes (proliferating)", "Cardiomyocytes (ttn.2) A", "Cardiomyocytes (ttn.2) V", "Cardiomyocytes A",
                             "Cardiomyocytes V", "Dead cells", "Endocardium (A)", "Endocardium (V)",
                             "Endocardium (frzb)", "Epicardium (Atrium)", "Epicardium (Ventricle)", "Fibroblast",
                             "Fibroblast (cfd)", "Fibroblast (col11a1a)", "Fibroblast (col12a1a)", "Fibroblast (cxcl12a)",
                             "Fibroblast (mpeg1.1)","Fibroblast (nppc)", "Fibroblast (proliferating)", "Fibroblast (spock3)",
                             "Fibroblast-like cells", "Macrophage (CM duplex)", "Macrophage (Endothelia duplex)", "Macrophage (Ery duplex)",
                             "Macrophage (Fibroblast duplex)", "Macrophage (apoeb)", "Macrophage (cd59)", "Macrophage (epdl)",
                             "Macrophage (il1b)", "Macrophage (proliferating)", "Macrophages", "Monocytes",
                             "Myelin cells", "Neuronal cells", "Neutrophils", "Perivascular cells", 
                             "Proliferating cells", "Smooth muscle cells", "T-cells", "T-cells (il4/13)",
                             "T-cells (proliferating)")
cell2loc_pred[cell2loc_pred < 1] <- 0
cell2loc_pred_prop <- cell2loc_pred / rowSums(cell2loc_pred)
cell2loc_pred_prop[is.na(cell2loc_pred_prop)] <- 0
cell2loc[["cell2loc"]] <- CreateAssayObject(t(cell2loc_pred_prop))
DefaultAssay(cell2loc) <- "cell2loc"

celltypes <- as_tibble(rownames(seurat_lognorm@assays$prediction.score.id@data))
colnames(celltypes) <- "cell_type"

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Comparison deconvolution methods"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          helpText("Which cell type do you want to plot?"),
          
          selectInput("celltype", label = "Choose a cellype:",
                      choices = celltypes$cell_type
                      #, selected = swiss$Fertility
                      ),
          
          width = 2
        ),
        

        # Show a plot of the generated distribution
        mainPanel(
          fluidRow(
            verticalLayout(
              splitLayout(cellWidths = c("50%", "50%"), 
                          plotOutput("seurat_lognorm"), 
                          plotOutput("seurat_sct")), 
              splitLayout(cellWidths = c("50%", "50%"), 
                          plotOutput("spotlight"), 
                          plotOutput("cell2loc"))
            )
           # plotOutput("seurat_lognorm"),
           # plotOutput("seurat_sct"),
           # plotOutput("spotlight")
        ),
        
        width = 10
        )
    )
)

# Define server logic ----
server <- function(input, output) {
  
  # Plot 
  output$seurat_lognorm <- renderPlot({
    SpatialFeaturePlot(seurat_lognorm, features = input$celltype, pt.size.factor = 2, crop = TRUE) +
      ggplot2::scale_fill_gradientn(colours = MyPalette, limits = c(0.0, 1.001), breaks = c(0.0, 0.5, 1.0)) +
      ggplot2::ggtitle("Seurat (LogNorm)") +
      ggplot2::theme(legend.position = "none")
  })
  
  output$seurat_sct <- renderPlot({
    SpatialFeaturePlot(seurat_sct, features = input$celltype, pt.size.factor = 2, crop = TRUE) +
      ggplot2::scale_fill_gradientn(colours = MyPalette, limits = c(0.0, 1.001), breaks = c(0.0, 0.5, 1.0)) +
      ggplot2::ggtitle("Seurat (SCT)") +
      ggplot2::theme(legend.position = "right")
  })
  
  output$spotlight <- renderPlot({
    SpatialFeaturePlot(spotlight, features = input$celltype, pt.size.factor = 2, crop = TRUE) +
      ggplot2::scale_fill_gradientn(colours = MyPalette, limits = c(0.0, 1.001), breaks = c(0.0, 0.5, 1.0)) +
      ggplot2::ggtitle("SPOTlight") +
      ggplot2::theme(legend.position = "none")
  })
  
  output$cell2loc <- renderPlot({
    SpatialFeaturePlot(cell2loc, features = input$celltype, pt.size.factor = 2, crop = TRUE) +
      ggplot2::scale_fill_gradientn(colours = MyPalette, limits = c(0.0, 1.001), breaks = c(0.0, 0.5, 1.0)) +
      ggplot2::ggtitle("Cell2location") +
      ggplot2::theme(legend.position = "none")
  })
  

  
}

# Run the app ----
shinyApp(ui = ui, server = server)