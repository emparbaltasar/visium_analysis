library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(ggpubr)

# Load and preprocess visium data
visium <- Load10X_Spatial("./data/SCD-VI-i001-004/spaceranger_output/SCD-VI-i003_full_area",
                          filename = "filtered_feature_bc_matrix.h5",
                          slice = "slice1",
                          filter.matrix = TRUE)

regeneration_state <- read.csv("data/SCD-VI-i001-004/spaceranger_output/SCD-VI-i003/outs/regeneration_state.csv", row.names = "Barcode") # This file has to be created with loupe browser
visium <- AddMetaData(visium, regeneration_state, col.name = "regeneration_state")

SpatialDimPlot(visium, group.by = "regeneration_state", pt.size.factor = 1) +
  theme(legend.position="right")


# SCTransform data
visium <- SCTransform(visium, assay = "Spatial", verbose = FALSE)

# Dimensionality reduction
visium <- RunPCA(visium, features = VariableFeatures(object = visium))
visium <- FindNeighbors(visium, reduction = "pca", dims = 1:30)
visium <- FindClusters(visium, verbose = TRUE)
visium <- RunUMAP(visium, reduction = "pca", dims = 1:30)

SpatialDimPlot(visium, label = TRUE)
dir.create("./output/processed_data/SCD-VI-i003_full_area")
saveRDS(visium, "./output/processed_data/SCD-VI-i003_full_area/SCD-VI-i003_full_area_sctransform_pca_clustered.rds")

# Plotting clustering some genes' expression
visium <- readRDS("./output/processed_data/SCD-VI-i003_full_area/SCD-VI-i003_full_area_sctransform_pca_clustered.rds")
DefaultAssay(visium) <- "SCT"


SpatialFeaturePlot(visium, features = c("nppa"), pt.size.factor = 1)
dir.create("./output/images/SCD-VI-i003_full_area_analysis")
ggsave("./output/images/SCD-VI-i003_full_area_analysis/SCT_spatial_nppa.png")

spatial_clustering <- SpatialDimPlot(visium, label = FALSE, pt.size.factor = 1)+
  theme(legend.position="none")
clustering <- DimPlot(visium, label = TRUE, label.size = 6)
wrap_plots(clustering, spatial_clustering, heights = c(1, 1), widths = c(1, 1))
ggsave("./output/images/SCD-VI-i003_full_area_analysis/SCD-VI-i003_full_area_SCT_clustering.png")
