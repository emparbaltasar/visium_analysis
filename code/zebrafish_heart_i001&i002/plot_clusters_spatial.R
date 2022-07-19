library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(ggpubr)

# Load and preprocess visium data
i001 <- readRDS("./output/processed_data/SCD-VI-i001/SCD-VI-i001_subset.rds")
i001 <- SCTransform(i001, assay = "Spatial", verbose = FALSE)

# Load and preprocess visium data
i002 <- Load10X_Spatial("./data/SCD-VI-i001-004/spaceranger_output/SCD-VI-i002/outs",
                          filename = "filtered_feature_bc_matrix.h5",
                          slice = "slice1",
                          filter.matrix = TRUE)
i002 <- SCTransform(i002, assay = "Spatial", verbose = FALSE)

cardmyo <- merge(i001,i002)


DefaultAssay(cardmyo) <- "SCT"
VariableFeatures(cardmyo) <- c(VariableFeatures(i001), VariableFeatures(i002))
cardmyo <- RunPCA(cardmyo, verbose = FALSE)
cardmyo <- FindNeighbors(cardmyo, dims = 1:30)
cardmyo <- FindClusters(cardmyo, verbose = FALSE)
cardmyo <- RunUMAP(cardmyo, dims = 1:30)

palettePaired <- c(
  "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c",
  "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#fff116",
  "#d8a283", "#b15928", "#b0c288", "#81b806", "#bcffde", "#00ff82",
  "#e7aee0", "#a40192", "#d77070", "#ff0052", "#c6bfbf", "#6a5b5b"
)

spatial_clustering <- SpatialDimPlot(cardmyo, label = FALSE, pt.size.factor = 2.5, cols = palettePaired) +
  theme(legend.position="none")


clustering <- DimPlot(cardmyo, label = TRUE, label.size = 6, cols = palettePaired)
nCount <- FeaturePlot(cardmyo, features = "nCount_Spatial") +
  labs(color = "UMI counts")
nFeature <- FeaturePlot(cardmyo, features = "nFeature_Spatial")


wrap_plots(clustering, spatial_clustering, nCount, nFeature, heights = c(1, 1), widths = c(0.5, 1))
dir.create("./output/images/SCD-VI-i001-i002_analysis/")
ggsave("./output/images/SCD-VI-i001-i002_analysis/SCD-VI-i001-i002_SCT_clustering.png")


SpatialDimPlot(cardmyo, label = FALSE, pt.size.factor = 2.5, cols = palettePaired) +
  guides(fill = guide_legend(title = "Cluster"))
ggsave("./output/images/SCD-VI-i001-i002_analysis/SCD-VI-i001-i002_SCT_clustering_2.png")

wrap_plots(clustering, nCount, heights = c(1, 1))




