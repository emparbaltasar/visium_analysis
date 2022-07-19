library(Seurat)
library(patchwork)
library(tidyverse)

# Load and preprocess visium data
visium <- Load10X_Spatial("./data/SCD-VI-i001-004/spaceranger_output/SCD-VI-i002/outs",
                          filename = "filtered_feature_bc_matrix.h5",
                          slice = "slice1",
                          filter.matrix = TRUE
)

# LogNormalize data
visium <- NormalizeData(visium, normalization.method = "LogNormalize", assay = "Spatial")
visium <- FindVariableFeatures(visium, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(visium)
visium <- ScaleData(visium, features = all.genes)

# Dimensionality reduction
visium <- RunPCA(visium, features = VariableFeatures(object = visium))
visium <- FindNeighbors(visium, reduction = "pca", dims = 1:30)
visium <- FindClusters(visium, verbose = FALSE)
visium <- RunUMAP(visium, reduction = "pca", dims = 1:30)

SpatialDimPlot(visium, label = TRUE)

dir.create("./output/processed_data/SCD-VI-i002")
saveRDS(visium, "./output/processed_data/SCD-VI-i002/SCD-VI-i002_lognorm_pca_clustered.rds")


# Plotting clustering some genes' expression

visium <- readRDS("./output/processed_data/SCD-VI-i002/SCD-VI-i002_lognorm_pca_clustered.rds")

spatial_clustering <- SpatialDimPlot(visium, label = FALSE, pt.size.factor = 2.5)+
  theme(legend.position="none")
clustering <- DimPlot(visium, label = TRUE, label.size = 6)
nCount <- FeaturePlot(visium, features = "nCount_Spatial")
nFeature <- FeaturePlot(visium, features = "nFeature_Spatial")


wrap_plots(clustering, spatial_clustering, nCount, nFeature, heights = c(1, 1), widths = c(1, 1))
