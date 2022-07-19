library(Seurat)

# Load and preprocess visium data
visium_subset <- readRDS("./output/processed_data/SCD-VI-i004/SCD-VI-i004_subset.rds")

# LogNormalize data
visium_subset <- NormalizeData(visium_subset, normalization.method = "LogNormalize", assay = "Spatial")
visium_subset <- FindVariableFeatures(visium_subset, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(visium_subset)
visium_subset <- ScaleData(visium_subset, features = all.genes)

# Dimensionality reduction
visium_subset <- RunPCA(visium_subset, features = VariableFeatures(object = visium_subset))
visium_subset <- FindNeighbors(visium_subset, reduction = "pca", dims = 1:30)
visium_subset <- FindClusters(visium_subset, verbose = FALSE)
visium_subset <- RunUMAP(visium_subset, reduction = "pca", dims = 1:30)

SpatialDimPlot(visium_subset, label = TRUE)

saveRDS(visium_subset, "./output/processed_data/SCD-VI-i004/SCD-VI-i004_subset_lognorm_pca_clustered.rds")


# Plotting clustering some genes' expression

visium_subset <- readRDS("./output/processed_data/SCD-VI-i004/SCD-VI-i004_subset_lognorm_pca_clustered.rds")

spatial_clustering <- SpatialDimPlot(visium_subset, label = FALSE, pt.size.factor = 2.5)+
  theme(legend.position="none")
clustering <- DimPlot(visium_subset, label = TRUE, label.size = 6)
nCount <- FeaturePlot(visium_subset, features = "nCount_Spatial")
nFeature <- FeaturePlot(visium_subset, features = "nFeature_Spatial")


wrap_plots(clustering, spatial_clustering, nCount, nFeature, heights = c(1, 1), widths = c(1, 1))
