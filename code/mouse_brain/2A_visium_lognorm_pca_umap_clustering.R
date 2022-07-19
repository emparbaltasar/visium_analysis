library(Seurat)

# Load and preprocess visium data
visium_subset <- readRDS("./output/processed_data/mouse_brain/mouse_brain_spatial_cortex.rds")

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

saveRDS(visium_subset, "./output/processed_data/mouse_brain/mouse_brain_spatial_cortex_lognorm_pca_clustered.rds")
