library(Seurat)
library(ggplot2)

# Load and preprocess visium data
visium_subset <- readRDS("./output/processed_data/mouse_brain/mouse_brain_spatial_cortex.rds")

# SCTransform data
visium_subset <- SCTransform(visium_subset, assay = "Spatial", verbose = FALSE)

# Dimensionality reduction
visium_subset <- RunPCA(visium_subset, features = VariableFeatures(object = visium_subset))
visium_subset <- FindNeighbors(visium_subset, reduction = "pca", dims = 1:30)
visium_subset <- FindClusters(visium_subset, verbose = FALSE)
visium_subset <- RunUMAP(visium_subset, reduction = "pca", dims = 1:30)

SpatialDimPlot(visium_subset, label = TRUE)

saveRDS(visium_subset, "./output/processed_data/mouse_brain/mouse_brain_spatial_cortex_sctransform_pca_clustered.rds")

# Plotting clustering some genes' expression
visium_subset <- readRDS("./output/processed_data/mouse_brain/mouse_brain_spatial_cortex_sctransform_pca_clustered.rds")
DefaultAssay(visium_subset) <- "SCT"

palettePaired <- c(
  "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c",
  "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#fff116",
  "#d8a283", "#b15928", "#b0c288", "#81b806", "#bcffde", "#00ff82",
  "#e7aee0", "#a40192", "#d77070", "#ff0052", "#c6bfbf", "#6a5b5b")

SpatialDimPlot(visium_subset, label = FALSE, cols = palettePaired) +
  guides(fill = guide_legend(title = "Cluster"))
ggsave("./output/images/mouse_brain_cortex_analysis/SCT_spatial_clustering.png")

DimPlot(visium_subset, cols = palettePaired)
ggsave("./output/images/mouse_brain_cortex_analysis/SCT_spatial_clustering_umap.png", height = 4, width = 4)


plot <- SpatialFeaturePlot(visium_subset, features = c("Hpca", "Ttr"))
png("./output/images/mouse_brain_cortex_analysis/SCT_spatial_Hpca_Ttr.png", height = 700, width = 1200)
print(plot)
dev.off()


spatial_clustering <- SpatialDimPlot(visium_subset, label = FALSE, pt.size.factor = 2)+
  theme(legend.position="none")
clustering <- DimPlot(visium_subset, label = TRUE, label.size = 6)
nCount <- FeaturePlot(visium_subset, features = "nCount_Spatial")
nFeature <- FeaturePlot(visium_subset, features = "nFeature_Spatial")
wrap_plots(clustering, spatial_clustering, nCount, nFeature, heights = c(1, 1), widths = c(1, 1))
ggsave("./output/images/mouse_brain_cortex_analysis/SCT_clustering.png")
