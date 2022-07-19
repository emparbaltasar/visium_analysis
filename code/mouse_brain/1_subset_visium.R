library(Seurat)
library(patchwork)

# Load and preprocess spatial dataset
visium <- Load10X_Spatial("./data/mouse_brain/brain_spatial",
  filename = "filtered_feature_bc_matrix.h5",
  slice = "anterior1",
  filter.matrix = TRUE
)

# Plot nCounts
SpatialFeaturePlot(visium, features = "nCount_Spatial") + 
  theme(legend.position = "right",
        plot.margin = unit(c(0,0,0,0), units = "cm"))

ggsave("./output/images/mouse_brain_cortex_analysis/nCounts.png")

SpatialFeaturePlot(visium, features = "nFeature_Spatial") + 
  theme(legend.position = "right",
        plot.margin = unit(c(0,0,0,0), units = "cm"))

ggsave("./output/images/mouse_brain_cortex_analysis/nFeatures.png")

visium[["percent_mt_spatial"]] <- PercentageFeatureSet(visium, pattern = "^mt-")


SpatialFeaturePlot(visium, features = "percent_mt_spatial") + 
  theme(legend.position = "right",
        plot.margin = unit(c(0,0,0,0), units = "cm"))
ggsave("./output/images/mouse_brain_cortex_analysis/pct_mitochondrial_genes.png")



visium <- NormalizeData(visium, normalization.method = "LogNormalize", assay = "Spatial")
visium <- FindVariableFeatures(visium, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(visium)
visium <- ScaleData(visium, features = all.genes)
visium <- RunPCA(visium, features = VariableFeatures(object = visium))
visium <- FindNeighbors(visium, reduction = "pca", dims = 1:30)
visium <- FindClusters(visium, verbose = FALSE)
visium <- RunUMAP(visium, reduction = "pca", dims = 1:30)

SpatialDimPlot(visium, label = TRUE)


# Subset cortex
visium_subset <- subset(visium, idents = c(1, 2, 3, 5, 12, 6, 8))
visium_subset <- subset(visium_subset, anterior1_imagerow > 400 | anterior1_imagecol < 150, invert = TRUE)
visium_subset <- subset(visium_subset, anterior1_imagerow > 285 & anterior1_imagecol > 350, invert = TRUE)
visium_subset <- subset(visium_subset, anterior1_imagerow > 250 & anterior1_imagecol > 440, invert = TRUE)
SpatialDimPlot(visium_subset, label = TRUE)
visium_subset <- subset(visium, cells = Cells(visium_subset))
saveRDS(visium_subset, "./output/processed_data/mouse_brain/mouse_brain_spatial_cortex.rds")

visium_subset <- readRDS("./output/processed_data/mouse_brain/mouse_brain_spatial_cortex.rds")

SpatialFeaturePlot(visium_subset, features = "nCount_Spatial") + 
  theme(legend.position = "right",
        plot.margin = unit(c(0,0,0,0), units = "cm"))
ggsave("./output/images/mouse_brain_cortex_analysis/cortex_nCounts.png")
