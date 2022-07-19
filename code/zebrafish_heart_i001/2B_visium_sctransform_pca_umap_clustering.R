library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(ggpubr)

# 2 states -----
# Load and preprocess visium data
visium_subset <- readRDS("./output/processed_data/SCD-VI-i001/SCD-VI-i001_subset.rds")


## nCounts #########

# Plot nCounts
plot1 <- VlnPlot(visium_subset, features = "nCount_Spatial", pt.size = 0.1) +
  NoLegend() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_blank(),
        plot.margin = unit(c(0,0,0,0), units = "cm")) +
  ylab("nCount")
  
plot2 <- SpatialFeaturePlot(visium_subset, features = "nCount_Spatial", pt.size.factor = 2.5) + 
  theme(legend.position = "right",
        plot.margin = unit(c(0,0,0,0), units = "cm"))


wrap_plots(plot1, plot2)+
  plot_annotation(title = "UMIs per spot",
                  theme = theme(plot.title = element_text(hjust = 0.5)))

ggsave("./output/images/SCD-VI-i001_analysis/SCD-VI-i001_ncount.png")




## nFeature #########

# Plot n genes per spot
plot1 <- VlnPlot(visium_subset, features = "nFeature_Spatial", pt.size = 0.1) +
  NoLegend() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_blank())+
  ylab("nFeature")

plot2 <- SpatialFeaturePlot(visium_subset, features = "nFeature_Spatial", pt.size.factor = 2.5) + 
  theme(legend.position = "right")


wrap_plots(plot1, plot2)+
  plot_annotation(title = "Genes per spot",
                  theme = theme(plot.title = element_text(hjust = 0.5)))

ggsave("./output/images/SCD-VI-i001_analysis/SCD-VI-i001_ngenes.png")


## % mitochondrial genes #########

visium_subset[["percent_mt_spatial"]] <- PercentageFeatureSet(visium_subset, pattern = "^mt-")

# Plot mitochondrial genes per spot
plot1 <- VlnPlot(visium_subset, features = "percent_mt_spatial", pt.size = 0.1) +
  NoLegend() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_blank()) +
  ylab("% Mitochondrial")

plot2 <- SpatialFeaturePlot(visium_subset, features = "percent_mt_spatial", pt.size.factor = 2.5) + 
  theme(legend.position = "right")


wrap_plots(plot1, plot2)+
  plot_annotation(title = "% mitochondrial genes per spot",
                  theme = theme(plot.title = element_text(hjust = 0.5)))

ggsave("./output/images/SCD-VI-i001_analysis/SCD-VI-i001_mitochondrial_genes.png")



# Transformation, PCA, clustering ----
# SCTransform data
visium_subset <- SCTransform(visium_subset, assay = "Spatial", verbose = FALSE)

# Dimensionality reduction
visium_subset <- RunPCA(visium_subset, features = VariableFeatures(object = visium_subset))
visium_subset <- FindNeighbors(visium_subset, reduction = "pca", dims = 1:30)
visium_subset <- FindClusters(visium_subset, verbose = FALSE)
visium_subset <- RunUMAP(visium_subset, reduction = "pca", dims = 1:30)

SpatialDimPlot(visium_subset, label = TRUE)

saveRDS(visium_subset, "./output/processed_data/SCD-VI-i001/SCD-VI-i001_subset_sctransform_pca_clustered.rds")

# Plotting clustering some genes' expression
visium_subset <- readRDS("./output/processed_data/SCD-VI-i001/SCD-VI-i001_subset_sctransform_pca_clustered.rds")
DefaultAssay(visium_subset) <- "SCT"

dir.create("./output/images/SCD-VI-i001_analysis")

SpatialFeaturePlot(visium_subset, features = c("nppa"))


spatial_clustering <- SpatialDimPlot(visium_subset, label = FALSE, pt.size.factor = 2.5)+
  theme(legend.position="none")
clustering <- DimPlot(visium_subset, label = TRUE, label.size = 6)
nCount <- FeaturePlot(visium_subset, features = "nCount_Spatial")
nFeature <- FeaturePlot(visium_subset, features = "nFeature_Spatial")


wrap_plots(clustering, spatial_clustering, nCount, nFeature, heights = c(1, 1), widths = c(1, 1))
ggsave("./output/images/SCD-VI-i001_analysis/SCD-VI-i001_SCT_clustering.png")

