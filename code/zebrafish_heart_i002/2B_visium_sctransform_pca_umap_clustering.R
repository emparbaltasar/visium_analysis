library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(ggpubr)

# Load and preprocess visium data
visium <- Load10X_Spatial("./data/SCD-VI-i001-004/spaceranger_output/SCD-VI-i002/outs",
                          filename = "filtered_feature_bc_matrix.h5",
                          slice = "slice1",
                          filter.matrix = TRUE
                          )

## nCounts #########
# Plot nCounts
plot1 <- VlnPlot(visium, features = "nCount_Spatial", pt.size = 0.1) +
  NoLegend() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_blank(),
        plot.margin = unit(c(0,0,0,0), units = "cm")) +
  ylab("nCount")
  
plot2 <- SpatialFeaturePlot(visium, features = "nCount_Spatial", pt.size.factor = 2.5) + 
  theme(legend.position = "right",
        plot.margin = unit(c(0,0,0,0), units = "cm"))


wrap_plots(plot1, plot2)+
  plot_annotation(title = "UMIs per spot",
                  theme = theme(plot.title = element_text(hjust = 0.5)))

dir.create("./output/images/SCD-VI-i002_analysis")
ggsave("./output/images/SCD-VI-i002_analysis/SCD-VI-i002_ncount.png")




## nFeature #########
# Plot n genes per spot
plot1 <- VlnPlot(visium, features = "nFeature_Spatial", pt.size = 0.1) +
  NoLegend() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_blank())+
  ylab("nFeature")

plot2 <- SpatialFeaturePlot(visium, features = "nFeature_Spatial", pt.size.factor = 2.5) + 
  theme(legend.position = "right")


wrap_plots(plot1, plot2)+
  plot_annotation(title = "Genes per spot",
                  theme = theme(plot.title = element_text(hjust = 0.5)))

ggsave("./output/images/SCD-VI-i002_analysis/SCD-VI-i002_ngenes.png")


## % mitochondrial genes #########

visium[["percent_mt_spatial"]] <- PercentageFeatureSet(visium, pattern = "^mt-")

# Plot mitochondrial genes per spot
plot1 <- VlnPlot(visium, features = "percent_mt_spatial", pt.size = 0.1) +
  NoLegend() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_blank()) +
  ylab("% Mitochondrial")

plot2 <- SpatialFeaturePlot(visium, features = "percent_mt_spatial", pt.size.factor = 2.5) + 
  theme(legend.position = "right")


wrap_plots(plot1, plot2)+
  plot_annotation(title = "% mitochondrial genes per spot",
                  theme = theme(plot.title = element_text(hjust = 0.5)))

ggsave("./output/images/SCD-VI-i002_analysis/SCD-VI-i002_mitochondrial_genes.png")



# Transformation, PCA, clustering ----
# SCTransform data
visium <- SCTransform(visium, assay = "Spatial", verbose = FALSE)

# Dimensionality reduction
visium <- RunPCA(visium, features = VariableFeatures(object = visium))
visium <- FindNeighbors(visium, reduction = "pca", dims = 1:30)
visium <- FindClusters(visium, verbose = FALSE)
visium <- RunUMAP(visium, reduction = "pca", dims = 1:30)

SpatialDimPlot(visium, label = TRUE)

saveRDS(visium, "./output/processed_data/SCD-VI-i002/SCD-VI-i002_sctransform_pca_clustered.rds")

# Plotting clustering some genes' expression
visium <- readRDS("./output/processed_data/SCD-VI-i002/SCD-VI-i002_sctransform_pca_clustered.rds")
DefaultAssay(visium) <- "SCT"
SpatialFeaturePlot(visium, features = c("nppa"))


spatial_clustering <- SpatialDimPlot(visium, label = FALSE, pt.size.factor = 2.5)+
  theme(legend.position="none")
clustering <- DimPlot(visium, label = TRUE, label.size = 6)
nCount <- FeaturePlot(visium, features = "nCount_Spatial")
nFeature <- FeaturePlot(visium, features = "nFeature_Spatial")


wrap_plots(clustering, spatial_clustering, nCount, nFeature, heights = c(1, 1), widths = c(1, 1))
ggsave("./output/images/SCD-VI-i002_analysis/SCD-VI-i002_SCT_clustering.png")

