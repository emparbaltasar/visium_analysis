library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(ggpubr)


# 2 states ----------------------------------------------------------------

# Load and preprocess visium data
visium_subset <- readRDS("./output/processed_data/SCD-VI-i003/SCD-VI-i003_subset.rds")
regeneration_state <- read.csv("data/SCD-VI-i001-004/spaceranger_output/SCD-VI-i003/outs/regeneration_state.csv", row.names = "Barcode") # This file has to be created with loupe browser
visium_subset <- AddMetaData(visium_subset, regeneration_state, col.name = "regeneration_state")

SpatialDimPlot(visium_subset, group.by = "regeneration_state", pt.size.factor = 2.5) +
  theme(legend.position="right")
dir.create("./output/images/SCD-VI-i003_analysis")
ggsave("./output/images/SCD-VI-i003_analysis/SCD-VI-i003_regeneration_state.png")



## nCounts #########

# Statistics
regeneration_state <- bind_cols(regeneration_state, visium_subset[["nCount_Spatial"]])

# Shapiro-Wilk normality test for regenerating nCount
with(regeneration_state, shapiro.test(nCount_Spatial[state == "regenerating"])) # p = 0.05729
# Shapiro-Wilk normality test for non-regenerating nCount
with(regeneration_state, shapiro.test(nCount_Spatial[state == "non-regenerating"])) # p = 2.2e-16
# --> we can't assume normality 
# --> Wilcoxon rank test

# Wilcoxon rank test
compare_means(nCount_Spatial ~ state, data = regeneration_state, method = "wilcox.test") # p = 2.3e-11 --> significantly different


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

plot3 <- VlnPlot(visium_subset, group.by = "regeneration_state", features = "nCount_Spatial", pt.size = 0.1) +
  NoLegend() +
  theme(axis.title.x = element_blank(),
        plot.title = element_blank(),
        plot.margin = unit(c(0,0,0,0), units = "cm"))+
  ylab("nCount") +
  stat_compare_means()

plot4 <- SpatialDimPlot(visium_subset, group.by = "regeneration_state", pt.size.factor = 2.5) +
  theme(legend.position="right")

wrap_plots(plot1, plot2, plot3, plot4, ncol = 2)+
  plot_annotation(title = "UMIs per spot - SCD-VI-i003",
                  theme = theme(plot.title = element_text(hjust = 0.5)))

ggsave("./output/images/SCD-VI-i003_analysis/SCD-VI-i003_ncount.png")





## nFeature #########

# Statistics
regeneration_state <- bind_cols(regeneration_state, visium_subset[["nFeature_Spatial"]])

# Shapiro-Wilk normality test for regenerating nCount
with(regeneration_state, shapiro.test(nFeature_Spatial[state == "regenerating"]))# p = 0.4957
# Shapiro-Wilk normality test for non-regenerating nCount
with(regeneration_state, shapiro.test(nFeature_Spatial[state == "non-regenerating"])) # p = 7.461e-13 
# --> we can't assume normality --> Wilcoxon rank test

# Wilcoxon rank test
compare_means(nFeature_Spatial ~ state, data = regeneration_state, method = "wilcox.test") # p = 5.5e-15 --> significantly different


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

plot3 <- VlnPlot(visium_subset, group.by = "regeneration_state", features = "nFeature_Spatial", pt.size = 0.1) +
  NoLegend() +
  theme(axis.title.x = element_blank(),
        plot.title = element_blank(),
        plot.margin = unit(c(0,0,0,0), units = "cm"))+
  ylab("nFeature") +
  stat_compare_means()

plot4 <- SpatialDimPlot(visium_subset, group.by = "regeneration_state", pt.size.factor = 2.5) +
  theme(legend.position="right")

wrap_plots(plot1, plot2, plot3, plot4, ncol = 2)+
  plot_annotation(title = "Genes per spot - SCD-VI-i003",
                  theme = theme(plot.title = element_text(hjust = 0.5)))

ggsave("./output/images/SCD-VI-i003_analysis/SCD-VI-i003_ngenes.png")


## % mitochondrial genes #########

visium_subset[["percent_mt_spatial"]] <- PercentageFeatureSet(visium_subset, pattern = "^mt-")

# Statistics
regeneration_state <- bind_cols(regeneration_state, visium_subset[["percent_mt_spatial"]])

# Shapiro-Wilk normality test for regenerating nCount
with(regeneration_state, shapiro.test(percent_mt_spatial[state == "regenerating"]))# p = 3.206e-06
# Shapiro-Wilk normality test for non-regenerating nCount
with(regeneration_state, shapiro.test(percent_mt_spatial[state == "non-regenerating"])) # p = 4.191e-06
# --> we can't assume normality 
# --> Wilcoxon rank test

# Wilcoxon rank test
compare_means(percent_mt_spatial ~ state, data = regeneration_state, method = "wilcox.test") # p = 7.9e-14 --> significantly different


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

plot3 <- VlnPlot(visium_subset, group.by = "regeneration_state", features = "percent_mt_spatial", pt.size = 0.1) +
  NoLegend() +
  theme(axis.title.x = element_blank(),
        plot.title = element_blank(),
        plot.margin = unit(c(0,0,0,0), units = "cm")) +
  ylab("% Mitochondrial") +
  stat_compare_means()

plot4 <- SpatialDimPlot(visium_subset, group.by = "regeneration_state", pt.size.factor = 2.5) +
  theme(legend.position="right")

wrap_plots(plot1, plot2, plot3, plot4, ncol = 2)+
  plot_annotation(title = "% mitochondrial genes per spot - SCD-VI-i003",
                  theme = theme(plot.title = element_text(hjust = 0.5)))

ggsave("./output/images/SCD-VI-i003_analysis/SCD-VI-i003_mitochondrial_genes.png")



# 3 states ----
# Load and preprocess visium data
visium_subset <- readRDS("./output/processed_data/SCD-VI-i003/SCD-VI-i003_subset.rds")
regeneration_state <- read.csv("data/SCD-VI-i001-004/spaceranger_output/SCD-VI-i003/outs/regeneration_state_2.csv", row.names = "Barcode") # This file has to be created with loupe browser
regeneration_state$state <- factor(regeneration_state$state,
                                   levels = c('normal','border','injury'),ordered = TRUE)
visium_subset <- AddMetaData(visium_subset, regeneration_state, col.name = "regeneration_state")

SpatialDimPlot(visium_subset, group.by = "regeneration_state", pt.size.factor = 2.5) +
  theme(legend.position="right")
ggsave("./output/images/SCD-VI-i003_analysis/SCD-VI-i003_regeneration_3states.png")

## nCounts #########

# Statistics
regeneration_state <- bind_cols(regeneration_state, visium_subset[["nCount_Spatial"]])

# Shapiro-Wilk normality test for border nCount
with(regeneration_state, shapiro.test(nCount_Spatial[state == "border"])) # p = 0.007
# Shapiro-Wilk normality test for injury nCount
with(regeneration_state, shapiro.test(nCount_Spatial[state == "injury"])) # p = 0.024
# Shapiro-Wilk normality test for normal nCount
with(regeneration_state, shapiro.test(nCount_Spatial[state == "normal"])) # p < 2.2e-16
# --> we can't assume normality 
# --> Wilcoxon rank test

# Wilcoxon rank test
compare_means(nCount_Spatial ~ state, data = regeneration_state, method = "wilcox.test")

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

plot3 <- VlnPlot(visium_subset, group.by = "regeneration_state", features = "nCount_Spatial", pt.size = 0.1) +
  NoLegend() +
  theme(axis.title.x = element_blank(),
        plot.title = element_blank(),
        plot.margin = unit(c(0,0,0,0), units = "cm"))+
  ylab("nCount") +
  ylim(0, 60000) +
  stat_compare_means(label = "p.signif",
                     comparisons = list(c("normal", "injury"), c("border", "injury"), c("normal", "border")))

plot4 <- SpatialDimPlot(visium_subset, group.by = "regeneration_state", pt.size.factor = 2.5) +
  theme(legend.position="right")

wrap_plots(plot1, plot2, plot3, plot4, ncol = 2)+
  plot_annotation(title = "UMIs per spot - SCD-VI-i003",
                  theme = theme(plot.title = element_text(hjust = 0.5)))

ggsave("./output/images/SCD-VI-i003_analysis/SCD-VI-i003_ncount_3states.png")






## nFeature #########

# Statistics
regeneration_state <- bind_cols(regeneration_state, visium_subset[["nFeature_Spatial"]])

# Shapiro-Wilk normality test for border nFeature
with(regeneration_state, shapiro.test(nFeature_Spatial[state == "border"])) # p = 0.4785
# Shapiro-Wilk normality test for injury nFeature
with(regeneration_state, shapiro.test(nFeature_Spatial[state == "injury"])) # p = 0.3423
# Shapiro-Wilk normality test for normal nFeature
with(regeneration_state, shapiro.test(nFeature_Spatial[state == "normal"])) # p < 7.936e-13
# --> we can't assume normality 
# --> Wilcoxon rank test


# Wilcoxon rank test
compare_means(nFeature_Spatial ~ state, data = regeneration_state, method = "wilcox.test") 

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

plot3 <- VlnPlot(visium_subset, group.by = "regeneration_state", features = "nFeature_Spatial", pt.size = 0.1) +
  NoLegend() +
  theme(axis.title.x = element_blank(),
        plot.title = element_blank(),
        plot.margin = unit(c(0,0,0,0), units = "cm"))+
  ylab("nFeature") +
  ylim(0, 8000) +
  stat_compare_means(label = "p.signif",
                     comparisons = list(c("normal", "injury"), c("border", "injury"), c("normal", "border")))


plot4 <- SpatialDimPlot(visium_subset, group.by = "regeneration_state", pt.size.factor = 2.5) +
  theme(legend.position="right")

wrap_plots(plot1, plot2, plot3, plot4, ncol = 2)+
  plot_annotation(title = "Genes per spot - SCD-VI-i003",
                  theme = theme(plot.title = element_text(hjust = 0.5)))

ggsave("./output/images/SCD-VI-i003_analysis/SCD-VI-i003_ngenes_3states.png")

## % mitochondrial genes #########

visium_subset[["percent_mt_spatial"]] <- PercentageFeatureSet(visium_subset, pattern = "^mt-")

# Statistics
regeneration_state <- bind_cols(regeneration_state, visium_subset[["percent_mt_spatial"]])

# Shapiro-Wilk normality test for border percent_mt_spatial
with(regeneration_state, shapiro.test(percent_mt_spatial[state == "border"])) # p = 0.0115
# Shapiro-Wilk normality test for injury percent_mt_spatial
with(regeneration_state, shapiro.test(percent_mt_spatial[state == "injury"])) # p = 0.0012
# Shapiro-Wilk normality test for normal percent_mt_spatial
with(regeneration_state, shapiro.test(percent_mt_spatial[state == "normal"])) # p < 4.356e-06
# --> we can't assume normality 
# --> Wilcoxon rank test


# Wilcoxon rank test
compare_means(percent_mt_spatial ~ state, data = regeneration_state, method = "wilcox.test") 

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

plot3 <- VlnPlot(visium_subset, group.by = "regeneration_state", features = "percent_mt_spatial", pt.size = 0.1) +
  NoLegend() +
  theme(axis.title.x = element_blank(),
        plot.title = element_blank(),
        plot.margin = unit(c(0,0,0,0), units = "cm")) +
  ylab("% Mitochondrial") +
  ylim(0, 50) +
  stat_compare_means(label = "p.signif",
                     comparisons = list(c("normal", "injury"), c("border", "injury"), c("normal", "border")))


plot4 <- SpatialDimPlot(visium_subset, group.by = "regeneration_state", pt.size.factor = 2.5) +
  theme(legend.position="right")

wrap_plots(plot1, plot2, plot3, plot4, ncol = 2)+
  plot_annotation(title = "% mitochondrial genes per spot - SCD-VI-i003",
                  theme = theme(plot.title = element_text(hjust = 0.5)))

ggsave("./output/images/SCD-VI-i003_analysis/SCD-VI-i003_mitochondrial_genes_3states.png")


# Transformation, PCA, clustering ----
# SCTransform data
visium_subset <- SCTransform(visium_subset, assay = "Spatial", verbose = FALSE)

# Dimensionality reduction
visium_subset <- RunPCA(visium_subset, features = VariableFeatures(object = visium_subset))
visium_subset <- FindNeighbors(visium_subset, reduction = "pca", dims = 1:30)
visium_subset <- FindClusters(visium_subset, verbose = TRUE)
visium_subset <- RunUMAP(visium_subset, reduction = "pca", dims = 1:30)

SpatialDimPlot(visium_subset, label = TRUE)

saveRDS(visium_subset, "./output/processed_data/SCD-VI-i003/SCD-VI-i003_subset_sctransform_pca_clustered.rds")

# Plotting clustering some genes' expression
visium_subset <- readRDS("./output/processed_data/SCD-VI-i003/SCD-VI-i003_subset_sctransform_pca_clustered.rds")
DefaultAssay(visium_subset) <- "SCT"


SpatialFeaturePlot(visium_subset, features = c("nppa"))
ggsave("./output/images/SCD-VI-i003_analysis/SCT_spatial_nppa.png")

spatial_clustering <- SpatialDimPlot(visium_subset, label = FALSE, pt.size.factor = 2, cols=c("#ABA300","#00B8E7","#F8766D","#C77CFF","#00C19A"))+
  theme(legend.position="none")
clustering <- DimPlot(visium_subset, label = TRUE, label.size = 6, cols = c("#ABA300","#00B8E7","#F8766D","#C77CFF","#00C19A")) 
nCount <- FeaturePlot(visium_subset, features = "nCount_Spatial")
nFeature <- FeaturePlot(visium_subset, features = "nFeature_Spatial")
wrap_plots(clustering, spatial_clustering, nCount, nFeature, heights = c(1, 1), widths = c(1, 1))
ggsave("./output/images/SCD-VI-i003_analysis/SCD-VI-i003_SCT_clustering.png")
