library(Seurat)
library(RColorBrewer)
library(scater)
library(SPOTlight)
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(patchwork)

MyPalette <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(n = 100)

visium_subset <- readRDS("./output/processed_data/SCD-VI-i004/SCD-VI-i004_sctransform_seurat_predictions.rds")
DefaultAssay(visium_subset) <- "prediction.score.id"

for (i in rownames(visium_subset@assays$prediction.score.id@data)) {
  p <- SpatialFeaturePlot(visium_subset, features = i, pt.size.factor = 3, crop = TRUE) +
    ggplot2::scale_fill_gradientn(colours = MyPalette, limits = c(0.0, 1.001), breaks = c(0.0, 0.5, 1.0))

  if (grepl("/", i)) {
    i <- sub("/", "-", i)
  }
  if (grepl(" ", i)) {
    i <- sub(" ", "_", i)
  }
  dir.create("./output/images/SCD-VI-i004_integration")
  ggsave(paste("./output/images/SCD-VI-i004_integration/seurat_sctransform_prediction_", i, ".png", sep = ""))
}

# Plot spatial scatterpie
seurat_prediction <- visium_subset@assays[["prediction.score.id"]]@data %>% t()

pal <- readRDS("./output/processed_data/SCD-VI-i004/spatial_scatterpie_palette.rds")

plotSpatialScatterpie(
  x = visium_subset,
  y = seurat_prediction,
  pie_scale = 1
) +
  scale_fill_manual(
    values = pal[names(pal) %in% names(colSums(seurat_prediction) != 0)[colSums(seurat_prediction) != 0]],
    breaks = names(pal)[names(pal) %in% names(colSums(seurat_prediction) != 0)[colSums(seurat_prediction) != 0]]
  ) +
  scale_y_reverse() +
  theme(legend.title = element_blank(), legend.text = element_text(size = 6)) +
  guides(fill = guide_legend(ncol = 1))

ggsave("./output/images/SCD-VI-i004_integration/seurat_sctransform_prediction_spatial_scatterpie.png")
ggsave("./output/images/SCD-VI-i004_integration/seurat_sctransform_prediction_spatial_scatterpie.pdf")

# Plot spatial scatterpie with 0.04 threshold (1/25)
seurat_prediction_filtered <- seurat_prediction
seurat_prediction_filtered[seurat_prediction_filtered < 0.04] <- 0

plotSpatialScatterpie(
  x = visium_subset,
  y = seurat_prediction_filtered,
  cell_types = colnames(seurat_prediction),
  pie_scale = 1.2
) +
  scale_fill_manual(
    values = pal[names(pal) %in% names(colSums(seurat_prediction) != 0)[colSums(seurat_prediction) != 0]],
    breaks = names(pal)[names(pal) %in% names(colSums(seurat_prediction) != 0)[colSums(seurat_prediction) != 0]]
  ) +
  scale_y_reverse() +
  theme(legend.title = element_blank(), legend.text = element_text(size = 6)) +
  guides(fill = guide_legend(ncol = 1))

ggsave("./output/images/SCD-VI-i004_integration/seurat_sctransform_prediction_filtered_spatial_scatterpie.png")
ggsave("./output/images/SCD-VI-i004_integration/seurat_sctransform_prediction_filtered_spatial_scatterpie.pdf")

# Compare proportions of prediction with proportions in scrnaseq dataset
scrnaseq <- readRDS("./output/processed_data/SCD-VI-i004/zebrafish_heart_scrnaseq_sctransform_pca.rds")
scRNAseq_proportions <- table(scrnaseq@meta.data[["plot.ident2"]]) / length(scrnaseq@meta.data[["plot.ident2"]])

seurat_prediction <- visium_subset@assays[["prediction.score.id"]]@data %>% t()

seurat_prediction <- colMeans(seurat_prediction)

proportions_comparison <- rbind(seurat_prediction, scRNAseq_proportions)
proportions_comparison <- as.data.frame(proportions_comparison) %>%
  t() %>%
  as.data.frame()


proportions_comparison %>%
  rownames_to_column("cell_type") %>%
  ggplot(aes(x = scRNAseq_proportions, y = seurat_prediction)) +
  scale_color_manual(values = pal) +
  theme_bw() +
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank(),
    aspect.ratio = 1 / 1, legend.text = element_text(size = 6)
  ) +
  geom_point(aes(color = cell_type)) +
  stat_cor(method = "pearson", label.x = 0.08, label.y = 0.45, size = 3) +
  stat_cor(method = "spearman", label.x = 0.08, label.y = 0.43, cor.coef.name = "rho", size = 3)


ggsave("./output/images/SCD-VI-i004_integration/seurat_sctransform_comparison_celltype_proportions_scatter.png", width = 10, height = 6)


# Plot cell type correlation heatmap:
seurat_prediction <- visium_subset@assays[["prediction.score.id"]]@data %>% t()
plotCorrelationMatrix(seurat_prediction)
ggsave("./output/images/SCD-VI-i004_integration/seurat_sctransform_celltype_colocalization.png")


# Plot Shannon index
seurat_prediction <- visium_subset@assays[["prediction.score.id"]]@data

seurat_prediction/colSums(seurat_prediction)
p = seurat_prediction
-colSums(p*log(p), na.rm = TRUE)

visium_subset[["shannon_index"]] <- -colSums(p*log(p), na.rm = TRUE)

SpatialFeaturePlot(visium_subset, features = "shannon_index", pt.size.factor = 2.5) + 
  theme(legend.position = "right",
        plot.margin = unit(c(0,0,0,0), units = "cm"))
ggsave("./output/images/SCD-VI-i004_integration/seurat_sctransform_shannon_index.png")


regeneration_state <- read.csv("data/SCD-VI-i001-004/spaceranger_output/SCD-VI-i004/outs/regeneration_state.csv", row.names = "Barcode")
visium_subset <- AddMetaData(visium_subset, regeneration_state, col.name = "regeneration_state")


# Statistics
regeneration_state <- bind_cols(regeneration_state, visium_subset[["shannon_index"]])

# Shapiro-Wilk normality test for regenerating nCount
with(regeneration_state, shapiro.test(shannon_index[state == "regenerating"]))# p = 0.3091
# Shapiro-Wilk normality test for non-regenerating nCount
with(regeneration_state, shapiro.test(shannon_index[state == "non-regenerating"])) # p = 1.253e-05
# --> we can't assume normality 
# --> Wilcoxon rank test

# Wilcoxon rank test
compare_means(shannon_index ~ state, data = regeneration_state, method = "wilcox.test") # p = 1,57e-9 --> significantly different


# Plot nCounts
Idents(visium_subset) <- "orig.ident"
  
plot1 <- VlnPlot(visium_subset, features = "shannon_index", pt.size = 0.1) +
  NoLegend() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_blank(),
        plot.margin = unit(c(0,0,0,0), units = "cm")) +
  ylab("Shannon_index")

plot2 <- SpatialFeaturePlot(visium_subset, features = "shannon_index", pt.size.factor = 2.5) + 
  theme(legend.position = "right",
        plot.margin = unit(c(0,0,0,0), units = "cm"))

plot3 <- VlnPlot(visium_subset, group.by = "regeneration_state", features = "shannon_index", pt.size = 0.1) +
  NoLegend() +
  theme(axis.title.x = element_blank(),
        plot.title = element_blank(),
        plot.margin = unit(c(0,0,0,0), units = "cm"))+
  ylab("Shanon index") +
  stat_compare_means()

plot4 <- SpatialDimPlot(visium_subset, group.by = "regeneration_state", pt.size.factor = 2.5) +
  theme(legend.position="right")

wrap_plots(plot1, plot2, plot3, plot4, ncol = 2)+
  plot_annotation(title = "Shanon index",
                  theme = theme(plot.title = element_text(hjust = 0.5)))

ggsave("./output/images/SCD-VI-i004_integration/seurat_sctransform_shannon_index.png")


# Plot Shannon index comparing 3 states (normal, border, injury)
seurat_prediction <- visium_subset@assays[["prediction.score.id"]]@data

seurat_prediction/colSums(seurat_prediction)
p = seurat_prediction
-colSums(p*log(p), na.rm = TRUE)

visium_subset[["shannon_index"]] <- -colSums(p*log(p), na.rm = TRUE)

SpatialFeaturePlot(visium_subset, features = "shannon_index", pt.size.factor = 2.5) + 
  theme(legend.position = "right",
        plot.margin = unit(c(0,0,0,0), units = "cm"))


regeneration_state <- read.csv("data/SCD-VI-i001-004/spaceranger_output/SCD-VI-i004/outs/regeneration_3states.csv", row.names = "Barcode")
visium_subset <- AddMetaData(visium_subset, regeneration_state, col.name = "regeneration_state")


# Statistics
regeneration_state <- bind_cols(regeneration_state, visium_subset[["shannon_index"]])

# Shapiro-Wilk normality test 
with(regeneration_state, shapiro.test(shannon_index[state == "normal"]))# p = 3.162e-05
with(regeneration_state, shapiro.test(shannon_index[state == "border"])) # p = 0.08785
with(regeneration_state, shapiro.test(shannon_index[state == "injury"])) # p = 0.07797
# --> we can't assume normality 
# --> Wilcoxon rank test

# Wilcoxon rank test
compare_means(shannon_index ~ state, data = regeneration_state, method = "wilcox.test") 

# Plot
Idents(visium_subset) <- "orig.ident"

plot1 <- VlnPlot(visium_subset, features = "shannon_index", pt.size = 0.1) +
  NoLegend() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_blank(),
        plot.margin = unit(c(0,0,0,0), units = "cm")) +
  ylab("Shannon_index")

plot2 <- SpatialFeaturePlot(visium_subset, features = "shannon_index", pt.size.factor = 2.5) + 
  theme(legend.position = "right",
        plot.margin = unit(c(0,0,0,0), units = "cm"))

plot3 <- VlnPlot(visium_subset, group.by = "regeneration_state", features = "shannon_index", pt.size = 0.1) +
  NoLegend() +
  theme(axis.title.x = element_blank(),
        plot.title = element_blank(),
        plot.margin = unit(c(0,0,0,0), units = "cm"))+
  ylab("Shanon index") +
  ylim(0.5,2.6) +
  stat_compare_means(label = "p.signif",
                     comparisons = list(c("normal", "injury"), c("border", "injury"), c("normal", "border")))

plot4 <- SpatialDimPlot(visium_subset, group.by = "regeneration_state", pt.size.factor = 2.5) +
  theme(legend.position="right")

wrap_plots(plot1, plot2, plot3, plot4, ncol = 2)+
  plot_annotation(title = "Shanon index",
                  theme = theme(plot.title = element_text(hjust = 0.5)))

ggsave("./output/images/SCD-VI-i004_integration/seurat_sctransform_shannon_index_3states.png")
