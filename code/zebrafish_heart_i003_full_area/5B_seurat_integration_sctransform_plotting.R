library(Seurat)
library(RColorBrewer)
library(scater)
library(SPOTlight)
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(patchwork)

MyPalette <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(n = 100)

visium_subset <- readRDS("./output/processed_data/SCD-VI-i003_full_area/SCD-VI-i003_full_area_sctransform_seurat_predictions.rds")
DefaultAssay(visium_subset) <- "prediction.score.id"

for (i in rownames(visium_subset@assays$prediction.score.id@data)) {
  p <- SpatialFeaturePlot(visium_subset, features = i, pt.size.factor = 1, crop = TRUE) +
    ggplot2::scale_fill_gradientn(colours = MyPalette, limits = c(0.0, 1.001), breaks = c(0.0, 0.5, 1.0))

  if (grepl("/", i)) {
    i <- sub("/", "-", i)
  }
  if (grepl(" ", i)) {
    i <- sub(" ", "_", i)
  }
  dir.create("./output/images/SCD-VI-i003_full_area_integration")
  ggsave(paste("./output/images/SCD-VI-i003_full_area_integration/seurat_sctransform_prediction_", i, ".png", sep = ""))
}

# Plot spatial scatterpie
seurat_prediction <- visium_subset@assays[["prediction.score.id"]]@data %>% t()

pal <- readRDS("./output/processed_data/SCD-VI-i004/spatial_scatterpie_palette.rds")

plotSpatialScatterpie(
  x = visium_subset,
  y = seurat_prediction,
  pie_scale = 0.3
) +
  scale_fill_manual(
    values = pal[names(pal) %in% names(colSums(seurat_prediction) != 0)[colSums(seurat_prediction) != 0]],
    breaks = names(pal)[names(pal) %in% names(colSums(seurat_prediction) != 0)[colSums(seurat_prediction) != 0]]
  ) +
  scale_y_reverse() +
  theme(legend.title = element_blank(), legend.text = element_text(size = 6)) +
  guides(fill = guide_legend(ncol = 1))

dir.create("./output/images/SCD-VI-i003_full_area_integration")
ggsave("./output/images/SCD-VI-i003_full_area_integration/seurat_sctransform_prediction_spatial_scatterpie.png")
ggsave("./output/images/SCD-VI-i003_full_area_integration/seurat_sctransform_prediction_spatial_scatterpie.pdf")


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
  stat_cor(method = "pearson", label.x = 0.08, label.y = 0.2, size = 4) +
  stat_cor(method = "spearman", label.x = 0.08, label.y = 0.19, cor.coef.name = "rho", size = 4)


ggsave("./output/images/SCD-VI-i003_full_area_integration/seurat_sctransform_comparison_celltype_proportions_scatter.png", width = 10, height = 6)


# Plot cell type correlation heatmap:
seurat_prediction <- visium_subset@assays[["prediction.score.id"]]@data %>% t()
plotCorrelationMatrix(seurat_prediction) +
  theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7))
ggsave("./output/images/SCD-VI-i003_full_area_integration/seurat_sctransform_celltype_colocalization.png")



