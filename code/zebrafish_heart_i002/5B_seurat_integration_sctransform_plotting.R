library(Seurat)
library(RColorBrewer)
library(scater)
library(SPOTlight)
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(patchwork)

MyPalette <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(n = 100)

visium <- readRDS("./output/processed_data/SCD-VI-i002/SCD-VI-i002_sctransform_seurat_predictions.rds")
DefaultAssay(visium) <- "prediction.score.id"

for (i in rownames(visium@assays$prediction.score.id@data)) {
  p <- SpatialFeaturePlot(visium, features = i, pt.size.factor = 2, crop = TRUE) +
    ggplot2::scale_fill_gradientn(colours = MyPalette, limits = c(0.0, 1.001), breaks = c(0.0, 0.5, 1.0))

  if (grepl("/", i)) {
    i <- sub("/", "-", i)
  }
  if (grepl(" ", i)) {
    i <- sub(" ", "_", i)
  }
  dir.create("./output/images/SCD-VI-i002_integration")
  ggsave(paste("./output/images/SCD-VI-i002_integration/seurat_sctransform_prediction_", i, ".png", sep = ""))
}

# Plot some cell types with different legend scale and palette to adjust it for low proportions
NewPalette <- rev(colorRampPalette(colors = rev(x = brewer.pal(n = 9, name = "Reds")))(n = 100))

SpatialFeaturePlot(visium, features = "Fibroblast (mpeg1.1)", pt.size.factor = 2, crop = TRUE) +
  ggplot2::scale_fill_gradientn(colours = NewPalette, limits = c(0, 0.5), breaks = c(0.0, 0.25, 0.5))
ggsave("./output/images/SCD-VI-i002_integration/seurat_sctransform_prediction_proportion_adjusted_fibroblast_(mpeg1.1).png")

SpatialFeaturePlot(visium, features = "Fibroblast (col12a1a)", pt.size.factor = 2, crop = TRUE) +
  ggplot2::scale_fill_gradientn(colours = NewPalette, limits = c(0, 0.5), breaks = c(0.0, 0.25, 0.5))
ggsave("./output/images/SCD-VI-i002_integration/seurat_sctransform_prediction_proportion_adjusted_fibroblast_(col12a1a).png")

SpatialFeaturePlot(visium, features = "Fibroblast (cxcl12a)", pt.size.factor = 2, crop = TRUE) +
  ggplot2::scale_fill_gradientn(colours = NewPalette, limits = c(0, 0.5), breaks = c(0.0, 0.25, 0.5))
ggsave("./output/images/SCD-VI-i002_integration/seurat_sctransform_prediction_proportion_adjusted_fibroblast_(cxcl12a).png")







# Plot spatial scatterpie
seurat_prediction <- visium@assays[["prediction.score.id"]]@data %>% t()

pal <- readRDS("./output/processed_data/SCD-VI-i004/spatial_scatterpie_palette.rds")

plotSpatialScatterpie(
  x = visium,
  y = seurat_prediction,
  pie_scale = 0.8
  ) +
  scale_fill_manual(
    values = pal[names(pal) %in% names(colSums(seurat_prediction) != 0)[colSums(seurat_prediction) != 0]],
    breaks = names(pal)[names(pal) %in% names(colSums(seurat_prediction) != 0)[colSums(seurat_prediction) != 0]]
  ) +
  scale_y_reverse() +
  theme(legend.title = element_blank(), legend.text = element_text(size = 6)) +
  guides(fill = guide_legend(ncol = 1))

ggsave("./output/images/SCD-VI-i002_integration/seurat_sctransform_prediction_spatial_scatterpie.png")
ggsave("./output/images/SCD-VI-i002_integration/seurat_sctransform_prediction_spatial_scatterpie.pdf")

# Plot spatial scatterpie with 0.04 threshold (1/25)
seurat_prediction_filtered <- seurat_prediction
seurat_prediction_filtered[seurat_prediction_filtered < 0.04] <- 0

plotSpatialScatterpie(
  x = visium,
  y = seurat_prediction_filtered,
  cell_types = colnames(seurat_prediction),
  pie_scale = 0.8
) +
  scale_fill_manual(
    values = pal[names(pal) %in% names(colSums(seurat_prediction) != 0)[colSums(seurat_prediction) != 0]],
    breaks = names(pal)[names(pal) %in% names(colSums(seurat_prediction) != 0)[colSums(seurat_prediction) != 0]]
  ) +
  scale_y_reverse() +
  theme(legend.title = element_blank(), legend.text = element_text(size = 6)) +
  guides(fill = guide_legend(ncol = 1))

ggsave("./output/images/SCD-VI-i002_integration/seurat_sctransform_prediction_filtered_spatial_scatterpie.png")
ggsave("./output/images/SCD-VI-i002_integration/seurat_sctransform_prediction_filtered_spatial_scatterpie.pdf")


# Plot cell type correlation heatmap:
seurat_prediction <- visium@assays[["prediction.score.id"]]@data %>% t()
plotCorrelationMatrix(seurat_prediction)
ggsave("./output/images/SCD-VI-i002_integration/seurat_sctransform_celltype_colocalization.png")


# Plot Shannon index
seurat_prediction <- visium@assays[["prediction.score.id"]]@data

seurat_prediction/colSums(seurat_prediction)
p = seurat_prediction
-colSums(p*log(p), na.rm = TRUE)

visium[["shannon_index"]] <- -colSums(p*log(p), na.rm = TRUE)

SpatialFeaturePlot(visium, features = "shannon_index", pt.size.factor = 2.5) + 
  theme(legend.position = "right",
        plot.margin = unit(c(0,0,0,0), units = "cm"))
ggsave("./output/images/SCD-VI-i002_integration/seurat_sctransform_shannon_index.png")


