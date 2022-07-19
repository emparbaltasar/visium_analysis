# Plotting cell2location results
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(scater)
library(SPOTlight)
library(tidyverse)
library(ggpubr)

MyPalette <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(n = 100)

palettePaired <- c(
  "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c",
  "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#fff116",
  "#d8a283", "#b15928", "#b0c288", "#81b806", "#bcffde", "#00ff82",
  "#e7aee0", "#a40192", "#d77070", "#ff0052", "#c6bfbf", "#6a5b5b"
)

visium_subset <- readRDS("./output/processed_data/mouse_brain/mouse_brain_spatial_cortex.rds")
cell2loc_pred <- read.csv("./output/cell2location/mouse_brain/cell2location_map/cell2location_prediction.csv",
                          row.names=1)

cell2loc_pred_subset <- cell2loc_pred[rownames(cell2loc_pred) %in% colnames(visium_subset),]
cell2loc_pred_subset_prop <- as_tibble(prop.table(as.matrix(cell2loc_pred_subset), 1), rownames = NA)
colnames(cell2loc_pred_subset_prop)[4] <- "L2/3 IT"
colnames(cell2loc_pred_subset_prop)[6] <- "L5 IT"
colnames(cell2loc_pred_subset_prop)[7] <- "L5 PT"
colnames(cell2loc_pred_subset_prop)[8] <- "L6 CT"
colnames(cell2loc_pred_subset_prop)[9] <- "L6 IT"

visium_subset[["cell2loc"]] <- CreateAssayObject(t(cell2loc_pred_subset_prop))

# Plot spatial scatterpie
pal <- readRDS("./output/processed_data/mouse_brain/spatial_scatterpie_palette.rds")
pal <- pal[names(pal) != "max"]

DefaultAssay(visium_subset) <- "cell2loc"

plotSpatialScatterpie(
  x = visium_subset,
  y = cell2loc_pred_subset_prop,
  cell_types = colnames(cell2loc_pred_subset_prop),
  img = FALSE,
  scatterpie_alpha = 1,
  pie_scale = 0.4
) +
  scale_fill_manual(
    values = pal,
    breaks = names(pal)
  ) +
  scale_y_reverse()

ggsave("./output/images/mouse_brain_integration/cell2location_prediction_spatial_scatterpie.png")
ggsave("./output/images/mouse_brain_integration/cell2location_prediction_spatial_scatterpie.pdf")


# Plot proportions
cell2loc_pred_subset[cell2loc_pred_subset < 1] <- 0
cell2loc_pred_subset_prop <- cell2loc_pred_subset / rowSums(cell2loc_pred_subset)
cell2loc_pred_subset_prop[is.na(cell2loc_pred_subset_prop)] <- 0
colnames(cell2loc_pred_subset_prop)[4] <- "L2/3 IT"
colnames(cell2loc_pred_subset_prop)[6] <- "L5 IT"
colnames(cell2loc_pred_subset_prop)[7] <- "L5 PT"
colnames(cell2loc_pred_subset_prop)[8] <- "L6 CT"
colnames(cell2loc_pred_subset_prop)[9] <- "L6 IT"

visium_subset[["cell2loc"]] <- CreateAssayObject(t(cell2loc_pred_subset_prop))

# Plot spatial scatterpie
pal <- readRDS("./output/processed_data/mouse_brain/spatial_scatterpie_palette.rds")

DefaultAssay(visium_subset) <- "cell2loc"

plotSpatialScatterpie(
  x = visium_subset,
  y = cell2loc_pred_subset_prop,
  cell_types = colnames(cell2loc_pred_subset_prop),
  img = FALSE,
  scatterpie_alpha = 1,
  pie_scale = 0.4
) +
  scale_fill_manual(
    values = pal[names(pal) %in% names(colSums(cell2loc_pred_subset_prop) != 0)[colSums(cell2loc_pred_subset_prop) != 0]],
    breaks = names(pal)[names(pal) %in% names(colSums(cell2loc_pred_subset_prop) != 0)[colSums(cell2loc_pred_subset_prop) != 0]]
  ) +
  scale_y_reverse() +
  theme(legend.title = element_blank(), legend.text = element_text(size = 6)) +
  guides(fill = guide_legend(ncol = 1))

ggsave("./output/images/mouse_brain_integration/cell2location_prediction_filtered_spatial_scatterpie.png")
ggsave("./output/images/mouse_brain_integration/cell2location_prediction_filtered_spatial_scatterpie.pdf")

ggsave("./output/images/SCD-VI-i004_integration/cell2location_prediction_spatial_scatterpie.png")
ggsave("./output/images/SCD-VI-i004_integration/cell2location_prediction_spatial_scatterpie.pdf")


# Plot per cell type

for (i in rownames(visium_subset@assays$cell2loc@data)) {
  p <- SpatialFeaturePlot(visium_subset, features = i, pt.size.factor = 1, crop = FALSE) +
    ggplot2::scale_fill_gradientn(colours = MyPalette, limits = c(0.0, 1.001), breaks = c(0.0, 0.5, 1.0))
  
  if (grepl("/", i)) {
    i <- sub("/", "-", i)
  }
  if (grepl(" ", i)) {
    i <- sub(" ", "_", i)
  }
  
  ggsave(paste("./output/images/mouse_brain_integration/cell2location_prediction_proportion_", i, ".png", sep = ""))
}

# Compare proportions of prediction with proportions in scrnaseq dataset
scrnaseq <- readRDS("./output/processed_data/mouse_brain/mouse_brain_scrnaseq_lognorm_pca.rds")
scRNAseq_proportions <- table(scrnaseq@meta.data[["subclass"]]) / length(scrnaseq@meta.data[["subclass"]])

cell2loc_prediction <- visium_subset@assays[["cell2loc"]]@data %>% t()

cell2loc_prediction <- colMeans(cell2loc_prediction)

proportions_comparison <- rbind(cell2loc_prediction, scRNAseq_proportions)
proportions_comparison <- as.data.frame(proportions_comparison) %>%
  t() %>%
  as.data.frame()


proportions_comparison %>%
  rownames_to_column("cell_type") %>%
  ggplot(aes(x = scRNAseq_proportions, y = cell2loc_prediction)) +
  scale_color_manual(values = pal) +
  theme_bw() +
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank(),
    aspect.ratio = 1 / 1, legend.text = element_text(size = 6)
  ) +
  geom_point(aes(color = cell_type)) +
  stat_cor(method = "pearson", label.x = 0.10, label.y = 0.20, size = 3) +
  stat_cor(method = "spearman", label.x = 0.10, label.y = 0.19, cor.coef.name = "rho", size = 3)


ggsave("./output/images/mouse_brain_integration/cell2location_comparison_celltype_proportions_scatter.png", width = 5, height = 4)
