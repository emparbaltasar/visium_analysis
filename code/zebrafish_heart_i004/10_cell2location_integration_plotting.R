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


visium_subset <- readRDS("./output/processed_data/SCD-VI-i004/SCD-VI-i004_subset.rds")
cell2loc_pred <- read.csv("./output/cell2location/zebrafish_heart/SCD-VI-i004/cell2location_map/cell2location_prediction.csv",
                          row.names=1)

#cell2loc_pred_subset <- cell2loc_pred[rownames(cell2loc_pred) %in% colnames(visium_subset),]
#cell2loc_pred_subset_prop <- as_tibble(prop.table(as.matrix(cell2loc_pred_subset), 1), rownames = NA)
colnames(cell2loc_pred) <- c("B-cells","Bl.ves.EC (apnln)","Bl.ves.EC (lyve1)","Bl.ves.EC (plvapb)",
                                         "Cardiomyocytes (proliferating)", "Cardiomyocytes (ttn.2) A", "Cardiomyocytes (ttn.2) V", "Cardiomyocytes A",
                                         "Cardiomyocytes V", "Dead cells", "Endocardium (A)", "Endocardium (V)",
                                         "Endocardium (frzb)", "Epicardium (Atrium)", "Epicardium (Ventricle)", "Fibroblast",
                                         "Fibroblast (cfd)", "Fibroblast (col11a1a)", "Fibroblast (col12a1a)", "Fibroblast (cxcl12a)",
                                         "Fibroblast (mpeg1.1)","Fibroblast (nppc)", "Fibroblast (proliferating)", "Fibroblast (spock3)",
                                         "Fibroblast-like cells", "Macrophage (CM duplex)", "Macrophage (Endothelia duplex)", "Macrophage (Ery duplex)",
                                         "Macrophage (Fibroblast duplex)", "Macrophage (apoeb)", "Macrophage (cd59)", "Macrophage (epdl)",
                                         "Macrophage (il1b)", "Macrophage (proliferating)", "Macrophages", "Monocytes",
                                         "Myelin cells", "Neuronal cells", "Neutrophils", "Perivascular cells", 
                                         "Proliferating cells", "Smooth muscle cells", "T-cells", "T-cells (il4/13)",
                                         "T-cells (proliferating)")
visium_subset$nCells <- cell2loc_pred %>% t() %>% colSums()
SpatialFeaturePlot(visium_subset, features = "nCells")

# Plot nCells per cell type
visium_subset[["cell2loc"]] <- CreateAssayObject(t(cell2loc_pred))
DefaultAssay(visium_subset) <- "cell2loc"

for (i in rownames(visium_subset@assays$cell2loc@data)) {
  p <- SpatialFeaturePlot(visium_subset, features = i, pt.size.factor = 3, crop = TRUE) +
    ggplot2::scale_fill_gradientn(colours = MyPalette)
  
  if (grepl("/", i)) {
    i <- sub("/", "-", i)
  }
  if (grepl(" ", i)) {
    i <- sub(" ", "_", i)
  }
  
  ggsave(paste("./output/images/SCD-VI-i004_integration/cell2location_prediction_abundance_", i, ".png", sep = ""))
}

# Plot proportions
cell2loc_pred[cell2loc_pred < 1] <- 0
cell2loc_pred_prop <- cell2loc_pred / rowSums(cell2loc_pred)
cell2loc_pred_prop[is.na(cell2loc_pred_prop)] <- 0
visium_subset[["cell2loc"]] <- CreateAssayObject(t(cell2loc_pred_prop))

# Plot spatial scatterpie
pal <- readRDS("./output/processed_data/SCD-VI-i004/spatial_scatterpie_palette.rds")

DefaultAssay(visium_subset) <- "cell2loc"

plotSpatialScatterpie(
  x = visium_subset,
  y = cell2loc_pred_prop,
  cell_types = colnames(cell2loc_pred_prop),
  img = FALSE,
  scatterpie_alpha = 1,
  pie_scale = 1.2
) +
  scale_fill_manual(
    values = pal[names(pal) %in% names(colSums(cell2loc_pred_prop) != 0)[colSums(cell2loc_pred_prop) != 0]],
    breaks = names(pal)[names(pal) %in% names(colSums(cell2loc_pred_prop) != 0)[colSums(cell2loc_pred_prop) != 0]]
  ) +
  scale_y_reverse() +
  theme(legend.title = element_blank(), legend.text = element_text(size = 6)) +
  guides(fill = guide_legend(ncol = 1))

ggsave("./output/images/SCD-VI-i004_integration/cell2location_prediction_spatial_scatterpie.png")
ggsave("./output/images/SCD-VI-i004_integration/cell2location_prediction_spatial_scatterpie.pdf")


# Plot proportion per cell type

for (i in rownames(visium_subset@assays$cell2loc@data)) {
  p <- SpatialFeaturePlot(visium_subset, features = i, pt.size.factor = 4, crop = TRUE) +
    ggplot2::scale_fill_gradientn(colours = MyPalette, limits = c(0.0, 1.001), breaks = c(0.0, 0.5, 1.0))
  
  if (grepl("/", i)) {
    i <- sub("/", "-", i)
  }
  if (grepl(" ", i)) {
    i <- sub(" ", "_", i)
  }
  
  ggsave(paste("./output/images/SCD-VI-i004_integration/cell2location_prediction_proportion_", i, ".png", sep = ""))
}

# Compare proportions of prediction with proportions in scrnaseq dataset
scrnaseq <- readRDS("./output/processed_data/SCD-VI-i004/zebrafish_heart_scrnaseq_lognorm_pca.rds")
scRNAseq_proportions <- table(scrnaseq@meta.data[["plot.ident2"]]) / length(scrnaseq@meta.data[["plot.ident2"]])

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
  stat_cor(method = "pearson", label.x = 0.08, label.y = 0.20, size = 3) +
  stat_cor(method = "spearman", label.x = 0.08, label.y = 0.19, cor.coef.name = "rho", size = 3)


ggsave("./output/images/SCD-VI-i004_integration/cell2location_comparison_celltype_proportions_scatter.png", width = 10, height = 6)
