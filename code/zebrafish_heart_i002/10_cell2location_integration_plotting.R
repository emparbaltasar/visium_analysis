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


visium <- Load10X_Spatial("./data/SCD-VI-i001-004/spaceranger_output/SCD-VI-i002/outs",
                          filename = "filtered_feature_bc_matrix.h5",
                          slice = "slice1",
                          filter.matrix = TRUE)

cell2loc_pred <- read.csv("./output/cell2location/zebrafish_heart/SCD-VI-i002/cell2location_map/cell2location_prediction.csv",
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
visium$nCells <- cell2loc_pred %>% t() %>% colSums()
SpatialFeaturePlot(visium, features = "nCells")

# Plot nCells per cell type
visium[["cell2loc"]] <- CreateAssayObject(t(cell2loc_pred))
DefaultAssay(visium) <- "cell2loc"

for (i in rownames(visium@assays$cell2loc@data)) {
  p <- SpatialFeaturePlot(visium, features = i, pt.size.factor = 2, crop = TRUE) +
    ggplot2::scale_fill_gradientn(colours = MyPalette)
  
  if (grepl("/", i)) {
    i <- sub("/", "-", i)
  }
  if (grepl(" ", i)) {
    i <- sub(" ", "_", i)
  }
  
  ggsave(paste("./output/images/SCD-VI-i002_integration/cell2location_prediction_abundance_", i, ".png", sep = ""))
}

# Plot proportions
cell2loc_pred[cell2loc_pred < 1] <- 0
cell2loc_pred_prop <- cell2loc_pred / rowSums(cell2loc_pred)
cell2loc_pred_prop[is.na(cell2loc_pred_prop)] <- 0
visium[["cell2loc_prop"]] <- CreateAssayObject(t(cell2loc_pred_prop))

# Plot spatial scatterpie
pal <- readRDS("./output/processed_data/SCD-VI-i004/spatial_scatterpie_palette.rds")

DefaultAssay(visium) <- "cell2loc_prop"

plotSpatialScatterpie(
  x = visium,
  y = cell2loc_pred_prop,
  cell_types = colnames(cell2loc_pred_prop),
  img = FALSE,
  scatterpie_alpha = 1,
  pie_scale = 0.8
) +
  scale_fill_manual(
    values = pal[names(pal) %in% names(colSums(cell2loc_pred_prop) != 0)[colSums(cell2loc_pred_prop) != 0]],
    breaks = names(pal)[names(pal) %in% names(colSums(cell2loc_pred_prop) != 0)[colSums(cell2loc_pred_prop) != 0]]
  ) +
  scale_y_reverse() +
  theme(legend.title = element_blank(), legend.text = element_text(size = 6)) +
  guides(fill = guide_legend(ncol = 1))

ggsave("./output/images/SCD-VI-i002_integration/cell2location_prediction_spatial_scatterpie.png")
ggsave("./output/images/SCD-VI-i002_integration/cell2location_prediction_spatial_scatterpie.pdf")

# Plot cell type correlation heatmap:
plotCorrelationMatrix(as.matrix(cell2loc_pred_prop),
                      insig = "pch",
                      pch.cex = 2)
ggsave("./output/images/SCD-VI-i002_integration/cell2location_correlation_matrix_proportions.png", width = 12, height = 12)
plotCorrelationMatrix(as.matrix(cell2loc_pred),
                      insig = "pch",
                      pch.cex = 2)
ggsave("./output/images/SCD-VI-i002_integration/cell2location_correlation_matrix_abundance.png", width = 12, height = 12)


#Cell type co-localization heatmap:
library(ggcorrplot)
cor_matrix = cell2loc_pred
cor_matrix[cor_matrix > 1] <- 1
cor_matrix[cor_matrix < 1] <- 0
cor_matrix <- cor_matrix[, colSums(cor_matrix) > 0] # Remove columns that are all 0
corr <- cor(cor_matrix, method = "pearson")
p.mat <- cor_pmat(corr)

ggcorrplot(
  corr = corr,
  p.mat = p.mat,
  hc.order = TRUE,
  insig = "pch",
  pch = 4,
  pch.cex = 2,
  lab = FALSE,
  colors = c("#6D9EC1", "white", "#E46726")) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 60, vjust = 1),
    axis.text = element_text(vjust = 0.5))
ggsave("./output/images/SCD-VI-i002_integration/cell2location_colocalization_heatmap.png", width = 12, height = 12)



# Keep only cell types with some significant value
sig <- rowSums(p.mat > 0.05) == ncol(p.mat) - 1
sig <- names(sig[sig == FALSE])

sig_corr <- corr[rownames(corr) %in% sig, colnames(corr) %in% sig]
sig_pmat <- p.mat[rownames(p.mat) %in% sig, colnames(p.mat) %in% sig]

ggcorrplot(
  corr = sig_corr,
  p.mat = sig_pmat,
  hc.order = TRUE,
  insig = "pch",
  pch = 4,
  pch.cex = 2,
  lab = FALSE,
  colors = c("#6D9EC1", "white", "#E46726")) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 60, vjust = 1),
    axis.text = element_text(vjust = 0.5),
    legend.position = "top")
ggsave("./output/images/SCD-VI-i002_integration/cell2location_colocalization_heatmap_significant.png", width = 7, height = 7)


# Jaccard similarity
jaccard = function(mat) {
  present = !is.na(mat) & mat > 0
  intersect = crossprod(present)
  union = nrow(mat) - crossprod(!present)
  J = intersect / union
  return(J)
}
jaccard <- jaccard(cor_matrix)

jaccard %>% 
  as.data.frame() %>%
  rownames_to_column("ct_1") %>%
  pivot_longer(-c(ct_1), names_to = "ct_2", values_to = "Jaccard") %>%
  ggplot(aes(x=ct_1, y=ct_2, fill=Jaccard)) + 
  geom_raster() + 
  scale_fill_distiller(palette = "Blues", direction = 1) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
    axis.text = element_text(vjust = 0.5)) +
  xlab("")+
  ylab("")

library(pheatmap)

png(file = "./output/images/SCD-VI-i002_integration/cell2location_colocalization_jaccard_heatmap.pdf", heigth = 820, width = 650)
pheatmap(jaccard,
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         angle_col = 45)
dev.off()


# Plot proportion per cell type

for (i in rownames(visium@assays$cell2loc@data)) {
  p <- SpatialFeaturePlot(visium, features = i, pt.size.factor = 2, crop = TRUE) +
    ggplot2::scale_fill_gradientn(colours = MyPalette, limits = c(0.0, 1.001), breaks = c(0.0, 0.5, 1.0))
  
  if (grepl("/", i)) {
    i <- sub("/", "-", i)
  }
  if (grepl(" ", i)) {
    i <- sub(" ", "_", i)
  }
  
  ggsave(paste("./output/images/SCD-VI-i002_integration/cell2location_prediction_proportion_", i, ".png", sep = ""))
}

# Plot some cell types with different legend scale and palette to adjust it for low proportions
NewPalette <- rev(colorRampPalette(colors = rev(x = brewer.pal(n = 9, name = "Reds")))(n = 100))

SpatialFeaturePlot(visium, features = "Fibroblast (mpeg1.1)", pt.size.factor = 2, crop = TRUE) +
  ggplot2::scale_fill_gradientn(colours = NewPalette, limits = c(0, 0.5), breaks = c(0.0, 0.25, 0.5))
ggsave("./output/images/SCD-VI-i002_integration/cell2location_prediction_proportion_adjusted_fibroblast_(mpeg1.1).png")

SpatialFeaturePlot(visium, features = "Fibroblast (col12a1a)", pt.size.factor = 2, crop = TRUE) +
  ggplot2::scale_fill_gradientn(colours = NewPalette, limits = c(0, 0.5), breaks = c(0.0, 0.25, 0.5))
ggsave("./output/images/SCD-VI-i002_integration/cell2location_prediction_proportion_adjusted_fibroblast_(col12a1a).png")

SpatialFeaturePlot(visium, features = "Fibroblast (cxcl12a)", pt.size.factor = 2, crop = TRUE) +
  ggplot2::scale_fill_gradientn(colours = NewPalette, limits = c(0, 0.5), breaks = c(0.0, 0.25, 0.5))
ggsave("./output/images/SCD-VI-i002_integration/cell2location_prediction_proportion_adjusted_fibroblast_(cxcl12a).png")



# Compare proportions of prediction with proportions in scrnaseq dataset
scrnaseq <- readRDS("./output/processed_data/SCD-VI-i001/zebrafish_heart_scrnaseq_lognorm_pca.rds")
scRNAseq_proportions <- table(scrnaseq@meta.data[["plot.ident2"]]) / length(scrnaseq@meta.data[["plot.ident2"]])

cell2loc_prediction <- visium@assays[["cell2loc"]]@data %>% t()

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


ggsave("./output/images/SCD-VI-i002_integration/cell2location_comparison_celltype_proportions_scatter.png", width = 10, height = 6)
