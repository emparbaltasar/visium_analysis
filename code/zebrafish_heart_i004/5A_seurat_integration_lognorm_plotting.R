library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(scater)
library(SPOTlight)
library(tidyverse)
library(ggpubr)

MyPalette <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(n = 100)

visium_subset <- readRDS("./output/processed_data/SCD-VI-i004/SCD-VI-i004_lognorm_seurat_predictions.rds")

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
  ggsave(paste("./output/images/SCD-VI-i004_integration/seurat_lognorm_prediction_", i, ".png", sep = ""))
}

# Plot spatial scatterpie
seurat_prediction <- visium_subset@assays[["prediction.score.id"]]@data %>% t()

# ct <- colnames(seurat_prediction)
# ct <- sort(ct)
# pal <- c("B-cells" = "grey30",
#          "Bl.ves.EC (apnln)" = "sienna4",
#          "Bl.ves.EC (lyve1)" = "sienna2",
#          "Bl.ves.EC (plvapb)" = "tan",
#          "Cardiomyocytes (proliferating)" = "#FF0000",
#          "Cardiomyocytes (ttn.2) A" = "#A50A0A",
#          "Cardiomyocytes (ttn.2) V" = "#FF9279",
#          "Cardiomyocytes A" = "#E21149",
#          "Cardiomyocytes V" = "#FD9090",
#          "Dead cells" = "#000000",
#          "Endocardium (A)" = "#1500FF",
#          "Endocardium (frzb)" = "#B1AAFD",
#          "Endocardium (V)" = "#D2E0FF",
#          "Epicardium (Atrium)" = "#5BA5C6",
#          [15] "Epicardium (Ventricle)"         "Fibroblast"                    
#          [17] "Fibroblast (cfd)"               "Fibroblast (col11a1a)"         
#          [19] "Fibroblast (col12a1a)"          "Fibroblast (cxcl12a)"          
#          [21] "Fibroblast (mpeg1.1)"           "Fibroblast (nppc)"             
#          [23] "Fibroblast (proliferating)"     "Fibroblast (spock3)"           
#          [25] "Fibroblast-like cells"          "Macrophage (apoeb)"            
#          [27] "Macrophage (cd59)"              "Macrophage (CM duplex)"        
#          [29] "Macrophage (Endothelia duplex)" "Macrophage (epdl)"             
#          [31] "Macrophage (Ery duplex)"        "Macrophage (Fibroblast duplex)"
#          [33] "Macrophage (il1b)"              "Macrophage (proliferating)"    
#          [35] "Macrophages"                    "Monocytes"                     
#          [37] "Myelin cells"                   "Neuronal cells"                
#          [39] "Neutrophils"                    "Perivascular cells"            
#          [41] "Proliferating cells"            "Smooth muscle cells"           
#          [43] "T-cells"                        "T-cells (il4/13)"              
#          [45] "T-cells (proliferating)"       
#          )

# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# pal <- sample(col_vector, length(ct))
# names(pal) <- ct
# saveRDS(pal, "./output/processed_data/SCD-VI-i004/spatial_scatterpie_palette.rds")
pal <- readRDS("./output/processed_data/SCD-VI-i004/spatial_scatterpie_palette.rds")

plotSpatialScatterpie(
  x = visium_subset,
  y = seurat_prediction,
  cell_types = colnames(seurat_prediction),
  pie_scale = 1
) +
  scale_fill_manual(
    values = pal[names(pal) %in% names(colSums(seurat_prediction) != 0)[colSums(seurat_prediction) != 0]],
    breaks = names(pal)[names(pal) %in% names(colSums(seurat_prediction) != 0)[colSums(seurat_prediction) != 0]]
  ) +
  scale_y_reverse() +
  theme(legend.title = element_blank(), legend.text = element_text(size = 6)) +
  guides(fill = guide_legend(ncol = 1))

ggsave("./output/images/SCD-VI-i004_integration/seurat_lognorm_prediction_spatial_scatterpie.png")
ggsave("./output/images/SCD-VI-i004_integration/seurat_lognorm_prediction_spatial_scatterpie.pdf")

# Plot spatial scatterpie with 0.04 threshold
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

ggsave("./output/images/SCD-VI-i004_integration/seurat_lognorm_prediction_filtered_spatial_scatterpie.png")
ggsave("./output/images/SCD-VI-i004_integration/seurat_lognorm_prediction_filtered_spatial_scatterpie.pdf")



# Compare proportions of prediction with proportions in scrnaseq dataset
scrnaseq <- readRDS("./output/processed_data/SCD-VI-i004/zebrafish_heart_scrnaseq_lognorm_pca.rds")
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
  stat_cor(method = "pearson", label.x = 0.06, label.y = 0.8, size = 3) +
  stat_cor(method = "spearman", label.x = 0.06, label.y = 0.75, cor.coef.name = "rho", size = 3)


ggsave("./output/images/SCD-VI-i004_integration/seurat_lognorm_comparison_celltype_proportions_scatter.png", width = 10, height = 6)
