library(ggplot2)
library(Seurat)

scrnaseq <- readRDS("./output/processed_data/SCD-VI-i001/zebrafish_heart_scrnaseq_sctransform_pca.rds")

# Plot UMAP
Idents(scrnaseq) <- "plot.ident2"
pal <- readRDS("./output/processed_data/SCD-VI-i004/spatial_scatterpie_palette.rds")

DimPlot(scrnaseq, cols = pal) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 8)) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  guides(colour=guide_legend(override.aes = list(size = 3), ncol=2))

ggsave("./output/images/SCD-VI-i001_analysis/scrnaseq_ctrl_7dpi_umap_celltypes.png")
ggsave("./output/images/SCD-VI-i001_analysis/scrnaseq_ctrl_7dpi_umap_celltypes.pdf")
