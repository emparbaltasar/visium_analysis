library(Seurat)
library(tidyverse)

scrnaseq <- readRDS("./output/processed_data/SCD-VI-i001/zebrafish_heart_scrnaseq_lognorm_pca.rds")

Idents(scrnaseq) <- "plot.ident2"
macrophages <- subset(scrnaseq, idents = str_subset(unique(Idents(scrnaseq)), "Macrophag."))

markers <- FindAllMarkers(macrophages)
write.csv(markers, "./output/processed_data/markergenes_macrophages.csv")

markers_pct0.5 <- FindAllMarkers(macrophages, min.pct = 0.25)

markers_apoeb <- markers[markers$cluster == "Macrophage (apoeb)",]
markers_il1b <- markers[markers$cluster == "Macrophage (il1b)",]
markers_cd59 <- markers[markers$cluster == "Macrophage (cd59)",]
markers_epdl <- markers[markers$cluster == "Macrophage (epdl)",]

markers_pct0.5_apoeb <- markers_pct0.5[markers_pct0.5$cluster == "Macrophage (apoeb)",]
markers_pct0.5_il1b <- markers_pct0.5[markers_pct0.5$cluster == "Macrophage (il1b)",]
markers_pct0.5_cd59 <- markers_pct0.5[markers_pct0.5$cluster == "Macrophage (cd59)",]
markers_pct0.5_epdl <- markers_pct0.5[markers_pct0.5$cluster == "Macrophage (epdl)",]

VlnPlot(macrophages, "apoeb", group.by = "plot.ident2")
pal <- readRDS("./output/processed_data/SCD-VI-i004/spatial_scatterpie_palette.rds")
DimPlot(scrnaseq,
                      cells.highlight = WhichCells(scrnaseq, idents = c("Macrophage (apoeb)"))) +
  scale_color_manual(labels = c("Unselected", "Macrophage (apoeb)"),
                     values = c("lightgrey", unname(pal[names(pal) == "Macrophage (apoeb)"]))) +
  theme(legend.position = "top", legend.text = element_text(size = 6))



fibroblasts <- subset(scrnaseq, idents = str_subset(unique(Idents(scrnaseq)), "^Fibroblas."))

markers_fibroblasts <- FindAllMarkers(fibroblasts)
markers_fib_cxcl12a <- markers_fibroblasts[markers_fibroblasts$cluster == "Fibroblast (cxcl12a)",]
markers_fib_col12a1a <- markers_fibroblasts[markers_fibroblasts$cluster == "Fibroblast (col12a1a)",]
markers_fib_col11a1a <- markers_fibroblasts[markers_fibroblasts$cluster == "Fibroblast (col11a1a)",]
