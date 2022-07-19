library(Seurat)
library(SingleCellExperiment)
library(tidyverse)
library(patchwork)

# Load and preprocess scrnaseq dataset
load("./data/SCD-VI-i001-004/scRNAseq/SCE_fullobject.Robj")
logcounts(final.all.hearts.sce) <- assay(final.all.hearts.sce, "counts")
scrnaseq <- as.Seurat(final.all.hearts.sce)
scrnaseq <- RenameAssays(scrnaseq, originalexp = "RNA")
remove(final.all.hearts.sce)

scrnaseq <- subset(scrnaseq, morphine == "Ctrl")
scrnaseq <- subset(scrnaseq, time == "7dpi")


# QC
scrnaseq[["percent.mt"]] <- PercentageFeatureSet(scrnaseq, pattern = "^mt-")

VlnPlot(scrnaseq, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # It seems like it has been already filtered

FeatureScatter(scrnaseq, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(scrnaseq, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


# Normalization
scrnaseq <- SCTransform(scrnaseq, ncells = 3000, verbose = FALSE)

scrnaseq <- RunPCA(scrnaseq, npcs = 30, verbose = FALSE)
scrnaseq <- RunUMAP(scrnaseq, reduction = "pca", dims = 1:30, verbose = FALSE)

saveRDS(scrnaseq, "./output/processed_data/SCD-VI-i004/zebrafish_heart_scrnaseq_sctransform_pca.rds")

scrnaseq <- readRDS("./output/processed_data/SCD-VI-i004/zebrafish_heart_scrnaseq_sctransform_pca.rds")

# Plot UMAP
Idents(scrnaseq) <- "plot.ident2"
pal <- readRDS("./output/processed_data/SCD-VI-i004/spatial_scatterpie_palette.rds")

DimPlot(scrnaseq, cols = pal) +
  theme(legend.position = "none")
ggsave("./output/images/SCD-VI-i004_analysis/scrnaseq_umap_celltypes.png")

fibroblast <- DimPlot(scrnaseq, 
                      cells.highlight = WhichCells(scrnaseq, idents = c("Fibroblast (spock3)"))) +
  scale_color_manual(labels = c("Unselected", "Fibroblast (spock3)"), 
                     values = c("lightgrey", unname(pal[names(pal) == "Fibroblast (spock3)"]))) +
  theme(legend.position = "top", legend.text = element_text(size = 6))

macrophage_end_dup <- DimPlot(scrnaseq, 
                              cells.highlight = WhichCells(scrnaseq, idents = c("Macrophage (Endothelia duplex)"))) +
  scale_color_manual(labels = c("Unselected", "Macrophage (Endothelia duplex)"), 
                     values = c("lightgrey", unname(pal[names(pal) == "Macrophage (Endothelia duplex)"]))) +
  theme(legend.position = "top", legend.text = element_text(size = 6))

macrophage_cm_dup <- DimPlot(scrnaseq, 
                              cells.highlight = WhichCells(scrnaseq, idents = c("Macrophage (CM duplex)"))) +
  scale_color_manual(labels = c("Unselected", "Macrophage (CM duplex)"), 
                     values = c("lightgrey", unname(pal[names(pal) == "Macrophage (CM duplex)"]))) +
  theme(legend.position = "top", legend.text = element_text(size = 6))

macrophage_fib_dup <- DimPlot(scrnaseq, 
                             cells.highlight = WhichCells(scrnaseq, idents = c("Macrophage (Fibroblast duplex)"))) +
  scale_color_manual(labels = c("Unselected", "Macrophage (Fibroblast duplex)"), 
                     values = c("lightgrey", unname(pal[names(pal) == "Macrophage (Fibroblast duplex)"]))) +
  theme(legend.position = "top", legend.text = element_text(size = 6))

cardiomyocytes_a <- DimPlot(scrnaseq, 
                              cells.highlight = WhichCells(scrnaseq, idents = c("Cardiomyocytes A"))) +
  scale_color_manual(labels = c("Unselected", "Cardiomyocytes A"), 
                     values = c("lightgrey", unname(pal[names(pal) == "Cardiomyocytes A"]))) +
  theme(legend.position = "top", legend.text = element_text(size = 6))

macrophages <- DimPlot(scrnaseq, 
                            cells.highlight = WhichCells(scrnaseq, idents = c("Macrophages"))) +
  scale_color_manual(labels = c("Unselected", "Macrophages"), 
                     values = c("lightgrey", unname(pal[names(pal) == "Macrophages"]))) +
  theme(legend.position = "top", legend.text = element_text(size = 6))

wrap_plots(macrophage_fib_dup, macrophage_end_dup, fibroblast, macrophage_cm_dup, cardiomyocytes_a, macrophages)
ggsave("./output/images/SCD-VI-i004_analysis/scrnaseq_umap_spotlight_problematic_celltypes.png")


# Simplify cell types
plot.ident2.simpl <- scrnaseq$plot.ident2
levels(plot.ident2.simpl)
plot.ident2.simpl <- recode_factor(plot.ident2.simpl,
  "Epicardium (Atrium)" = "Epicardium",
  "Epicardium (Ventricle)" = "Epicardium",
  "Macrophage (il1b)" = "Macrophages",
  "Macrophage (apoeb)" = "Macrophages",
  "Macrophage (proliferating)" = "Macrophages",
  "Macrophage (CM duplex)" = "Macrophages",
  "Macrophage (Endothelia duplex)" = "Macrophages",
  "Macrophage (epdl)" = "Macrophages",
  "Macrophage (cd59)" = "Macrophages",
  "Macrophage (Ery duplex)" = "Macrophages",
  "Macrophage (Fibroblast duplex)" = "Macrophages",
  "T-cells (il4/13)" = "T-cells",
  "T-cells (proliferating)" = "T-cells",
  "Fibroblast" = "Fibroblasts",
  "Fibroblast (cxcl12a)" = "Fibroblasts",
  "Fibroblast (col11a1a)" = "Fibroblasts",
  "Fibroblast (cfd)" = "Fibroblasts",
  "Fibroblast-like cells" = "Fibroblasts",
  "Fibroblast (col12a1a)" = "Fibroblasts",
  "Fibroblast (nppc)" = "Fibroblasts",
  "Fibroblast (spock3)" = "Fibroblasts",
  "Fibroblast (mpeg1.1)" = "Fibroblasts",
  "Fibroblast (proliferating)" = "Fibroblasts",
  "Endocardium (A)" = "Endocardium",
  "Endocardium (V)" = "Endocardium",
  "Endocardium (frzb)" = "Endocardium",
  "Cardiomyocytes (ttn.2) A" = "Cardiomyocytes",
  "Cardiomyocytes V" = "Cardiomyocytes",
  "Cardiomyocytes A" = "Cardiomyocytes",
  "Cardiomyocytes (ttn.2) V" = "Cardiomyocytes",
  "Cardiomyocytes (proliferating)" = "Cardiomyocytes",
  "Bl.ves.EC (apnln)" = "Bl.ves.EC",
  "Bl.ves.EC (plvapb)" = "Bl.ves.EC",
  "Bl.ves.EC (lyve1)" = "Bl.ves.EC",
)
scrnaseq$plot.ident2.simpl <- plot.ident2.simpl
saveRDS(scrnaseq, "./output/processed_data/SCD-VI-i004/zebrafish_heart_scrnaseq_sctransform_pca.rds")
