library(Seurat)
library(SingleCellExperiment)
library(stringr)
library(dplyr)

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
scrnaseq <- NormalizeData(scrnaseq, normalization.method = "LogNormalize")
scrnaseq <- FindVariableFeatures(scrnaseq, selection.method = "vst", nfeatures = 3000)
scrnaseq <- ScaleData(scrnaseq)

scrnaseq <- RunPCA(scrnaseq, npcs = 30, verbose = TRUE)
scrnaseq <- RunUMAP(scrnaseq, reduction = "pca", dims = 1:30, verbose = TRUE)

saveRDS(scrnaseq, "./output/processed_data/SCD-VI-i004/zebrafish_heart_scrnaseq_lognorm_pca.rds")
scrnaseq <- readRDS("./output/processed_data/SCD-VI-i004/zebrafish_heart_scrnaseq_lognorm_pca.rds")

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
saveRDS(scrnaseq, "./output/processed_data/SCD-VI-i004/zebrafish_heart_scrnaseq_lognorm_pca.rds")
