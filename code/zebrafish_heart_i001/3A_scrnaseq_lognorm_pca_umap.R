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
scrnaseq <- subset(scrnaseq, time == "Ctrl" | time == "7dpi")


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

saveRDS(scrnaseq, "./output/processed_data/SCD-VI-i001/zebrafish_heart_scrnaseq_lognorm_pca.rds")

