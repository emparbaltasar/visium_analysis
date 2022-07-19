library(Seurat)

# Load and preprocess scrnaseq dataset
scrnaseq <- readRDS("./data/mouse_brain/allen_cortex.rds")

scrnaseq <- NormalizeData(scrnaseq, normalization.method = "LogNormalize")
scrnaseq <- FindVariableFeatures(scrnaseq, selection.method = "vst", nfeatures = 3000)
scrnaseq <- ScaleData(scrnaseq)

scrnaseq <- RunPCA(scrnaseq, npcs = 30, verbose = FALSE)
scrnaseq <- RunUMAP(scrnaseq, reduction = "pca", dims = 1:30, verbose = FALSE)

saveRDS(scrnaseq, "./output/processed_data/mouse_brain/mouse_brain_scrnaseq_lognorm_pca.rds")
