library(Seurat)
library(SeuratDisk)

# Visium data
rds <- readRDS("./output/processed_data/SCD-VI-i001/SCD-VI-i001_subset.rds")
rds <- UpdateSeuratObject(rds)
SaveH5Seurat(rds, filename = "./output/processed_data/SCD-VI-i001/SCD-VI-i001_subset.h5Seurat")

Convert("./output/processed_data/SCD-VI-i001/SCD-VI-i001_subset.h5Seurat", dest = "h5ad")


# scrnaseq data
load("./data/SCD-VI-i001-004/scRNAseq/SCE_fullobject.Robj")
logcounts(final.all.hearts.sce) <- assay(final.all.hearts.sce, "counts")
scrnaseq <- as.Seurat(final.all.hearts.sce)
scrnaseq <- RenameAssays(scrnaseq, originalexp = "RNA")
remove(final.all.hearts.sce)

scrnaseq <- subset(scrnaseq, morphine == "Ctrl")
scrnaseq <- subset(scrnaseq, time == "Ctrl" | time == "7dpi")

scrnaseq@meta.data[["plot.ident2"]] <- as.character(scrnaseq@meta.data[["plot.ident2"]])
scrnaseq@meta.data[["time"]] <- as.character(scrnaseq@meta.data[["time"]])

rds <- UpdateSeuratObject(scrnaseq)
SaveH5Seurat(rds, filename = "./output/processed_data/SCD-VI-i001/zebrafish_heart_scrnaseq_subset.h5Seurat")
Convert("./output/processed_data/SCD-VI-i001/zebrafish_heart_scrnaseq_subset.h5Seurat", dest = "h5ad")
