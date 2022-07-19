library(Seurat)
library(SeuratDisk)

scrnaseq <- readRDS("./data/mouse_brain/allen_cortex.rds")
scrnaseq <- UpdateSeuratObject(scrnaseq)
SaveH5Seurat(scrnaseq, filename = "./data/mouse_brain/allen_cortex_2.h5Seurat")

Convert("./data/mouse_brain/allen_cortex_2.h5Seurat", dest = "h5ad")
