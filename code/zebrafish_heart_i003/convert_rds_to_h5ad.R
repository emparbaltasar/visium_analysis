library(Seurat)
library(SeuratDisk)

# Visium data
rds <- readRDS("./output/processed_data/SCD-VI-i003/SCD-VI-i003_subset.rds")
rds <- UpdateSeuratObject(rds)
SaveH5Seurat(rds, filename = "./output/processed_data/SCD-VI-i003/SCD-VI-i003_subset.h5Seurat")

Convert("./output/processed_data/SCD-VI-i003/SCD-VI-i003_subset.h5Seurat", dest = "h5ad")
