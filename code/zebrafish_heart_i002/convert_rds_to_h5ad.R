library(Seurat)
library(SeuratDisk)

# Visium data
visium <- Load10X_Spatial("./data/SCD-VI-i001-004/spaceranger_output/SCD-VI-i002/outs",
                          filename = "filtered_feature_bc_matrix.h5",
                          slice = "slice1",
                          filter.matrix = TRUE)

SaveH5Seurat(visium, filename = "./output/processed_data/SCD-VI-i002/SCD-VI-i002.h5Seurat")

Convert("./output/processed_data/SCD-VI-i002/SCD-VI-i002.h5Seurat", dest = "h5ad")

