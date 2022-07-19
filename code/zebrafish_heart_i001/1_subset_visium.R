library(Seurat)

# Load and preprocess spatial dataset
visium <- Load10X_Spatial("./data/SCD-VI-i001-004/spaceranger_output/SCD-VI-i001/outs",
  filename = "filtered_feature_bc_matrix.h5",
  slice = "slice1",
  filter.matrix = TRUE
)

SpatialDimPlot(visium, cells.highlight = WhichCells(visium, expression = slice1_imagerow < 500 & slice1_imagecol > 260 & slice1_imagecol < 475))
visium_subset <- subset(visium, slice1_imagerow < 500 & slice1_imagecol > 260 & slice1_imagecol < 475, invert = FALSE)
SpatialDimPlot(visium_subset, label = FALSE)


saveRDS(visium_subset, "./output/processed_data/SCD-VI-i001/SCD-VI-i001_subset.rds")
