library(Seurat)

# Load and preprocess spatial dataset
visium <- Load10X_Spatial("./data/SCD-VI-i001-004/spaceranger_output/SCD-VI-i003/outs",
  filename = "filtered_feature_bc_matrix.h5",
  slice = "slice1",
  filter.matrix = TRUE
)

SpatialDimPlot(visium, cells.highlight = WhichCells(visium, expression = slice1_imagerow < 220 & slice1_imagecol > 240))
visium_subset <- subset(visium, slice1_imagerow < 220 & slice1_imagecol > 240, invert = TRUE)
SpatialDimPlot(visium_subset, label = FALSE)

dir.create("./output/processed_data/SCD-VI-i003")
saveRDS(visium_subset, "./output/processed_data/SCD-VI-i003/SCD-VI-i003_subset.rds")
