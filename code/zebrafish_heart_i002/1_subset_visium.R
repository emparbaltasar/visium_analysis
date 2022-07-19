library(Seurat)

# Load and preprocess spatial dataset
visium <- Load10X_Spatial("./data/SCD-VI-i001-004/spaceranger_output/SCD-VI-i002/outs",
  filename = "filtered_feature_bc_matrix.h5",
  slice = "slice1",
  filter.matrix = TRUE
)

SpatialDimPlot(visium)

# No need to subset