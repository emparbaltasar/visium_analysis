library(Seurat)
library(ggplot2)
library(RColorBrewer)

visium_subset <- readRDS("./output/processed_data/mouse_brain/mouse_brain_spatial_cortex_lognorm_seurat_predictions.rds")

DefaultAssay(visium_subset) <- "predictions"

visium_subset <- FindSpatiallyVariableFeatures(visium_subset,
  assay = "predictions", selection.method = "markvariogram",
  features = rownames(visium_subset), r.metric = 5, slot = "data"
)



MyPalette <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(n = 100)

SpatialPlot(object = visium_subset, features = SpatiallyVariableFeatures(visium_subset), ncol = 5, crop = FALSE) +
  ggplot2::scale_fill_gradientn(colours = MyPalette, limits = c(0.0, 1.0), breaks = c(0.0, 0.5, 1.0))
