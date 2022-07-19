library(Seurat)

scrnaseq <- readRDS("./output/processed_data/mouse_brain/mouse_brain_scrnaseq_lognorm_pca.rds")
visium_subset <- readRDS("./output/processed_data/mouse_brain/mouse_brain_spatial_cortex_lognorm_pca_clustered.rds")

# Anchor-based integration
anchors <- FindTransferAnchors(
  reference = scrnaseq, query = visium_subset,
  normalization.method = "LogNormalize"
)

visium_subset <- TransferData(
  anchorset = anchors, refdata = scrnaseq$subclass, query = visium_subset,
  prediction.assay = TRUE, weight.reduction = visium_subset[["pca"]], dims = 1:30
)

saveRDS(visium_subset, "./output/processed_data/mouse_brain/mouse_brain_spatial_cortex_lognorm_seurat_predictions.rds")
