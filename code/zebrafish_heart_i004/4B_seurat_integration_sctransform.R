library(Seurat)

scrnaseq <- readRDS("./output/processed_data/SCD-VI-i004/zebrafish_heart_scrnaseq_sctransform_pca.rds")

visium_subset <- readRDS("./output/processed_data/SCD-VI-i004/SCD-VI-i004_subset_sctransform_pca_clustered.rds")

# Anchor-based integration
anchors <- FindTransferAnchors(
  reference = scrnaseq, query = visium_subset,
  normalization.method = "SCT"
)

visium_subset <- TransferData(
  anchorset = anchors, refdata = scrnaseq$plot.ident2, query = visium_subset,
  prediction.assay = TRUE, weight.reduction = visium_subset[["pca"]], dims = 1:30
)

saveRDS(visium_subset, "./output/processed_data/SCD-VI-i004/SCD-VI-i004_sctransform_seurat_predictions.rds")
