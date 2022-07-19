library(Seurat)

scrnaseq <- readRDS("./output/processed_data/SCD-VI-i004/zebrafish_heart_scrnaseq_lognorm_pca.rds")
visium_subset <- readRDS("./output/processed_data/SCD-VI-i003/SCD-VI-i003_subset_lognorm_pca_clustered.rds")

# Anchor-based integration
anchors <- FindTransferAnchors(
  reference = scrnaseq, query = visium_subset, k.filter = NA,
  normalization.method = "LogNormalize"
)

visium_subset <- TransferData(
  anchorset = anchors, refdata = scrnaseq$plot.ident2, query = visium_subset,
  prediction.assay = TRUE, weight.reduction = visium_subset[["pca"]], dims = 1:30
)

saveRDS(visium_subset, "./output/processed_data/SCD-VI-i003/SCD-VI-i003_lognorm_seurat_predictions.rds")
