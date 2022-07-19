library(Seurat)

scrnaseq <- readRDS("./output/processed_data/SCD-VI-i004/zebrafish_heart_scrnaseq_sctransform_pca.rds")

visium <- readRDS("./output/processed_data/SCD-VI-i003_full_area/SCD-VI-i003_full_area_sctransform_pca_clustered.rds")

# Anchor-based integration
anchors <- FindTransferAnchors(
  reference = scrnaseq, query = visium,
  normalization.method = "SCT"
)

visium <- TransferData(
  anchorset = anchors, refdata = scrnaseq$plot.ident2, query = visium,
  prediction.assay = TRUE, weight.reduction = visium[["pca"]], dims = 1:30
)

saveRDS(visium, "./output/processed_data/SCD-VI-i003_full_area/SCD-VI-i003_full_area_sctransform_seurat_predictions.rds")
