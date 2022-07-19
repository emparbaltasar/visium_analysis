library(Seurat)

scrnaseq <- readRDS("./output/processed_data/SCD-VI-i001/zebrafish_heart_scrnaseq_lognorm_pca.rds")

DefaultAssay(scrnaseq) <- "RNA"

# Feature selection
hvg <- VariableFeatures(scrnaseq)

# Get vector indicating which genes are neither ribosomal or mitochondrial
genes <- !grepl(pattern = "^rp[l|s]|mt", x = rownames(scrnaseq))

# Compute marker genes
Idents(scrnaseq) <- "plot.ident2"
mgs <- FindAllMarkers(scrnaseq,
  assay = "RNA", slot = "data", test.use = "roc", features = genes,
  logfc.threshold = 0.4, only.pos = TRUE
)
write.csv(mgs, "./output/processed_data/SCD-VI-i001/scrnaseq_lognorm_markers_auc.csv")
