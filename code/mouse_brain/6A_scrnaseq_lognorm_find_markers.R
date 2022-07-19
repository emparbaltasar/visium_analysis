library(Seurat)

scrnaseq <- readRDS("./output/processed_data/mouse_brain/mouse_brain_scrnaseq_lognorm_pca.rds")

DefaultAssay(scrnaseq) <- "RNA"

# Feature selection
hvg <- VariableFeatures(scrnaseq)

# Get vector indicating which genes are neither ribosomal or mitochondrial
genes <- !grepl(pattern = "^Rp[l|s]|Mt", x = rownames(scrnaseq))

# Compute marker genes
Idents(scrnaseq) <- "subclass"
mgs <- FindAllMarkers(scrnaseq,
  assay = "RNA", slot = "data", test.use = "roc", features = genes,
  logfc.threshold = 0.4, only.pos = TRUE
)
write.csv(mgs, "./output/processed_data/mouse_brain/scrnaseq_lognorm_markers_auc.csv")
