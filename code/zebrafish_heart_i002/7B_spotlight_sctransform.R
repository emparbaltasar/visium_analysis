library(SPOTlight)
library(Seurat)

# Load datasets
scrnaseq <- readRDS("./output/processed_data/SCD-VI-i001/zebrafish_heart_scrnaseq_sctransform_pca.rds")
visium_subset <- readRDS("./output/processed_data/SCD-VI-i002/SCD-VI-i002_sctransform_pca_clustered.rds")
mgs <- read.csv("./output/processed_data/SCD-VI-i001/scrnaseq_sctransform_markers_auc.csv", row.names = 1)

# Feature selection
hvg <- VariableFeatures(scrnaseq)

# keep only those genes that are relevant for each cell identity
mgs_df <- mgs[mgs$myAUC > 0.7, ]


## Cell downsampling
# split cell indices by identity
idx.seurat <- split(seq(ncol(scrnaseq)), scrnaseq$plot.ident2)
# downsample to at most 20 per identity & subset
n_cells <- 100
cs_keep <- lapply(idx.seurat, function(i) {
  n <- length(i)
  if (n < n_cells) {
    n_cells <- n
  }
  sample(i, n_cells)
})
scrnaseq.reduced <- scrnaseq[, unlist(cs_keep)]

## Deconvolution
res <- SPOTlight(
  x = GetAssayData(scrnaseq.reduced, slot = "data", assay = "SCT"),
  y = GetAssayData(visium_subset, slot = "data", assay = "SCT"),
  groups = scrnaseq.reduced$plot.ident2,
  mgs = mgs_df,
  hvg = hvg,
  weight_id = "myAUC",
  group_id = "cluster",
  gene_id = "gene",
  verbose = TRUE
)

saveRDS(res, file = "./output/processed_data/SCD-VI-i002/spotlight_sctransform-hvg_rawdata.rds")
