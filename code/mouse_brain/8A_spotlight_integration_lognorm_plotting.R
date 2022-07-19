library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(scater)
library(SPOTlight)
library(tidyverse)
library(ggpubr)

MyPalette <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(n = 100)

# Load datasets
res <- readRDS("./output/processed_data/mouse_brain/spotlight_lognorm.rds")
visium_subset <- readRDS("./output/processed_data/mouse_brain/mouse_brain_spatial_cortex_lognorm_pca_clustered.rds")
scrnaseq <- readRDS("./output/processed_data/mouse_brain/mouse_brain_scrnaseq_lognorm_pca.rds")

## Cell downsampling
# split cell indices by identity
idx.seurat <- split(seq(ncol(scrnaseq)), scrnaseq$subclass)
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



DefaultAssay(visium_subset)

# Extract deconvolution matrix
head(mat <- res$mat)[, seq_len(3)]
mat_raw <- mat

# Extract NMF model fit
mod <- res$NMF


# Plotting
plotTopicProfiles(
  x = mod,
  y = scrnaseq.reduced$subclass,
  facet = FALSE,
  min_prop = 0.01,
  ncol = 1
) +
  theme(aspect.ratio = 1)
ggsave("./output/images/mouse_brain_integration/spotlight_lognorm_topic_profiles.png")


plotTopicProfiles(
  x = mod,
  y = scrnaseq.reduced$subclass,
  facet = TRUE,
  min_prop = 0.01,
  ncol = 8
)
ggsave("./output/images/mouse_brain_integration/spotlight_lognorm_topic_profiles_facet.png", height = 9.5, width = 11.5)


library(NMF)
sign <- basis(mod)
colnames(sign) <- paste0("Topic", seq_len(ncol(sign)))
head(sign)

library(ggcorrplot)
plotCorrelationMatrix(mat)
ggsave("./output/images/mouse_brain_integration/spotlight_lognorm_correlation_matrix.png")


plotInteractions(mat, "heatmap")
ggsave("./output/images/mouse_brain_integration/spotlight_lognorm_interactions_heatmap.png")


# Plot spatial scatterpie
pal <- readRDS("./output/processed_data/mouse_brain/spatial_scatterpie_palette.rds")
pal <- pal[names(pal) != "max"]

plotSpatialScatterpie(
  x = visium_subset,
  y = mat,
  cell_types = colnames(y),
  img = FALSE,
  scatterpie_alpha = 1,
  pie_scale = 0.4
) +
  scale_fill_manual(
    values = pal,
    breaks = names(pal)
  ) +
  scale_y_reverse()

ggsave("./output/images/mouse_brain_integration/spotlight_lognorm_prediction_spatial_scatterpie.png")
ggsave("./output/images/mouse_brain_integration/spotlight_lognorm_prediction_spatial_scatterpie.pdf")


# Plot per cell type
visium_subset[["predictionsspotlight"]] <- CreateAssayObject(t(mat))

DefaultAssay(visium_subset) <- "predictionsspotlight"

for (i in rownames(visium_subset@assays$predictions@data)) {
  p <- SpatialFeaturePlot(visium_subset, features = i, pt.size.factor = 1, crop = FALSE) +
    ggplot2::scale_fill_gradientn(colours = MyPalette, limits = c(0.0, 1.001), breaks = c(0.0, 0.5, 1.0))

  if (grepl("/", i)) {
    i <- sub("/", "-", i)
  }
  if (grepl(" ", i)) {
    i <- sub(" ", "_", i)
  }

  ggsave(paste("./output/images/mouse_brain_integration/spotlight_lognorm_prediction_", i, ".png", sep = ""))
}


# Compare proportions of prediction with proportions in scrnaseq dataset
spotlight_prediction <- colMeans(mat_raw)
scRNAseq_proportions <- table(scrnaseq@meta.data[["subclass"]]) / length(scrnaseq@meta.data[["subclass"]])

proportions_comparison <- rbind(spotlight_prediction, scRNAseq_proportions)
proportions_comparison <- as.data.frame(proportions_comparison) %>%
  t() %>%
  as.data.frame()

proportions_comparison %>%
  rownames_to_column("cell_type") %>%
  ggplot(aes(x = scRNAseq_proportions, y = spotlight_prediction)) +
  scale_color_manual(values = pal) +
  theme_bw() +
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank(),
    aspect.ratio = 1 / 1, legend.text = element_text(size = 6)
  ) +
  geom_point(aes(color = cell_type)) +
  stat_cor(method = "pearson", label.x = 0.07, label.y = 0.45, size = 3) +
  stat_cor(method = "spearman", label.x = 0.07, label.y = 0.42, cor.coef.name = "rho", size = 3)

ggsave("./output/images/mouse_brain_integration/spotlight_lognorm_comparison_celltype_proportions_scatter.png", width = 5, height = 4)
