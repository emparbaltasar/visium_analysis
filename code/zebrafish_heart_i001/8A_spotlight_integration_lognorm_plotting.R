library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(scater)
library(SPOTlight)
library(tidyverse)
library(ggpubr)
library(NMF)
library(ggcorrplot)

MyPalette <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(n = 100)

# Load datasets
res <- readRDS("./output/processed_data/SCD-VI-i001/spotlight_lognorm-hvg_rawdata.rds")
visium_subset <- readRDS("./output/processed_data/SCD-VI-i001/SCD-VI-i001_subset_lognorm_pca_clustered.rds")
scrnaseq <- readRDS("./output/processed_data/SCD-VI-i001/zebrafish_heart_scrnaseq_lognorm_pca.rds")

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



DefaultAssay(visium_subset)

# Extract deconvolution matrix
head(mat <- res$mat)[, seq_len(3)]
mat_raw <- mat

# Extract NMF model fit
mod <- res$NMF


# Plotting
plotTopicProfiles_test <- function(
    x,
    y,
    facet = FALSE,
    min_prop = 0.01,
    ncol = NULL) {
  
  # Convert y to character
  y <- as.character(y)
  
  # check validity of input arguments
  stopifnot(
    is(x, "NMFfit"),
    is.character(y),
    length(y) == ncol(coef(x)),
    setequal(
      colnames(basis(x)), unique(y) ############## modified from original code: 
      # colnames(basis(x)), paste0("topic_", seq_len(length(unique(y)))
      # to make it work
    ),
    is.logical(facet), length(facet) == 1,
    is.numeric(min_prop), length(min_prop) == 1)
  
  # get group proportions
  mat <- prop.table(t(coef(x)), 1)
  
  if (facet) {
    # stretch for plotting
    df <- data.frame(
      id = seq_len(nrow(mat)),
      weight = c(mat),
      group = rep(y, ncol(mat)),
      topic = rep(seq_len(ncol(mat)), each = nrow(mat)))
    
    # drop cells with 'weight < min_prop'
    df <- df[df$weight >= min_prop, ]
    
    # set aesthetics
    x <- "id"
    f <- facet_wrap(~group, ncol = ncol, scales = "free_x")
  } else {
    # get topic medians
    df <- aggregate(mat, list(y), median)[, -1]
    rownames(df) <- unique(y)
    
    # stretch for plotting
    df <- data.frame(
      weight = unlist(df),
      group = rep(rownames(df), each = nrow(df)),
      topic = rep(seq_len(nrow(df)), ncol(df)))
    
    # set aesthetics
    x <- "group"
    f <- NULL
  }
  # fix topic order
  df$topic <- factor(df$topic, seq_along(unique(y)))
  
  # render plot
  ggplot(df, aes_string(x, "topic",
                        col = "weight", size = "weight")) +
    f + geom_point() +
    guides(col = guide_legend(override.aes = list(size = 2))) +
    scale_size_continuous(range = c(0, 3)) +
    scale_color_continuous(low = "lightgrey", high = "#3d2bff") +
    xlab(if (facet) x) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.key.size = unit(0.5, "lines"),
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1))
}

plotTopicProfiles_test(
  x = mod,
  y = scrnaseq.reduced$plot.ident2,
  facet = FALSE,
  min_prop = 0.01,
  ncol = 1
) +
  theme(aspect.ratio = 1)
ggsave("./output/images/SCD-VI-i001_integration/spotlight_lognorm-hvg_rawdata_topic_profiles.png")


plotTopicProfiles_test(
  x = mod,
  y = scrnaseq.reduced$plot.ident2,
  facet = TRUE,
  min_prop = 0.01,
  ncol = 12
)
ggsave("./output/images/SCD-VI-i001_integration/spotlight_lognorm-hvg_rawdata_topic_profiles_facet.png", height = 20, width = 25)



sign <- basis(mod)
colnames(sign) <- paste0("Topic", seq_len(ncol(sign)))
head(sign)


plotCorrelationMatrix(mat)
ggsave("./output/images/SCD-VI-i001_integration/spotlight_lognorm-hvg_rawdata_correlation_matrix.png")

plotInteractions(mat, "heatmap")
ggsave("./output/images/SCD-VI-i001_integration/spotlight_lognorm-hvg_rawdata_interactions_heatmap.png")



# Plot spatial scatterpie
pal <- readRDS("./output/processed_data/SCD-VI-i004/spatial_scatterpie_palette.rds")

plotSpatialScatterpie(
  x = visium_subset,
  y = mat,
  pie_scale = 0.8
) +
  scale_fill_manual(
    values = pal[names(pal) %in% names(colSums(mat) != 0)[colSums(mat) != 0]],
    breaks = names(pal)[names(pal) %in% names(colSums(mat) != 0)[colSums(mat) != 0]]
  ) +
  scale_y_reverse() +
  theme(legend.title = element_blank(), legend.text = element_text(size = 6)) +
  guides(fill = guide_legend(ncol = 1))

ggsave("./output/images/SCD-VI-i001_integration/spotlight_lognorm-hvg_rawdata_prediction_spatial_scatterpie.png")
ggsave("./output/images/SCD-VI-i001_integration/spotlight_lognorm-hvg_rawdata_prediction_spatial_scatterpie.pdf")

# Plot spatial scatterpie with 0.04 threshold (1/25)
mat_filtered <- mat
mat_filtered[mat_filtered < 0.04] <- 0

plotSpatialScatterpie(
  x = visium_subset,
  y = mat_filtered,
  pie_scale = 0.8
) +
  scale_fill_manual(
    values = pal[names(pal) %in% names(colSums(mat) != 0)[colSums(mat) != 0]],
    breaks = names(pal)[names(pal) %in% names(colSums(mat) != 0)[colSums(mat) != 0]]
  ) +
  scale_y_reverse() +
  theme(legend.title = element_blank(), legend.text = element_text(size = 6)) +
  guides(fill = guide_legend(ncol = 1))

ggsave("./output/images/SCD-VI-i001_integration/spotlight_lognorm-hvg_rawdata_prediction_filtered_spatial_scatterpie.png")
ggsave("./output/images/SCD-VI-i001_integration/spotlight_lognorm-hvg_rawdata_prediction_filtered_spatial_scatterpie.pdf")


# Plot per cell type
visium_subset[["predictionsspotlight"]] <- CreateAssayObject(t(mat))

DefaultAssay(visium_subset) <- "predictionsspotlight"

for (i in rownames(visium_subset@assays$predictions@data)) {
  p <- SpatialFeaturePlot(visium_subset, features = i, pt.size.factor = 2, crop = TRUE) +
    ggplot2::scale_fill_gradientn(colours = MyPalette, limits = c(0.0, 1.001), breaks = c(0.0, 0.5, 1.0))

  if (grepl("/", i)) {
    i <- sub("/", "-", i)
  }
  if (grepl(" ", i)) {
    i <- sub(" ", "_", i)
  }

  ggsave(paste("./output/images/SCD-VI-i001_integration/spotlight_lognorm-hvg_rawdata_prediction_", i, ".png", sep = ""))
}

