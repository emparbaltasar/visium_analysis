library(dplyr)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)

MyPalette <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(n = 100)

# Violin plots for location of celltypes ----
# cell2location ----
visium_subset_cell2loc <- readRDS("./output/processed_data/mouse_brain/mouse_brain_spatial_cortex.rds")
cell2loc_pred <- read.csv("./output/cell2location/mouse_brain/cell2location_map/cell2location_prediction.csv",
                          row.names=1)

cell2loc_pred_subset <- cell2loc_pred[rownames(cell2loc_pred) %in% colnames(visium_subset_cell2loc),]
cell2loc_pred_subset_prop <- as_tibble(prop.table(as.matrix(cell2loc_pred_subset), 1), rownames = NA)
colnames(cell2loc_pred_subset_prop)[4] <- "L2/3 IT"
colnames(cell2loc_pred_subset_prop)[6] <- "L5 IT"
colnames(cell2loc_pred_subset_prop)[7] <- "L5 PT"
colnames(cell2loc_pred_subset_prop)[8] <- "L6 CT"
colnames(cell2loc_pred_subset_prop)[9] <- "L6 IT"

visium_subset_cell2loc[["prediction.score.id"]] <- CreateAssayObject(t(cell2loc_pred_subset_prop))

DefaultAssay(visium_subset_cell2loc) <- "prediction.score.id"

prediction_cell2loc <- visium_subset_cell2loc@assays[["prediction.score.id"]]@data
SpatialFeaturePlot(visium_subset_cell2loc, features = "Astro", pt.size.factor = 1, crop = FALSE) +
  ggplot2::scale_fill_gradientn(colours = MyPalette, limits = c(0.0, 1.001), breaks = c(0.0, 0.5, 1.0))

SpatialFeaturePlot(visium_subset_cell2loc, features = "L6b", pt.size.factor = 1, crop = FALSE) +
  ggplot2::scale_fill_gradientn(colours = MyPalette, limits = c(0.0, 1.001), breaks = c(0.0, 0.5, 1.0))


# Divide barcodes into squares
# SpatialDimPlot(subset(visium_subset, anterior1_imagerow > 330 | anterior1_imagecol < 200, invert = TRUE))

sq_1_cell2loc <- subset(visium_subset_cell2loc, anterior1_imagerow > 330 | anterior1_imagecol < 200, invert = TRUE)

# SpatialDimPlot(sq_1)
SpatialFeaturePlot(sq_1_cell2loc, features = "Astro", pt.size.factor = 1, crop = FALSE) +
  ggplot2::scale_fill_gradientn(colours = MyPalette, limits = c(0.0, 1.001), breaks = c(0.0, 0.5, 1.0))


# Transform barcodes into xy
coordinates_cell2loc <- sq_1_cell2loc@images[["anterior1"]]@coordinates[,c("col","row")] %>% 
  rownames_to_column(var = "barcodes")
prediction_cell2loc <- sq_1_cell2loc@assays[["prediction.score.id"]]@data %>% as_tibble(rownames = NA) %>% t() %>% as_tibble(rownames = NA) %>% 
  rownames_to_column(var = "barcodes")
prediction_coord_cell2loc <- right_join(prediction_cell2loc, coordinates_cell2loc, by = "barcodes")


cortex_surface_seurat_cell2loc <- prediction_coord_cell2loc[prediction_coord_cell2loc$Astro > 0.3,c("col","row")]
ggplot(cortex_surface_seurat_cell2loc, mapping = aes(x = col, y = row)) +
  geom_point()

cortex_surface_seurat_cell2loc <- cortex_surface_seurat_cell2loc[cortex_surface_seurat_cell2loc$col < 100,c("col","row")]
ggplot(cortex_surface_seurat_cell2loc, mapping = aes(x = col, y = row)) +
  geom_point()

line <- lm(row ~ col, data = cortex_surface_seurat_cell2loc)

ggplot(cortex_surface_seurat_cell2loc, mapping = aes(x = col, y = row)) +
  geom_point() 

#Perpendicular distance from point 'a' to a line with 'slope' and 'intercept'
dist_point_line <- function(a, slope, intercept) {
  b = c(1, intercept+slope)
  c = c(-intercept/slope,0)       
  v1 <- b - c
  v2 <- a - b
  m <- cbind(v1,v2)
  return(abs(det(m))/sqrt(sum(v1*v1)))
}
prediction_coord_cell2loc$dist_to_cortex <- apply(prediction_coord_cell2loc, 1, function(x) dist_point_line(as.numeric(x[25:26]), slope = line[["coefficients"]][["col"]], intercept = line[["coefficients"]][["(Intercept)"]]) )

ggplot(prediction_coord_cell2loc, mapping = aes(x = col, y = row, color=dist_to_cortex)) +
  geom_point() 

# Normalize, 0 being surface and 1 being the median distance of L6b to cortex surface
SpatialFeaturePlot(sq_1_cell2loc, features = "L6b", pt.size.factor = 1, crop = FALSE) +
  ggplot2::scale_fill_gradientn(colours = MyPalette, limits = c(0.0, 1.001), breaks = c(0.0, 0.5, 1.0))

cortex_l6b_cell2loc <- prediction_coord_cell2loc[prediction_coord_cell2loc$L6b > 0.2,c("col","row")]
ggplot(cortex_l6b_cell2loc, mapping = aes(x = col, y = row)) +
  geom_point()

cortex_l6b_cell2loc <- cortex_l6b_cell2loc[cortex_l6b_cell2loc$col < 100,c("col","row")]
ggplot(cortex_l6b_cell2loc, mapping = aes(x = col, y = row)) +
  geom_point()

cortex_l6b_cell2loc <- prediction_coord_cell2loc[prediction_coord_cell2loc$L6b > 0.2,c("col","row","dist_to_cortex")]
cortex_l6b_cell2loc <- cortex_l6b_cell2loc[cortex_l6b_cell2loc$col < 100,c("dist_to_cortex")]
median(cortex_l6b_cell2loc$dist_to_cortex)

prediction_coord_cell2loc$dist_to_cortex = prediction_coord_cell2loc$dist_to_cortex / median(cortex_l6b_cell2loc$dist_to_cortex)

ggplot(prediction_coord_cell2loc, mapping = aes(x = col, y = row, color=dist_to_cortex)) +
  geom_point() 


# Transform matrix into 0 and 1
long_format_cell2loc <- pivot_longer(prediction_coord_cell2loc, cols = 2:24, names_to = "celltype", values_to = "prob")
long_format_cell2loc <- long_format_cell2loc[long_format_cell2loc$prob > 0.05,]

## Multiplicate rows with higher probabilities:
# How many replicates you want of each row
duptimes <- c(long_format_cell2loc$prob %/% 0.01)

# Create an index of the rows you want with duplications
idx <- rep(1:nrow(long_format_cell2loc), duptimes)

# Use that index to generate your new data frame
long_format_cell2loc <- long_format_cell2loc[idx,]

long_format_cell2loc$celltype <- as.factor(long_format_cell2loc$celltype)
long_format_cell2loc$method <- "Cell2location"
ggplot(long_format_cell2loc, aes(x=celltype, y=dist_to_cortex)) +
  geom_violin() +
  scale_y_reverse()


# Seurat lognorm ----
visium_subset_seurat_lognorm <- readRDS("./output/processed_data/mouse_brain/mouse_brain_spatial_cortex_lognorm_seurat_predictions.rds")
DefaultAssay(visium_subset_seurat_lognorm) <- "prediction.score.id"

# SpatialFeaturePlot(visium_subset_seurat_lognorm, features = "Astro", pt.size.factor = 1, crop = FALSE) +
#   ggplot2::scale_fill_gradientn(colours = MyPalette, limits = c(0.0, 1.001), breaks = c(0.0, 0.5, 1.0))



# Divide barcodes into squares
# SpatialDimPlot(subset(visium_subset_seurat_lognorm, anterior1_imagerow > 330 | anterior1_imagecol < 200, invert = TRUE))

sq_1_seurat_lognorm <- subset(visium_subset_seurat_lognorm, anterior1_imagerow > 330 | anterior1_imagecol < 200, invert = TRUE)

# SpatialDimPlot(sq_1_seurat_lognorm)
# SpatialFeaturePlot(sq_1_seurat_lognorm, features = "Astro", pt.size.factor = 1, crop = FALSE) +
  # ggplot2::scale_fill_gradientn(colours = MyPalette, limits = c(0.0, 1.001), breaks = c(0.0, 0.5, 1.0))



# Transform barcodes into xy
coordinates_seurat_lognorm <- sq_1_seurat_lognorm@images[["anterior1"]]@coordinates[,c("col","row")] %>% 
  rownames_to_column(var = "barcodes")
prediction_seurat_lognorm <- sq_1_seurat_lognorm@assays[["prediction.score.id"]]@data %>% t() %>% as_tibble(rownames = NA) %>% 
  rownames_to_column(var = "barcodes")
prediction_coord_seurat_lognorm <- right_join(prediction_seurat_lognorm, 
                               coordinates_seurat_lognorm, 
                               by = "barcodes")


prediction_coord_seurat_lognorm$dist_to_cortex <- apply(prediction_coord_seurat_lognorm, 1, function(x) dist_point_line(as.numeric(x[25:26]), slope = line[["coefficients"]][["col"]], intercept = line[["coefficients"]][["(Intercept)"]]) )

ggplot(prediction_coord_seurat_lognorm, mapping = aes(x = col, y = row, color=dist_to_cortex)) +
  geom_point() 

prediction_coord_seurat_lognorm$dist_to_cortex = prediction_coord_seurat_lognorm$dist_to_cortex / median(cortex_l6b_cell2loc$dist_to_cortex)

ggplot(prediction_coord_seurat_lognorm, mapping = aes(x = col, y = row, color=dist_to_cortex)) +
  geom_point() 


# Transform matrix into 0 and 1
long_format_lognorm <- pivot_longer(prediction_coord_seurat_lognorm, cols = 2:24, names_to = "celltype", values_to = "prob")
long_format_lognorm <- long_format_lognorm[long_format_lognorm$prob > 0.1,]

## Multiplicate rows with higher probabilities:
# How many replicates you want of each row
duptimes <- c(long_format_lognorm$prob %/% 0.1)
 
# Create an index of the rows you want with duplications
idx <- rep(1:nrow(long_format_lognorm), duptimes)
 
# Use that index to generate your new data frame
long_format_lognorm <- long_format_lognorm[idx,]
 


long_format_lognorm$celltype <- as.factor(long_format_lognorm$celltype)
long_format_lognorm$method <- "Seurat (LogNorm)"
ggplot(long_format_lognorm, aes(x=celltype, y=dist_to_cortex)) +
  geom_violin() +
  scale_y_reverse()


# seurat sctransform ----
visium_subset_seurat_sctransform <- readRDS("./output/processed_data/mouse_brain/mouse_brain_spatial_cortex_sctransform_seurat_predictions.rds")
DefaultAssay(visium_subset_seurat_sctransform) <- "prediction.score.id"

# SpatialFeaturePlot(visium_subset_seurat_sctransform, features = "Astro", pt.size.factor = 1, crop = FALSE) +
#   ggplot2::scale_fill_gradientn(colours = MyPalette, limits = c(0.0, 1.001), breaks = c(0.0, 0.5, 1.0))


# Divide barcodes into squares
# SpatialDimPlot(subset(visium_subset, anterior1_imagerow > 330 | anterior1_imagecol < 200, invert = TRUE))

sq_1_seurat_sctransform <- subset(visium_subset_seurat_sctransform, anterior1_imagerow > 330 | anterior1_imagecol < 200, invert = TRUE)

# SpatialDimPlot(sq_1)
# SpatialFeaturePlot(sq_1, features = "Astro", pt.size.factor = 1, crop = FALSE) +
  # ggplot2::scale_fill_gradientn(colours = MyPalette, limits = c(0.0, 1.001), breaks = c(0.0, 0.5, 1.0))



# Transform barcodes into xy
coordinates_seurat_sctransform <- sq_1_seurat_sctransform@images[["anterior1"]]@coordinates[,c("col","row")] %>% 
  rownames_to_column(var = "barcodes")
prediction_seurat_sctransform <- sq_1_seurat_sctransform@assays[["prediction.score.id"]]@data %>% t() %>% as_tibble(rownames = NA) %>% 
  rownames_to_column(var = "barcodes")
prediction_coord_seurat_sctransform <- right_join(prediction_seurat_sctransform, coordinates_seurat_sctransform, by = "barcodes")

prediction_coord_seurat_sctransform$dist_to_cortex <- apply(prediction_coord_seurat_sctransform, 1, function(x) dist_point_line(as.numeric(x[25:26]), slope = line[["coefficients"]][["col"]], intercept = line[["coefficients"]][["(Intercept)"]]) )

ggplot(prediction_coord_seurat_sctransform, mapping = aes(x = col, y = row, color=dist_to_cortex)) +
  geom_point() 

prediction_coord_seurat_sctransform$dist_to_cortex = prediction_coord_seurat_sctransform$dist_to_cortex / median(cortex_l6b_cell2loc$dist_to_cortex)

ggplot(prediction_coord_seurat_sctransform, mapping = aes(x = col, y = row, color=dist_to_cortex)) +
  geom_point() 


# Transform matrix into 0 and 1

long_format_sctransform <- pivot_longer(prediction_coord_seurat_sctransform, cols = 2:24, names_to = "celltype", values_to = "prob")
long_format_sctransform <- long_format_sctransform[long_format_sctransform$prob > 0.1,]
## Multiplicate rows with higher probabilities:
# How many replicates you want of each row
duptimes <- c(long_format_sctransform$prob %/% 0.1)

# Create an index of the rows you want with duplications
idx <- rep(1:nrow(long_format_sctransform), duptimes)

# Use that index to generate your new data frame
long_format_sctransform <- long_format_sctransform[idx,]

long_format_sctransform$celltype <- as.factor(long_format_sctransform$celltype)
long_format_sctransform$method <- "Seurat (SCTransform)"
ggplot(long_format_sctransform, aes(x=celltype, y=dist_to_cortex)) +
  geom_violin() +
  scale_y_reverse()


# spotlight ----
visium_subset_spotlight <- readRDS("./output/processed_data/mouse_brain/mouse_brain_spatial_cortex_lognorm_pca_clustered.rds")
res <- readRDS("./output/processed_data/mouse_brain/spotlight_with_raw_data/spotlight_lognorm.rds")
mat <- res$mat
visium_subset_spotlight[["spotlight"]] <- CreateAssayObject(t(mat))

DefaultAssay(visium_subset_spotlight) <- "spotlight"


# SpatialFeaturePlot(visium_subset, features = "Astro", pt.size.factor = 1, crop = FALSE) +
  # ggplot2::scale_fill_gradientn(colours = MyPalette, limits = c(0.0, 1.001), breaks = c(0.0, 0.5, 1.0))


# Divide barcodes into squares
# SpatialDimPlot(subset(visium_subset, anterior1_imagerow > 330 | anterior1_imagecol < 200, invert = TRUE))

sq_1_spotlight <- subset(visium_subset_spotlight, anterior1_imagerow > 330 | anterior1_imagecol < 200, invert = TRUE)

# SpatialDimPlot(sq_1)
# SpatialFeaturePlot(sq_1, features = "Astro", pt.size.factor = 1, crop = FALSE) +
#   ggplot2::scale_fill_gradientn(colours = MyPalette, limits = c(0.0, 1.001), breaks = c(0.0, 0.5, 1.0))



# Transform barcodes into xy
coordinates_spotlight <- sq_1_spotlight@images[["anterior1"]]@coordinates[,c("col","row")] %>% 
  rownames_to_column(var = "barcodes")
prediction_spotlight <- sq_1_spotlight@assays[["spotlight"]]@data %>% as_tibble(rownames=NA) %>% t() %>% as_tibble(rownames=NA) %>% rownames_to_column(var = "barcodes")
prediction_coord_spotlight <- right_join(prediction_spotlight, coordinates_spotlight, by = "barcodes")

prediction_coord_spotlight$dist_to_cortex <- apply(prediction_coord_spotlight, 1, function(x) dist_point_line(as.numeric(x[25:26]), slope = line[["coefficients"]][["col"]], intercept = line[["coefficients"]][["(Intercept)"]]) )

ggplot(prediction_coord_spotlight, mapping = aes(x = col, y = row, color=dist_to_cortex)) +
  geom_point() 

prediction_coord_spotlight$dist_to_cortex = prediction_coord_spotlight$dist_to_cortex / median(cortex_l6b_cell2loc$dist_to_cortex)

ggplot(prediction_coord_spotlight, mapping = aes(x = col, y = row, color=dist_to_cortex)) +
  geom_point() 


# Transform matrix into 0 and 1

long_format_spotlight <- pivot_longer(prediction_coord_spotlight, cols = 2:24, names_to = "celltype", values_to = "prob")
long_format_spotlight <- long_format_spotlight[long_format_spotlight$prob > 0.1,]

## Multiplicate rows with higher probabilities:
# How many replicates you want of each row
duptimes <- c(long_format_spotlight$prob %/% 0.1)

# Create an index of the rows you want with duplications
idx <- rep(1:nrow(long_format_spotlight), duptimes)

# Use that index to generate your new data frame
long_format_spotlight <- long_format_spotlight[idx,]

long_format_spotlight$celltype <- as.factor(long_format_spotlight$celltype)
long_format_spotlight$method <- "SPOTlight"
ggplot(long_format_spotlight, aes(x=celltype, y=dist_to_cortex)) +
  geom_violin() +
  scale_y_reverse()





# Ground truth MERFISH ----
gt <- read.csv("./data/mouse_brain/MERFISH_cell_depth.csv")
gt$method <- "MERFISH"
gt$subclass[gt$subclass == "L5 ET"] <- "L5 PT"
gt$subclass[gt$subclass == "L4/5 IT"] <- "L4"
gt$subclass[gt$subclass == "PVM"] <- "Macrophage"
gt$subclass[gt$subclass == "L5/6 NP"] <- "NP"

gt <- gt[gt$subclass %in% c("L2/3 IT","L4","L5 IT", "L5 PT",
                            "L6 IT", "L6 CT", "L6b", "Lamp5",
                            "Sncg","Vip","Sst","Pvalb",
                            "Astro","Oligo","Endo","VLMC",
                            "SMC","Peri","Macrophage","NP",
                            "Serpinf1","CR","Meis2"),]
ggplot(gt, aes(x=subclass, y=normalized_depth)) +
  geom_violin() +
  scale_y_reverse()


long_format <- bind_rows(long_format_lognorm, long_format_sctransform, long_format_spotlight, long_format_cell2loc)
pal <- readRDS("./output/processed_data/mouse_brain/spatial_scatterpie_palette.rds")

ggplot(long_format, aes(x=method, y=dist_to_cortex, fill = celltype)) +
  geom_violin() +
  facet_wrap('celltype', ncol = 6) +
  scale_y_reverse() +
  scale_fill_manual(
    values = pal,
    breaks = names(pal)
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5)) +
  scale_x_discrete(limits = c("Seurat (LogNorm)", "Seurat (SCTransform)", "SPOTlight", "Cell2location")) +
  ylab("Cortical depth (normalized)") +
  xlab("") +
  geom_hline(yintercept=0, color = "gray", linetype="dashed", size=0.5) +
  geom_hline(yintercept=1, color = "gray", linetype="dashed", size=0.5) +
  theme(legend.position = "none") 

ggsave("./output/images/mouse_brain_integration/violin_plots_prob01cutoff.png", width = 7, height = 5)



colnames(long_format)[colnames(long_format) == "dist_to_cortex"] <- "normalized_depth"
colnames(gt)[colnames(gt) == "subclass"] <- "celltype"
long_format <- bind_rows(long_format, gt)
long_format$celltype = factor(long_format$celltype, 
                              levels=c("L2/3 IT","L4","L5 IT", "L5 PT",
                                       "L6 IT", "L6 CT", "L6b", "Lamp5",
                                       "Sncg","Vip","Sst","Pvalb","Astro",
                                       "Oligo","Endo","VLMC","SMC","Peri",
                                       "Macrophage","NP","Serpinf1","CR","Meis2"))

ggplot(long_format, aes(x=method, y=normalized_depth, fill = celltype)) +
  geom_violin() +
  facet_wrap('celltype', ncol = 4) +
  scale_y_reverse() +
  scale_fill_manual(
    values = pal,
    breaks = names(pal)
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5)) +
  scale_x_discrete(limits = c("MERFISH","Seurat (LogNorm)", "Seurat (SCTransform)", "SPOTlight", "Cell2location")) +
  ylab("Cortical depth (normalized)") +
  xlab("") +
  geom_hline(yintercept=0, color = "gray", linetype="dashed", size=0.5) +
  geom_hline(yintercept=1, color = "gray", linetype="dashed", size=0.5) +
  theme(legend.position = "none") 

ggsave("./output/images/mouse_brain_integration/violin_plots_prob01cutoff_merfish.png", width = 5, height = 8)

ggplot(long_format, aes(x=method, y=normalized_depth, fill = celltype)) +
  geom_violin() +
  facet_wrap('celltype', ncol = 6) +
  scale_y_reverse() +
  scale_fill_manual(
    values = pal,
    breaks = names(pal)
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5)) +
  scale_x_discrete(limits = c("MERFISH","Seurat (LogNorm)", "Seurat (SCTransform)", "SPOTlight", "Cell2location")) +
  ylab("Cortical depth (normalized)") +
  xlab("") +
  geom_hline(yintercept=0, color = "gray", linetype="dashed", size=0.5) +
  geom_hline(yintercept=1, color = "gray", linetype="dashed", size=0.5) +
  theme(legend.position = "none") 
ggsave("./output/images/mouse_brain_integration/violin_plots_prob01cutoff_merfish_6col.png", width = 9, height = 6)


# Plot cell type counts
df2 <- gt %>% group_by(celltype) %>% mutate(count_name_occurr = n())

ggplot(data=df2, aes(x=reorder(celltype,-count_name_occurr), fill = celltype)) +
  geom_bar(stat="count") +
  scale_fill_manual(
    values = pal,
    breaks = names(pal)
  ) +
  ylab("Counts") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        legend.position = "none")

ggsave("./output/images/mouse_brain_integration/merfish_celltype_counts.png", width = 4, height = 3)


