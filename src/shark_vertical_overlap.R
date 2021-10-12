
###############################################################################
##' @title Shark vertical niche overlap
##' @author Maurice Goodman
##' @email goodman.maurice@gmail.com
##' @date July 2021
###############################################################################

## Load packages and data -----------------------------------------------------

library(tidyverse)
library(cluster)
library(factoextra)
library(ggdendro)
library(proxy)
library(PNWColors)
library(ggtext)

## Read in individual tag metadata
metadata <- read.csv("./data/Individual_tag_metadata.csv")

## Read in species metadata
species_metadata <- read.csv("./data/species_metadata.csv")

## list time series files
ts_files <- list.files("./data/time series data", full.names = TRUE)

## empty list to store proportion time at depth for each individual
tag_depths_binned <- vector("list", length(ts_files))

## Read in tag data, store proportion time at depths
for (i in seq_along(ts_files)) {
  
  print(paste0(i, " (", round(i/length(ts_files), 3) * 100, "%)"))
  
  tag_data <- read.csv(ts_files[i])
  
  tag_data$Depth[tag_data$Depth < 0] <- 0
  
  tag_depths_binned[[i]] <- tag_data %>%
    mutate(Depth_bin = floor(Depth / 10) * 10) %>% 
    group_by(Species, Individual_ID, Depth_bin) %>% 
    summarize(p = n()/nrow(tag_data), .groups = "drop")
}

## Bind individual proportion time at depth data together
tag_depths_binned <- do.call("rbind", tag_depths_binned)

write.csv(tag_depths_binned, "./output/individual_proportion_time_at_depth.csv")

## Format and summarize data ------------------------------------------

## define depth breaks
depth_max <- max(tag_depths_binned$Depth_bin)
depth_breaks <- seq(0, depth_max, by = 10)
tag_depths_binned$Depth_bin <- factor(tag_depths_binned$Depth_bin, levels = depth_breaks, ordered = TRUE)

## Summarize number of individuals by species, select those with 5 or more individuals
species_counts <- tag_depths_binned %>% group_by(Species) %>% summarize(n = length(unique(Individual_ID)))
species_sub <- species_counts$Species[species_counts$n >= 5]
n_sp <- length(species_sub)

## Summarize depth distribution for each species, averaging across individuals
depth_binned <- tag_depths_binned %>%
  complete(nesting(Species, Individual_ID), Depth_bin, fill = list(p = 0)) %>% 
  filter(Species %in% species_sub) %>% 
  group_by(Species, Depth_bin) %>%
  summarize(p = mean(p))
  
write.csv(depth_binned, "./output/depth_binned.csv")

## Pivot to wide format for calculating distance matrices
depth_binned_wide <- depth_binned %>% 
  pivot_wider(names_from = Depth_bin, values_from = p, values_fill = 0) %>% 
  column_to_rownames("Species")

## Similarity and Dissimilarity matrices --------------------------------------

## Bhattacharyya's coefficient - dissimilarity function
## Measures similarity between two probability distributions
## If dist = TRUE, returns a dissimilarity
bhattacharyya <- function(x, y, dist = TRUE) {
  if (dist) 1 - sum(sqrt(x * y)) else sum(sqrt(x * y))
} 

## Schoener's D - dissimilarity function
## Measures how equally two species use space relative to availability
schoener <- function(x, y, dist = TRUE) {
  if (dist) 1 - (1 - 0.5 * sum(abs(x - y))) else 1 - 0.5 * sum(abs(x - y))
}

## Create distance matrices for each metric
dists <- list(
  bhattacharyya = proxy::dist(depth_binned_wide, method = bhattacharyya),
  schoener = proxy::dist(depth_binned_wide, method = schoener), 
  euclidian = proxy::dist(depth_binned_wide)
)

saveRDS(dists, "./output/distance_matrices.rds")

## Create Bhattacharyya and Schoener's D similarity matrices
simils <- list(
  bhattacharyya = proxy::simil(depth_binned_wide, method = bhattacharyya, dist = FALSE),
  schoener = proxy::simil(depth_binned_wide, method = schoener, dist = FALSE)
)

saveRDS(simils, "./output/similarity_matrices.rds")

## Format Bhattacharyya similarity matrix for plotting
bhattacharyya_long <- as.matrix(simils$bhattacharyya)
bhattacharyya_long[upper.tri(bhattacharyya_long)] <- NA
bhattacharyya_long <- data.frame(
  species1 = colnames(bhattacharyya_long)[col(bhattacharyya_long)], 
  species2 = rownames(bhattacharyya_long)[row(bhattacharyya_long)], 
  bhattacharyya = c(bhattacharyya_long)
)
bhattacharyya_long <- na.omit(bhattacharyya_long)

## Plot Bhattacharyya similarity matrix
bhattacharyya_plot <- bhattacharyya_long %>% 
  mutate(species2 = forcats::fct_rev(species2)) %>% 
  ggplot(aes(species1, species2, fill = bhattacharyya)) + 
  geom_tile(color = "black", size = 0.6) + 
  scale_fill_viridis_c(option = "magma", breaks = seq(0, 1, 0.25), limits = c(0, 1), 
                       labels = c("0", "0.25", "0.5", "0.75", "1")) + 
  theme_minimal() + 
  theme(axis.title = element_blank(),
        axis.text.x = element_text(color = "black", face = "bold", angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(color = "black", face = "bold"),
        legend.position = c(0.75, 0.76), 
        legend.direction = "horizontal", 
        plot.background = element_rect(fill = "White")) + 
  guides(fill = guide_colorbar("Bhattacharyya's coefficient", title.position = "top", 
                               title.hjust = 0.5, barwidth = 10, barheight = 0.7, 
                               ticks.linewidth = 1, ticks.colour = "black", 
                               frame.colour = "black", frame.linewidth = 1)) + 
  coord_fixed()


ggsave("./output/bhattacharyya_matrix.png", bhattacharyya_plot, height = 7, 
       width = 7, units = "in", dpi = 500)


## Hierarchical clustering ----------------------------------------------------

## Hierarchical clustering for each metric
clusters <- lapply(dists, hclust, method = "ward.D2")

## Determine optimal number of clusters using Bhattacharyya dissimilarity
fviz_nbclust(
  depth_binned_wide, diss = dists$bhattacharyya, 
  FUNcluster = hcut, method = "wss", hc_method = "ward.D2"
) + ggtitle("Bhattacharyya dissimilarity")

ggsave("./output/cluster_number_bhattacharyya.png", height = 3, width = 5, units = "in")

## Determine optimal number of clusters using Schoener dissimilarity
fviz_nbclust(
  depth_binned_wide, diss = dists$schoener, 
  FUNcluster = hcut, method = "wss", hc_method = "ward.D2"
) + ggtitle("Schoener dissimilarity")

ggsave("./output/cluster_number_schoener.png", height = 3, width = 5, units = "in")

## Determine optimal number of clusters using Euclidian distance
fviz_nbclust(
  depth_binned_wide, diss = dists$euclidian, 
  FUNcluster = hcut, method = "wss", hc_method = "ward.D2"
) + ggtitle("Euclidian distance")

ggsave("./output/cluster_number_euclidian.png", height = 3, width = 5, units = "in")


## Plot heat maps and clusters ------------------------------------------------

## Habitats corresponding to each species
species_metadata$habitat_short <- ifelse(
  trimws(species_metadata$Habitat) == "Oceanic", "o", 
  ifelse(trimws(species_metadata$Habitat) == "Coastal", "c", "t")
)
species_habitats <- setNames(
  species_metadata$habitat_short, 
  trimws(species_metadata$Species..common.)
)

## Function to convert hierarchical clustering output to data frame
## With column indicating the cluster each line segment belongs to
## Function is from atrebas.github.io/post/2019-06-08-lightweight-dendrograms/
dendro_data_k <- function(hc, k) {
  
  require("ggdendro")
  
  hcdata    <-  ggdendro::dendro_data(hc, type = "rectangle")
  seg       <-  hcdata$segments
  labclust  <-  cutree(hc, k)[hc$order]
  segclust  <-  rep(0L, nrow(seg))
  heights   <-  sort(hc$height, decreasing = TRUE)
  height    <-  mean(c(heights[k], heights[k - 1L]), na.rm = TRUE)
  
  for (i in 1:k) {
    xi      <-  hcdata$labels$x[labclust == i]
    idx1    <-  seg$x    >= min(xi) & seg$x    <= max(xi)
    idx2    <-  seg$xend >= min(xi) & seg$xend <= max(xi)
    idx3    <-  seg$yend < height
    idx     <-  idx1 & idx2 & idx3
    segclust[idx] <- i
  }
  
  idx                    <-  which(segclust == 0L)
  segclust[idx]          <-  segclust[idx + 1L]
  hcdata$segments$clust  <-  segclust
  hcdata$segments$line   <-  as.integer(segclust < 1L)
  hcdata$labels$clust    <-  labclust
  
  hcdata
}

## Function to plot heatmap of species depth distributions up to `max_depth`
## With dendrogram plot of cluster results, rescaled to `tree_depth`
## Pass depth_binned as `depth_data`
## Past results of hierarchical clustering (i.e. clusters$...) as cluster_data
cluster_heatmap <- function(depth_data, cluster_data, k, habitats = species_habitats, 
                            max_depth = 100, tree_depth = 80, tile_border = 0.6, 
                            tree_spacing = 1) {
  
  require("tidyverse")
  require("ggtext")
  
  ## reorder species according to hierarchical clustering results
  clust_order <- cluster_data$labels[cluster_data$order]
  habitats <- habitats[clust_order]
  
  ## Colors for clusters on dendrogram
  if (k <= 4) {
    d_colors <- c("grey50", pnw_palette("Bay", n = k + 2)[-c(1, k + 2)])
  } else {
    d_colors <- c("grey50", pnw_palette("Bay", n = k))
  }
  
  ### Convert hierarchical clustering output to data frame
  cluster_ggdata <- dendro_data_k(cluster_data, k)
  ## Shift and rescale dendrogram for plotting beside heatmap
  cluster_ggdata <- cluster_ggdata$segments %>% 
    mutate(y = tree_spacing + max_depth + ((y/max(yend)) * tree_depth),
           yend = tree_spacing + max_depth + ((yend/max(yend)) * tree_depth))
  
  ## Select vertices at which to place cluster labels
  cluster_labels <- cluster_ggdata %>% filter(clust != 0) %>% 
    group_by(clust) %>% filter(yend == max(yend)) %>% ungroup() %>% 
    mutate(rank = order(xend))
  
  ## Order and re-number clusters by vertical position in plot
  cluster_ggdata$clust <- as.numeric(factor(
    cluster_ggdata$clust, 
    levels = cluster_labels$clust[cluster_labels$rank]
  ))
  cluster_ggdata$clust[is.na(cluster_ggdata$clust)] <- 0
  cluster_ggdata$clust <- factor(cluster_ggdata$clust)
  cluster_labels$rank <- factor(cluster_labels$rank)
  
  ## Change depth bins to numeric variable so ticks can be placed at 10m intervals
  depth_data <- depth_data %>% 
    mutate(Species = factor(Species, levels = clust_order), 
           y = as.numeric(Species)) %>%
    mutate(Depth = (as.numeric(Depth_bin) - 0.5) * 10)
  depth_data$Species_label <- factor(paste0(
    "**", depth_data$Species, "** *", habitats[depth_data$Species], "*"
  ), levels = paste0("**", clust_order, "** *", habitats, "*")
  )
  
  ## Plot heatmap and dendrogram
  depth_data %>% 
    filter(Depth <= max_depth) %>% 
    ggplot(aes(Depth, Species_label, fill = p)) + 
    geom_tile(color = "black", size = tile_border) + 
    geom_segment(aes(x = y, y = x, xend = yend, yend = xend, color = clust), 
                 data = cluster_ggdata, size = 0.8, inherit.aes = FALSE, 
                 show.legend = FALSE, lineend = "round") +
    geom_text(aes(x = yend, y = xend, label = rank, color = rank), nudge_x = 0.038 * tree_depth,
              data = cluster_labels, inherit.aes = FALSE, show.legend = FALSE, 
              fontface = "bold") + 
    scale_fill_viridis_c(option = "inferno") + 
    coord_cartesian(expand = FALSE, clip = "off") + 
    labs(x = "Depth (m)") + 
    scale_x_continuous(breaks = seq(0, max_depth, 10)) + 
    scale_color_manual(values = d_colors) +  
    guides(fill = guide_colorbar("proportion time at depth", title.position = "top", 
                                 title.hjust = 0.5, barwidth = 15, barheight = 0.5, 
                                 ticks.linewidth = 1, ticks.colour = "black", 
                                 frame.colour = "black", frame.linewidth = 1)) + 
    theme(legend.position = "top", 
          axis.title.y = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_line(color = "black"), 
          axis.text.x = element_text(color = "black"),
          axis.text.y = element_markdown(color = "black"),
          legend.justification = c(0, 0), 
          legend.text = element_text(size = 8), 
          legend.title = element_text(size = 10))
  
}

## Plot: Clustered by Bhattacharyya (dis)similarity; depth displayed up to 100 m
bhattacharyya_heatmap <- cluster_heatmap(depth_binned, clusters$bhattacharyya, k = 4) + 
  viridis::scale_fill_viridis(option = "rocket")
  scale_fill_viridis_c(option = "cividis", limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1))

ggsave("./output/bhattacharyya_heatmap.png", bhattacharyya_heatmap, height = 5, width = 8, 
       units = "in", dpi = 500)


## Plot: Clustered by Bhattacharyya (dis)similarity; all depth bins displayed
bhattacharyya_heatmap_full <- cluster_heatmap(depth_binned, clusters$bhattacharyya, k = 4, 
                                               max_depth = depth_max, tree_depth = 1000, 
                                               tile_border = 0, tree_spacing = 10) + 
  scale_x_continuous(breaks = seq(0, depth_max, 100)) + 
  theme(axis.text.x = element_text(color = "black", angle = 90, hjust = 1, vjust = 0.5))

ggsave("./output/bhattacharyya_heatmap_full.png", bhattacharyya_heatmap_full, height = 6, 
       width = 9, units = "in", dpi = 500)


## Plot: Clustered by Schoener's D; depth displayed up to 100m
schoener_heatmap <- cluster_heatmap(depth_binned, clusters$schoener, k = 4)

ggsave("./output/schoener_heatmap.png", schoener_heatmap, height = 5, width = 8, 
       units = "in", dpi = 500)


## Plot: Clustered by Schoener's D; all depth bins displayed
schoener_heatmap_full <- cluster_heatmap(depth_binned, clusters$schoener, k = 4, 
                                          max_depth = depth_max, tree_depth = 1000, 
                                          tile_border = 0, tree_spacing = 10) + 
  scale_x_continuous(breaks = seq(0, depth_max, 100)) + 
  theme(axis.text.x = element_text(color = "black", angle = 90, hjust = 1, vjust = 0.5))

ggsave("./output/schoener_heatmap_full.png", schoener_heatmap_full, height = 5, width = 8, 
       units = "in", dpi = 500)


## Plot: Clustered by Euclidian distance; depth displayed up to 100m
euclidian_heatmap <- cluster_heatmap(depth_binned, clusters$euclidian, k = 4)

ggsave("./output/euclidian_heatmap.png", euclidian_heatmap, height = 5, width = 8, 
       units = "in", dpi = 500)

## Plot: Clustered by Euclidian distance; all depth bins diplayed
euclidian_heatmap_full <- cluster_heatmap(depth_binned, clusters$euclidian, k = 4, 
                                           max_depth = depth_max, tree_depth = 1000, 
                                           tile_border = 0, tree_spacing = 10) + 
  scale_x_continuous(breaks = seq(0, depth_max, 100)) + 
  theme(axis.text.x = element_text(color = "black", angle = 90, hjust = 1, vjust = 0.5))

ggsave("./output/euclidian_heatmap_full.png", euclidian_heatmap_full, height = 5, width = 8, 
       units = "in", dpi = 500)
