
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

## empty list to store median depth data frames
med_depths <- vector("list", length(ts_files))

## Read in tag data, store daily median depths
for (i in seq_along(ts_files)) {
  
  print(paste0(i, " (", round(i/length(ts_files), 3) * 100, "%)"))
  
  tag_data <- read.csv(ts_files[i])
  
  med_depths[[i]] <- tag_data %>%
    group_by(Species, Individual_ID, Date) %>% 
    summarize(Depth = median(Depth), .groups = "drop")
  
}

## Bind median depth summaries together 
med_depths <- do.call("rbind", med_depths)
med_depths <- na.omit(med_depths) ## 3 NA date values

write.csv(med_depths, "./output/median_depths.csv")

## Format and summarize data ------------------------------------------

## define depth breaks
ceiling10 <- function(x) if (x %% 10 == 0) return(x) else return((10 - x %% 10) + x)
depth_max <- ceiling10(max(med_depths$Depth))
depth_breaks <- seq(0, depth_max, by = 10)

## bin depth values
med_depths$Depth[med_depths$Depth < 0] <- 0
med_depths$Depth_bin <- cut(med_depths$Depth, depth_breaks, right = FALSE)

## Summarize number of individuals by species, select those with 5 or more individuals
species_counts <- med_depths %>% group_by(Species) %>% summarize(n = length(unique(Individual_ID)))
species_sub <- species_counts$Species[species_counts$n >= 5]
n_sp <- length(species_sub)

## Summarize depth distribution for each species as proportion by depth bin
depth_binned <- med_depths %>% 
  filter(Species %in% species_sub) %>% 
  group_by(Species) %>% 
  mutate(n_obs = n()) %>% 
  group_by(Species, Depth_bin, n_obs) %>% 
  tally() %>% 
  mutate(p = n/n_obs) %>% select(-n_obs, -n) %>%
  ungroup() %>% 
  complete(Species, Depth_bin, fill = list(p = 0))

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
  scale_fill_viridis_c(option = "magma") + 
  theme_minimal() + 
  theme(axis.title = element_blank(),
        axis.text.x = element_text(color = "black", angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(color = "black"),
        legend.position = c(0.75, 0.75), 
        legend.direction = "horizontal") + 
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

ggsave("./output/cluster_number_bhattacharyya.png")

## Determine optimal number of clusters using Schoener dissimilarity
fviz_nbclust(
  depth_binned_wide, diss = dists$schoener, 
  FUNcluster = hcut, method = "wss", hc_method = "ward.D2"
) + ggtitle("Schoener dissimilarity")

ggsave("./output/cluster_number_schoener.png")

## Determine optimal number of clusters using Euclidian distance
fviz_nbclust(
  depth_binned_wide, diss = dists$euclidian, 
  FUNcluster = hcut, method = "wss", hc_method = "ward.D2"
) + ggtitle("Euclidian distance")

ggsave("./output/cluster_number_euclidian.png")

## Compute Bhattacharyya's coefficient and save output
bhattacharyya_matrix <- proxy::dist(
  depth_binned_wide, 
  method = function(x, y) bhattacharyya
)

bhattacharyya_matrix <- proxy::dist(
  depth_binned_wide, bhattacharyya, convert_similarities = FALSE
)


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
  d_colors <- c("grey40", pnw_palette("Bay", n = k))
  
  ### Convert hierarchical clustering output to data frame
  cluster_ggdata <- dendro_data_k(cluster_data, k)
  ## Shift and rescale dendrogram for plotting beside heatmap
  cluster_ggdata <- cluster_ggdata$segments %>% 
    mutate(y = tree_spacing + max_depth + ((y/max(yend)) * tree_depth),
           yend = tree_spacing + max_depth + ((yend/max(yend)) * tree_depth), 
           line = factor(clust == 0))
  
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
    geom_segment(aes(x = y, y = x, xend = yend, yend = xend, linetype = line, color = clust), 
                 data = cluster_ggdata, size = 0.8, inherit.aes = FALSE, 
                 show.legend = FALSE) +
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
bhattacharyya_heatmap <- cluster_heatmap(depth_binned, clusters$bhattacharyya, k = 6)

ggsave("./output/bhattacharyya_heatmap.png", bhattacharyya_heatmap, height = 5, width = 8, 
       units = "in", dpi = 500)


## Plot: Clustered by Bhattacharyya (dis)similarity; all depth bins displayed
bhattacharyya_heatmap_1500m <- cluster_heatmap(depth_binned, clusters$bhattacharyya, k = 6, 
                                               max_depth = 1510, tree_depth = 1000, 
                                               tile_border = 0, tree_spacing = 10) + 
  scale_x_continuous(breaks = seq(0, 1500, 100)) + 
  theme(axis.text.x = element_text(color = "black", angle = 90, hjust = 1, vjust = 0.5))

ggsave("./output/bhattacharyya_heatmap_1500m.png", bhattacharyya_heatmap_1500m, height = 6, 
       width = 9, units = "in", dpi = 500)


## Plot: Clustered by Schoener's D; depth displayed up to 100m
schoener_heatmap <- cluster_heatmap(depth_binned, clusters$schoener, k = 6)

ggsave("./output/schoener_heatmap.png", schoener_heatmap, height = 5, width = 8, 
       units = "in", dpi = 500)


## Plot: Clustered by Schoener's D; all depth bins displayed
schoener_heatmap_1500m <- cluster_heatmap(depth_binned, clusters$schoener, k = 6, 
                                          max_depth = 1510, tree_depth = 1000, 
                                          tile_border = 0, tree_spacing = 10) + 
  scale_x_continuous(breaks = seq(0, 1500, 100)) + 
  theme(axis.text.x = element_text(color = "black", angle = 90, hjust = 1, vjust = 0.5))

ggsave("./output/schoener_heatmap_1500m.png", schoener_heatmap_1500m, height = 5, width = 8, 
       units = "in", dpi = 500)


## Plot: Clustered by Euclidian distance; depth displayed up to 100m
euclidian_heatmap <- cluster_heatmap(depth_binned, clusters$euclidian, k = 6)

ggsave("./output/euclidian_heatmap.png", euclidian_heatmap, height = 5, width = 8, 
       units = "in", dpi = 500)

## Plot: Clustered by Euclidian distance; all depth bins diplayed
euclidian_heatmap_1500m <- cluster_heatmap(depth_binned, clusters$euclidian, k = 6, 
                                           max_depth = 1510, tree_depth = 1000, 
                                           tile_border = 0, tree_spacing = 10) + 
  scale_x_continuous(breaks = seq(0, 1500, 100)) + 
  theme(axis.text.x = element_text(color = "black", angle = 90, hjust = 1, vjust = 0.5))

ggsave("./output/euclidian_heatmap_1500m.png", euclidian_heatmap_1500m, height = 5, width = 8, 
       units = "in", dpi = 500)
