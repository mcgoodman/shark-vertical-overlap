
###############################################################################
##' @title Shark vertical niche overlap
##' @author Maurice Goodman
##' @email goodman.maurice@gmail.com
##' @date July 2021
###############################################################################

## Load packages and data -----------------------------------------------------

library(tidyverse)
library(cluster)
library(proxy)

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

## Pivot to wide format for calculating distance matrices
depth_binned_wide <- depth_binned %>% 
  pivot_wider(names_from = Depth_bin, values_from = p, values_fill = 0) %>% 
  column_to_rownames("Species")

## Dissimilarity matrices & Hierarchical clustering ---------------------------

## Bhattacharya's coefficient - dissimilarity function
## Measures similarity between two probability distributions
bhattacharya <- function(x, y) 1 - sum(sqrt(x * y))

## Schoener's D - dissimilarity function
## Measures how equally two species use space relative to availability
schoener <- function(x, y) 1 - (1 - 0.5 * sum(abs(x - y)))

## Create distance matrices for each metric
dists <- list(
  bhattacharya = dist(depth_binned_wide, method = bhattacharya),
  schoener = dist(depth_binned_wide, method = schoener), 
  euclidian = dist(depth_binned_wide)
)

saveRDS(dists, "./output/distance_matrices.rds")

## Hierarchical clustering for each metric
clusters <- lapply(dists, hclust, method = "ward.D2")

## Plot Heat maps and clusters ------------------------------------------------

## Function to plot heatmap of species depth distributions up to `max_depth`
## With dendrogram plot of cluster results, rescaled to `tree_depth`
## Pass depth_binned as `depth_data`
## Past results of hierarchical clustering (i.e. clusters$...) as cluster_data
cluster_heatmap <- function(depth_data, cluster_data, species_habitats, habitat_colors,
                            max_depth = 100, tree_depth = 80) {
  
  require("tidyverse")
  require("ggdendro")
  
  clust_order <- cluster_data$labels[cluster_data$order]
  species_habitats <- species_habitats[clust_order]
  
  cluster_ggdata <- ggdendro::dendro_data(cluster_data)
  cluster_ggdata <- cluster_ggdata$segments %>% 
    mutate(y = 1 + max_depth + ((y/max(yend)) * tree_depth),
           yend = 1 + max_depth + ((yend/max(yend)) * tree_depth))
  
  depth_data <- depth_data %>% 
    mutate(Species = factor(Species, levels = clust_order), 
           y = as.numeric(Species)) %>%
    mutate(Depth = (as.numeric(Depth_bin) - 0.5) * 10)
  
  depth_data %>% 
    filter(Depth <= max_depth) %>% 
    ggplot(aes(Depth, Species, fill = p)) + 
    geom_tile(color = "black", size = 0.6) + 
    geom_segment(aes(x = y, y = x, xend = yend, yend = xend), 
                 data = cluster_ggdata, size = 0.8, inherit.aes = FALSE) +
    scale_fill_viridis_c(option = "inferno") + 
    coord_cartesian(expand = FALSE, clip = "off") + 
    labs(x = "Depth (m)") + 
    scale_x_continuous(breaks = seq(0, max_depth, 10)) + 
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
          axis.text = element_text(color = "black"), 
          legend.justification = c(0, 0), 
          legend.text = element_text(size = 8), 
          legend.title = element_text(size = 10))
  
}

## Labels for species nodes
species_habitats <- setNames(
  trimws(species_metadata$Habitat), 
  trimws(species_metadata$Species..common.)
)

## Cluster by Bhattacharya (dis)similarity
bhattacharya_heatmap <- cluster_heatmap(depth_binned, clusters$bhattacharya)

ggsave("./output/bhattacharya_heatmap.png", height = 5, width = 8, 
       units = "in", dpi = 500)


## Cluster by Schoener's D
schoener_heatmap <- cluster_heatmap(depth_binned, clusters$schoener)

ggsave("./output/schoener_heatmap.png", height = 5, width = 8, 
       units = "in", dpi = 500)


## Cluster by Euclidian distance
euclidian_heatmap <- cluster_heatmap(depth_binned, clusters$euclidian)

ggsave("./output/euclidian_heatmap.png", height = 5, width = 8, 
       units = "in", dpi = 500)
