
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

## Read in metadata
metadata <- read.csv("./data/Individual_tag_metadata.csv")

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
clusters <- lapply(dists, hclust, method = "complete")

## Plot Heat maps and clusters ------------------------------------------------

depth_binned %>% 
  mutate(Species = factor(Species, levels = x$Species)) %>% 
  filter(as.numeric(Depth_bin) < 11) %>% 
  ggplot(aes(Depth_bin, Species, fill = p)) + geom_tile() + 
  scale_fill_viridis_c(option = "inferno") + 
  theme(legend.position = "top", 
        axis.title.y = element_blank()) + 
  coord_cartesian(expand = FALSE) + 
  labs(x = "Depth bin")
