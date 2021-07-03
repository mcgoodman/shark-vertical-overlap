
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

## Define overlap functions ---------------------------------------------------

#' Bhattacharya's coefficient
#' Measures similarity between two probability distributions
#' @param x A vector of depths for the first species
#' @param y A vector of depths for the second species
#' @param bin Logical; use if data are not already binned by depth
#' @param breaks The breaks used to bin depths for both species, if necessary
#' @citation Bhattacharyya, A. (1943). Bulletin of the Calcutta Mathematical Society, 35, 99–110
bhattacharya_depth <- function(x, y, bin = FALSE, breaks) {
  if (bin) {
    x <- cut(x, breaks)
    y <- cut(y, breaks)
  } 
  xdens <- as.numeric(table(x)/length(x))
  ydens <- as.numeric(table(y)/length(y))
  sum(sqrt(xdens * ydens))
}

#' Schoener's D
#' Measures how equally two species use space relative to availability
#' @param x A vector of depths for the first species
#' @param y A vector of depths for the second species
#' @param bin Logical; use if data are not already binned by depth
#' @param breaks The breaks used to bin depths for both species, if necessary
#' @citation Schoener, T. W. (1970). Ecology, 51, 408–418
schoener_depth <- function(x, y, bin = FALSE, breaks) {
  if (bin) {
    x <- cut(x, breaks)
    y <- cut(y, breaks)
  } 
  xdens <- as.numeric(table(x)/length(x))
  ydens <- as.numeric(table(y)/length(y))
  1 - 0.5 * sum(abs(xdens - ydens))
}

#' Euclidian Distance
#' Distance between proportion depth utilization in n-bin space
#' Provided for comparability with Madigan et al. (2020)
#' @param x A vector of depths for the first species
#' @param y A vector of depths for the second species
#' @param bin Logical; use if data are not already binned by depth
#' @param breaks The breaks used to bin depths for both species, if necessary
#' @citation Madigan et al. (2020). ICES Journal of Marine Science
euclidian_depth <- function(x, y, bin = FALSE, breaks) {
  if (bin) {
    x <- cut(x, breaks)
    y <- cut(y, breaks)
  } 
  xdens <- as.numeric(table(x)/length(x))
  ydens <- as.numeric(table(y)/length(y))
  sqrt(sum((xdens - ydens)^2))
}

## Calculate overlap between species ------------------------------------------

## quick function to round number up by 10
ceiling10 <- function(x) if (x %% 10 == 0) return(x) else return((10 - x %% 10) + x)

## define depth breaks
depth_max <- ceiling10(max(med_depths$Depth))
depth_breaks <- seq(0, depth_max, by = 10)

## bin depth values
med_depths$Depth[med_depths$Depth < 0] <- 0
med_depths$Depth_bin <- cut(med_depths$Depth, depth_breaks, right = FALSE)

## Summarize number of individuals by species, select those with 5 or more individuals
species_counts <- med_depths %>% group_by(Species) %>% summarize(n = length(unique(Individual_ID)))
species_sub <- species_counts$Species[species_counts$n >= 5]
n_sp <- length(species_sub)

## Combinations of species for which to calculate overlap metrics
species_combn <- expand.grid(species1 = species_sub, species2 = species_sub)
species_combn$bhattacharya <- species_combn$schoener <- species_combn$euclidian <- NA

## Compute overlap metrics for each pair of species
for (i in 1:nrow(species_combn)) {
  x <- med_depths$Depth_bin[med_depths$Species == species_combn$species1[i]]
  y <- med_depths$Depth_bin[med_depths$Species == species_combn$species2[i]]
  species_combn$bhattacharya[i] <- bhattacharya_depth(x, y)
  species_combn$schoener[i] <- schoener_depth(x, y)
  species_combn$euclidian[i] <- euclidian_depth(x, y)
}

write.csv(species_combn, "./output/overlap_metrics.csv")

## Hierarchical clustering ----------------------------------------------------

## Create distance matrices for each metric
## Subtract Bhattacharya and Schoener's D from 1 to change from similarity metric to dissimilarity
dists <- list(
  bhattacharya = as.dist(1 - matrix(species_combn[,"bhattacharya"], n_sp, n_sp)),
  schoener = as.dist(1 - matrix(species_combn[,"schoener"], n_sp, n_sp)),
  euclidian = as.dist(matrix(species_combn[,"euclidian"], n_sp, n_sp))
)

## Hierarchical clustering for each metric
clusters <- lapply(dists, hclust, method = "complete")
