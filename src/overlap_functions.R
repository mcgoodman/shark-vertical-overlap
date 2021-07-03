
###############################################################################
##' @title Functions for calculating overlap
##' @author Maurice Goodman
##' @email goodman.maurice@gmail.com
##' @date July 2021
###############################################################################

## These functions are for using with raw depths vectors for pairs of species
## The functions used in the main script are simplified because they are applied
## to data pre-summarized into proportions by depth bin

#' Bhattacharya's coefficient
#' Measures similarity between two probability distributions
#' @param x A vector of depths for the first species
#' @param y A vector of depths for the second species
#' @param bin Logical; use if data are not already binned by depth
#' @param breaks The breaks used to bin depths for both species, if necessary
#' @param dissimilarity If TRUE, returns dissimilarity instead of similarity
#' @citation Bhattacharyya, A. (1943). Bulletin of the Calcutta Mathematical Society, 35, 99–110
bhattacharya_depth <- function(x, y, bin = FALSE, breaks, dissimilarity = FALSE) {
  if (bin) {
    x <- cut(x, breaks)
    y <- cut(y, breaks)
  } 
  xdens <- as.numeric(table(x)/length(x))
  ydens <- as.numeric(table(y)/length(y))
  if (dissimilarity) {
    return(1 - sum(sqrt(xdens * ydens)))
  } else {
    return(sum(sqrt(xdens * ydens)))
  }
}

#' Schoener's D
#' Measures how equally two species use space relative to availability
#' @param x A vector of depths for the first species
#' @param y A vector of depths for the second species
#' @param bin Logical; use if data are not already binned by depth
#' @param breaks The breaks used to bin depths for both species, if necessary
#' @param dissimilarity If TRUE, returns dissimilarity instead of similarity
#' @citation Schoener, T. W. (1970). Ecology, 51, 408–418
schoener_depth <- function(x, y, bin = FALSE, breaks, dissimilarity = FALSE) {
  if (bin) {
    x <- cut(x, breaks)
    y <- cut(y, breaks)
  } 
  xdens <- as.numeric(table(x)/length(x))
  ydens <- as.numeric(table(y)/length(y))
  if (dissimilarity) {
    return(1 - (1 - 0.5 * sum(abs(xdens - ydens))))
  } else {
    return(1 - 0.5 * sum(abs(xdens - ydens)))
  }
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
