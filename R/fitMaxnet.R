#' \pkg{fitMaxnet}
#'
#' A package of functions to provide methods for preparing data files and
#' fitting MaxEnt models using the function maxnet in the R-package maxnet.
#'
#' @name fitMaxnet
#' @docType package
#' @import Rcpp
#' @import maxnet
#' @import ggplot2
#' @import dplyr
#' @importFrom stats prcomp
#' @importFrom grDevices chull png
#' @importFrom raster xyFromCell cellFromXY extract plot stack values mask dropLayer calc writeRaster
#' @importFrom utils read.csv write.csv
#' @importFrom sp spTransform coordinates CRS spDistsN1
#' @importFrom sf st_geometry st_transform st_as_sf as_Spatial st_intersection
#' @importFrom rgeos gBuffer
NULL
