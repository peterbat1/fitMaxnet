#' \pkg{fitMaxnet}
#'
#' A package of functions to provide methods for preparing data files and
#' fitting MaxEnt models using the function maxnet in the R-package maxnet.
#'
#' @name fitMaxnet
#' @import maxnet
#' @import ggplot2
#' @import dplyr
#' @importFrom stats prcomp binomial cor ecdf glm na.omit predict
#' @importFrom grDevices chull png dev.off
#' @importFrom graphics points
#' @importFrom terra cellFromXY xyFromCell subset extract nlyr values writeRaster
#' @importFrom geodist geodist
#' @importFrom utils read.csv write.csv
#' @importFrom sf st_geometry st_transform st_as_sf as_Spatial st_intersection st_crs st_buffer st_union st_is
#' @importFrom ggcorrplot ggcorrplot
#' @importFrom energy eqdist.etest
#' @importFrom splines ns
#' @importFrom ggpubr ggarrange ggexport ggbarplot
#' @importFrom dismo mess
#' @importFrom PRROC pr.curve roc.curve
#' @importFrom Rfast rowsums
#' @importFrom ggsci pal_npg
#' @importFrom stringr str_pad
NULL
