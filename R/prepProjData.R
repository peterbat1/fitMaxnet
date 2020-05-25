#' Prepare environmental projection data
#'
#' Load a set of environmental covariate to be used to spatially project a fitted maxnet model
#'
#' @details {
#' Loads a raster stack of environmental data layers to be used for projecting a maxnet model, and creates a matrix of data values and an array of indices pointing to rows with no missing values (\emph{goodCellInd}).
#'
#' The often \emph{very} large data objects are placed in the global environment for direct access by functions as an efficiency measure.
#'
#'  A future version is planned wherein very large output will be stored in Rd files for fast reloading.
#'  }
#'
#' @param dataPath String giving the file system path to the folder containing the geoTIFF rasters representing the environmental layers to be used to spatially project a fitted maxnet model
#' @param quiet Logical. Default of TRUE supresses progress messages.
#' @return Invisibly returns NULL but has the side effect of creating several objects in the R global environment.
#' @export
#'
#' @examples
prepProjData <- function(dataPath, quiet = TRUE)
{
  if (!quiet) cat("Preparing projection data: ")
  projStack <<- raster::stack(list.files(dataPath, "*.tif", full.names = TRUE))
  projData <<- raster::values(projStack)
  goodCellInd <<- which(!is.na(rowSums(projData)))
  if (!quiet) cat(length(projStack), " layers loaded\n")
  invisible(NULL)
}
