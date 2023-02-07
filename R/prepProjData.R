#' Prepare environmental projection data
#'
#' Load a set of environmental covariates to be used to spatially project a fitted maxnet model
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
#' @param quiet Logical. Default of TRUE suppresses progress messages.
#' @return Invisibly returns NULL but has the side effect of creating several objects in the R global environment.
#' @export
#'
#' @examples
#' \dontrun{}
prepProjData <- function(dataPath, quiet = TRUE)
{
  if (!dir.exists(dataPath)) stop("dataPath not found")

  if (!quiet) cat("Preparing projection data: ")
  theFilePaths <- list.files(dataPath, "*.tif", full.names = TRUE)
  projStack <- terra::rast(theFilePaths)
  names(projStack) <- gsub(".tif$", "", basename(theFilePaths))

  # "Export" data as a large matrix where rows are raster cells and columns are environmental variables
  projData <<- terra::values(projStack)
  #goodCellInd <<- which(!is.na(rowSums(projData)))

  # "Export" a template raster layer
  rasTemplate <<- projStack[[1]]
  if (!quiet) cat(length(projStack), " layers loaded\n")

  invisible(NULL)
}
