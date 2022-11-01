



#' Select background cells
#'
#' Select a set of background cells from a base raster layer
#'
#' @param occData a data frame or matrix with at least two numeric columns which can be associated with longitude (or X) and latitude (or Y)
#' @param baseRaster A raster layer object representing the geometry of the environmental data layers to be used for modelling
#' @param boundsPolygon An sf or SpatialPolygons* object defining the boundary within which background points must be selected
#' @param nBkgSamples Numeric. The number of samples to constitute the background set with the default = maxBkgSamples
#' @param maxBkgSamples Numeric. Maximum number of background samples to be generated (default = 10,000)
#' @param trace Logical. Emits additional diagnostic messages when TRUE; default is FALSE
#'
#' @details
#' {The number of background samples may need to be adjusted to suit particular situations.
#'
#' The default value for maxBkgSamples matches the value used in the original MaxEnt Java implementation (see Phillips, S. J., and M. Dudík. 2008. Modeling of species distributions with Maxent: new extensions and a comprehensive evaluation. Ecography 31:161–175. for the reason it was selected).
#'
#' For many taxa, sampling a constrained area means that less than 10,000 cells will be available. When this occurs, the user is warned and the number of background samples returned is the maximum physically available within the specified buffer.
#'
#' }
#'
#' @return A integer array of indices to selected background cells on the base raster
#' @export
#'
#' @examples
#' \dontrun{}
sampleBackground <- function(occData, baseRaster, boundsPolygon, nBkgSamples = maxBkgSamples, maxBkgSamples = 10000, trace = FALSE)
{
  if (!inherits(baseRaster, c("RasterLayer", "SpatRaster"))) stop("'baseRaster' must be a sp::RasterLayer or terra::SpatRaster object")

  # Convert to class terra::SpatRaster
  if (inherits(baseRaster, "RasterLayer")) baseRaster <- terra::rast(baseRaster)

  if (any(grepl("sf|sfc", class(boundsPolygon))))
  {
    boundsPolygon <- sf::as_Spatial(boundsPolygon)
  }

  if (!grep("SpatialPolygons", class(boundsPolygon)))
    stop("'boundsPolygon' must be a sf, SpatialPolygons, or SpatialPolygonsDataFrame object")

  # Try to identify longitude and latitude, or X & Y, columns in occData
  ind <- grep("LONG|X$", toupper(colnames(occData)))
  if (length(ind) >= 1)
    X_ind <- ind #colnames(occData)[ind[1]] <- "longitude"
  else
    stop("Cannot identify the 'longitude' or 'X' column in 'occData'")

  ind <- grep("LAT|Y$", toupper(colnames(occData)))
  if (length(ind) >= 1)
    Y_ind <- ind #colnames(occData)[ind[1]] <- "latitude"
  else
    stop("Cannot identify the 'latitude' or 'Y' column in 'occData'")

  occCells <- terra::cellFromXY(baseRaster, as.matrix(occData[, c(X_ind, Y_ind)]))

  if (trace)
  {
    cat(length(occCells), " occCells found\n")
    cat("----------------------------------------------------\n")
  }

  terra::values(baseRaster)[which(!is.na(terra::values(baseRaster)))] <- 1

  # Remove occupied cells from the set available for selection
  terra::values(baseRaster)[occCells] <- NA

  activeArea <- terra::mask(baseRaster, terra::vect(boundsPolygon))
  availableCells <- which(terra::values(activeArea) == 1)

  if (trace)
  {
    cat("Recoding of baseRaster completed\n")
    cat("----------------------------------------------------\n")
  }

  if (length(availableCells) < nBkgSamples)
  {
    warning("Number of available raster cells for background selection < nBkgSamples; all available cells returned")
    return(availableCells)
  }

  #availableCells <- which(values(baseRaster) == 1)

  # selectedCells <- NULL
  # for (i in 1:nBkgSamples)
  # {
  #   newInd <- sample(1:length(availableCells), 1)
  #   selectedCells <- c(selectedCells, availableCells[newInd])
  #   availableCells <- availableCells[-newInd]
  # }

  selectedCells <- sample(availableCells, nBkgSamples)

  return(selectedCells)
}
