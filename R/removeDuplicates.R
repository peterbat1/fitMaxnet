# Find and remove duplicate occurrence records.
#
#
# Peter D. Wilson 9 December 2013
#
# Some bug fixes: 13 March 2014
# Remove NA lat-longs before testing for duplicates; repect "quiet" paramater: 5 November 2014
# Converted to use raster package functions: 3 January 2016

#' Remove duplicate occurrence records
#'
#' Remove duplicate occurrence records by exact coordinate matching or by ensuring only one record per grid cell of a raster. Coordinates are assumed be geographic coordinates (i.e. latitude and longitude) on the WGS84 datum
#'
#' @param data A data.frame or matrix with at least a latitude and longitude column.
#' @param baseGrid A RasterLayer or SpatRaster object representing the grid for enforcing one record per cell when \emph{byGrid} = TRUE.
#' @param byGrid Logical. If TRUE then use the rasterLayer passed in \emph{baseGrid} to select opne record per cell. If FALSE (default) then ignore anything passed in \emph{baseGrid} and remove exact duplicates based on there coordinates.
#' @param quiet Logical. If TRUE then emit some progress message. If FALSE (default) procede quietly.
#'
#' @return A version of \emph{data} with duplicates removed
#' @export
#'
#' @examples
#' \dontrun{}
removeDuplicates <- function(data = NULL, baseGrid = NULL, byGrid = FALSE, quiet = FALSE)
{
  if (is.null(data) || (nrow(data) < 2) || (ncol(data) < 2))
  {
    stop("removeDuplicates: Degenerate data matrix - no processing performed.")
  }

  # Try to identify longitude (X) and latitude (Y) columns
  ind <- grep("LONG|X", toupper(colnames(data)))
  if (length(ind) >= 1)
    X_ind <- ind[1]
  else
    stop("removeDuplicates: Cannot identify the 'longitude' or 'X' column in 'data'")

  ind <- grep("LAT|Y", toupper(colnames(data)))
  if (length(ind) >= 1)
    Y_ind <- ind[1]
  else
    stop("removeDuplicates: Cannot identify the 'latitude' or 'Y' column in 'data'")

  # Remove NA lat-longs
  na.latlong <- union(which(is.na(data[, Y_ind])), which(is.na(data[, X_ind])))
  if (length(na.latlong) > 0)
  {
    data <- data[-na.latlong,]
  }

  if (nrow(data) == 0)
  {
    stop("removeDuplicates: Degenerate data matrix: no valid lat-long values after removing NAs - no processing performed.")
  }

  if (byGrid)
  {
    if (!inherits(baseGrid, c("RasterLayer", "SpatRaster")))
      stop("removeDuplicates: byGrid = TRUE so parameter baseGrid must of class RasterLayer or class SpatRaster.")

    # Convert to class terra::SpatRaster
    if (inherits(baseGrid, "RasterLayer")) baseRaster <- terra::rast(baseGrid)

    cellInd <- terra::cellFromXY(baseGrid, cbind(data[, X_ind], data[, Y_ind]))
    dataInd <- 1:nrow(data)

    # Test to see if all lat and all long are out of bounds - it does happen sometimes!
    if (all(is.na(cellInd)))
    {
      if (!quiet) message("removeDuplicates: Degenerate data matrix: all lat or all long values are NAs - no processing performed.")
      return(NULL)
    }

    # Remove out-of-bounds lat-longs
    na.latlong <- which(is.na(cellInd))
    if (length(na.latlong) > 0)
    {
      data <- data[-na.latlong,]
      cellInd <- cellInd[-na.latlong]
    }

    duplicates <- which(duplicated(cellInd))
  }
  else
  {
    sortIndex <- order(data[, Y_ind], data[, X_ind])
    data <- data[sortIndex, ]
    xDuplicates <- which(duplicated(data[, X_ind]))
    yDuplicates <- which(duplicated(data[, Y_ind]))

    duplicates <- intersect(xDuplicates, yDuplicates)
  }

  #print(duplicates)
  # Use the index values stored in duplicates (produced by one of two methods) to trim data:
  if (length(duplicates) > 0)
  {
    if (!quiet) message(paste0("removeDuplicates: ",length(duplicates)," removed."))
    return(data[-duplicates,])
  }
  else
  {
    if (!quiet) message("removeDuplicates: No duplicates found.")
    return(data)
  }
}

