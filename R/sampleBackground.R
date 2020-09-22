



###############
#' Selected background cells from disks of given radius
#'
#' Perform background cell selection from a union of disks of specified radius (in km) around each occurrence cell. Based on VanDerWal et al. (2009).
#'
#' @param baseMap A RasterLayer object representing the spatial grid to be sampled
#' @param seedCells Cell indices on the \emph{baseMap} grid representing occupied cells to act as seed cells
#' @param diskRadius Radius in kilometres of the disk to be used fro background cell sleection
#'
#' @return An integer array of selected grid cells on the \emph{baseMap} raster
#' @export
#'
#' @examples
unionOfDisks <- function(baseMap, seedCells, diskRadius = 200)
{
  ##cat("   Hello from UnionOfDisks\n")

  if (class(baseMap)=="RasterLayer")
  {
    #print(extent(baseMap))
    ndCells <- which(is.na(values(baseMap)))
    dataCells <- which(!is.na(values(baseMap)))
    #print(length(data.cells))
    ######seed.cells <- cellFromXY(baseMap,cbind(x=latlon.pts[,2],y=latlon.pts[,1]))
    #print(length(seed.cells))

    values(baseMap)[dataCells] <- 0
    values(baseMap)[seedCells] <- 1

    tmp <- doDilation2(length(seedCells), as.matrix(baseMap), nrow(baseMap),
                       ncol(baseMap), extent(baseMap)@ymax, extent(baseMap)@xmin,
                       xres(baseMap), rowFromCell(baseMap, seedCells),
                       colFromCell(baseMap, seedCells), diskRadius)

    values(baseMap) <- as.vector(t(tmp))
    ####plot(baseMap)

    #print(sum(values(baseMap)==1,na.rm=TRUE))

    values(baseMap)[ndCells] <- NA

    # Make sure occurrence cells are not included in the background sample
    values(baseMap)[seedCells] <- NA

    #print(sum(values(baseMap)==1,na.rm=TRUE))

    return(baseMap)
  }
  else
  {
    stop("unionOfDisks: Parameter basemap must be of class RasterLayer")
  }
}


###################
#' Perform background cell selection
#'
#' Given a set of occupied grid cells on a specified raster, perform background cell selection using one of two methods: unconstrained, or 'union of disks' based on VanDerWal et al. (2009).
#'
#' @param occupiedCells Gridcell indices of occupied cells on the raster defined in \emph{baseMap}.
#' @param nSamples Number fo samples to be drawn to represent the background.
#' @param baseMap A RasterLayer object defining the grid to be sampled.
#' @param method Character string selecting which of two methods to use to draw the background sample.
#' @param diskRadius Radius of disk in kilometers around each occurrence cell for "Union of Disks' method. Ignored otherwise.
#'
#' @return An integer array giving the indices of selected grid cells on the raster defined in \emph{baseMap}
#' @export
#'
#' @examples
sampleBackgroundOld <- function(occupiedCells, nSamples = 1000, baseMap, method = "unconstrained", diskRadius = 200)
{  # Samples background points using a sampling strategy determined by parameter method.
  # Sampling is currently done without replacement.
  #
  # As of 12 May 2014 changes: Parameter occupiedCells must be a numeric matrix of integers with 2
  # columns labelled with lat and long
  # And, method must be one of "unconstrained", "unionOfDisks",...???
  #
  # Returns a matrix of row, col pairs for grid cells included in the background sample

  ##cat("   Hello from sampleBackground\n")

  if (class(baseMap) == "RasterLayer")
  {
    if (is.matrix(occupiedCells))
    {
      if (!is.null(hasLatLongCols(occupiedCells)))
      {
        if (method %in% validBkgMethods)
        {
          # Simplify the task ahead by removing duplicates in cells (added 2 Jan 2016):
          cellInd <- cellFromXY(baseMap,cbind(occupiedCells[,2],occupiedCells[,1]))
          goodRows <- 1:nrow(occupiedCells)

          naCells <- which(is.na(cellInd))

          if (length(naCells) > 0)
          {
            cellInd <- cellInd[-naCells]
            goodRows <- goodRows[-naCells]
          }

          duplCells <- which(duplicated(cellInd))
          if (length(duplCells) > 0)
          {
            cellInd <- cellInd[-duplCells]
            goodRows <- goodRows[-duplCells]
          }

          if (method=="unconstrained")
          {
            # "Switch off" occupied cells by setting them to NA - this will prevent them being sampled
            # for the background:
            values(baseMap)[cellInd] <- NA

            # Make a matrix of row, col values for non-NA cells
            available <- which(!is.na(values(baseMap)))
          }
          else
          {
            baseMap <- unionOfDisks(baseMap, cellInd, diskRadius)
            available <- which(values(baseMap) == 1)
          }

          if (length(available) == 0)
          {
            warning("sampleBackground: There are no available grid cells to form a background sample!")
            return(NULL)
          }
          else
          {
            if (length(available) <= nSamples)
            {
              warning("sampleBackground: Number of available grid cells for backgound is less than or equal to requested number of samples. All available cells returned.")
              return(available)
            }
            else
            {
              selectedCells <- sample(1:length(available),nSamples,replace=F)
              return(available[selectedCells])
            }
          }
        }
        else
        {
          warning("sampleBackground: Value passed in parameter method is not a recognised sampling strategy")
          return(NULL)
        }
      }
      else
      {
        warning("sampleBackground: Parameter occupiedCells does not have lat/long columns")
        return(NULL)
      }
    }
    else
    {
      warning("sampleBackground: Parameter occupiedCells must be a matrix")
      return(NULL)
    }
  }
  else
  {
    warning("sampleBackground: Parameter baseMap must of class RasterLayer")
    return(NULL)
  }
}


################################################################################


#' Select background cells
#'
#' Select a set of background cells from a base raster layer
#'
#' @param occData a data frame or matrix with at least two numeric columns which can be associated with longitude and latitude
#' @param baseRaster A raster layer object representing the geometry of the environmental data layers to be used for modelling
#' @param boundsPolygon An sf or SpatialPolygons* object defning the boundary within which background points must be selected
#' @param nBkgSamples Numeric. The number of samples to constitute the background set with the default = maxBkgSamples
#' @param maxBkgSamples Numeric. Maximum number of background samples to be generated (default = 10,000)
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
sampleBackground <- function(occData, baseRaster, boundsPolygon, nBkgSamples = maxBkgSamples, maxBkgSamples = 10000)
{
  if (class(baseRaster) != "RasterLayer") stop("'baseRaster' must be a RasterLayer object")

  if (any(grepl("sf|sfc", class(boundsPolygon))))
  {
    boundsPolygon <- sf::as_Spatial(boundsPolygon)
  }

  if (!grep("SpatialPolygons", class(boundsPolygon)))
    stop("'boundsPolygon' must be a sf, SpatialPolygons, or SpatialPolygonsDataFrame object")

  # Try to identify longitude and latitude, or X & Y, columns in occData
  ind <- grep("LONG|X", toupper(colnames(occData)))
  if (length(ind) >= 1)
    X_ind <- ind #colnames(occData)[ind[1]] <- "longitude"
  else
    stop("Cannot identify the 'longitude' or 'X' column in 'occData'")

  ind <- grep("LAT|Y", toupper(colnames(occData)))
  if (length(ind) >= 1)
    Y_ind <- ind #colnames(occData)[ind[1]] <- "latitude"
  else
    stop("Cannot identify the 'latitude' or 'Y' column in 'occData'")

  occCells <- raster::cellFromXY(baseRaster, occData[, c(X_ind, Y_ind)])

  raster::values(baseRaster)[which(!is.na(raster::values(baseRaster)))] <- 1

  # Remove occupied cells from the set available for selection
  raster::values(baseRaster)[occCells] <- NA

  activeArea <- raster::mask(baseRaster, boundsPolygon)
  availableCells <- which(raster::values(activeArea) == 1)

  if (length(availableCells) < nBkgSamples)
  {
    warning("Number of available raster cells for background selection < nBkgSamples; all available cells returned")
    return(availableCells)
  }

  #availableCells <- which(values(baseRaster) == 1)

  selectedCells <- NULL
  for (i in 1:nBkgSamples)
  {
    newInd <- sample(1:length(availableCells), 1)
    selectedCells <- c(selectedCells, availableCells[newInd])
    availableCells <- availableCells[-newInd]
  }

  return(selectedCells)
}
