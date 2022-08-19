#' Prepare SWD files of occurrence and background environmental data
#'
#' Load a set of environmental covariate to be used to spatially project a fitted maxnet model and build pairs of SWD-format files
#'
#' @details Load raster stack of environmental data layers to be used for fitting a maxnet model and, using coordinate values of occurrence records supplied in \emph{occData} and a background-constraining spatial polygon object, output samples-with-data (SWD) files defined in the original Java implementation of MaxEnt. Two files are produced: an occurrence SWD file, and a background points SWD file which are diferentiated by the file name assigned.
#'
#' NOTE: No checks are made that the projection of occurrence lat/longs or X/Y values matches the raster layers used for environmental data. You must ensure they match before calling this function.
#'
#'
#'
#' @param taxonName String (character object) giving the taxonomic name of the entity for which data is being supplied.
#' @param occData A numeric matrix or data.frame with at least two columns holding latitude and longitude values for occurrence records. Column names are matched against 'lat' (or 'Y") and 'long' (or 'X') in a effort to automagically identify these columns.
#' @param boundsPoly An sf or SpatialPolygons* object defining the region within which background points will be selected.
#' @param excludedVars A character vector of variable names to be excluded from output files. Default is NULL.
#' @param envDataPath String giving the file system path to the folder containing the geoTIFF rasters representing the environmental layers to be used to spatially project a fitted maxnet model
#' @param outputPath String giving the file system path to the folder into which occurrence and background SWD files will be written. See \emph{Details} for information about the file naming convention.
#' @param appendDate Logical. Should the current system date be appended to output file names? Default is FALSE.
#' @param trace Logical. Emits additional diagnostic messages when TRUE; default is FALSE
#' @return Invisibly returns NULL but has side-effect of outputting data files
#' @export
#'
#' @examples
#' \dontrun{}
prepData <- function(taxonName, occData, excludedVars = NULL, boundsPoly,
                     envDataPath, outputPath, appendDate = FALSE, trace = FALSE)
{
  taxon_Name <- gsub(" ", "_", taxonName, fixed = TRUE)

  cat("Preparing SWD data for", taxonName, "\n")

  if (dir.exists(envDataPath))
  {
    dataFileSet <- list.files(envDataPath, "*.tif", full.names = TRUE)
    if (length(dataFileSet) == 0)
      stop("Cannot find any geoTIFF files in the folder given in 'envDataPath'")
    else
    {
      envStack <- terra::rast(dataFileSet)
    }

    # Trim envStack using excludedVars
    if (!is.null(excludedVars))
    {
      if (trace)
      {
        cat("Excluded vars:\n")
        print(excludedVars)
        cat("----------------------------------------------------\n")
      }
      # Are all vars names in excludedVars present in the stack?
      if (all(excludedVars %in% names(envStack)))
      {
        dropInd <- which(names(envStack) %in% excludedVars)

        if (trace)
          {
          cat("\nTotal indices in stack =", terra::nlyr(envStack), "\n")
          cat("Drop var indices:\n")
          print(dropInd)
          cat("\nTotal indices in stack =", terra::nlyr(envStack), "\n")
          cat("----------------------------------------------------\n")
        }

         envStack <- terra::subset(envStack, dropInd, negate = TRUE)
      }
      else
      {
        stop("Some variable names in 'excludedVars' not found in environmental data layers")
      }
    }

    # Automagically try to identify longitude and latitude columns:
    Xind <- grep("LONG|^X$", toupper(colnames(occData)))[1]
    if (length(Xind) == 0) stop("Cannot identify a 'longitude' or 'X' column in occurrence data file")
    if (length(Xind) > 1)
    {
      warning("Cannot identify a unique 'longitude' column in occurrence data file; using first hit come-what-may...")
    }

    Yind <- grep("LAT|^Y$", toupper(colnames(occData)))[1]
    if (length(Yind) == 0) stop("Cannot identify a 'latitude' or 'Y' column in occurrence data file")
    if (length(Yind) > 1)
    {
      warning("Cannot identify a unique 'latitude' or 'Y' column in occurrence data file; using first hit come-what-may...")
    }

    if (trace)
    {
      cat("Xind =", Xind, " : Yind = ", Yind, "\n")
      cat("----------------------------------------------------\n")
    }

    # Set a flag to be used when determining column names of SWD files
    isLatLong <- ifelse(any(grepl("LAT|LONG", toupper(colnames(occData)))), TRUE, FALSE)

    if (appendDate)
    {
      outputSWDname <- paste0(outputPath, "/", taxon_Name, "_SWD_", Sys.Date(), ".csv")
      outputSWDBKGname <- paste0(outputPath, "/", taxon_Name, "_SWD_BKG_", Sys.Date(), ".csv")
    }
    else
    {
      outputSWDname <- paste0(outputPath, "/", taxon_Name, "_SWD.csv")
      outputSWDBKGname <- paste0(outputPath, "/", taxon_Name, "_SWD_BKG.csv")
    }

    if (trace)
    {
      cat("Output filenames:\n")
      cat("SWD filename:", outputSWDname, "\n")
      cat("BKG SWD filename:", outputSWDBKGname, "\n")
      cat("----------------------------------------------------\n")
    }

    cat("   Extracting background samples...")

    if (trace)
    {
      cat("\nFilename of first envStack layer:", envStack[[1]]@file@name, "\n")
      cat("----------------------------------------------------\n")
    }

    bkg.cells <- sampleBackground(as.matrix(occData[, c(Xind, Yind)]),
                                  nBkgSamples = 10000,
                                  baseRaster = envStack[[1]],
                                  boundsPolygon = boundsPoly,
                                  trace = trace)

    if (trace)
    {
      cat("Returned from sampleBackground()\n")
      cat("----------------------------------------------------\n")
    }

    bkg.xy <- terra::xyFromCell(envStack[[1]], bkg.cells)        #print(names(bkg.latlong))

    if (trace)
    {
      cat("Setting bkg.xy values completed\n")
      cat("----------------------------------------------------\n")
    }

    if (isLatLong)
      bkgData <- data.frame(species = rep("background", length(bkg.cells)),
                            longitude = bkg.xy[, 1],
                            latitude = bkg.xy[, 2],
                            terra::extract(envStack, bkg.cells),
                            stringsAsFactors = FALSE)
    else
      bkgData <- data.frame(species = rep("background", length(bkg.cells)),
                            x = bkg.xy[, 1],
                            y = bkg.xy[, 2],
                            terra::extract(envStack, bkg.cells),
                            stringsAsFactors = FALSE)

    if (trace)
    {
      cat("Dataframe 'bkgData' created\n")
      cat("----------------------------------------------------\n")
    }


    badRows <- unique(which(is.na(bkgData), arr.ind = TRUE)[ , 1]) #nBackground <- length(bkg.cells)
    if (length(badRows) > 0) bkgData <- bkgData[-badRows, ]
    cat("done.\n")

    cat("   Extracting data for occupied grid cells...")
    # Find grid cell row, col coordinates of occupied grid cells:
    occ.cells <- terra::cellFromXY(envStack[[1]], cbind(occData[, Xind], occData[, Yind]))
    #nPresence <- length(occ.cells)

    if (isLatLong)
      occData <- data.frame(species = rep(taxon_Name, length(occ.cells)),
                            longitude = occData[, Xind],
                            latitude = occData[, Yind],
                            terra::extract(envStack, occ.cells),
                            stringsAsFactors = FALSE)
    else
      occData <- data.frame(species = rep(taxon_Name, length(occ.cells)),
                            x = occData[, Xind],
                            y = occData[, Yind],
                            terra::extract(envStack, occ.cells),
                            stringsAsFactors = FALSE)

    if (trace)
    {
      cat("Dataframe 'occData' created\n")
      cat("----------------------------------------------------\n")
    }

    badRows <- unique(which(is.na(occData), arr.ind = TRUE)[ , 1]) #nBackground <- length(bkg.cells)
    if (length(badRows) > 0) occData <- occData[-badRows, ]

    cat("done.\nSaving files...")
    write.csv(bkgData, outputSWDBKGname, row.names = FALSE)
    write.csv(occData, outputSWDname, row.names = FALSE)
    cat("done.\n")
  }
  else
  {
    stop("Cannot find data path for environmental data")
  }

  invisible(NULL)
}
