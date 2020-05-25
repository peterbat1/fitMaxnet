#' Prepare SWD files of occurrence and background environmental data
#'
#' Load a set of environmental covariate to be used to spatially project a fitted maxnet model and build pairs of SWD-format files
#'
#' @details Load raster stack of environmental data layers to be used for fitting a maxnet model and, using lat/long values of occurrence records supplied in \emph{occData} and a background-constraining spatial polygon object, output samples-with-data (SWD) files defined in the original Java implementation of MaxEnt. Two files are produced: an occurence SWD file, and a background points SWD file which are diferentiated by the file name assigned.
#'
#'
#'
#' @param taxonName String (character object) giving the taxonomic name of the entity for which data is being supplied.
#' @param occData A numeric matrix or data.frame with at least two columns holding latitude and longitude values for occurrence records. Column names are matched against 'lat' and 'long' in a effort to automagically identify these columns.
#' @param boundsPoly A SpatialPoylgons* object defining the region within which background points will be selected.
#' @param envDataPath String giving the file system path to the folder containing the geoTIFF rasters representing the environmental layers to be used to spatially project a fitted maxnet model
#' @param outputPath String giving the file system path to the folder into which occurrence and background SWD files will be written. See \emph{Details} for information about the file naming convention.
#' @param appendDate Logical. Should the current system date be appended to output file names? Default is FALSE.
#' @return Invisibly returns NULL but has side-effect of outputing data files
#' @export
#'
#' @examples
prepData <- function(taxonName, occData, excludedVars = NULL, boundsPoly, envDataPath, outputPath, appendDate = FALSE)
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
      envStack <- raster::stack(dataFileSet)
    }

    # Trim envStack using excludedVars
    if (!is.null(excludedVars))
    {
      # Are all vars names in excludedVars present in the stack?
      if (all(excludedVars %in% names(envStack)))
      {
        dropInd <- which(names(envStack) %in% excludedVars)
        envStack <- raster::dropLayer(envStack, dropInd)
      }
      else
      {
        stop("Some variable names in 'excludedVars' not found in environmental data layers")
      }
    }

    # Identify lat/long columns in occData
    latColInd <- grep("LAT", toupper(colnames(occData)))
    longColInd <- grep("LONG", toupper(colnames(occData)))

    if ((length(latColInd) == 0) || (length(longColInd) == 0)) stop("Cannot find a lat or long column in occData")

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

    cat("   Extracting background samples...")
    bkg.cells <- sampleBackground(as.matrix(occData[, c(latColInd, longColInd)]) , nBkgSamples = 10000, baseRaster = envStack[[1]], boundsPolygon = boundsPoly)
    bkg.latlong <- xyFromCell(envStack[[1]], bkg.cells)        #print(names(bkg.latlong))

    bkgData <- data.frame(species = rep("background", length(bkg.cells)),
                          longitude = bkg.latlong[, 1],
                          latitude = bkg.latlong[, 2],
                          extract(envStack, bkg.cells),
                          stringsAsFactors = FALSE)
    badRows <- unique(which(is.na(bkgData), arr.ind = TRUE)[ , 1]) #nBackground <- length(bkg.cells)
    if (length(badRows) > 0) bkgData <- bkgData[-badRows, ]
    cat("done.\n")

    cat("   Extracting data for occupied grid cells...")
    # Find grid cell row, col coordinates of occupied grid cells:
    occ.cells <- cellFromXY(envStack[[1]], cbind(occData[, longColInd], occData[, latColInd]))
    #nPresence <- length(occ.cells)

    occData <- data.frame(species = rep(taxon_Name, length(occ.cells)),
                          longitude = occData[, longColInd],
                          latitude = occData[, latColInd],
                          extract(envStack, occ.cells),
                          stringsAsFactors = FALSE)
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
