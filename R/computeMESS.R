






#' Compute MESS rasters
#'
#' Compute set of MESS rasters for a maxnet model
#'
#' @param thisModel A maxnet model object
#' @param envData raster stack. The environmental predictor raster layers used to fit the maxnet model
#' @param swdData data.frame. SWD data used in model fitting and here used as tyhe reference data set to compute the MESS raster
#' @param outPath Character. Path to folder into which output raster layers will be written
#' @param MESSonly Logical. Compute only the MESS raster layer? Default is TRUE. If FALSE, then MESS component rasters for each variable plus the MESS raster layer are computed
#' @param ... Optional parameters passed to the raster function writeRaster
#' @details
#'
#' A MESS computation is performed using the function mess in the package dismo. This function is a wrapper which orchestrates a call to the mess function using information stored in the maxnet model object.
#'
#' To facilitate the processing of very large, high-resolution rasters, only those variables which are used in the final model are included in the computation of MESS output.
#'
#'
#' @return Stuff like what is produced by dismo::mess()
#' @export
#'
#' @examples
#' \dontrun{}
#'
computeMESS <- function(thisModel = NULL, envData = "", swdData = "", outPath = "", MESSonly = TRUE, ...)
{
  if (!("maxnet" %in% class(thisModel)))
    stop("'thisModel' is not a maxnet model object")

  # if (!dir.exists(envDataPath))
  #   stop("Cannot find folder given in 'envDataPath'")
  #
  # if (!file.exists(swdFilename))
  #   stop("Cannot find SWD file referenced in 'swdFilename'")

  # Which variables where used in the model fit?
  featureBetas <- thisModel$betas

  # Trim to non-zero coefficient values
  if (any(featureBetas == 0))
    featureBetas <- featureBetas[-which(featureBetas == 0)]

  tmpNames <- gsub("I(", "", names(featureBetas), fixed = TRUE)
  ii <- grep("^2)", tmpNames, fixed = TRUE)
  tmpNames[ii] <- gsub("^2)", "", tmpNames[ii], fixed = TRUE)
  tmpNames[ii] <- paste0(tmpNames[ii],":", tmpNames[ii])

  varNames <- sort(unique(unlist(strsplit(tmpNames, ":", fixed = TRUE))))

  # Load as a raster stack only those variables in the model
  #envDataFiles <- list.files(envDataPath, pattern = "tif", full.names = TRUE)

  #keepInd <- unlist(lapply(varNames, function(el) { grep(el, basename(envDataFiles)) }))
  keepInd <- which(varNames %in% names(envData))

  #envData <- raster::stack(envDataFiles[keepInd])

  envData <- envData[[keepInd]]

  # Load SWD file and do variable check
  #refData <- read.csv(swdFilename, stringsAsFactors = FALSE)[, -c(1:3)]
  refData <- swdData

  if (!all(varNames %in% colnames(refData)))
    stop("Variable names in model do not match variables names in the SWD file")
  else
  {
    refData <- refData[, varNames]
  }

  ans <- dismo::mess(x = envData, v = refData, full = MESSonly)
}
