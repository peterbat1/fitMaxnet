#' Compute MESS rasters
#'
#' Compute set of MESS rasters for a maxnet model
#'
#' @param thisModel A maxnet model object
#' @param envPath Full specified path to the environmental predictor raster layers used to fit the maxnet model
#' @param occSWD numeric matrix. Reference data; typically created by concatenating occurrence SWD and background SWD files used to fit the model
#' @param bkgSWD numeric matrix
#' @param outPath Character. Path to folder into which output raster layers will be written
#' @param varImpThreshold Numeric. Value used to filter variables by their importance scores; percentage importance therefore ranges from 0 to 100.
#' @param MESSonly Logical. Compute only the MESS raster layer? Default is TRUE. If FALSE, then MESS component rasters for each variable plus the MESS raster layer are computed
#' @param ... Optional parameters passed to the raster function writeRaster
#' @details
#'
#' A MESS computation is performed blah blah.
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
computeMESS <- function(thisModel = NULL, envPath = NULL, occSWD, bkgSWD, varImpThreshold = 0, outPath = "", MESSonly = TRUE, ...)
{
  if (!("maxnet" %in% class(thisModel)))
    stop("'thisModel' is not a maxnet model object")

  envFileSet <- list.files(envPath, "*.tif", full.names = TRUE)
  envStack <- terra::rast(envFileSet)

  varImp <- fitMaxnet::varImportance(thisModel, occSWD, bkgSWD, "cloglog")

  goodVars <- names(varImp)[varImp > varImpThresold]

  refMat <- as.matrix(rbind(occSWD[, 4:ncol(occSWD)],
                            bkgSWD[, 4:ncol(bkgSWD)]))
  refMat <- refMat[, goodVars]

  varMin <- Rfast::colMins(refMat[, goodVars], value = TRUE)
  names(varMin) <- goodVars

  varMax <- Rfast::colMaxs(refMat[, goodVars], value = TRUE)
  names(varMax) <- goodVars

  varRange <- varMax - varMin
  names(varRange) <- goodVars

  rasTemplate <- envStack[[1]]

  envMat <- terra::as.matrix(envStack)
  #print(colnames(envMat))
  envMat <- envMat[, goodVars]

  #remove(envStack)
  rowStatus <- Rfast::rowAll(!is.na(envMat))

  goodRows <- which(rowStatus)

  N <- length(goodRows)
  nRefVars <- length(goodVars)
  nRefPts <- nrow(refMat)

  varMESS <- matrix(0, N, nRefVars)
  colnames(varMESS) <- goodVars

  f <- matrix(0, N, nRefVars)
  colnames(f) <- goodVars

  p <- envMat[goodRows, goodVars]
  #colnames(p) <- goodVars

  # Borrow a clever idea from ecospat MESS implementation:
  varECDF <- vector("list", nRefVars)
  names(varECDF) <- goodVars

  for (thisVar in goodVars)
  {
    varECDF[[thisVar]] <- ecdf(refMat[, thisVar])
  }

  for (thisVar in goodVars)
  {
    f[, thisVar] <- 100*varECDF[[thisVar]](p[, thisVar])
  }

  for (thisVar in goodVars)
  {
    firstInd <- na.omit(which(f[, thisVar] == 0))
    if (length(firstInd) > 0)
      varMESS[firstInd, thisVar] <- 100*(p[firstInd, thisVar] - varMin[thisVar])/varRange[thisVar]

    secondInd <- which((f[, thisVar] > 0) & (f[, thisVar] <= 50))
    if (length(secondInd) > 0)
      varMESS[secondInd, thisVar] <- 2 * f[secondInd, thisVar]

    thirdInd <- which((f[, thisVar] > 50) & (f[, thisVar] < 100))
    if (length(thirdInd) > 0)
      varMESS[thirdInd, thisVar] <- 2 * (100 - f[thirdInd, thisVar])

    fourthInd <- which(f[, thisVar] == 100)
    if (length(fourthInd) > 0)
      varMESS[fourthInd, thisVar] <- 100*(varMax[thisVar] - p[fourthInd, thisVar])/varRange[thisVar]
  }

  MESS <- Rfast::rowMins(varMESS, value = TRUE)

  if (outPath != "")
  {
    rasTemplate[goodRows] <- MESS
    if (!dir.exists(outPath)) dir.create(outPath, recursive = TRUE)
    outFilename <- file.path(outPath, "MESS.tif")
    terra::writeRaster(rasTemplate, overwrite = TRUE, options = "COMPRESS=DEFLATE")
  }

  return(MESS)
}
