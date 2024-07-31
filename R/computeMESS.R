#' Compute MESS rasters
#'
#' Compute set of MESS rasters for a maxnet model
#'
#' @param thisModel A maxnet model object
#' @param envPath String (character). Full specified path to the environmental predictor raster layers used to fit the maxnet model
#' @param occSWD Data.frame. Reference environmental data at the occurrence locations used to fit the model in MaxEnt SWD format
#' @param bkgSWD Data.frame. Reference environmental data for locations used as background points used to fit the model in MaxEnt SWD format
#' @param outPath Character. Path to folder into which output raster layers will be written
#' @param varImpThreshold Numeric. Value used to filter variables by their importance scores; percentage importance therefore ranges from 0 to 100.
#' @param MESSonly Logical. Compute only the MESS raster layer? Default is TRUE. If FALSE, then MESS component rasters for each variable plus the MESS raster layer are computed
#' @details
#'
#' Multivariate Environmental Similarity Surfaces (MESS) were introduced by Elith, Kearney & Phillips (2010; The art of modelling range-shifting species. Methods in Ecology and Evolution 1: 330-342) to assist in identifying the degree of similarity between the ranges of environmental predictors presented to model training (the 'reference' set)  and the values of those predictors across other areas. It's purpose is to identify areas of extrapolation. In the present implementation, the reference set is extracted from the tables of variable maximum and minimumthose in areas or times to which the model is projected.
#'
#' A MESS computation is performed using the supplied maxnet model and set of environmental rasters.
#'
#' To facilitate the processing of very large, high-resolution rasters, only those variables which pass the variable importance threshold are used in the final model are included in the computation of MESS output.
#'
#' Please note that this function is tailored to work with data structures and workflows associated with the primary model fitting function in this package, namely \link{fit_maxnet}. It can be used for MESS computations in other settings, but users will have to transform data into the required structure for each parameter. An alternate and more generic MEWS function in found in the R-package \emph{modEVA}.
#'
#' If output rasters are to be written, note that the following options are set as defaults for terra::writeRaster(): 'overwrite = TRUE', and 'gdal = "COMPRESS=DEFLATE"'
#'
#' @return A named list:
#' \describe{
#' \item{MESS}{Matrix of OVERALL MESS scores}
#' \item{varMESS}{Matrix of MESS scores for each included variable}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{}
#'
computeMESS <- function(thisModel = NULL,
                        envPath = "",
                        occSWD,
                        bkgSWD,
                        varImpThreshold = 0,
                        outPath = "",
                        MESSonly = TRUE) #,
#                        ...)
{
  if (!inherits(thisModel, "maxnet"))
    stop("'thisModel' is not a maxnet model object")

  if (envPath == "")
    stop("'envPath' is empty. Please supply a path to a folder to storing environmental rasters")

  if (!dir.exists(envPath))
    stop("Path given in 'envPath' cannot be found")

  if (!inherits(occSWD, "data.frame"))
    stop("occSWD must be a data.frame in MaxEnt SWD format")

  if (!inherits(bkgSWD, "data.frame"))
    stop("bkgSWD must be a data.frame in MaxEnt SWD format")

  if (!is.numeric(varImpThreshold))
    stop("'varImpThreshold' must be a numeric value between 0 and 100")

  if ((varImpThreshold < 0) | (varImpThreshold > 100))
    stop("'varImpThreshold' must be a numeric value between 0 and 100")

  if (outPath == "")
    stop("'outPath' is empty. Please supply a path to a folder to store output rasters")

  if (!dir.exists(outPath))
    stop("Path given in 'outPath' cannot be found")

  envFileSet <- list.files(envPath, "*.tif", full.names = TRUE)
  envStack <- terra::rast(envFileSet)
  envVarNames <- terra::names(envStack)

  varImp <- fitMaxnet::varImportance(thisModel, occSWD, bkgSWD, "cloglog")

  if (!all(envVarNames %in% colnames(occSWD)[4:ncol(occSWD)]))
    stop("Environmental variable names do not match variables in 'occSWD'")

  if (!all(envVarNames %in% colnames(bkgSWD)[4:ncol(bkgSWD)]))
    stop("Environmental variable names do not match variables in 'bkgSWD'")

  if (!all(envVarNames %in% names(varImp)))
    stop("Environmental variable names do not match variables usewd in model")

  goodVars <- names(varImp)[varImp > varImpThreshold]

  varMin <- thisModel$varmin
  varMax <- thisModel$varmax
  varRange <- varMax - varMin

  # refMat <- as.matrix(rbind(occSWD[, 4:ncol(occSWD)],
  #                           bkgSWD[, 4:ncol(bkgSWD)]))
  # refMat <- refMat[, goodVars]
  #
  # varMin <- Rfast::colMins(refMat[, goodVars], value = TRUE)
  # names(varMin) <- goodVars
  #
  # varMax <- Rfast::colMaxs(refMat[, goodVars], value = TRUE)
  # names(varMax) <- goodVars

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
  nRefPts <- thisModel$nobs

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
    # Overall MESS raster is ALWAYS produced...
    rasTemplate[goodRows] <- MESS
    if (!dir.exists(outPath)) dir.create(outPath, recursive = TRUE)
    outFilename <- file.path(outPath, "MESS.tif")
    terra::writeRaster(rasTemplate, outFilename, overwrite = TRUE, gdal = "COMPRESS=DEFLATE", ...)

    # Optionally, produce MESS raster for each included variable
    if (!MESSonly)
    {
      for (thisVar in goodVars)
      {
        rasTemplate[goodRows] <- varMESS[, thisVar]
        outFilename <- file.path(outPath, paste0("MESS_", thisVar, ".tif"))
        terra::writeRaster(rasTemplate, outFilename, overwrite = TRUE, gdal = "COMPRESS=DEFLATE")
      }
    }
  }

  return(list(MESS = MESS, varMESS = varMESS))
}
