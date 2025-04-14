
#' Mask a projection raster removing extrapolation areas
#'
#' Compute and apply a mask of areas within which strict extrapolation has been detected
#'
#' @param maxnetModel Character. Full path to a maxnet model fitted by \link{fit_maxnet} and saved as an .Rd file.
#' @param projRasFile Character. Path to a raster of the projected maxnet model to which the computed extrapolation mask will be applied.
#' @param maskOutpath Character. Full path to the folder into which the mask raster and projection raster with mask applied will be written.
#' @param fileLabel Character. Label to be used to distinguish the filename of the saved mask raster.
#' @param makePlots Logical. Make basic plots of raster objects for interpretation, review and quality control? Default is TRUE.
#' @param saveMask Logical. Should the mask raster (as distinct from the the masked version of projRas) be saved? Default is TRUE.
#' @param silent Logical. If TRUE (default), no progress messages are written to the console.
#'
#' @details {
#' The nominated MaxEnt model object produced by \link{fit_maxnet} is interrogated to find the list of variables used in the model. The function then proceeds to:
#' \enumerate{
#' \item Load a raster stack of these variables
#' \item Code each non-NA cell in a raster layer 1 if it is within the range of the variable recorded in the model object, 0 otherwise
#' \item Compute an output raster with only cells set to 1 when all values in the stack are 1, 0 otherwise and NA cells retained
#' \item Save the output raster
#'}
#' The file name for the output raster defaults to 'extrapolationMask.tif'. If \emph{fileLabel} is not NULL, then the output file name is 'extrapolationMask_' + '\emph{fileLabel}' + '.tif'.
#' }
#'
#' @return Nothing but has side-effect of writing up to two raster files: mask raster (if saveMask == TRUE) and projection raster with mask applied plus, if \emph{makePlots} == TRUE, a set of PNG graphics files.
#' @export
#'
#' @examples
#' {
#' \dontrun{
#' ## Very basic use case:
#' maskExtrapolation("furryCreature.Rd", "currentClimate/envVarFolder",
#'                   "furryCreature_run1.tif", makePlots = FALSE)
#' }
#' }

maskExtrapolation <- function(maxnetModel,
                              projRasFile,
                              maskOutpath = dirname(projRasFile),
                              fileLabel = NULL,
                              makePlots = TRUE,
                              saveMask = TRUE,
                              silent = TRUE)
{
  if (!exists("projData")) stop("Global object 'projData' not found. Please run 'prepProjData'")

  if (file.exists(projRasFile))
    ras <- terra::rast(projRasFile)
  else
    stop("Cannot find projection raster specified in 'projRas'")

  if (!all(as.vector(terra::ext(rasTemplate)) == as.vector(terra::ext(ras))))
    stop("Raster parameters of global object 'rasTemplate' are not same as parameters for 'projRas'")

  # Fetch a fitted glmnet/maxnet model
  if (!silent)
  {
    cat("  Compute and apply extrapolation mask to ", basename(maxnetModel),":\n", sep = "")
    cat("     Loading maxnet model object\n")
  }

  if (file.exists(maxnetModel))
    load(maxnetModel)
  else
    stop("File referenced in 'maxnetModel' does not exist: check path")

  if (length(maxnet_model$betas) > 0)
  {
    modelFeatures <- names(maxnet_model$betas)

    mm <- gsub("I(", "", modelFeatures, fixed = TRUE)
    mm <- gsub("^2)", "", mm, fixed = TRUE)
    modelVars <- sort(unique(unlist(strsplit(mm, ":", fixed = TRUE))))

    minVarVals <- maxnet_model$varmin[modelVars]
    maxVarVals <- maxnet_model$varmax[modelVars]

    # Filter envFiles to remove rasters not need to evaluate the model
    envLayers <- base::match(modelVars, colnames(projData))
    localProjData <- projData[, envLayers]

    # Cell indices of NA values in the raster template
    naInd <- which(is.na(rasTemplate[]))

    if (!silent) cat("     Processing raster layers and coding extrapolated cells\n")

    for (i in 1:length(modelVars))
    {
      rasVals <- localProjData[, i]
      localProjData[, i] <- ifelse((rasVals >= minVarVals[modelVars[i]]) & (rasVals <= maxVarVals[modelVars[i]]), 1, 0)
    }

    if (!silent) cat("     Making mask raster layer\n")

    mask_ras <- rasTemplate
    mask_ras[] <- 1
    tmpVals <- Rfast::rowsums(localProjData)
    mask_ras[which(tmpVals != length(modelVars))] <- 0
    mask_ras[naInd] <- NA

    if (makePlots)
    {
      outStack <- terra::rast()

      layer_ras <- rasTemplate

      localProjData[naInd, i] <- NA

      terra::values(layer_ras) <- localProjData[, i]

      grDevices::png(paste0(maskOutpath, "/layer_mask_", modelVars[i], ".png"))
      terra::plot(layer_ras,
                  main = modelVars[i],
                  col = c("orange", "grey40"),
                  type = "classes",
                  levels = c("Extrapolation", "OK"))
      dev.off()

      grDevices::png(paste0(maskOutpath, "/mask_layer.png"))
      terra::plot(mask_ras,
                  main = "Mask layer",
                  col = c("orange", "grey40"),
                  type = "classes",
                  levels = c("Extrapolation", "OK"))
      dev.off()
    }

    if (saveMask)
    {
      if (is.null(fileLabel))
        maskFile <- paste0(maskOutpath, "/extrapolationMask.tif")
      else
        maskFile <- paste0(maskOutpath, "/extrapolationMask_", fileLabel, ".tif")

      terra::writeRaster(mask_ras, maskFile, overwrite = TRUE)
    }

    if (!silent) cat("     Applying mask raster to projected model raster\n")
    offInd <- which(terra::values(mask_ras) != 1)
    terra::values(ras)[offInd] <- 0

    if (makePlots)
    {
      png(paste0(maskOutpath, "/Masked_projected_raster.png"))
      plot(ras, main = "Masked projected raster")
      dev.off()
    }

    ##########################################
    if (!silent) cat("     Saving masked projected raster\n")
    if (is.null(fileLabel))
      outFilename <- gsub(".tif", "_masked.tif", projRasFile, fixed = TRUE)
    else
      outFilename <- gsub(".tif", paste0("_masked_", fileLabel, ".tif"), projRasFile, fixed = TRUE)

    terra::writeRaster(ras, outFilename, filetype = "GTiff", overwrite = TRUE)
    if (!silent) cat("     End of masking operation.\n\n")
  }
  else
  {
    cat("     Maxnet model is degenerate: no masking operation possible.\n\n")
  }
}
