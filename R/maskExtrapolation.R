# Explore methods for marking rasters with areas exhibiting strict extrapolation
# after fitting a model using the maxnet package
#
#
# Peter D. Wilson
#
# 2019-10-15

#' Mask a projection raster removing extrapolation areas
#'
#' Compute and apply a mask of areas within which strict extrapolation has been detected
#'
#' @param maxnetModel Character. Full path to a maxnet model fitted by \link{fit_maxnet} and saved as an .Rd file.
#' @param envStack Raster stack object. Stack of environmental rasters onto which the model has been projected.
#' @param projRas RasterLayer. Raster of the projected maxnet model to which the computed extrapolation mask will be applied.
#' @param maskOutpath Character. Full path to the folder into which the mask raster and projection raster with mask applied will be written
#' @param fileLabel Character. Label to be used to distinguish the filename of the saved mask raster.
#' @param makePlots Logical. Make basic plots of raster objects for interpretation, review and quality control? Default is TRUE.
#' @param silent Logical. If TRUE (default), no progress messages are written to the console.
#'
#' @details {
#' The nominated MaxEnt model object produced by \link{fit_maxnet} is interrogated to find the list of variables used in the model. The function then proceeds to:
#' \enumerate{
#' \item Load a raster stack of these variables
#' \item Code each non-NA cell a raster layer 1 if it is within the range of the variable recorded in the model object, 0 otherwise
#' \item Compute an output raster with only cells set to 1 when all values in the stack are 1, 0 otherwise and NA cells retained
#' \item Save the output raster
#'}
#' The file name for the output raster defaults to 'extrapolationMask.tif'. If \emph{fileLabel} is not NULL, then the output file name is 'extrapolationMask_' + '\emph{fileLabel}' + '.tif'.
#' }
#'
#' @return Nothing but has side-effect of writing two raster files: mask raster and projection raster with mask applied plus, if \emph{makePlots} == TRUE, a set of PNG graphics files.
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
                              envStack = NULL,
                              projRas,
                              maskOutpath = dirname(projRas),
                              fileLabel = NULL,
                              makePlots = TRUE)
{
  # Valid raster stack is available?
  if (is.null(envStack)) stop("Parameter 'envStack' has not be given a value")
  if (class(envStack) != "RasterStack") stop("Parameter 'envStack' must be class 'RasterStack'")

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

    mm <- gsub("I(","",modelFeatures, fixed = TRUE)
    mm <- gsub("^2)", "", mm, fixed = TRUE)
    modelVars <- sort(unique(unlist(strsplit(mm, ":", fixed = TRUE))))

    minVarVals <- maxnet_model$varmin[modelVars]
    maxVarVals <- maxnet_model$varmax[modelVars]

    # Filter envFiles to remove rasters not need to evaluate the model

    envLayers <- base::match(modelVars, names(envStack))
    localEnvStack <- raster::subset(envStack, envLayers)

    outStack <- localEnvStack

    # Trim outStack by removing vars not found in the list of vars i.e. mmm
    # dropInd <- which(!(names(outStack) %in% modelVars))
    # outStack <- raster::dropLayer(outStack, dropInd)

    if (!silent) cat("     Processing raster layers and coding extrapolated cells\n")

    for (i in 1:length(modelVars))
    {
      aRas <- localEnvStack[[modelVars[i]]]
      rasVals <- raster::getValues(aRas)
      naInd <- which(is.na(rasVals))
      #### ?.bincode
      rasVals[(rasVals >= minVarVals[modelVars[i]]) & (rasVals <= maxVarVals[modelVars[i]])] <- 1
      rasVals[rasVals != 1] <- 0
      rasVals[naInd] <- NA
      raster::values(aRas) <- rasVals
      outStack[[modelVars[i]]] <- aRas
    }

    if (!silent) cat("     Making mask raster layer\n")
    ans <- raster::calc(outStack, min)
    if (is.null(fileLabel))
      maskFile <- paste0(maskOutpath, "/extrapolationMask.tif")
    else
      maskFile <- paste0(maskOutpath, "/extrapolationMask_", fileLabel, ".tif")

    raster::writeRaster(ans, maskFile, overwrite = TRUE)

    if (makePlots)
    {
      png(paste0(maskOutpath, "/layer_masks.png"))
      plot(outStack)
      dev.off()
      png(paste0(maskOutpath, "/mask_layer.png"))
      plot(ans, main = "Mask layer")
      dev.off()
    }

    if (!silent) cat("     Applying mask raster to projected model raster\n")
    ras <- raster::raster(projRas)

    offInd <- which(values(ans) != 1)
    raster::values(ras)[offInd] <- 0

    if (makePlots)
    {
      png(paste0(maskOutpath, "/Masked_projected_raster.png"))
      plot(ras, main = "Masked projected raster")
      dev.off()
    }

    if (!silent) cat("     Saving masked projected raster\n")
    outFilename <- gsub(".tif", "_masked.tif", projRas, fixed = TRUE)
    raster::writeRaster(ras, outFilename, format = "GTiff", overwrite = TRUE)
    if (!silent) cat("     End of masking operation.\n\n")
  }
  else
  {
    cat("     Maxnet model is degenerate: no masking operation possible.\n\n")
  }
}
