# Explore methods for marking rasters with areas exhibiting strict extrapolation
# after fitting a model using the maxnet package
#
#
# Peter D. Wilson
#
# 2019-10-15

#library(raster)


#' Mask a projection raster removing extrapolation areas
#'
#' Compute and apply a mask of areas within which strict extrapolation has been detected
#'
#' @param maxnetModel Full path to a maxnet model fitted by \link{fit_maxnet} and saved as an .Rd file.
#' @param envPath Path to a set of environmental rasters showing the environment onto which the model has been projected.
#' @param projRas Raster of the projected maxnet model.
#' @param maskOutpath Full path to the folder into which the mask raster and projection raster with mask applied will be written
#' @param fileLabel Label to be used to distinguish the filename of the saved mask raster.
#' @param makePlots. Logical. Make basic plots of raster objects for review and quality control. Default is TRUE.
#'
#' @details {
#' Some info on method + details of coding used in the output raster...
#'
#' Note method for forming an output name for the mask raster
#' }
#'
#' @return Nothing but has side-effect of writing two raster files: mask raster and projection raster with mask applied.
#' @export
#'
#' @examples
maskExtrapolation <- function(maxnetModel, envPath, projRas, maskOutpath = dirname(projRas), fileLabel, makePlots = TRUE)
{
  # Fetch a fitted glmnet/maxnet model
  cat("Compute and apply extrapolation mask to ", basename(maxnetModel),":\n", sep = "")
  cat("  Loading maxnet model object\n")
  load(maxnetModel)

  if (length(maxnet_model$betas) > 0)
  {
    modelFeatures <- names(maxnet_model$betas)

    mm <- gsub("I(","",modelFeatures, fixed = TRUE)
    mm <- gsub("^2)", "", mm, fixed = TRUE)
    modelVars <- sort(unique(unlist(strsplit(mm, ":", fixed = TRUE))))

    minVarVals <- maxnet_model$varmin[modelVars]
    maxVarVals <- maxnet_model$varmax[modelVars]

    cat("  Loading environmental data layers as a raster stack\n")
    envFiles <- list.files(envPath, "*.tif", full.names = TRUE)

    envStack <- raster::stack(envFiles)

    outStack <- envStack

    # Trim outStack by removing vars not found in the list of vars i.e. mmm
    dropInd <- which(!(names(outStack) %in% modelVars))
    outStack <- raster::dropLayer(outStack, dropInd)

    cat("  Processing raster layers and coding extrapolated cells\n")

    for (i in 1:length(modelVars))
    {
      aRas <- envStack[[modelVars[i]]]
      rasVals <- raster::getValues(aRas)
      naInd <- which(is.na(rasVals))
      rasVals[(rasVals >= minVarVals[modelVars[i]]) & (rasVals <= maxVarVals[modelVars[i]])] <- 1
      rasVals[rasVals != 1] <- 0
      rasVals[naInd] <- NA
      raster::values(aRas) <- rasVals
      outStack[[modelVars[i]]] <- aRas
    }

    cat("  Making mask raster layer\n")
    ans <- raster::calc(outStack, min)
    maskFile <- paste0(maskOutpath, "/extrapolationMask_", fileLabel, ".tif")

    raster::writeRaster(ans, maskFile, overwrite = TRUE)

    if (makePlots)
    {
      plot(outStack)
      plot(ans, main = "Mask layer")
    }

    cat("  Applying mask raster to projected model raster\n")
    ras <- raster::raster(projRas)

    offInd <- which(values(ans) != 1)
    raster::values(ras)[offInd] <- 0

    if (makePlots) plot(ras, main = "Masked projection raster")

    cat("  Saving masked projection raster\n")
    outFilename <- gsub(".tif", "_masked.tif", projRas, fixed = TRUE)
    raster::writeRaster(ras, outFilename, format = "GTiff", overwrite = TRUE)
    cat("  End of masking operation.\n\n")
  }
  else
  {
    cat("  Maxnet model is degenerate: no masking operation possible.\n\n")
  }
}
