# Organise data and fit a set of MaxEnt models using the maxnet package
#' Fit a set of MaxEnt models using the maxnet package
#'
#' A set of MaxEnt models will be fitted using the maxnet function in the package maxnet. A model will be fitted for each combination of replTags and multSet values.
#'
#' @param taxonName The taxonomic name associated with a fitted model.
#' @param occPath File system path to the file holding occurrence records and associated environmental data.
#' @param backgroundPath File system path to a file holding a matrix of environmental data representing the background environment for the supplied occurrence data.
#' @param nBkgPoints Number (integer) giving the number of background points to sub-sample from data in the file pointed to by \emph{backgroundPath}.
#' @param baseOutputPath File system path to the base folder to be used for model output. A sub-folder with the taxonomic name supplied in \emph{taxonName}. It will be created if it doesn't already exist to receive output.
#' @param replTags A character array of tags representing replicates to be fitted.
#' @param multSet A numeric array of values to be used for overall regularisation multipliers. At least one value must be present which may be the default value of 1.
#'
#' @return Nothing
#' @export
#'
#' @examples
fitModels <- function(taxonName = NULL,
                     occPath = NULL,
                     backgroundPath = NULL,
                     nBkgPoints = 5000,
                     baseOutputPath = NULL,
                     replTags = "",
                     multSet = 1)
{
  if (is.null(taxonName)) stop("taxonName must have value")
  if (is.null(occPath)) stop("occPath must have a value")
  if (is.null(backgroundPath)) stop("backgroundPath must have a value")
  if (is.null(baseOutputPath)) stop("baseOutputPath must have a value")
  if (class(multSet) != "numeric") stop("multSet must be an array of one or more numeric values")

  if (!file.exists(occPath)) stop("Cannot find file referenced in parameter 'occPath'")
  if (!file.exists(backgroundPath)) stop("Cannot find file referenced in parameter 'backgroundPath'")

  if (!exists("projData")) stop("Global object 'projData' is missing. Did you run prepData()?")

  # Load occurrence data and clean by removing rows with missing values
  occ <- read.csv(occPath, stringsAsFactors = FALSE)
  badInd <- which(is.na(rowSums(occ[, 4:ncol(occ)])))
  if (length(badInd) > 0) occ <- occ[-badInd, ]

  # Load background data and clean by removing rows with missing values
  bkg <- read.csv(backgroundPath, stringsAsFactors = FALSE)
  badInd <- which(is.na(rowSums(bkg[, 4:ncol(bkg)])))
  if (length(badInd) > 0) bkg <- bkg[-badInd, ]

  if (nrow(bkg) < nBkgPoints)
  {
    nBkgPoints <- nrow(bkg)
  }

  # Prep data set for model fitting call
  envData <- rbind(occ[, 4:ncol(occ)],
                   bkg[sample(1:nrow(bkg), nBkgPoints), 4:ncol(bkg)])

  predVar <- c(rep(1, nrow(occ)), rep(0, nBkgPoints))

  cat("  Start model fit:\n")
  for (replTag in replTags)
  {
    for (thisRegVal in multSet)
    {
      #cat(thisRegVal, "\n")

      thisRegStr <- paste0("reg_", as.character(sub(".", "_", thisRegVal, fixed = TRUE)))

      if (replTag == "")
        thisTag <- NULL
      else
        thisTag <- replTag

      ans <- fit_maxnet(taxonName,
                        replTag = thisTag,
                        baseOutputPath = baseOutputPath,
                        predVar = predVar,
                        envData = envData,
                        featureTypes = "lpq",
                        regMult = thisRegVal)

      cat("    Projecting model\n")
      projMod <- predict(ans, projData, type = "cloglog")

      cat("    Preparing and saving projection raster\n")
      projRas <- projStack[[1]]
      raster::values(projRas) <- NA
      raster::values(projRas)[goodCellInd] <- projMod[,1]

      outputPath <- paste0(baseOutputPath, "/", taxonName)
      if (!is.null(thisTag)) outputPath <- paste0(outputPath, "/", thisTag)
      #outputPath <- paste0(outputPath, "/", thisTag)

      regStr <- paste0("reg_", gsub(".","_", thisRegVal, fixed = TRUE))

      outputPath <- paste0(outputPath, "/", regStr, "/", paste0(gsub(" ", "_", taxonName, fixed = TRUE), ifelse(is.null(thisTag), "", paste0("_", thisTag)), "_", thisRegStr, ".tif"))
      raster::writeRaster(projRas, outputPath, format = "GTiff", overwrite = TRUE)
      #plot(projRas, main = paste0("Regularization = ", thisRegVal))
      cat("  End model fit\n\n")
    }
  }
}
