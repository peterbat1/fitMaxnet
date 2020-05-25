library(maxnet)
library(raster)

# A wrapper around the function maxnet in the maxnet package to suit the
# workflow I prefer for fitting MaxEnt models to a set of taxa and tuning models
# across a spectrum of spatial thinning and regularisation parameter values
fit_maxnet <- function(taxonName = NULL,
                       replTag = NULL,
                       regTag = paste0("reg_", gsub(".","_", regMult, fixed = TRUE)),
                       baseOutputPath = NULL,
                       predVar,
                       envData,
                       featureTypes = "lpq",
                       regMult = 1)
{
  if (is.null(taxonName)) stop("taxonNamne cannot be NULL")

  if (is.null(baseOutputPath)) stop("baseOutputPath cannot be NULL")

  if (length(predVar) == 0) stop("Array predVar contains no values")

  if (class(predVar) != "numeric") stop("Array predVar must be numeric")

  if ((nrow(envData) == 0) || (is.null(envData))) stop("Matrix envData contains no values")

  if (is.null(replTag))
    cat("    Fitting model for", taxonName, "with regularisation multiplier =", regMult, "\n")
  else
    cat("    Fitting model for", taxonName, "with replTag = ", replTag, "and regularisation multiplier =", regMult, "\n")

  ### try() ??
  maxnet_model <- maxnet(predVar, envData, f = maxnet.formula(predVar, envData, classes = featureTypes), regmult = regMult)

  # Save object
  outputFolder <- paste0(baseOutputPath, "/", taxonName)
  if (!is.null(replTag)) outputFolder <- paste0(outputFolder, "/", replTag)
  outputFolder <- paste0(outputFolder, "/", regTag)
  #print(outputFolder)
  if (!dir.exists(outputFolder)) dir.create(outputFolder, recursive = TRUE)

  fileName <- paste0(outputFolder, "/", gsub(" ", "_", taxonName, fixed = TRUE), ifelse(is.null(replTag), "", paste0("_", replTag)), "_", regTag, ".Rd")

  cat("    Saving model object:", fileName, "\n")
  save(maxnet_model, file = fileName)
  return(maxnet_model)
}


# Load raster stack of env data layers to be used for projecting a maxnet model.
# Matrix of data values and array of indices pointing to rows with no missing
# values are placed in the global environment.
prepData <- function(dataPath)
{
  cat("Preparing projection data: ")
  projStack <<- stack(list.files(dataPath, "*.tif", full.names = TRUE))
  projData <<- values(projStack)
  goodCellInd <<- which(!is.na(rowSums(projData)))
  cat(length(projStack), " layers loaded\n")
  invisible(NULL)
}


# Organise data and fit a set of MaxEnt models using the maxnet package
fitModel <- function(taxonName = NULL,
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
      values(projRas) <- NA
      values(projRas)[goodCellInd] <- projMod[,1]

      outputPath <- paste0(baseOutputPath, "/", taxonName)
      if (!is.null(thisTag)) outputPath <- paste0(outputPath, "/", thisTag)
      #outputPath <- paste0(outputPath, "/", thisTag)

      regStr <- paste0("reg_", gsub(".","_", thisRegVal, fixed = TRUE))

      outputPath <- paste0(outputPath, "/", regStr, "/", paste0(gsub(" ", "_", taxonName, fixed = TRUE), ifelse(is.null(thisTag), "", paste0("_", thisTag)), "_reg_", thisRegVal, ".tif"))
      writeRaster(projRas, outputPath, format = "GTiff", overwrite = TRUE)
      #plot(projRas, main = paste0("Regularization = ", thisRegVal))
      cat("  End model fit\n")
    }
  }
}


