#' Fit a maxnet model
#'
#' A wrapper around the function maxnet in the maxnet package to suit the workflow I prefer for fitting MaxEnt models to a set of taxa and tuning models across a spectrum of spatial thinning and regularisation parameter values
#'
#' @param taxonName Taxonomic name associated with the model
#' @param replTag String of characters identifying a replicate model fit
#' @param regTag String of characters representing the identifiers for regularisation values to be used. Default is prefix "reg_" added to the numeric value(s) of the regularisation values passed in the parameter regMult.
#' @param baseOutputPath Character string indicating the file system path used as a base for model fit output. The taxonomic name will be used to create a sub-folder under this base path when \emph{createFolders} is TRUE and used to prefix the output filename always.
#' @param predVar An array of integers representing the predicted variable, and corresponding to rows in envData, showing to which class each row belongs: 1 indicates presence or occupied cell, and 0 shows unoccupied background cell.
#' @param envData A numeric matrix of environmental covariates or predictor variables to be considered in the model fit
#' @param featureTypes Character value indicating the list of feature types to be used of the model fit. Default is 'lpq' meaning that linear, product and quadratic features will be used. Other options are available (see maxnet and MaxEnt documentation/literature) but not encouraged for ENMs as they lead to serious overfitting.
#' @param regMult A numeric array of regularisation values to be used. A single numeric value maybe passed leading to just one fitted model at that regularisation value. More than one value will gerenate a sequence of fitted models, one for each regularisation value in the array.
#' @param outputType Character. Output scaling of fitted model: "link" is raw linear predictor scores, while "exponential", "logistic" and "cloglog" generate a non-linear re-scaling of the raw linear predictor scores. See \link[maxnet]{maxnet} for further information.
#' @param createFolders Logical. Should the function create sub-folders below \emph{baseOutPath} using taxon name, replicate number and regularisation values? If FALSE (default) then replTag and regTag are ignored
#' @param quiet Logical. Should the function proceed without emitting messages?
#'
#' @return A maxnet object
#' @export
#'
#' @examples
fit_maxnet <- function(taxonName = NULL,
                       replTag = NULL,
                       regTag = paste0("reg_", gsub(".","_", regMult, fixed = TRUE)),
                       baseOutputPath = NULL,
                       predVar,
                       envData,
                       featureTypes = "lpq",
                       regMult = 1,
                       outputType = "link",
                       createFolders = FALSE,
                       quiet = TRUE)
{
  if (is.null(taxonName)) stop("taxonName cannot be NULL")

  if (is.null(baseOutputPath)) stop("baseOutputPath cannot be NULL")

  if (length(predVar) == 0) stop("Array predVar contains no values")

  if (class(predVar) != "numeric") stop("Array predVar must be numeric")

  if ((nrow(envData) == 0) || (is.null(envData))) stop("Matrix envData contains no values")

  if (!(outputType %in% c("link", "exponential", "logistic", "cloglog"))) stop("Unkown value in 'outputType'")

  if (!quiet)
  {
    if (is.null(replTag))
      cat("    Fitting model for", taxonName, "with regularisation multiplier =", regMult, "\n")
    else
      cat("    Fitting model for", taxonName, "with replTag = ", replTag, "and regularisation multiplier =", regMult, "\n")
  }

  ### try() ??
  maxnet_model <- maxnet::maxnet(predVar,
                                 envData,
                                 f = maxnet::maxnet.formula(predVar, envData, classes = featureTypes),
                                 regmult = regMult,
                                 type = outputType)

  # Save object
  if (createFolders)
  {
    outputFolder <- paste0(baseOutputPath, "/", taxonName)
    if (!is.null(replTag)) outputFolder <- paste0(outputFolder, "/", replTag)
    outputFolder <- paste0(outputFolder, "/", regTag)
    #print(outputFolder)
    if (!dir.exists(outputFolder)) dir.create(outputFolder, recursive = TRUE)
  }
  else
    outputFolder <- baseOutputPath

  fileName <- paste0(outputFolder, "/", gsub(" ", "_", taxonName, fixed = TRUE), ifelse(is.null(replTag), "", paste0("_", replTag)), "_", regTag, ".Rd")

  if (!quiet) cat("    Saving model object:", fileName, "\n")
  save(list(maxnet_model, type = outputType, runDate = Sys.Date()), file = fileName)
  return(maxnet_model)
}
