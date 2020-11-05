
#' Project maxnet model
#'
#' Project a maxnet model fitted using fitModel
#'
#' @details {
#' Function \link{prepProjData} \strong{MUST} be run before calling this function to prepare global data objects needed for the projection.}
#'
#' @param taxonName String. Taxonomic name associated with this model; used to construct a file name for the output raster file.
#' @param maxnetModel String. \emph{Full} path name to the .Rd file storing a fitted maxnet model produced by the companion functions \link{fit_maxnet} and \link{fitModels}.
#' @param baseOutputPath String. Full path to the output folder to receive the produced raster.
#' @param fileLabel String. An identifying tag to be included in the output filename.
#' @param makeTaxonFolder Logical. Should a sub-folder on \emph{baseOutputPath} be created using \emph{taxonName}?
#'
#' @return Nothing
#' @export
#'
#' @examples
#' \dontrun{}
projectMaxnet <- function(taxonName = NULL,
                          maxnetModel,
                          type = "exponential",
                          baseOutputPath,
                          fileLabel = NULL,
                          makeTaxonFolder = FALSE)
{
  # Check params...
  if(!(type %in% c('link', 'exponential', 'cloglog', 'logistic')))
    stop("'type' must be one of 'link', 'exponential', 'cloglog', 'logistic'")

  # Check for existence of projStack, etc in the global environment...

  #featureSet_names <- names(maxnet_model$betas)
  #varNames <- sort(unique(unlist(strsplit(gsub("^2)","",gsub("I(","",featureSet_names, fixed = TRUE), fixed = TRUE), ":", fixed = TRUE))))
  cat("Maxnet model projection:\n    Loading model object\n")

  load(maxnetModel)

  goodRows <- which(!is.na(rowSums(projData)))

  cat("    Projecting model\n")
  projMod <- predict(maxnet_model, projData[goodRows, ], type = type)

  cat("    Preparing and saving projection raster\n")
  projRas <- rasTemplate #projStack[[1]]
  raster::values(projRas) <- NA
  raster::values(projRas)[goodRows] <- projMod[,1]

  if (makeTaxonFolder)
    outputPath <- paste0(baseOutputPath, "/", taxonName)
  else
    outputPath <- baseOutputPath

  if (!dir.exists(outputPath)) dir.create(outputPath, recursive = TRUE)

  if (is.null(fileLabel))
    outputPath <- paste0(outputPath, "/", paste0(gsub(" ", "_", taxonName, fixed = TRUE)), "_projection.tif")
  else
    outputPath <- paste0(outputPath, "/", paste0(gsub(" ", "_", taxonName, fixed = TRUE)), "_projection_",fileLabel,".tif")

  raster::writeRaster(projRas, outputPath, format = "GTiff", overwrite = TRUE)
  #plot(projRas, main = paste0("Regularization = ", thisRegVal))
  cat("  End model projection\n\n")
}
