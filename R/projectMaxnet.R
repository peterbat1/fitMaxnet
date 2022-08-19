
#' Project maxnet model
#'
#' Project a maxnet model fitted using fitModel
#'
#' @details {
#' Function \link{prepProjData} \strong{MUST} be run before calling this function to prepare global data objects needed for the projection.
#'
#' \strong{NOTE:} The resulting raster is written to the specified output folder with the file name composed by concatenating the taxon name spaces replaced by underscores), "projection" and, if supplied, the character object \emph{fileTag}.
#'
#' \strong{If a file with this name already exists in the output folder, it will be overwritten.}
#'
#' If dataType = "INT4S" is given, output is multiplied by 1,000 and truncated to an integer. Therefore, real values (to 4 dec. places) must be regenerated after loading the TIFF file by dividing the raster values by 1,000.
#' }
#'
#' @param taxonName String. Taxonomic name associated with this model; used to construct a file name for the output raster file.
#' @param maxnetModel String. \emph{Full} path name to the .Rd file storing a fitted maxnet model produced by the companion functions \link{fit_maxnet} and \link{fitModels}.
#' @param type String. The type of scaling applied to predicted model values.
#' @param doClamp Logical. Should values of predictors (covariates) be clamped to those seen during model fitting? Default is FASLE, no clamping.
#' @param dataType String. Data format to be used. Default is "FLT4S". Only alternative at moment is "INT4S".
#' @param baseOutputPath String. Full path to the output folder to receive the produced raster.
#' @param fileLabel String. An identifying tag to be included in the output filename.
#' @param makeTaxonFolder Logical. Should a sub-folder on \emph{baseOutputPath} be created using \emph{taxonName}?
#' @param quiet Logical. Proceed without emitting progress messages?
#'
#' @return Nothing
#' @export
#'
#' @examples
#' \dontrun{}
projectMaxnet <- function(taxonName = NULL,
                          maxnetModel,
                          type = "exponential",
                          doClamp = FALSE,
                          dataType = "FLT4S",
                          baseOutputPath,
                          fileLabel = NULL,
                          makeTaxonFolder = FALSE,
                          quiet = TRUE)
{
  # Check params...
  if(!(type %in% c('link', 'exponential', 'cloglog', 'logistic')))
    stop("'type' must be one of 'link', 'exponential', 'cloglog', 'logistic'")

  # Check for existence of projData, etc in the global environment...
  if (!exists("projData")) stop("Object projData not found in global environment. Please run 'prepProjData'")
  if (!exists("rasTemplate")) stop("Object rasTemplate not found in global environment. Please run 'prepProjData'")

  if (!quiet) cat("Maxnet model projection:\n    Loading model object\n")

  load(maxnetModel)

  goodRows <- which(!is.na(Rfast::rowsums(projData)))

  if (!quiet) cat("    Projecting model\n")
  projMod <- predict(maxnet_model, projData[goodRows, ], type = type, clamp = doClamp)

  if (!quiet) cat("    Preparing and saving projection raster\n")
  projRas <- rasTemplate
  terra::values(projRas) <- NA
  terra::values(projRas)[goodRows] <- round(projMod[,1], 4)

  if (makeTaxonFolder)
    outputPath <- paste0(baseOutputPath, "/", taxonName)
  else
    outputPath <- baseOutputPath

  if (!dir.exists(outputPath)) dir.create(outputPath, recursive = TRUE)

  if (is.null(fileLabel))
    outputPath <- paste0(outputPath, "/", paste0(gsub(" ", "_", taxonName, fixed = TRUE)), "_projection.tif")
  else
    outputPath <- paste0(outputPath, "/", paste0(gsub(" ", "_", taxonName, fixed = TRUE)), "_projection_", fileLabel, ".tif")

  terra::writeRaster(projRas,
                     outputPath,
                     filetype = "GTiff",
                     overwrite = TRUE,
                     datatype = dataType,
                     gdal = "COMPRESS=DEFLATE")

  if (!quiet) cat("  End model projection\n\n")
}
