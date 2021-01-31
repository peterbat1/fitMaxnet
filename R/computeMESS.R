






#' Compute MESS rasters
#'
#' Compute set of MESS rasters for a maxnet model
#'
#' @param thisModel maxnet object. A maxnet model object.
#' @param envDataPath Character. Path to folder storing the environmental predictor raster layers used to fit the maxnet model.
#' @param outPath Character. Path to folder into which output layers will be written.
#' @param MESSonly Logical. Compute only the MESS raster layer? Default is TRUE. If FALSE, then MESS component rasters for each variable plus the MESS raster layer are computed
#' @param ... Optional parameters passed to the raster function writeRaster.
#' @details {
#' A MESS computation is performed using the function mess in the package dismo. This function is a wrapper which orchestrates a call to the mess function using inforamtion stored in the maxnet model object.
#'
#' To facilitate the processing of very large, high-resolution rasters, only those variables which are used in the final model are included in the computation of MESS output.
#'
#' }
#'
#' @return Nothing
#' @export
#'
#' @examples
#' \dontrun {
#' }
computeMESS <- function(thisModel, envDataPath, outPath, MESSonly, ...)
{

  # Which variables where used in the model fit?


  # Load only the required variables into a raster stack:

ans <- dismo::mess()
}
