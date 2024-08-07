% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/maskExtrapolation.R
\name{maskExtrapolation}
\alias{maskExtrapolation}
\title{Mask a projection raster removing extrapolation areas}
\usage{
maskExtrapolation(
  maxnetModel,
  projRasFile,
  maskOutpath = dirname(projRasFile),
  fileLabel = NULL,
  makePlots = TRUE,
  saveMask = TRUE,
  silent = TRUE
)
}
\arguments{
\item{maxnetModel}{Character. Full path to a maxnet model fitted by \link{fit_maxnet} and saved as an .Rd file.}

\item{projRasFile}{Character. Path to a raster of the projected maxnet model to which the computed extrapolation mask will be applied.}

\item{maskOutpath}{Character. Full path to the folder into which the mask raster and projection raster with mask applied will be written.}

\item{fileLabel}{Character. Label to be used to distinguish the filename of the saved mask raster.}

\item{makePlots}{Logical. Make basic plots of raster objects for interpretation, review and quality control? Default is TRUE.}

\item{saveMask}{Logical. Should the mask raster (as distinct from the the masked version of projRas) be saved? Default is TRUE.}

\item{silent}{Logical. If TRUE (default), no progress messages are written to the console.}
}
\value{
Nothing but has side-effect of writing up to two raster files: mask raster (if saveMask == TRUE) and projection raster with mask applied plus, if \emph{makePlots} == TRUE, a set of PNG graphics files.
}
\description{
Compute and apply a mask of areas within which strict extrapolation has been detected
}
\details{
{
The nominated MaxEnt model object produced by \link{fit_maxnet} is interrogated to find the list of variables used in the model. The function then proceeds to:
\enumerate{
\item Load a raster stack of these variables
\item Code each non-NA cell in a raster layer 1 if it is within the range of the variable recorded in the model object, 0 otherwise
\item Compute an output raster with only cells set to 1 when all values in the stack are 1, 0 otherwise and NA cells retained
\item Save the output raster
}
The file name for the output raster defaults to 'extrapolationMask.tif'. If \emph{fileLabel} is not NULL, then the output file name is 'extrapolationMask_' + '\emph{fileLabel}' + '.tif'.
}
}
\examples{
{
\dontrun{
## Very basic use case:
maskExtrapolation("furryCreature.Rd", "currentClimate/envVarFolder",
                  "furryCreature_run1.tif", makePlots = FALSE)
}
}
}
