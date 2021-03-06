% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/computeMESS.R
\name{computeMESS}
\alias{computeMESS}
\title{Compute MESS rasters}
\usage{
computeMESS(thisModel, envDataPath, outPath, MESSonly, ...)
}
\arguments{
\item{thisModel}{A maxnet model object.}

\item{envDataPath}{Character. Path to folder storing the environmental predictor raster layers used to fit the maxnet model.}

\item{outPath}{Character. Path to folder into which output layers will be written.}

\item{MESSonly}{Logical. Compute only the MESS raster layer? Default is TRUE. If FALSE, then MESS component rasters for each variable plus the MESS raster layer are computed}

\item{...}{Optional parameters passed to the raster function writeRaster.}
}
\value{
Nothing
}
\description{
Compute set of MESS rasters for a maxnet model
}
\details{
A MESS computation is performed using the function mess in the package dismo. This function is a wrapper which orchestrates a call to the mess function using information stored in the maxnet model object.

To facilitate the processing of very large, high-resolution rasters, only those variables which are used in the final model are included in the computation of MESS output.
}
\examples{
\dontrun{}

}
