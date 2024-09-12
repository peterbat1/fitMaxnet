
#' Occurrence thinning
#' @param occ Data.frame. Data matrix to be thinned with at least two columns storing longitude (X) and latitude (Y) coordinates.
#' @param xCol Integer or Character. The column index or name of the longitude or X coordinate.
#' @param yCol Integer or Character. The column index or name of the latitude or Y coordinate.
#' @param thinDist Numeric. The distance in kilometres used to filter points.
#' @param isLatLong Logical. Are the coordinates latitude and longitude? Default is TRUE. If set to FALSE, the points are assumed to be on a projection with X & Y coordinates in metres relative an origin.
#' @param quiet Logical. Should progress messages be emitted? Default is FALSE.
#'
#' @details {
#'     This function is based on source code for the function \emph{ecospat.occ.dessagregation()} written by Olivier Broennimann and included in the R-package \emph{ecospat}. I think that the original algorithm is extremely clever and very efficient; it out-performs the \emph{thin()} function in package \emph{spThin} by at least an order of magnitude, and is faster than my "fast" interpretation of the \emph{thin()} algorithm by a factor of at least 5. However, the original code for this function had a number of quirks which made it tricky to use in a "production environment" ie for bulk processing hundreds or even thousands of species occurrence files. The following changes and improvements where made:
#' \describe{
#'  \item{Plotting}{The original source code could generate plots. I decided to leave that out allowing users to make their own plots.}
#'  \item{Input data.frame}{The original code used a rather odd and convoluted way of identifying and using columns representing the x- and y-coordinates which ultimately meant that the object returned was not the full input dataframe with 'bad' rows removed. This function does return the original dataframe with bad rows removed.}
#'  \item{Parameters}{Simplified the suite of parameters (or arguments) and giving them more meaningful names.}
#'  \item{X & Y columns}{The original code used a rather odd and convoluted way of identifying and using columns representing the x- and y-coordinates which ultimately meant that the object returned was not the full input data.frame with 'bad' rows removed. This function does return the original but thinned data.frame.}
#'  }
#' }
#'
#' @return  A data.frame with exactly the same column structure as passed in parameter \emph{occ}, but with rows for occurrence records less than \emph{thinDist} from nearest neighbours removed.
#'
#' @export
occThin <- function(occ = NA, xCol = NULL, yCol = NULL, thinDist = 0, isLatLong = TRUE, quiet = TRUE)
{
  #if (is.na(occ)) stop("No data supplied in paramater 'occ'")

  if ((is.null(xCol) || (is.null(yCol)))) stop("xCol and yCol must both have values")

  if (thinDist <= 0) stop("'thinDist' cannot be 0 or negative")

  if (any(is.na(occ[, c(xCol, yCol)])))
  {
    stop ("NA values in argument 'occ'.")
  }

  train <- occ
  keep <- NULL

  repeat
  {
    #if (!quiet) cat("Before nrow(train): ", nrow(train), "\n")

    i <- sample(1:nrow(train), 1)

    #distVals <- sp::spDistsN1(as.matrix(train[, c(xCol, yCol)]), pt = unlist(train[i, c(xCol, yCol)]), longlat = isLatLong)
    # Returned distance will be in metres, so convert to kilometres:
    if (isLatLong)
    distVals <- geodist::geodist(train[, c(xCol, yCol)], train[i, c(xCol, yCol)], measure = "geodesic")[, 1]/1000
    else
      distVals <- geodist::geodist(train[, c(xCol, yCol)], train[i, c(xCol, yCol)], measure = "cheap")[, 1]/1000

    # Test if the only distance in distVals <= thinDist is the distance 0.0
    # representing the distance between the test point and itself:
    if ((sum(distVals <= thinDist) == 1))
    {
      keep <- c(keep, row.names(train)[i])
      #if (!quiet) cat("Keeping: ", i, "\n")
    }

    # Trim the training set
    train <- train[-i,]

    #if (!quiet) cat("After nrow(train): ", nrow(train), "\n")

    if (nrow(train) == 0) break
  }

  return(occ[keep, ])
}

