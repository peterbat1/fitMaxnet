###################################################################################
# NOTE: Because of frustrating limitations in the way roxygen markup works, I have
# made the help file for this function manually - edit the Rd file directly to
# make changes.
#
# Peter D. Wilson 11 October 2019
###################################################################################
#' @export
occThin <- function(occ = NA, xCol = NULL, yCol = NULL, thinDist = 0, isLatLong = TRUE)
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

  while (nrow(train) > 0)
  {
    i <- sample(1:nrow(train), 1)

    distVals <- sp::spDistsN1(as.matrix(train[, c(xCol, yCol)]), pt = unlist(train[i, c(xCol, yCol)]), longlat = isLatLong)

    # If isLatLong == FALSE, returned distance will be in metres, so convert to kilometres:
    if (!isLatLong) distVals <- distVals/1000

    # Test if the only distance in distVals <= thinDist is the distance 0.0
    # representing the distance between the test point and itself:
    if ((sum(distVals <= thinDist) == 1))
    {
      keep <- c(keep, row.names(train)[i])
    }

    # Trim the training set
    train <- train[-i,]
  }

  return(occ[keep, ])
}

