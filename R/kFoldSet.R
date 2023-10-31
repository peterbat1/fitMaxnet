
#' Make k sets of indices to samples for k-fold cross-validation
#'
#' The function randomly allocates sample indices to k numeric vectors. The length of each vector is determined as the nearest integer to the number of samples times the test fraction.
#'
#'
#' @param numSamples Integer. The number of samples to be partitioned into k folds.
#' @param k Integer. The number of folds.
#' @param testSplit Numeric. The fraction of samples to be used in making the test set
#'
#' @return A list of numeric vectors of length k, where the elements of each vector are the indices of the samples forming the TEST set for that fold.
#' @export
#'
#' @examples
#' \dontrun{ ## Default testSplit
#' these_folds <- fitMaxnet::makeFolds(324, 5)
#' ## User-selected testSplit
#' these_folds <- fitMaxnet::makeFolds(261, 5, 0.3) }
makeFolds <- function(numSamples, k, testSplit = 0.2)
{
  # Do checks...
  if (numSamples < 0) stop("'numSamples' must be a positive integer")
  k <- trunc(k)
  if (k < 1) stop("'k' must be a positive integer or truncate to one")
  if ((testSplit <= 0) | (testSplit > 1)) stop("'testSplit' must be a decimal fraction (ie between 0 and 1)")

  # Enforce a reset of the RNG to default settings so that full stochastic
  # selection is made on each iteration and for each call of the function
  set.seed(NULL)

  sampleSize <- trunc(testSplit * numSamples)
  sampleInd <- 1:numSamples
  foldMembers <- vector("list", k)

  for (i in 1:k)
  {
    foldMembers[[i]] <- sample(sampleInd, sampleSize, replace = FALSE)
  }

  return(foldMembers)
}
