
#' Make k sets of indices to samples for k-fold cross-validation
#'
#' The function randomly allocates sample indices to k numeric vectors. The length of each vector is determined as the nearest integer to the number of samples times the test fraction.
#'
#'
#' @param numSamples Integer. The number of samples to be partitioned into k folds.
#' @param k Integer. The number of folds.
#' @param trainSplit Numeric. The fraction of samples to be used in
#'
#' @return A list of numeric vectors of length k, where the elements of each vector are the indices of the samples forming the TEST set for that fold.
#' @export
#'
#' @examples
#' \dontrun{ }
makeFolds <- function(numSamples, k, testSplit = 0.2)
{
  # Do checks...

  set.seed(1953)

  sampleSize <- trunc(testSplit * numSamples)
  sampleInd <- 1:numSamples
  foldMembers <- vector("list", k)

  for (i in 1:k)
  {
    foldMembers[[i]] <- sample(sampleInd, sampleSize, replace = FALSE)
  }

  return(foldMembers)
}
