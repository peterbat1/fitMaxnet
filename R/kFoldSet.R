
#' Make k sets of indices to samples for k-fold cross-validation
#'
#' The function randomly allocates sample indices to k numeric vectors. The length of each vector is set initially to be the nearest whole number which, when multiplied by the number of folds, comes as close as possible to the number of samples. Any residual of un-allocated sample indices remaining after the initial allocation are distributed randomly to the k groups. Each selected group receives only 1 additional sample index.
#'
#' @param numSamples Integer. The number of samples to be partitioned into k folds.
#' @param k Integer. The number of folds.
#'
#' @return A list of numeric vectors of length k, where the elements of each vector are the indices of the samples to be included on that fold.
#' @export
#'
#' @examples
#' \dontrun{ }
makeFolds <- function(numSamples, k)
{
  # Do checks...


  set.seed(1953)
  targetGroupSize <- numSamples %/% k
  sampleInd <- 1:numSamples

  grpSizes <- rep(targetGroupSize, k)

  numResidualSamples <- numSamples %% k

  # Select numResidualSamples among the groups to receive one more sample
  extraInd <- sample(1:k, numResidualSamples, replace = FALSE)

  grpSizes[extraInd] <- grpSizes[extraInd] + 1

  grpMembers <- vector("list", length(grpSizes))

  for (i in 1:length(grpSizes))
  {
    grpMembers[[i]] <- sample(sampleInd, grpSizes[i], replace = FALSE)
    sampleInd <- sampleInd[-grpMembers[[i]]]
  }

  return(grpMembers)
}
