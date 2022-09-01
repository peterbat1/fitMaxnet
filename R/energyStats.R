# Energy statistics for difference between two sampoles form a multivariate distribution
#
# Peter D. Wilson
# 2020-11-26

#' Compute difference between two multivariate samples
#'
#' @param sample1 Numeric. Matrix or data.frame with values sampled from a multivariate distribution.
#' @param sample2 Numeric. Matrix or data.frame with a second example of values sampled from a multivariate distribution.
#' @param edist_method Character. Selects the method used in a call to function \link{eqdist.etest} in the package \emph{energy}.
#' @param edist_replicates Numeric (integer). Number of replicates used in energy statistics calculation. Default is 199 which appears to give reliable results without significant computational cost.
#'
#' @details {
#' This function is a wrapper around the function \emph{eqdist.etest} in the package \emph{energy} designed to facilitate its orderly and robust application to testing the difference between a base sample and a second sample (typically a sub-sample) of points in a multivariate environmental space.
#'
#'
#' }
#'
#' @return An object of class \emph{htest} as returned by function \link{eqdist.etest}.
#' @export
#'
#' @examples
#' \dontrun{
#' m1 <- matrix(runif(100), 10, 10)
#' df1 <- data.frame(m1)
#' m2 <- matrix(runif(100), 10, 10)
#' df2 <- data.frame(m2)
#' ans1 <- energyStats(m1, m2)
#' print(ans1)
#'
#'
#' ans2 <- energyStats(df1, df2)
#' print(ans2)
#'
#' m2 <- matrix(runif(120), 10, 12)
#' df2 <- data.frame(m2)
#' ans <- energyStats(df1, df2)
#' print(ans)
#' }
energyStats <- function(sample1 = NULL,
                        sample2 = NULL,
                        edist_method = "original",
                        edist_replicates = 199)

{
  if (is.null(sample1)) stop("'sample1' requires a value")
  if (is.null(sample2)) stop("'sample2' requires a value")
  #if (!dir.exists(envDataPath)) stop("Path given in 'envDataPath' cannot be found")

  if ("matrix" %in% class(sample1)) sample1 <- data.frame(sample1)
  if ("matrix" %in% class(sample2)) sample2 <- data.frame(sample2)

  if (!("data.frame" %in% class(sample1))) stop("'sample1' must be a data frame or matrix")
  if (!("data.frame" %in% class(sample2))) stop("'sample2' must be a data frame or matrix")

  if (!(any(colnames(sample2) %in% colnames(sample1))))
    stop("'sample1' and 'sample2' do not share any variables")

  if (length(colnames(sample1)) < length(colnames(sample2)))
  {
    missingVars <- colnames(sample2)[which(!(colnames(sample2) %in% colnames(sample1)))]
    sample2 <- sample2[, colnames(sample1)]
    cat("These variables are present in 'sample2' but not in 'sample1':\n")
    cat(paste(missingVars, collapse = ", "), "\n")
    cat("They will be ignored in energy statistic calculations\n")
  }

  if (length(colnames(sample2)) < length(colnames(sample1)))
  {
    missingVars <- colnames(sample1)[which(!(colnames(sample1) %in% colnames(sample2)))]
    sample1 <- sample1[, colnames(sample2)]
    cat("These variables are present in 'sample1' but not in 'sample2':\n")
    cat(paste(missingVars, collapse = ", "), "\n")
    cat("They will be ignored in energy statistic calculations\n")
  }

  if (!(edist_method %in% c("original", "discoB", "discoF")))
    stop("Unrecognised value given for 'edist_method'")

  if (!is.numeric(edist_replicates))
    stop("'edist_replicates' must be an integer value")

  comboData <- rbind(sample1, sample2)

  return(eqdist.etest(comboData,
                      sizes = c(nrow(sample1), nrow(sample2)),
                      method = edist_method,
                      R = edist_replicates))
}
