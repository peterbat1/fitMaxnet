

#' Variable importance for a maxnet model
#'
#' Computes variable importance scores for all variables present in a maxnet model object. Scores are indicative and should be interpreted with caution as variable interaction is not considered in the computation of the scores.
#'
#' @param theModel Object of class maxnet
#' @param occSWD Data.frame. Environmental data at occurrence locations in SWD format.
#' @param bkgSWD Data.frame. Background environmental data in SWD format.
#' @param responseType Character. A MaxEnt response scale; one of "link", "exponential", "logistic" or "cloglog".
#' @param numReplicates Integer. Number of permutations performed to compute importance.
#'
#' @details {
#' The method used to compute variable importance follows that used in R packages \pkg{biomod2} and \pkg{ecospat}. First, model predictions are made for each row of the combined environmental data table formed by stacking \emph{occSWD} and \emph{bkgSWD}. This is the reference or full-model result.
#'
#' For each variable in the maxnet model object, values for the variable are permuted between rows and a model prediction made for each row using the permuted or shuffled data table. The permutation is performed \emph{numReplicates} times for each variable.
#'
#' At each permutation, a Pearson correlation is computed between reference predictions and the predicted values from the shuffled table. The importance score is 1 - correlation coefficient.
#'
#' A vector of mean scores for each variable expressed as a percentage of the sum of all mean scores is returned.
#' }
#'
#' @return A named vector of percent importance scores for the variables present in the maxnet model object sorted from highest to lowest.
#' @export
#'
#' @examples
#' \dontrun{}
#'
varImportance <- function(theModel,
                          occSWD = NULL,
                          bkgSWD = NULL,
                          responseType = c("link","exponential","cloglog","logistic"),
                          numReplicates = 5)
{
  if (!("maxnet" %in% class(theModel))) stop("Parameter 'theModel' is a naxnet object")
  if (is.null(occSWD)) stop("occSWD cannot be NULL")
  if (is.null(bkgSWD)) stop("bkgSWD cannot be NULL")

  responseType <- match.arg(responseType, c("link","exponential","cloglog","logistic"))
  if (!(responseType %in% c("link","exponential","cloglog","logistic")))
    stop("Parameter 'responseType' must be one of link, exponential, cloglog or logistic (may be abbreviated)")

  varList <- names(theModel$samplemeans)

  importance <- vector("numeric", length(varList))
  names(importance) <- varList

  envData <- rbind(occSWD[, -c(1:3)], bkgSWD[, -c(1:3)])

  fullModelVals <- predict(theModel, envData, type = responseType)

  for (thisVar in varList)
  {
    correls <- vector("numeric", numReplicates)
    origVals <- envData[, thisVar]

    for (thisRep in 1:numReplicates)
    {
      permInd <- sample(1:nrow(envData), nrow(envData))
      envData[, thisVar] <- origVals[permInd]
      newModelVals <- predict(theModel, envData)
      correls[thisRep] <- cor(fullModelVals, newModelVals)
    }

    # Re-instate original values for the variable
    envData[, thisVar] <- origVals

    # Compute the mean correlation for this variable
    importance[thisVar] <- mean(correls)
  }

  importance <- 1 - importance

  return(round(100*importance/sum(importance), 2))
}
