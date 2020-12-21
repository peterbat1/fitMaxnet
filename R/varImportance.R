

#' Variable importance for a maxnet model
#'
#' @param theModel Object of class maxnet
#' @param occSWD Data.frame. Environmental data at occurrence locations in SWD format.
#' @param bkgSWD Data.frame. Background environmental data in SWD format.
#' @param responseType Character. A MaxEnt response scale; one of "link", "exponential", "logistic" or "cloglog".
#' @param numReplicates Integer. Number of permutations performed to compute importance.
#'
#' @return A named vector of percent importance scores for all variables present in the maxnet model object.
#' @export
#'
#' @examples
varImportance <- function(theModel,
                          occSWD,
                          bkgSWD,
                          responseType = c("link","exponential","cloglog","logistic"),
                          numReplicates = 5)
{
  varList <- names(theModel$samplemeans)

  importance <- vector("numeric", length(varList))
  names(importance) <- varList

  envData <- rbind(occSWD[, -c(1:3)], bkgSWD[, -c(1:3)])

  fullModelVals <- predict(theModel, envData)

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

    envData[, thisVar] <- origVals
    importance[thisVar] <- mean(correls)
  }

  importance <- 1 - importance

  return(round(100*importance/sum(importance), 2))
}
