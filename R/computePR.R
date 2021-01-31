


#' Compute Precision-Rrecall curve
#'
#' Compute a Precision-Recall (PR) curve and AUC value for a maxnet model
#'
#' @param thisModel maxnet object. A maxnet object representing the fitted model to be used
#' @param pres data.frame. Env data for pres records
#' @param bkg  data.frame. Env data for bkg records
#'
#' @return
#'  A PR curve object produced by the function \emph{pr.curve} in the package PRROC
#'
#' @export
#'
#' @examples
computePRcurve <- function(thisModel, pres, bkg)
{


  projMod <- predict(maxnet_model, projData[goodRows, ], type = type)

}



#' Title
#'
#' @param PRobj
#' @param plotFilename
#'
#' @return
#' @export
#'
#' @examples
plotPRcurve <- function(PRobj, plotFilename = NULL)
{

}
