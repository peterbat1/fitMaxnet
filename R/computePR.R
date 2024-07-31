


#' Compute Precision-Recall curve
#'
#' Compute a Precision-Recall (PR) curve and AUC value for a maxnet model
#'
#' @param thisModel maxnet object. A maxnet object representing the fitted model
#' @param occSWD data.frame. Environmental data for presence or occurrence records in MaxEnt SWD format
#' @param bkgSWD  data.frame. Environmental data for background locations in MaxEnt SWD format
#' @param responseType Character. A MaxEnt response scale; one of "link", "exponential", "logistic" or "cloglog". Default is "cloglog"
#'
#' @return
#'  A PR curve object produced by the function \emph{pr.curve} in the package \pkg{PRROC}
#'
#' @export
#'
#' @examples
computePRcurve <- function(thisModel, occSWD, bkgSWD, responseType = "cloglog")
{

  if (!inherits(thisModel, "maxnet"))
    stop("'thisModel' must be a maxnet model object")

  if (!inherits(occSWD, "data.frame"))
    stop("'occSWD' must be a data.frame")

  if (!inherits(bkgSWD, "data.frame"))
    stop("'occSWD' must be a data.frame")

  if (!(responseType %in% c("link","exponential","cloglog","logistic")))
    stop("'responseType' must be one of 'link', 'exponential', 'logistic' or 'cloglog'")

  modelVarNames <- names(thisModel$varmax)

  if (!all(modelVarNames %in% colnames(occSWD)[4:ncol(occSWD)]))
    stop("MaxEnt model variable names do not match variables in 'occSWD'")

  if (!all(modelVarNames %in% colnames(bkgSWD)[4:ncol(bkgSWD)]))
    stop("MaxEnt model variable names do not match variables in 'bkgSWD'")

  occ_scores <- predict(thisModel, occSWD[, -c(1:3)], type = responseType)
  bkg_scores <- predict(thisModel, bkgSWD[, -c(1:3)], type = responseType)

  return(PRROC::pr.curve(occ_scores, bkg_scores, curve = TRUE))
}



#' Plot a Precision-Recall (PR) Curve
#'
#' @param PRobj pr.curve object. A PR ROC object with type = "PR" produced by the function \link{pr.curve} in the package \pkg{PRROC} run with parameter curve = TRUE
#' @param plotFilename Character (string). Optional file name to save the plot output (as PNG formatted graphics file)
#'
#' @return
#' A ggplot2 object
#'
#' @export
#'
#' @examples
plotPRcurve <- function(PRobj, plotFilename = NULL)
{

  if (!inherits(PRobj, "PRROC"))
    stop("'PRobj' does not contain an 'PRROC' object")

  if (PRobj$type != "PR")
    stop("PRobj must be type 'PR'")

}
