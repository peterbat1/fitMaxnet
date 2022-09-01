

#' Response plots for maxnet-produced MaxEnt models
#'
#' Response plots are produced for some or all variables in the supplied maxnet model object.
#'
#' @param thisModel A maxnet model object.
#' @param variable Character. Names of one or more variables present in the model object. The default of 'all' causes all variables to included.
#' @param responseType Character. Response type (scaling of model output) to be used.
#' @param ylimits Numeric vector. Limits to be used on the plot y-axis.
#' @param ylabel Character. Label to be used on the \emph{y}-axis replacing the default value formed from the prefix "Response" and the value of \emph{responseType}.
#' @param filename Character. Filename (with full path) into which graphics will be plotted.
#' @param occSWD Data.frame. Environmental data at occurrence locations in SWD format; used to compute variable importance.
#' @param bkgSWD Data.frame. Environmental data at background locations in SWD format; used to compute variable importance.
#'
#' @details {
#' This function is closely based on the source code for the corresponding function in the R package \pkg{maxnet} and draws inspiration from the response plots produced for Boosted Regression Tree (BRT) models produced by the package \pkg{gbm}.
#'
#' Features of this function include:
#' \itemize{
#' \item Returns a list of ggplot objects allowing individual plots to be extracted and used for other purposes
#' \item Following the approach in package \pkg{gbm}, incorporation of variable importance scores supplied by the companion function \link{varImportance}
#' \item Ability to produce a plot for \emph{all} variables in the maxnet model or a selection
#' }
#' }
#'
#' @return A named list of ggplot graphics objects ordered in decreasing variable importance, plus the side-effect of a pdf of plots.
#' @export
#'
#' @examples
responsePlot <- function(thisModel,
                         variable = "all",
                         responseType = "",
                         ylimits = c(0, 1),
                         ylabel = NULL,
                         filename = NULL,
                         occSWD = NULL,
                         bkgSWD = NULL)
{
  if (!("maxnet" %in% class(thisModel))) stop("Paramater 'thisModel' must class 'maxnet'")

  responseInd <- pmatch(responseType, c("link","exponential","cloglog","logistic"))

  if ((is.na(responseInd)) | (length(responseInd) > 1))
    stop("Cannot understand value passed in parameter 'responseType'")
  else
    responseType <- c("link","exponential","cloglog","logistic")[responseInd]

  if (variable == "all")
  {
    varList <- names(thisModel$samplemeans)
  }
  else
    varList <- variable

  if (!all(varList %in% names(thisModel$samplemeans)))
  {
    missingVars <- names(thisModel$samplemeans)[which(!(varList %in% names(thisModel$samplemeans)))]
    stop(paste0("The following variables names in paramater 'variable' are not in the maxnet model: ", paste(missingVars, collapse = ", ")))
  }

  cat("Producing a response plot for these variables:\n")
  cat(paste0(varList, collapse = ", "), "\n")

  # Compute variable importance scores
  varImp <- varImportance(thisModel, occSWD, bkgSWD, responseType = responseType)
  names(varImp) <- varList
  varImpOrder <- order(varImp, decreasing = TRUE)

  # Reorder varList so plots, etc are sorted in descending order of variable
  # importance or contribution to the model fit
  varList <- varList[varImpOrder]
  varImp <- varImp[varImpOrder]

  plotyBits <- vector("list", length(varList))
  names(plotyBits) <- varList

  for (thisVar in varList)
  {
    # For the current variable, make a mean value matrix and then replace the
    # variable column with a sequence of values spanning the range of values for
    # the variable encountered during model fitting
    hasLevels <- !is.null(unlist(thisModel$levels[thisVar]))

    if (hasLevels)
      numRows <- unlist(thisModel$levels[thisVar])
    else
      numRows <- 100

    d <- data.frame(matrix(unlist(thisModel$samplemeans), numRows, length(thisModel$samplemeans), byrow = TRUE))
    colnames(d) <- names(thisModel$samplemeans)

    varMax <- thisModel$varmax[thisVar]
    varMin <- thisModel$varmin[thisVar]

    if (hasLevels) d[, thisVar] <- unlist(thisModel$levels[thisVar])     else
      d[, thisVar] <- seq(varMin - 0.1 * (varMax - varMin), varMax + 0.1 * (varMax - varMin), length = 100)

    if (is.null(ylabel))
      ylabel <- paste0("Response (", responseType, ")")
    predVals <- predict(thisModel, d, type = responseType)

    plotData <- data.frame(prediction = predVals[,1], d)

    if (hasLevels)
    {
      plotyBits[[thisVar]] <- ggpubr::ggbarplot(plotData, x = thisVar, y = "prediction", palette = ggsci::pal_npg, xlab = thisVar, ylab = ylab, ylim = ylim)
    }
    else
    {
      plotyBits[[thisVar]] <- ggplot2::ggplot(plotData, aes(x = plotData[, thisVar], y = prediction)) +
        geom_line(colour = "blue", size = 1) + ylab(ylabel) + xlab(thisVar) + ylim(ylimits) +
        theme(axis.title.x = element_text(size = 8),
              axis.title.y = element_text(size = 8),
              axis.text.x = element_text(size = 6),
              axis.text.y = element_text(size = 6),
              plot.title = element_text(size = 8)) +
        ggtitle(paste0("Importance: ", varImp[thisVar], "%"))
    }
  }

  if (!is.null(filename))
  {
    ggpubr::ggarrange(plotlist = plotyBits, ncol = 4, nrow = 3) %>%
      ggpubr::ggexport(filename = filename)
  }

  return(plotyBits)
}
