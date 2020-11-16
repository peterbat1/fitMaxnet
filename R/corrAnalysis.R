
#' Exploratory Correlation Analysis of Environmental Data
#'
#' Explore correlation relationships between environmental predictor values at occurrence locations
#'
#' @param taxon Character. The name of the taxon whose occurrence records are being analysed.
#' @param titleText Character. A title to be used in graphical output.
#' @param envDataPath Character. Path to the environmental data layers to be  used in the analysis.
#' @param occData Data.frame or matrix. At least two columns must be present to provide longitude/X and latitude/Y coordinates of occurrence locations.
#' @param xVar Character. Name of a variable in \emph{occData} which is interpreted as the x-coordinate. If NULL (default) a search is made for nearest match to 'longitude' or 'X'.
#' @param yVar Character. Name of a variable in \emph{occData} which is interpreted as the y-coordinate. If NULL (default) a search is made for nearest match to 'latitude' or 'Y'.
#' @param threshold Numeric. A correlation value (ie between 0 and 1) used to determine which variables in \emph{envData} will be recommended for removal. Correlations greater than or equal to \emph{threshold} will be listed.
#' @param outFile Character. A non-NULL value is used as a file name to save the graphical output as a PNG file. By default, the output is plotted to the default graphics device.
#' @param outPath Character. Path used by \link[ggplot2]{ggsave} in combination with \emph{outFile} to save the plot.
#' @return A character matrix listing the names of variables with absolute value of correlations greater than \emph{threshold} which may be candidates for removal, and the number of threshold-exceeding correlations in which a listed variable has been found.
#' @export
#'
#' @examples
#' \dontrun{}
envCorrAnalysis <- function(taxon = "",
                            titleText = NULL,
                            envDataPath,
                            occData,
                            xVar = NULL,
                            yVar = NULL,
                            threshold = 0.7,
                            outFile = NULL,
                            outPath = NULL)
{
  if (!is.null(xVar))
  {
    if (xVar %in% colnames(occData))
      xCoords <- occData[, xVar]
    else
      stop(paste0("Named x-coord variable, ", xVar, " not found in occData"))
  }
  else
  {
    # Try to find a "longitude" match in colnames of occData
    xInd <- grep("LONG|X", toupper(colnames(occData)))
    #print(xInd)
    if (length(xInd) == 0)
      stop("No longitude/X variable found in occData")
    else
    {
      if (length(xInd) > 1)
        stop("Multiple options for longitude/X found in occData: please select one")
      else
        xCoords <- occData[, xInd]
    }
  }


  if (!is.null(yVar))
  {
    if (yVar %in% colnames(occData))
      yCoords <- occData[, yVar]
    else
      stop(paste0("Named y-coord variable, ", yVar, " not found in occData"))
  }
  else
  {
    # Try to find a "longitude" match in colnames of occData
    yInd <- grep("LAT|Y", toupper(colnames(occData)))
    #print(yInd)
    if (length(yInd) == 0)
      stop("No latitude/Y variable found in occData")
    else
    {
      if (length(yInd) > 1)
        stop("Multiple options for latitude/Y found in occData: please select one")
      else
        yCoords <- occData[, yInd]
    }
  }

  rasStack <- raster::stack(list.files(envDataPath, "*.*", full.names = TRUE))
  envData <- raster::extract(rasStack, cbind(x = xCoords, y = yCoords))

  badInd <- which(is.na(envData), arr.ind = TRUE)
  if (nrow(badInd) > 0) envData <- envData[-badInd[, 1], ]
  corr <- stats::cor(envData)

  ind1 <- which(lower.tri(corr), arr.ind = TRUE)
  ind2 <- which(abs(corr[ind1]) >= threshold, arr.ind = TRUE)

  ind3 <- ind1[ind2, ]

  if (length(ind2) > 0)
  {
    for (i in 1:nrow(ind3))
    {
      ind3[i,] <- ind3[i, order(ind3[i, ])]
    }

    ans <- as.matrix(table(colnames(envData)[ind3[,1]]))
    colnames(ans) <- "Frequency"
  }
  else
    ans <- NULL

  if (!is.null(outFile))
    if (dir.exists(outPath))
    {
      ggplot2::ggsave(ggcorrplot::ggcorrplot(corr, hc.order = TRUE, method = "circle"),
                      file = outFile,
                      device = "png",
                      width = 100,
                      height = 100,
                      units = "mm")
    }

  return(ans)
}
