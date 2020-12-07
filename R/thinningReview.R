# Explore changes in coverage of individual env vars with application of
# thinning: can this help explain odd behaviour for some test taxa such as
# Hardenbergia violacea
#
#
# Peter D. Wilson
#
# 2019-12-10

#' Review impact of a sequence of thinning distances
#'
#' The function applies a standard sequence of thinning distances (in kilometers) to occupancy records and determines the maximum distance for thinning by examining the impact of each prospective thinning distance on coverage of environmental conditions.
#'
#' @param taxon Character. A taxon name to label output files and plots
#' @param occData Data frame of occurrence locations, or a character object (string) giving the path to this file. An attempt is made to identify the columns storing longitude (or X) and latitude (or Y) values by partial matching on "longitude" and "latitude"
#' @param envDataPath Character. Full path to the raster layers of environmental data to be used to assess environmental coverage. Any file format accepted by the raster package can be used, and all files in the set must have the same gridcell size, origin and resolution
#' @param outPath Character. A path to a folder into which out will be written
#' @param thinDistSet Numeric vector. Sequence of thinning distances (in km) to be processed
#' @param pThreshold Numeric. Threshold p-value used to determine optimal thinning distance.
#' @param numReplicates Integer. How many repeated thinning runs should be performed at each thinning distance?
#' @param isLatLong Logical. Are geographic coordinates (longitude/latitude) used in occData and environmental data layers?
#' @param doPlots Logical. Should a plot of results be written to the folder named in 'outPath'?
#' @param writeResults Logical. Should a results summary table be saved to the folder named in 'outPath'?
#' @param quiet Logical. Should progress messages be written to the console
#'
#' @details {
#' The sequence of thinning distances (in kilometres) should be chosen in relation to the grid cell size of the environmental data layers.
#'
#' An approach which seems to be useful is to start with half the width of a grid cell and set a sequence of increments which span successive sets of grid cells.
#'
#' \emph{pThreshold} default of 0.9 was chosen by trails which suggested that fitted ENMs began to show poor performance when thinning caused greater loss of environmental coverage.
#' }
#'
#' @return The largest distance in the standard sequence for which environmental coverage is the closest to the threshold approached from above (i.e. is equal to or just greater than the \emph{pThreshold})
#' @export
#'
#' @examples
#' \dontrun{}
thinningReview <- function(taxon = "",
                           occData = NULL,
                           envDataPath = NULL,
                           outPath = "",
                           thinDistSet = c(0.5, 1, 2, 3, 4, 5, 6, 7, 10),
                           pThreshold = 0.9,
                           numReplicates = 5,
                           isLatLong = TRUE,
                           doPlots = FALSE,
                           writeResults = FALSE,
                           quiet = TRUE)
{
  #plotColours <- c("blue", "darkorange", "magenta1")

  if (taxon == "") stop("Please supply a taxon name in parameter 'taxon'")

  if (!quiet) cat(taxon, "\n")

  if (is.null(occData))
    stop("Please supply a data.frame or path to csv file containing occurrence data in parameter 'occData'")

  if (is.data.frame(occData))
  {
    theseOccData <- occData
  }
  else
  {
    if (is.character(occData))
    {
      if (file.exists(occData))
        theseOccData <- read.csv(occData, stringsAsFactors = FALSE)
      else
        stop("Expected a data.frame or path to a csv file in parameter 'occData' but got neither")
    }
    else
    {
      if (is.matrix(occData))
      {
        theseOccData <- data.frame(occData)
      }
      else
        stop("Cannot convert 'occData' to a data.frame")
    }
  }

  if (is.null(envDataPath))
    stop("Please supply a path to a folder of environmental data layers in parameter 'envDataPath'")

  if (!dir.exists(envDataPath))
    stop("Cannot find environmental data folder passed in parameter 'envDataPath'")

  if (outPath == "")
    stop("Please supply a path to an output folder in parameter 'outPath'")
  else
    if (!dir.exists(outPath)) dir.create(outPath)

  if (!dir.exists(outPath)) dir.create(outPath)

  # Automagically try to identify longitude and latitude columns:
  xColInd <- grep("LONG|X", toupper(colnames(theseOccData)))
  if (length(xColInd) == 0) stop("Cannot identify a 'longitude' or 'x' column in occurrence data file")
  if (length(xColInd) > 1)
  {
    warning("Cannot identify a unique 'longitude' or 'x' column in occurrence data file; using first hit come-what-may...")
    xColInd <- xColInd[1]
  }

  yColInd <- grep("LAT|Y", toupper(colnames(theseOccData)))
  if (length(yColInd) == 0) stop("Cannot identify a 'latitude' or 'y' column in occurrence data file")
  if (length(yColInd) > 1)
  {
    warning("Cannot identify a unique 'latitude' or 'y' column in occurrence data file; using first hit come-what-may...")
    yColInd <- yColInd[1]
  }

  # Set this flag which is needed for the call to occThin()
  if (any(grepl("LONG|LAT", toupper(colnames(theseOccData)))))
    isLatLong <- TRUE
  else
    isLatLong <- FALSE

  if (!quiet) cat("  loading env data\n")
  envStack <- raster::stack(list.files(envDataPath, "*.tif", full.names = TRUE))

  if (!quiet) cat("  checking for duplicated occupied cells and removing any that are found\n")
  duplInd <- which(duplicated(cellFromXY(envStack[[1]], theseOccData[, c(xColInd, yColInd)])))
  if (length(duplInd) > 0)
  {
    if (!quiet) cat("    (removed", length(duplInd), "duplicates)\n")
    theseOccData <- theseOccData[-duplInd, ]
  }

  if (!quiet) cat("  extracting env data at occ locations\n")
  envData_orig <- raster::extract(envStack, theseOccData[, c(xColInd, yColInd)])

  badRows <- which(is.na(rowSums(envData_orig)))
  if (length(badRows) > 0)
  {
    if (!quiet) cat("  removed", length(badRows), "rows from occData and environmental data matrix with missing environmental data\n")
    envData_orig <- envData_orig[-badRows, ]
    theseOccData <- theseOccData[-badRows, ]
  }


  numDist <- length(thinDistSet)

  accumulResults <- NULL

  if (!quiet) cat("  start replicate sampling along thinning distance sequence\n")
  for (thisRepl in 1:numReplicates)
  {
    newResults <- data.frame(repl = rep(thisRepl, numDist),
                             thinningDist = thinDistSet,
                             energy_statistic = rep(0, numDist),
                             energy_p_value = rep(0, numDist),
                             numOrig = rep(nrow(theseOccData), numDist),
                             numThinned = rep(0, numDist),
                             percRetained = rep(0, numDist))

    for (thisDist in thinDistSet)
    {
      rowInd <- which(thinDistSet == thisDist)

      cat(" >>>> call occThin\n")
      theseOccData_thin <- occThin(theseOccData, xColInd, yColInd, thisDist, isLatLong, quiet)
      envData_thin <- raster::extract(envStack, theseOccData_thin[, c(xColInd, yColInd)])

      cat(" >>>> call energyStats\n")
      ans <- energyStats(envData_orig, envData_thin)

      if (!quiet)
      {
        cat("Thinning distance =", thisDist, "km\n")
        cat("Energy statistic =", ans$statistic, "\n")
        cat("Energy p-value =", round(ans$p.value, 3), "\n")
        cat("Number removed =", nrow(theseOccData) - nrow(theseOccData_thin), "\n")
        cat("Fraction retained =", round(nrow(theseOccData_thin)/nrow(theseOccData), 2), "\n")
        cat("========================================================\n\n")
      }

      newResults[rowInd, "energy_statistic"] <- ans$statistic
      newResults[rowInd, "energy_p_value"] <- ans$p.value
      newResults[rowInd, "numRemoved"] <- nrow(theseOccData) - nrow(theseOccData_thin)
      newResults[rowInd, "fracRetained"] <- round(nrow(theseOccData_thin)/nrow(theseOccData), 2)
    }

    accumulResults <- rbind(accumulResults, newResults)
  }

  if (writeResults)
  {
    write.csv(accumulResults, paste0(outPath, "/", taxon, "_spatialThinningExploration.csv"), row.names = FALSE)
  }

  if (doPlots)
  {
    plotData <- data.frame(thinDist = accumulResults[, "thinningDist"],
                           probability = accumulResults[, "energy_p_value"],
                           fracRetained = accumulResults[, "fracRetained"])

    grDevices::png(paste0(outPath, "/", taxon, "_thining_distanceExporation_resultSummary.png"))
    p <- ggplot2::ggplot(plotData, aes(x = thinDist, y = probability)) +
      ggplot2::geom_hline(yintercept = pThreshold, colour = "grey70", linetype = 3) +
      ggplot2::ylim(0, 1) +
      ggplot2::geom_point(colour = "slateblue1") +
      ggplot2::ylab("Probability & Fraction retained") + xlab("Thinning distance (km)") +
      ggplot2::geom_point(aes(x = thinDist, y = fracRetained), data = plotData, colour = "darkorange") +
      ggplot2::ggtitle(taxon) +
      theme(plot.title = element_text(face = "bold.italic", size = 16))
    print(p)
    dev.off()
  }

  # Produce a summary table and locate the best thinning distance using the default method
  summaryTable <- accumulResults %>% dplyr::group_by(thinningDist) %>% dplyr::summarise(mean_p_value = mean(p.value))
  xx <- summaryTable$mean_p_value - pThreshold
  ii <- which(xx == min(xx[xx > 0]))
  bestDist = summaryTable$thinningDist[ii]

  if (!quiet)
  {
    cat("\n-----------------------------------------------------\n")
    cat("  Best distance =", bestDist, "km\n")
    cat("\n  Results table:\n")
    print(summaryTable)
  }

  return(bestDist)
}


