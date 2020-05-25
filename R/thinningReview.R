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
#' The function applies a standard sequence of thinning ditances (in kilometers) to occupancy records and determines the maximum distance for thinning by examining the imapct of each prospective tinnning distance on coverage of environmental conditions.
#'
#' @param taxon Character. A taxon name to label output files and plots
#' @param occData Data frame of occurrence locations. Location coordinates are expected to be in decimal degrees. An attempt is made to identify the columns storing latitude and longitude values by partila matching on "longitude" and "latitude"
#' @param envDataPath Characater. Full path to the raster layers of environmental data to nbe used to acssess environmental coverage. Any file format accepted by the raster package can be used, and alll files in the set must have the same gridcell size, origin and resolution
#' @param outPath Character. A path to a folder into which out will be written
#' @param threshold Numeric. A value between 0 and 1 representing the fraction of environmental space which must be covered by the thinning operation
#' @param numReplicates Integer. How many repeated thinning runs should be performed at each thinning distance?
#' @param doPlots Logical. Should plots of results for each replicate and distance combination be produced in additional to a final summary plot?
#' @param writeResults Logical. Should a results summary table be saved?
#' @param quiet Logical. Should progress messages be written to the console
#'
#' @return The largest distance in the standard sequence for which environmental coverage is >= threshold
#' @export
#'
#' @examples
thinningReview <- function(taxon = "",
                           occData = NULL,
                           envDataPath = NULL,
                           outPath = "",
                           threshold = 0.9,
                           numReplicates = 5,
                           doPlots = FALSE,
                           writeResults = FALSE,
                           quiet = TRUE)
{
  # Automagically try to identify longitude and latitude columns:
  longColInd <- grep("LONG", toupper(colnames(occData)))
  latColInd <- grep("LAT", toupper(colnames(occData)))

  if (any(c(length(longColInd), length(latColInd)) != 1))
      stop("Cannot identify longitude and latitide cols in occData")

  envStack <- raster::stack(list.files(envDataPath, "*.tif", full.names = TRUE))

  if (!quiet) cat(taxon, "\n")
  envData_orig <- raster::extract(envStack, occData[, c(longColInd, latColInd)])

  badRows <- which(is.na(rowSums(envData_orig)))
  if (length(badRows) > 0)
  {
    envData_orig <- envData_orig[-badRows, ]
    occData <- occData[-badRows, ]
  }

  plotColours <- c("blue", "darkorange", "magenta1")

  basePCA <- stats::prcomp(envData_orig, center = TRUE, scale. = TRUE)
  hull_orig <- grDevices::chull(data.frame(basePCA$x[, c("PC1", "PC2")]))

  # Area of orig polygon in PC1-PC2 plane: may serve as an index of niche volume
  # Compute using an application of Green's Theorem
  pp <- rev(c(hull_orig, hull_orig[1]))

  pp_xy <- data.frame(basePCA$x[pp, c("PC1", "PC2")])
  ii <- 1:length(hull_orig)

  orig_area <- 0
  for (i in ii)
  {
    orig_area <- orig_area + (pp_xy$PC1[i + 1] + pp_xy$PC1[i]) * (pp_xy$PC2[i + 1] - pp_xy$PC2[i])/2
  }

  cellInd <- cellFromXY(envStack[[1]], occData[, c(longColInd, latColInd)])
  duplInd <- which(duplicated(cellInd))
  if (length(duplInd) > 0) cellInd <- cellInd[-duplInd]

  thinDistSet <- c(2, 5, 10, 15, 20, 25, 30)
  numDist <- length(thinDistSet)

  accumulResults <- NULL

  for (thisRepl in 1:numReplicates)
  {
    newResults <- data.frame(repl = rep(thisRepl, numDist),
                             thinningDist = thinDistSet,
                             origArea = rep(round(orig_area, 2), numDist),
                             thinnedArea = rep(0, numDist),
                             propArea = rep(0, numDist),
                             numOrig = rep(nrow(occData), numDist),
                             numThinned = rep(0, numDist),
                             percThinned = rep(0, numDist))

    #thisDist <- 1
    for (thisDist in thinDistSet)
    {
      rowInd <- which(thinDistSet == thisDist)
      occData_thin <- occThin(occData, longColInd, latColInd, thisDist)
      envData_thin <- raster::extract(envStack, occData_thin[, c(longColInd, latColInd)])

      thinPCA_proj <- predict(basePCA, envData_thin)

      hull_thin <- grDevices::chull(data.frame(thinPCA_proj[, c("PC1", "PC2")]))

      if (doPlots)
      {
        grDevices::png(paste0(outPath, "/", taxon, "_thining_distanceExporation_PCA_dist_", thisDist,"_repl_", thisRepl,".png"))
        p2 <- ggplot2::ggplot(data.frame(basePCA$x), aes(x = PC1, y = PC2)) +
          ggplot2::geom_point(data = data.frame(basePCA$x), shape = 16, colour = plotColours[1]) +
          ggplot2::geom_polygon(data = data.frame(basePCA$x[hull_orig, ]), aes(x = PC1, y = PC2), colour = plotColours[1], fill = plotColours[1], alpha = 0.2) +
          ggplot2::geom_point(data = data.frame(thinPCA_proj), shape = 16, colour = plotColours[2]) +
          ggplot2::geom_polygon(data = data.frame(thinPCA_proj[hull_thin, ]), aes(x = PC1, y = PC2), colour = plotColours[2], fill = plotColours[2], alpha = 0.2) +
          ggplot2::annotate("text", x = 3, y = 8, label = paste("Thinning distance =", thisDist,"km"))
        print(p2)
        dev.off()
      }

      # Compute area of the thinned polygon
      pp <- rev(c(hull_thin, hull_thin[1]))

      pp_xy <- data.frame(thinPCA_proj[pp, c("PC1", "PC2")])
      ii <- 1:length(hull_thin)

      thin_area <- 0
      for (i in ii)
      {
        thin_area <- thin_area + (pp_xy$PC1[i + 1] + pp_xy$PC1[i]) * (pp_xy$PC2[i + 1] - pp_xy$PC2[i])/2
      }

      if (!quiet)
      {
        cat("Thinning distance =", thisDist, "km\n")
        cat("Area original convex hull =", round(orig_area, 2), "\n")
        cat("Area of thinned convex hull =",round(thin_area, 2), "\n")
        cat("propArea =", round(thin_area/orig_area, 2), "\n")
        cat("numThinned =", nrow(occData_thin), "\n")
        cat("percThinned =", round(100*nrow(occData_thin)/nrow(occData), 2), "\n")
        cat("========================================================\n\n")
      }

      newResults[rowInd, "thinnedArea"] <- round(thin_area, 2)
      newResults[rowInd, "propArea"] <- round(thin_area/orig_area, 2)

      if (newResults[rowInd, "propArea"] >= threshold) bestDist <- thisDist

      newResults[rowInd, "numThinned"] <- nrow(occData_thin)
      newResults[rowInd, "percThinned"] <- round(100*nrow(occData_thin)/nrow(occData), 2)
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
                           percNicheVol = 100*accumulResults[, "propArea"],
                           percNumPoints = accumulResults[, "percThinned"])

    grDevices::png(paste0(outPath, "/", taxon, "_thining_distanceExporation_resultSummary.png"))
    p4 <- ggplot2::ggplot(plotData, aes(x = thinDist, y = percNicheVol)) +
      ggplot2::geom_hline(yintercept = 100, colour = "grey70", linetype = 3) +
      ggplot2::geom_point(colour = "slateblue1") +
      #geom_line(aes(x = thinDist, y = percNicheVol), data = plotData, colour = "slateblue1") +
      ggplot2::ylab("Percent") + xlab("Thinning distance (km)") +
      ggplot2::geom_point(aes(x = thinDist, y = percNumPoints), data = plotData, colour = "darkorange") +
      #geom_line(aes(x = thinDist, y = percNumPoints), data = plotData, colour = "darkorange") +
      #annotate("text", x = 20, y = 75, hjust = 0, label = "Niche volume") +
      #annotate("text", x = 10, y = 37.5, hjust = 0, label = "Number of points") +
      #geom_point(aes(x = 1, y = 100*thin_area_unique/orig_area), colour = "purple1") +
      ggplot2::geom_point(aes(x = 1, y = 100*length(cellInd)/nrow(envData_orig)), shape = 17, size = 2, colour = "magenta1") +
      ggplot2::ggtitle(taxon)
    print(p4)
    dev.off()
  }

  summaryTable <- accumulResults %>% dplyr::group_by(thinningDist) %>% dplyr::summarise(minPropArea = min(propArea))
  bestDist <- summaryTable$thinningDist[max(which(summaryTable$minPropArea >= threshold))]
  if (!quiet) cat("Best distance =", bestDist, "km\n")
  return(bestDist)
}



##### TEST DRIVER
# thisTaxon <- "Correa alba"
# bestDist <- thinningReview(thisTaxon,
#                            read.csv(paste0("/home/peterw/Data_and_Projects/RBG Projects/Restore and Renew/RandR-webtool-maintenance/New species staging/ENM fitting/", thisTaxon, "/", gsub(" ", "_", thisTaxon, fixed = TRUE), "_herbariumRecords_filtered_cleaned.csv")),
#                            "/home/peterw/Data_and_Projects/ENM_env_data/eastOZ-dataset/Current_climate/CHELSA",
#                            "/home/peterw/Data_and_Projects/Personal Projects/Spatial thinning/testing",
#                            doPlots = TRUE)


# thisTaxon <- "Banksia spinulosa"
# bestDist <- thinningReview(thisTaxon,
#                            read.csv("/home/peterw/Data_and_Projects/RBG Projects/Banksia modelling/B_spinulosa/Occurrence data from Trevor/28_4_2020_B_spinulosa_TOTAL.csv"),
#                            "/home/peterw/Data_and_Projects/ENM_env_data/eastOZ-dataset/Current_climate/CHELSA",
#                            "/home/peterw/Data_and_Projects/RBG Projects/Banksia modelling/B_spinulosa/thinningReview",
#                            doPlots = TRUE)



