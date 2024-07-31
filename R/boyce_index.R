# fitMaxnet: An adaptation of the function ecospat.boyce in the package ecospat
# which computes the continuous form of Boyce's index as developed by Hirzel et
# al. (2006)


#### boyce index

#### functions calculating Boyce index (Hirzel et al. 2006) by Blaise Petitpierre & Frank Breiner(28.06.2013)
# fit: A vector or Raster-Layer containing the predicted suitability values
#obs: A vector containing the predicted suitability values or xy-coordinates (if fit is a Raster-Layer) of the validation points (presence records)

#### internal function calculating predicted-to-expected ratio for each class-interval
# boycei <- function(interval, obs, fit) {
#   tic("boycei call")
#   fit.bin <- fit
#   obs.bin <- obs
#   fit.bin[fit[] >= interval[1] & fit[] <= interval[2]] <- "i"
#   fit.bin[fit.bin != "i"] <- 0
#   ecospat_fit.bin <<- fit.bin
#   obs.bin[obs[] >= interval[1] & obs[] <= interval[2]] <- "i"
#   obs.bin[obs.bin != "i"] <- 0
#   print(obs.bin)
#   pi <- length(which(obs.bin == "i"))/length(obs)
#   cat("pi =", pi, "\n")
#
#   ei <- length(which(fit.bin == "i"))/length(fit.bin)
#   ecospat_ei <<- ei
#   cat("length(which(fit.bin == 'i')) =", length(which(fit.bin == "i")), "\n")
#   cat("ei =", ei, "\n")
#   cat("length(fit.bin) =", length(fit.bin), "\n")
#   fi <- pi/ei
#   cat("fi =", fi, "\n")
#   #print(length(fi))
#   toc()
#   return(fi)
# }

fast_boycei <- function(interval, obs, fit)
{
  pi <- sum(.bincode(obs, interval, FALSE), na.rm = TRUE)/length(obs)
  #cat("pi =", pi, "\n")
  ei <- sum(.bincode(fit, interval, FALSE), na.rm = TRUE)/length(fit)
  #cat("ei =", ei, "\n")
  fi <- pi/ei
  return(fi)
}


#### Calculating Boyce index as in Hirzel et al. 2006
# fit: A vector or Raster-Layer containing the predicted suitability values
# obs: A vector containing the predicted suitability values or xy-coordinates (if fit is a Raster-Layer) of the validation points (presence records)
# nclass : number of classes or vector with classes threshold. If nclass=0, Boyce index is calculated with a moving window (see next parameters)
# windows.w : width of the moving window (by default 1/10 of the suitability range)
# res : resolution of the moving window (by default 100 focals)
# PEplot : if True, plot the predicted to expected ratio along the suitability class




#' Boyce's continuous index
#'
#' @param fit A numeric vector or RasterLayer or SpatRaster containing the predicted suitability values
#' @param obs A numeric vector containing the predicted suitability values or xy-coordinates (if fit is a Raster-Layer) of the validation points (presence records)
#' @param nclass Numeric. Number of classes or vector with classes threshold. If nclass=0, Boyce index is calculated with a moving window (see next parameters)
#' @param window.w Numeric. Width of the moving window (by default 1/10 of the suitability range)
#' @param res Numeric. Resolution of the moving window (by default 100 focals)
#' @param PEplot Logical. If TRUE, plot the predicted to expected ratio along the suitability class
#'
#' @details {
#' This is a daft bit of text.
#'
#' Hirzel et al. 2006. Evaluating the ability of habitat suitability models to predict species presences. Ecological Modelling 199:142-152.
#' }
#'
#' @return A named list...
#' @export
#'
#' @examples
#' \dontrun{}
boyce <- function(fit, obs, nclass = 0, window.w = "default", res = 100, PEplot = FALSE)
{

  if (!inherits(fit, c("RasterLayer", "SpatRaster")))
    stop("removeDuplicates: byGrid = TRUE so parameter baseGrid must of class RasterLayer or class SpatRaster.")
  else
  {
    # Convert to class terra::SpatRaster
    if (inherits(fit, "RasterLayer")) baseRaster <- terra::rast(fit)

    if (inherits(obs, "data.frame") | inherits(obs, "matrix"))
    {
      obs <- terra::extract(fit, obs)
    }

    fit <- fit[]
    fit <- fit[!is.na(fit)]
  }

  if (window.w == "default")
  {
    window.w <- (max(fit) - min(fit))/10
  }

  interval <- range(fit)
  mini <- interval[1]
  maxi <- interval[2]

  if (nclass == 0)
  {
    vec.mov <- seq(from = mini, to = maxi - window.w, by = (maxi - mini - window.w)/res)
    vec.mov[res + 1] <- vec.mov[res + 1] + 1  #Trick to avoid error with closed interval in R
    interval <- cbind(vec.mov, vec.mov + window.w)
  }
  else
    if (length(nclass) > 1)
    {
      vec.mov <- c(mini, nclass)
      interval <- cbind(vec.mov, c(vec.mov[-1], maxi))
    }
  else
    if (nclass > 0 & length(nclass) < 2)
    {
      vec.mov <- seq(from = mini, to = maxi, by = (maxi - mini)/nclass)
    }

  f <- apply(interval, 1, fast_boycei, obs, fit)

  to.keep <- which(f != "NaN")  # index to keep no NaN data
  f <- f[to.keep]

  if (length(f) < 2)
  {
    b <- NA  #at least two points are necessary to draw a correlation
  }
  else
  {
    r <- c(1:length(f))[f != c(f[-1], FALSE)]  #index to remove successive duplicates
    b <- cor(f[r], vec.mov[to.keep][r], method = "spearman")  # calculation of the spearman correlation (i.e. Boyce index) after removing successive duplicated values
  }

  HS <- apply(interval, 1, sum)/2  # mean habitat suitability in the moving window
  HS[length(HS)] <- HS[length(HS)] - 1  #Correction of the 'trick' to deal with closed interval
  HS <- HS[to.keep]  #exlude the NaN

  if (PEplot)
  {
    plot(HS, f, xlab = "Habitat suitability", ylab = "Predicted/Expected ratio", col = "grey", cex = 0.75)
    points(HS[r], f[r], pch = 19, cex = 0.75)
  }
  #results <- list(F.ratio = f, Spearman.cor = round(b, 3), HS = HS)

  return(list(F.ratio = f, Spearman.cor = round(b, 3), HS = HS))
}
