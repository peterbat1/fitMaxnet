# Produce a spatial polygon object representing the background constraint
# boundary for ENM fitting
#
# Peter D. Wilson
# Biodiversity Analyst
# Evolutionary Ecology Research Section
# Science and Conservation Branch
# Royal Botanic Garden, Sydney
#
# 2019-09-26

# library(sp)
# library(sf)
# library(rgeos)
# library(rgdal)


#' Make buffer around occurrence points
#'
#' Produce a spatial polygon object representing the boundary within which background points will be sampled for ENM fitting
#'
#' @param occ_pts A data frame or matrix with at least two columns representing latitude and longitude of occurrence records
#' @param bufferDist_km Size of buffer to be built around the occurrence points in kilometres. Default value of 200 km follows the advice of VanDerWal et al. (2009). See details below.
#'
#' @details {
#' \emph{bufferPoints()} implements the widely used simple circular buffering approach introduced to niche modelling by VanDerWal et al. 2009. 'Selecting pseudo-absence data for presence-only distribution modeling: How far should you stray from what you know?' \emph{Ecological Modelling} 220:589–594.
#' The default value of 200 km for the buffer distance is the value recommended by VanDerWal et al. (2009) for a range of taxa which they considered in their paper. Naturally, the value to be used is a matter for consideration by the user. If there is a clear indication from available information regarding a taxon's mode of dispersal, or there is guidance from population genetics data, then another value might be supplied for this parameter.
#' This implementation uses the spatial processing resources provided by the packages \emph{sp}, \emph{sf} and \emph{rgeos} to achieve good performance and accurate handling of spatial data.
#' }
#'@note {
#' Currently, this function can only process \bold{occurrence data from the Australian continent}. Providing global coverage is feasible and will be implemented shortly for the OzWeeds project.
#'
#' The function only outputs a spatial polygon object projected on the WGS84 datum. If you need to use the buffer polygon on other projections, please reproject using, say, \emph{spTransform()} from the package \emph{sp}.
#'}
#' @return A spatial polygon object from the package \emph{sp}
#' @export
#'
#' @examples
#' \dontrun{}
bufferPoints <- function(occ_pts, bufferDist_km = 200, trace = FALSE)
{
  # Convert buffer distance from kilometres to metres
  bufferDist <- 1000 * bufferDist_km

  if (any(grepl("LONG|LAT", toupper(colnames(occ_pts)))))
    pts_crs <- 4326
  else
    pts_crs <- 3577

  # Try to identify longitude and latitude columns and rename them to 'longitude' and 'latitude'
  ind <- grep("LONG|X", toupper(colnames(occ_pts)))
  if (length(ind) >= 1)
    X_ind <- ind #colnames(occ_pts)[ind[1]] <- "longitude"
  else
    stop("Cannot identify the 'longitude' column in 'occ_pts'")

  ind <- grep("LAT|Y", toupper(colnames(occ_pts)))
  if (length(ind) >= 1)
    Y_ind <- ind #colnames(occ_pts)[ind[1]] <- "latitude"
  else
    stop("Cannot identify the 'latitude' column in 'occ_pts'")

  occ_pts_sf <- sf::st_as_sf(occ_pts, coords = c(X_ind, Y_ind), crs = pts_crs)
  #sf::st_crs(occ_pts_sf) <- pts_crs

  if (trace) cat("step 1\n")
  if (pts_crs == 4326)
    occ_pts_albers <- sf::st_transform(occ_pts_sf, crs = 3577)
  else
    occ_pts_albers <- occ_pts_sf

  if (trace) cat("step 2\n")
  ptsBuffer <- sf::st_union(sf::st_buffer(occ_pts_albers, dist = bufferDist))

  if (trace)
  {
    cat("step 3\n")
    cat("\nCRS(ptsBuffer:\n")
    print(st_crs(ptsBuffer))

    cat("\nCRS(ozPolygon):\n")
    print(st_crs(ozPolygon))

    cat("\nst_crs(ptsBuffer) == st_crs(ozPolygon) is", st_crs(ptsBuffer) == st_crs(ozPolygon), "\n")
  }

  clippedBuffer <- sf::st_intersection(ptsBuffer, sf::st_transform(ozPolygon, 3577))

  if (trace) cat("step 4\n")
  if (pts_crs == 4326) clippedBuffer <- sf::st_transform(clippedBuffer, crs = 4326)
  return(clippedBuffer)
}
