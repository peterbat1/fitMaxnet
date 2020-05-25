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
bufferPoints <- function(occ_pts, bufferDist_km = 200)
{
  # Convert buffer distance from kilometres to metres
  bufferDist <- 1000 * bufferDist_km

  # Try to identify longitude and latitude columns and rename them to 'longitude' and 'latitude'
  ind <- grep("LONG", toupper(colnames(occ_pts)))
  if (length(ind) == 1)
    colnames(occ_pts)[ind] <- "longitude"
  else
    stop("Cannot identify the 'longitude' column in 'occ_pts'")

  ind <- grep("LAT", toupper(colnames(occ_pts)))
  if (length(ind) == 1)
    colnames(occ_pts)[ind] <- "latitude"
  else
    stop("Cannot identify the 'latitude' column in 'occ_pts'")

  sp::coordinates(occ_pts) <- ~longitude + latitude

  sp::proj4string(occ_pts) <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

  occ_pts <- sp::spTransform(occ_pts, "+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs ")

  ptsBuffer <- rgeos::gBuffer(occ_pts, width = bufferDist)

  #st <- Sys.time()
  clippedBuffer <- sf::st_intersection(st_geometry(st_as_sf(ptsBuffer)), ozPolygon)
  #print(Sys.time() - st)
  clippedBuffer <- sp:: spTransform(sf::as_Spatial(clippedBuffer), sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
  return(clippedBuffer)
}

