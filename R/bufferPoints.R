
#' Make buffer around occurrence points
#'
#' Produce an sf polygon object representing the boundary within which background points will be sampled for ENM fitting
#'
#' @param occ_pts A data frame or matrix with at least two columns representing longitude and latitude (X and Y) of occurrence records
#' @param bufferDist_km Numeric. Size of buffer to be built around the occurrence points in kilometres. Default value of 200 km follows the advice of VanDerWal et al. (2009). See details below.
#' @param trace Logical. Should messages be emitted to assist progress tracking or debugging
#'
#' @details {
#' \emph{bufferPoints()} implements the widely used simple circular buffering approach introduced to niche modelling by VanDerWal et al. 2009. 'Selecting pseudo-absence data for presence-only distribution modeling: How far should you stray from what you know?' \emph{Ecological Modelling} 220:589–594.
#' The default value of 200 km for the buffer distance is the value recommended by VanDerWal et al. (2009) for a range of taxa which they considered in their paper. Naturally, the value to be used is a matter for consideration by the user. If there is a clear indication from available information regarding a taxon's mode of dispersal, or there is guidance from population genetics data, then another value might be supplied for this parameter.
#' This implementation uses the spatial processing resources provided by the packages \emph{sp}, \emph{sf} and \emph{rgeos} to achieve good performance and accurate handling of spatial data.
#' }
#'@note {
#' Currently, this function can only process \bold{occurrence data from the Australian continent}. Providing global coverage is feasible and will be implemented shortly for the OzWeeds project.
#'
#' The function outputs a simple features ('sf') polygon object on the same projection which was inferred when \emph{occ_pts} was parsed. That is, either EPSG:4326 (WGS84) if 'longitude' and 'latitude' column names were found, or EPSG:3577 (Australian ALber's Equal Area) if 'X' and 'Y' column names were found. If you need to use the buffer polygon on other projections, please re-project the returned object using, say, \emph{st_transform()} from the package \emph{sf}.
#'}
#' @return An sf polygon object from the package \emph{sf}
#' @export
#'
#' @examples
#' \dontrun{}
bufferPoints <- function(occ_pts, bufferDist_km = 200, trace = FALSE)
{
  # Convert buffer distance from kilometres to metres
  bufferDist <- 1000 * bufferDist_km

  if ("sf" %in% class(occ_pts))
  {
    pts_crs <- gsub("EPSG:", "", st_crs(occ_pts)$input, fixed = TRUE)
    if (pts_crs != 3577)
      occ_pts_albers <- st_transform(occ_pts, 3577)
    else
      occ_pts_albers <- occ_pts
  }
  else
  {
    stop("occ_pts must be class 'sf'")
  }


  if (trace) cat("make buffer polygon\n")
  ptsBuffer <- sf::st_union(sf::st_buffer(occ_pts_albers, dist = bufferDist))

  # Clip the buffer polygon to the Australian coastline
  if (trace) cat("clip buffer polygon to OZ coastline\n")
  clippedBuffer <- st_sf(sf::st_union(sf::st_intersection(ptsBuffer, sf::st_transform(ozPolygon, 3577))))


  if (trace) cat("transform clipped buffer to CRS of occ_pts\n")
  if (pts_crs != 3577) clippedBuffer <- sf::st_transform(clippedBuffer, crs = pts_crs)
  return(clippedBuffer)
}
