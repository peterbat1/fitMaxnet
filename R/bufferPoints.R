
#' Make buffer around occurrence points
#'
#' Produce a spatial polygon object representing the boundary within which background points will be sampled for ENM fitting
#'
#' @param occ_pts A data frame or matrix with at least two columns representing longitude and latitude (X and Y) of occurrence records
#' @param bufferDist_km Size of buffer to be built around the occurrence points in kilometres. Default value of 200 km follows the advice of VanDerWal et al. (2009). See details below.
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
#' @return A spatial polygon object from the package \emph{sf}
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

  # Try to identify longitude (X) and latitude (Y) columns
  ind <- grep("LONG|X", toupper(colnames(occ_pts)))
  if (length(ind) >= 1)
    X_ind <- ind[1] #colnames(occ_pts)[ind[1]] <- "longitude"
  else
    stop("Cannot identify the 'longitude' or 'X' column in 'occ_pts'")

  ind <- grep("LAT|Y", toupper(colnames(occ_pts)))
  if (length(ind) >= 1)
    Y_ind <- ind[1] #colnames(occ_pts)[ind[1]] <- "latitude"
  else
    stop("Cannot identify the 'latitude' or 'Y' column in 'occ_pts'")

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

  # Clip the buffer polygon to the Australian coastline
  clippedBuffer <- sf::st_union(sf::st_intersection(ptsBuffer, sf::st_transform(ozPolygon, 3577)))

  if (trace) cat("step 4\n")
  if (pts_crs == 4326) clippedBuffer <- sf::st_transform(clippedBuffer, crs = 4326)
  return(clippedBuffer)
}
