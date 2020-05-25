#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]

NumericVector endLatLongC(double dist_km, double bearing, double lat1, double long1, bool isRadians)
{
  double lat2, long2;
  double degtorad = M_PI/180;
  double radtodeg = 180/M_PI;
  double Earth_radius = 6371.1;
  double dR = dist_km/Earth_radius;

  if (!isRadians)
  {
    lat1 = degtorad*lat1;
    long1 = degtorad*long1;
    bearing = degtorad*bearing;
  }

  NumericVector ans(2);
  dR = dist_km/Earth_radius;
  // lat2
  ans(0) = radtodeg*(asin(sin(lat1)*cos(dR) + cos(lat1)*sin(dR)*cos(bearing)));
  // long2
  ans(1) = radtodeg*(long1 + atan2(sin(bearing)*sin(dR)*cos(lat1),cos(dR)-sin(lat1)*sin(lat2)));

  return ans;
}
