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
NumericMatrix doDilation2(int npts, NumericMatrix d, int nRow, int nCol,
                                       double yul, double xll, double cellSize,
                                       NumericVector ptRows, NumericVector ptCols,
                                       double dist_km)
{
  int i, r, c, leftcol, rightcol, toprow, botrow;
  double endLat, endLong, lat2, long2, GCD;
  double Earth_radius = 6371.1;
  double degtorad = M_PI/180;
  double radtodeg = 180/M_PI;
  double lat1[npts];
  double long1[npts];
  double yulCentre = yul - cellSize/2;
  double xllCentre = xll + cellSize/2;
  double rowCentreLat[nRow];
  double colCentreLong[nCol];
  //double endPt[2];
  double dR = dist_km/Earth_radius;

  for (i = 0; i < npts; i++)
  {
    lat1[i] = degtorad*(yulCentre - (ptRows[i] - 1)*cellSize);
    long1[i] = degtorad*(xllCentre + (ptCols[i] - 1)*cellSize);
  }

  for (i = 0; i < nRow; i++)
  {
    rowCentreLat[i] = degtorad*(yulCentre - i*cellSize);
  }

  for (i = 0; i < nCol; i++)
  {
    colCentreLong[i] = degtorad*(xllCentre + i*cellSize);
  }

  for (i = 0; i < npts; i++)
  {
    endLat =  asin(sin(lat1[i])*cos(dR) + cos(lat1[i])*sin(dR)*cos(3*M_PI/2));
    endLong = long1[i] + atan2(sin(3*M_PI/2)*sin(dR)*cos(lat1[i]),cos(dR)-sin(lat1[i])*sin(endLat));
    leftcol = trunc((radtodeg*endLong - xllCentre)/cellSize);
    if (leftcol < 0) { leftcol = 0; }
    if (leftcol > nCol - 1) { leftcol = nCol - 1; }

    endLat =  asin(sin(lat1[i])*cos(dR) + cos(lat1[i])*sin(dR)*cos(M_PI/2));
    endLong = long1[i] + atan2(sin(M_PI/2)*sin(dR)*cos(lat1[i]),cos(dR)-sin(lat1[i])*sin(endLat));
    rightcol = trunc((radtodeg*endLong - xllCentre)/cellSize);
    if (rightcol < 0) { rightcol = 0; }
    if (rightcol > nCol - 1) { rightcol = nCol - 1; }

    endLat =  asin(sin(lat1[i])*cos(dR) + cos(lat1[i])*sin(dR)*cos(0.0));
    toprow = trunc((yulCentre - radtodeg*endLat)/cellSize);
    if (toprow < 0) { toprow = 0; }
    if (toprow > nRow) { toprow = nRow - 1; }

    endLat =  asin(sin(lat1[i])*cos(dR) + cos(lat1[i])*sin(dR)*cos(M_PI));
    botrow = trunc((yulCentre - radtodeg*endLat)/cellSize);
    if (botrow < 0) { botrow = 0; }
    if (botrow > nRow) { botrow = nRow - 1; }

    // toprow, botrow, leftcol, rightcol are assumed to already be zero-adjusted
    for (r = toprow; r < botrow; r++)
    {
      for (c = leftcol; c < rightcol; c++)
      {
        //if ((d(r,c) != 1) && (!(d(r,c) == NA_REAL)))
        //{
        lat2 = rowCentreLat[r];
        long2 = colCentreLong[c];
        GCD = Earth_radius*acos(sin(lat1[i])*sin(lat2)+cos(lat1[i])*cos(lat2)*cos(long2-long1[i]));
        if (GCD <= dist_km)
        {
          d(r,c) = 1;
        }
        //}
      }
    }
  }
  return d;
}

