#include <stdio.h>
#include <math.h>
#include "convcoord.h"

#define REQ 6378137.0  /* meters */
//#define EPS (1.0/298.257223)
// #define ECC1 0.9966471893289236e+00 // =(1.0-EPS)
// #define ECC2 0.6694380002756624e-02 // =EPS*(2.0-EPS)

// now based on 298.257223563 (WGS84 Wikipedia 2013-06-07)
#define ECC1 0.99664718933525251928015447138148
#define ECC2 0.00669437999014131699613723354004
#define CB(X) (X*X*X)

#define J2000 51544.5  /* J2000.0 in mean julian days */
#define A 18.69737458
#define B 24.06570982441908

/*
 *  Computes the matrix that transforms a vector in GEI coordinates (RA/Dec)
 *    into local terrestrial coordinates.  This is the matrix product:
 *
 *            | 1  0     0   |   | C(h) 0 -S(h) |   | 0  1  0 |
 *            | 0 C(l) -S(l) | x |  0   1   0   | x | 0  0  1 |
 *            | 0 S(l)  C(l) |   | S(h) 0  C(h) |   | 1  0  0 |
 *
 *    where 'l' is the local terrestrial latitude, and 'h' is the local hour 
 *    in Greenwich Mean sidereal time (GMST) plus the local terrestrial 
 *    longitude.  The first rotation matrix is a simple change of axes.
 */
void ctmatrix(double mjday, double lat, double lon, double ZA[3][3]) {
  double hour = ( A + B*(mjday-J2000) ) * M_PI/12.0 + lon;

  double sh = sin(hour);
  double ch = cos(hour);
  double sl = sin(lat);
  double cl = cos(lat);

  ZA[0][0] = -sh;
  ZA[0][1] = ch;
  ZA[0][2] = 0.00;

  ZA[1][0] = -ch * sl;
  ZA[1][1] = -sh * sl;
  ZA[1][2] = cl;

  ZA[2][0] = ch * cl;
  ZA[2][1] = sh * cl;
  ZA[2][2] = sl;
}

void ccs2tcs(double mjday, double lat, double lon, double v1[3], double v2[3]) {
  double ZA[3][3], vt[3];
  ctmatrix(mjday, lat, lon, ZA);

  /* in case v1 == v2 */
  vt[0] = v1[0];
  vt[1] = v1[1];
  vt[2] = v1[2];

  v2[0] = ZA[0][0]*vt[0] + ZA[0][1]*vt[1] + ZA[0][2]*vt[2];
  v2[1] = ZA[1][0]*vt[0] + ZA[1][1]*vt[1] + ZA[1][2]*vt[2];
  v2[2] = ZA[2][0]*vt[0] + ZA[2][1]*vt[1] + ZA[2][2]*vt[2];

}

void tcs2ccs(double mjday, double lat, double lon, double v1[3], double v2[3]) {
  double ZA[3][3], vt[3];
  ctmatrix(mjday, lat, lon, ZA);

  /* in case v1 == v2 */
  vt[0] = v1[0];
  vt[1] = v1[1];
  vt[2] = v1[2];

  /* ZA is unitary -- inverse ZA is same as transpose ZA */
  v2[0] = ZA[0][0]*vt[0] + ZA[1][0]*vt[1] + ZA[2][0]*vt[2];
  v2[1] = ZA[0][1]*vt[0] + ZA[1][1]*vt[1] + ZA[2][1]*vt[2];
  v2[2] = ZA[0][2]*vt[0] + ZA[1][2]*vt[1] + ZA[2][2]*vt[2];

}

#undef J2000
#undef A
#undef B

/*
 *  Gives the latitude (radians), longitude (radians), and altitude (meters) 
 *    of a coordinate (x, y, z) (meters) wrt. the center of the Earth.
 *
 *  Input:
 *    double x --- x-coordinate (Axis pointing from 180W, 0N to 0E, 0N)
 *    double y --- y-coordinate (Axis pointing from 90W, 0N to 90E, 0N)
 *    double z --- z-coordinate (Axis pointing from south pole to north pole)
 *
 *  Output:
 *    double *lat --- the latitude cooresponding to the given point in radians
 *    double *lon --- the longitude in radians
 *    double *alt --- the altitude above mean sea level in meters.
 *
 *  Returns
 *    void
 */
void xyz2lla(double x, double y, double z, double *lat, double *lon, double *alt) {
  double rho, sinq, cosq, num, den, norm;

  rho = ECC1 * sqrt( x*x + y*y );
  norm = sqrt( rho*rho + z*z );
  sinq = z / norm;
  cosq = rho / norm;

  num = z + ECC2*REQ/ECC1 * CB(sinq);
  den = rho/ECC1 - ECC2*REQ * CB(cosq);
  norm = sqrt( num*num + den*den );

  *lat = atan(num/den);
  *lon = atan2(y, x);

  sinq = num / norm;
  cosq = den / norm;

  *alt = rho/(ECC1*cosq) - REQ/sqrt(1.0 - ECC2*sinq*sinq);

  return;
}

/*
 *  Gives the coordinate (x, y, z) (meters) relative to the center of the 
 *    Earth for the given latitude (radians), longitude (radians), and 
 *    sea-level altitude (meters).
 *
 *  Input:
 *    double lat --- a latitude (radians)
 *    double lon --- a longitude (radians)
 *    double alt --- an altitude (meters above mean sea level)
 *
 *  Output:
 *    double *x --- x-coordinate (Axis pointing from China to Africa)
 *    double *y --- y-coordinate (Axis pointing from S. America to India)
 *    double *z --- z-coordinate (Axis pointing form S Pole to N pole)
 *
 *  Returns:
 *    void
 */
void lla2xyz(double lat, double lon, double alt, double *x, double *y, double *z) {
  double cosq, sinq, r0, as;

  cosq = cos(lat);
  sinq = sin(lat);

  r0 = REQ/sqrt(cosq*cosq+ECC1*ECC1*sinq*sinq);
  as = ECC1 * ECC1 * r0;

  *x = ( r0 + alt ) * cosq * cos(lon);
  *y = ( r0 + alt ) * cosq * sin(lon);
  *z = ( as + alt ) * sinq;

  return;
}

/*
 *  Convenience function for getting the latitude, longitude, and altitude 
 *    for a vector r=(x, y, z) in Earth-center coordinates.
 */
void r2lla(double *r, double *lat, double *lon, double *alt) {
  xyz2lla(r[0], r[1], r[2], lat, lon, alt);
  return;
}

/*
 *  Convenience function for getting the vector r=(x, y, z) (meters, Earth-
 *    center coords) from a latitude, longitude, and altitude.
 */
void lla2r(double lat, double lon, double alt, double *r) {
  lla2xyz(lat, lon, alt, &r[0], &r[1], &r[2]);
  return;
}

#undef REQ
#undef EPS
#undef ECC1
#undef ECC2
#undef CB
