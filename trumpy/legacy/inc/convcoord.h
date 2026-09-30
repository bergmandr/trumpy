#ifndef _CONVCOORD_H_
#define _CONVCOORD_H_

/*
 *  A set of optimized functions for converting between a coordinate in 
 *    Earth-center coordinates (units of meters) and latitude, longitude, 
 *    & altitude (radians, radians, meters above mean sea level).
 */
void xyz_to_latlonalt(double x, double y, double z, double *lat, double *lon, double *alt);
void R_to_latlonalt(double *r, double *lat, double *lon, double *alt);
void xyz2lla(double x, double y, double z, double *lat, double *lon, double *alt);
void lla2xyz(double lat, double lon, double alt, double *x, double *y, double *z);
void r2lla(double *r, double *lat, double *lon, double *alt);
void lla2r(double lat, double lon, double alt, double *r);

/*
 *  A set of functions for converting between celestial and terrestrial 
 *    coordinates.
 */

void ctmatrix(double mjday, double lat, double lon, double ZA[3][3]);
void ccs2tcs(double mjday, double lat, double lon, double v1[3], double v2[3]);
void tcs2ccs(double mjday, double lat, double lon, double v1[3], double v2[3]);

#endif
