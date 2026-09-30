/*
 * Convert geodectic Latitude, Longitude, Altitude to ECEF X, Y, Z
 * Altitude is geodectic altitude (km)
 * Latitude and Longitude in degrees!
 */

#include <math.h>
#include "earth.h"
#include "toolbox.h"

void latlonalt_to_xyz (double Latitude, double Longitude, double Altitude,
                       double *X, double *Y, double *Z) {
  double r0, as;
  double cL = cos(Latitude);
  double sL = sin(Latitude);
  double iFLAT = 1.0 - FLAT;

  r0 = R_EQ/sqrt(cL*cL + iFLAT*iFLAT*sL*sL);
  as = iFLAT * iFLAT * r0;

  *X = ( r0 + Altitude ) * cos(Latitude);
  *X *= cos(Longitude);
  *Y = ( r0 + Altitude ) * cos(Latitude);
  *Y *= sin(Longitude);
  *Z = ( as + Altitude ) * sin(Latitude);
}

void latlonalt_to_R(double A, double B, double C, double *R) {
  latlonalt_to_xyz(A,B,C,R,R+1,R+2);
}
