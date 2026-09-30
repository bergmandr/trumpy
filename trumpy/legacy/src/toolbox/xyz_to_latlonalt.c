/* Convert ECEF X, Y, Z to geodetic Latitude, Longitude, Altitude (km) */
#include <math.h>
#include "earth.h"
#include "toolbox.h"

void xyz_to_latlonalt (double X, double Y, double Z,
                       double *Latitude, double *Longitude, double *Altitude) {

  double ecc_prime2, R_Polar, theta;
  double p, top, bottom;

  *Longitude = atan2(Y,X);

  R_Polar = R_EQ*(1-FLAT);
  p = sqrt(X*X + Y*Y);
  theta = atan(Z*R_EQ/(p*R_Polar));
  ecc_prime2 = (R_EQ*R_EQ - (R_Polar*R_Polar))/(R_Polar*R_Polar);

  double st = sin(theta);
  double ct = cos(theta);
  top = Z + ecc_prime2 * R_Polar * st*st*st;
  bottom = p - ECC2 * R_EQ * ct*ct*ct;

  *Latitude = atan(top/bottom);

  double sl = sin(*Latitude);
  *Altitude = p/cos( *Latitude ) - 
    R_EQ/sqrt( 1 - ECC2 * sl * sl );

}

void R_to_latlonalt(double *R, double *A, double *B, double *C) {
  xyz_to_latlonalt(R[0],R[1],R[2],A,B,C);
}
