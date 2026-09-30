#include <stdio.h>
#include <math.h>

#include "event.h"
#include "control.h"
#include "toolbox.h"

#define  CORR  1.343e-06

/*
 *  From WGS84 data.
 */
#define REQ       6378137.0     // R_Earth @ equator (m)
#define INV_FLAT  298.257223    // Inverse flatness
#define RPO  (REQ-REQ/INV_FLAT) // R_Earth @ poles (m) (=REQ*(1-1/INV_FLAT))

/*
 *  Returns the distance from vector 'rr_ec' along unit vector 'uv_ec' to get 
 *    to an altitude of 'z' above sea level.  The vector to that point in 
 *    Earth's center coordinates is stored in 'rz_ec'.  Arguments are also in 
 *    Earth's center coordinates.
 */
double getDistanceToAltitude_ECC(double rr_ec[3], double uv_ec[3], double z, double rz_ec[3]) {

  double a, b, c, d, dp, ld, lu;
  double rr_sec[3], uv_sec[3];
  double rp[3], ru[3], rd[3], lv[3];

  /*
   *  Apply correction to 'z' to compensate for machine precision errors.
   */
  z *= 1.0 / ( 1.0 - CORR );

  /*
   *  Transform coordinates to 'stretched' Earth coorinates, where z-axis is 
   *    stretched so that altitude 'z' forms a sphere.
   */
  a = REQ + z;
  b = RPO + z;
  cpyvec(rr_ec, rr_sec);
  cpyvec(uv_ec, uv_sec);

  rr_sec[2] *= a / b;
  uv_sec[2] *= a / b;

  /*
   *  Renormalize the track direction vector.
   */
  c = magvec(uv_sec); 
  uv_sec[0] /= c; 
  uv_sec[1] /= c; 
  uv_sec[2] /= c;

  /*
   *  Get the impact parameter vector, then find the intersection of the EAS 
   *    track and a sphere of radius 'a' (takes advantage of spherical 
   *    symmetry, which is why we need a 'stretched' Earth conversion).
   */
  dp = dotprod(rr_sec, uv_sec);
  subvec(rr_sec, dp*uv_sec, rp);

  c = magvec(rp);
  d = sqrt( a*a - c*c );

  addvec(rp, d*uv_sec, rd);
  subvec(rp, d*uv_sec, ru);

  /*
   *  Convert 'rd' and 'ru' to real Earth coordinates by scaling back their 
   *    z-components.
   */
  rd[2] *= b / a;
  ru[2] *= b / a;

  /*
   *  Select the point nearest to the original vector
   */
  subvec(rr_ec, rd, lv);
  ld = dotprod(lv, uv_ec);
  double mld = fabs(ld);

  subvec(rr_ec, ru, lv);
  lu = dotprod(lv, uv_ec);
  double mlu = fabs(lu);

  if ( mlu < mld ) {
    cpyvec(ru, rz_ec);
    return lu;
  }
  else {
    cpyvec(rd, rz_ec);
    return ld;
  }
}

/*
 *  Returns the distance from 'rr' along unit vector 'uv' to an altitude 'z' 
 *    above sea level.  The vector to that point in CLF coordinates is stored 
 *    in 'rz'.  Arguments are all in CLF coordinates.
 */
double getDistanceToAltitude(double rr[3], double uv[3], double z, double rz[3]) {
  static double tv[3]={0.00}, rm[3][3]={{0.00}}, trm[3][3]={{0.00}};
  double rr_ec[3], uv_ec[3], ll, rz_ec[3];

  if ( fabs(tv[0]) < DOUBLE_FP_PRECISION ) {
    lla2r(CLF_LATITUDE, CLF_LONGITUDE, CLF_ALTITUDE, tv);
    trm[2][0] = cos(CLF_LATITUDE);
    trm[2][0] *= cos(CLF_LONGITUDE);
    trm[2][1] = cos(CLF_LATITUDE);
    trm[2][1] *= sin(CLF_LONGITUDE);
    trm[2][2] = sin(CLF_LATITUDE);
    zmatrix(trm[2], trm);
    matrixTranspose(trm, rm);

    pout("tv = ( %e, %e, %e )\n", tv[0], tv[1], tv[2]);
    pout("     | %9lf  %9lf  %9lf |\n", rm[0][0], rm[0][1], rm[0][2]);
    pout("rm = | %9lf  %9lf  %9lf |\n", rm[1][0], rm[1][1], rm[1][2]);
    pout("     | %9lf  %9lf  %9lf |\n", rm[2][0], rm[2][1], rm[2][2]);
  }

  /*
   *  Convert vectors to Earth-Center coordinates.
   */
  applyRotation(rm, rr, rr_ec);
  addvec(tv, rr_ec, rr_ec);

  applyRotation(rm, uv, uv_ec);

  ll = getDistanceToAltitude_ECC(rr_ec, uv_ec, z, rz_ec);

  /*
   *  Convert return vector to CLF coordinates.
   */
  subvec(rz_ec, tv, rz_ec);
  applyRotation(trm, rz_ec, rz);

  return ll;
}

//#define DISTTOALT_TEST_MODE
#ifdef DISTTOALT_TEST_MODE

int main(void) {
  static double T[3] = {
    -1.9243628455057752e+06,
    -4.5536647143797372e+06,
     4.0187567843789104e+06 
  };

  static double R[3][3] = {
    {  0.9211259691767480,  0.2465366275835532, -0.3012418300434471 },
    { -0.3892646258115379,  0.5833853758147893, -0.7128355731704474 },
    {  0.0000000000000000,  0.7738741464509364,  0.6333394077860895 } 
  };

  double rr[3] = { 0.00, 0.00, 0.00 };
  //double uv[3] = { 0.00, sqrt(3.0)/2.0, -0.5 };
  double uv[3] = { 0.00, 0.00, -1.0 };
  double rz[3];

  double q = 0.00;
  double z = 10000.0;
  while ( q <= 60.0 ) {
    uv[1] = sin( q * M_PI / 180.0 );
    uv[0] = 0.00;
    uv[2] = cos( q * M_PI / 180.0 );

    double ll = getDistanceToAltitude(rr, uv, z, rz);

    applyRotation(R, rz, rz);
    addvec(rz, T, rz);

    double lat, lon, alt;
    r2lla(rz, &lat, &lon, &alt);

    printf("%e  %e  %e  %e\n", q, alt, (alt-z)/z, ll);

    q += 1.0;
  }

  return 0;
}

#endif
