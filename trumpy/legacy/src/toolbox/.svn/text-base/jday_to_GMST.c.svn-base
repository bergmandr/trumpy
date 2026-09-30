#include <math.h>
#include "toolbox.h"

/* Returns Greenich Mean Standard Time (GMST) in radians */
double jday_to_GMST (double jday) {

  double T0, T;
  double gmst;

  /* Getting Julian days since Jan 1 2000 at noon */
  T0 = jday - 2451545.0;

  /* Getting T as julian centuries since Jan 1 2000 at noon */
  T = T0 / 36525.0;

  /* GMST in degrees */
  gmst = 280.46061837 + 360.98564736629 * T0
                      + 0.000387933 * T * T
                      - (1/38710000.0) * T * T * T;

  /* Getting gmst between 0 and 360 */
  gmst = range (gmst*M_PI/180., 2.*M_PI);

  return gmst;
}
