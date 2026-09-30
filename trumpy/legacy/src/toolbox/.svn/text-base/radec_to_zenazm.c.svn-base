#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "constants.h"
#include "control.h"
#include "toolbox.h"

/*
 * Provide: Julian day, Latitude, Longitude, Altitude, RA, and Dec (radians)
 * Returns: Zenith and Azimuth angles (radians)
 */

void radec_to_zenazm (double jday, double lat, double lon,
		      double ra, double dec, double *zen, double *azm) {

  double lmst, ha;

  double sinha, cosha;
  double sindec, cosdec;
  double sinlat, coslat;

  lmst = jday_to_GMST (jday) + lon;
  ha = lmst - ra;

  sinha  = sin(ha);
  cosha  = cos(ha);
  sindec = sin(dec);
  cosdec = cos(dec);
  sinlat = sin(lat);
  coslat = cos(lat);

  *zen = M_PI/2. - asin( cosha*cosdec*coslat + sindec*sinlat);

  *azm = atan2(sinha*cosdec, sinlat*cosha*cosdec - sindec*coslat);
  *azm = range (*azm, (2.*M_PI));

  /* Correction... */
  *azm -= M_PI;

  /* Converting from North 0 to East 0 */
  *azm = 5.*M_PI/2. - *azm;
  *azm = range (*azm, (2.*M_PI));
}

void radec_to_zenazm_ecc (double jday, double lat, double lon, double alt,
			  double ra, double dec, double *zen, double *azm) {

  double lmst, ha;
  double v[3], modlat;

  double sinha, cosha;
  double sindec, cosdec;
  double sinlat, coslat;

  lmst = jday_to_GMST (jday) + lon;
  ha = lmst - ra;

  lla2xyz(lat, lon, alt, &v[0], &v[1], &v[2]);
  unitVector (v, v);
  modlat = asin (v[2]);

  sinha  = sin(ha);
  cosha  = cos(ha);
  sindec = sin(dec);
  cosdec = cos(dec);
  sinlat = sin(modlat);
  coslat = cos(modlat);

  *zen = M_PI/2. - asin( cosha*cosdec*coslat + sindec*sinlat);

  *azm = atan2(sinha*cosdec, sinlat*cosha*cosdec - sindec*coslat);
  *azm = range (*azm, (2.*M_PI));

  /* Correction... */
  *azm -= M_PI;

  /* Converting from North 0 to East 0 */
  *azm = 5.*M_PI/2. - *azm;
  *azm = range (*azm, (2.*M_PI));
}
