#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "constants.h"
#include "control.h"
#include "toolbox.h"

/*
 * Provide: Julian day, Latitude, Longitude, Altitude, Zenith and 
 *          Azimuth (radians)
 * Returns: RA and Dec (radians)
 */

void zenazm_to_radec (double jday, double lat, double lon,
		      double zen, double azm, double *ra, double *dec) {

  double lmst, ha;

  double sinzen, coszen;
  double sinazm, cosazm;
  double sinlat, coslat;

  lmst = range (jday_to_GMST (jday) + lon, (2.*M_PI));

  sinzen = sin(zen);
  coszen = cos(zen);
  sinazm = sin(azm);
  cosazm = cos(azm);
  sinlat = sin(lat);
  coslat = cos(lat);

  *dec = asin ( coszen*sinlat + sinzen*sinazm*coslat );

  ha = atan2 ( -sinzen*cosazm, (coszen*coslat - sinzen*sinazm*sinlat ) );
  ha = range (ha, (2.*M_PI));

  *ra = range (lmst - ha, (2.*M_PI));
}

void zenazm_to_radec_ecc (double jday, double lat, double lon, double alt, 
			  double zen, double azm, double *ra, double *dec) {

  double lmst, ha;
  double v[3], modlat;

  double sinzen, coszen;
  double sinazm, cosazm;
  double sinlat, coslat;

  lmst = range (jday_to_GMST (jday) + lon, (2.*M_PI));

  lla2xyz (lat, lon, alt, &v[0], &v[1], &v[2]);
  unitVector (v, v);
  modlat = asin (v[2]);

  sinzen = sin(zen);
  coszen = cos(zen);
  sinazm = sin(azm);
  cosazm = cos(azm);
  sinlat = sin(modlat);
  coslat = cos(modlat);

  *dec = asin ( coszen*sinlat + sinzen*sinazm*coslat );

  ha = atan2 ( -sinzen*cosazm, (coszen*coslat - sinzen*sinazm*sinlat ) );
  ha = range (ha, (2.*M_PI));

  *ra = range (lmst - ha, (2.*M_PI));
}
