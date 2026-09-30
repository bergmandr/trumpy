#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "constants.h"
#include "control.h"
#include "event.h"
#include "airshower.h"
#include "showerlib.h"
#include "track.h"
#include "fdconstants.h"
#include "fdsite.h"

/*
 * loadFDSiteGeometry(char *, geofd_dst_common *)
 *
 * Author: Sean R. Stratton
 *         Rutgers University, Dept. of Physics & Astronomy
 *
 * Initializes FDSiteGeometry structure using the specified DST file.
 *
 * note:  no error checking is done yet.
 */
int loadFDSiteGeometry(char *path, geofd_dst_common *fdsg) {
  int rc, filedes=FDGEOM_FILEDES, mode=MODE_READ_DST;
  int w, g, s=100, e;

  pout("loading FD Site Geometry\n");

  w = newBankList(s);
  g = newBankList(s);

  rc = dstOpenUnit(filedes, path, mode);
  if ( rc != 0 ) {
    perr("unable to open '%s' for reading. (%d)\n", path, rc);
    return rc;
  }

  rc = eventRead(filedes, w, g, &e);

  if ( tstBankList(g, GEOBR_BANKID) != 0 )
    geofd_ = geobr_;
  else if ( tstBankList(g, GEOLR_BANKID) != 0 )
    geofd_ = geolr_;
  else if ( tstBankList(g, GEOMD_BANKID) != 0 )
    geofd_ = geomd_;
  else if ( tstBankList(g, GEOTL_BANKID) != 0 )
    geofd_ = geotl_;
  else {
    perr("DST file does not contain valid bank.\n");
    return 1;
  }

  dstCloseUnit(filedes);

  *fdsg = geofd_;

  return 0;
}

#if 0
/*
 *  setFDSiteGeometry(RuntimeParameters, geofd_dst_common)
 *
 *  Author: Sean R. Stratton
 *          Rutgers University, Dept. of Physics & Astronomy
 *
 *  Initializes the FDSiteGeometry structure using the given coordinates.
 *
 *  Input:
 *    par -- set of runtime parameters.
 *
 *  Output:
 *    fdsg -- pointer to geofd dst bank.
 *
 *  Returns:
 *    '0' upon successful completion (always).
 *
 */
void setFDSiteGeometry(const RuntimeParameters *par, geofd_dst_common *fdsg) {
  double clfv[3], clfm[3][3], clfim[3][3];
  double fdsv[3], fdsm[3][3], fdsim[3][3];

  double lat, lon, alt;
  if ( par->siteid == 0 ) {
    lat = BR_LATITUDE;
    lon = BR_LONGITUDE;
    alt = BR_ALTITUDE;
  }
  else if ( par->siteid == 1 ) {
    lat = LR_LATITUDE;
    lon = LR_LONGITUDE;
    alt = LR_ALTITUDE;
  }
  else if ( par->siteid == 2 ) {
    lat = MD_LATITUDE;
    lon = MD_LONGITUDE;
    alt = MD_ALTITUDE;
  }
  else {
    lat = CLF_LATITUDE;
    lon = CLF_LONGITUDE;
    alt = CLF_ALTITUDE;
  }

  /* step 1: get vectors to CLF and the given location in Cartesian coords */
  lla2xyz(CLF_LATITUDE, CLF_LONGITUDE, CLF_ALTITUDE, 
		   &clfv[0], &clfv[1], &clfv[2]);
  lla2xyz(lat, lon, alt, &fdsv[0], &fdsv[1], &fdsv[2]);

  /* step 2: get the vector that translates the origin from the CLF to the 
   *           given coordinate (in CLF coordinates). */
  //  subvec(fdsv, clfv, fdsg->origin);

  /* step 3: find the rotation matrix that transforms the CLF c.sys to the 
   *           given location's */
  clfm[2][0] = cos(CLF_LONGITUDE) * cos(CLF_LATITUDE);
  clfm[2][1] = sin(CLF_LONGITUDE) * cos(CLF_LATITUDE);
  clfm[2][2] = sin(CLF_LATITUDE);

  fdsm[2][0] = cos(lon) * cos(lat);
  fdsm[2][1] = sin(lon) * cos(lat);
  fdsm[2][2] = sin(lat);

  zmatrix(clfm[2], clfm);
  zmatrix(fdsm[2], fdsm);

  matrixInverse(clfm, clfim);
  //  TODO:  fix this function or remove it altogether
  //  matrixMultiply(fdsm, clfim, fdsg->rmatrx);

  //  applyRotation(clfm, fdsg->origin, fdsg->origin);
}

#endif
