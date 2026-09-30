#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "control.h"
#include "event.h"
#include "fdconstants.h"
#include "toolbox.h"

#include "calibration.h"
#include "raytrace.h"

void snellShift (TGeom *tgm, double n2, double thick, double g[], double s[]) {
  double sinthe1 = sqrt(1. - g[2]*g[2]);
  double sinthe2 = tgm->n_air * sinthe1 / n2;

  /* Distance that vector is translated after passing through medium */
  double dist = thick * (sinthe1 - sinthe2);

  /* Break this distance into x and y components */
  double phi = atan2 (g[1], g[0]);

  /* translation vector */
  s[0] = -dist * cos (phi);
  s[1] = -dist * sin (phi);
  s[2] = 0.;
}

int pmtDetect (geofd_dst_common *geo, TGeom *tgm, RayTrace *ray) {
  int i;
  int tube;
  double dist, mindist;
  double segdist[3];

  tube = MISSED_TUBE;
  mindist = 1000.;

  for (i=0; i<GEOFD_MIRTUBE; i++) {
    dist = sqrt ( (ray->xcam - geo->xtube[i]) * (ray->xcam - geo->xtube[i]) +
                  (ray->ycam - geo->ytube[i]) * (ray->ycam - geo->ytube[i]) );

    if ( dist < (geo->pmt_flat2flat / 2.) ) {
      tube = i;
      mindist = dist;
      break;
    }

    if (dist < mindist) {
      tube = i;
      mindist = dist;
    }
  }

  if (mindist > (geo->pmt_point2point / 2.) )
    return MISSED_TUBE;

  /* Get vectors from each point to impact point
   * Dot hexpt and segdist - throw out if the angle
   * between them is greater than 60 degrees.
   */
  for (i=0; i<6; i++) {
    segdist[0] = (geo->pmt_point2point/2) * tgm->hexpt[i][0] +
      ray->xcam - geo->xtube[tube];
    segdist[1] = (geo->pmt_point2point/2) * tgm->hexpt[i][1] +
      ray->ycam - geo->ytube[tube];
    segdist[2] = (geo->pmt_point2point/2) * tgm->hexpt[i][2];

    /* Dot hexmirpt and segdist - throw out if the angle
     * between them is greater than 60 degrees (dot product less than 0.5)
     */
    unitVector (segdist, segdist);

    if ( dotProduct (tgm->hexpt[i], segdist) < 0.5 )
      return MISSED_TUBE;
  }

  return tube;
}

double pmtQE (geofd_dst_common *geo, TGeom *tgm, RayTrace *ray, int tube) {
  int xbin, ybin;
  double x, y;

  /* get x and y in millimeters */
  x = 1000.0 * (ray->xcam - geo->xtube[tube]);
  y = 1000.0 * (ray->ycam - geo->ytube[tube]);

  if ( fabs(x) > 4000. || fabs(y) > 4000. )
    return 0.;

  
  // Warning! in ray trace coordinates, x increases TO THE LEFT.
  if (x < 0.) // ray is to the RIGHT of PMT center; need to sample 
    xbin = PMT_NUMDIV/2 + (int)(fabs(x) + 0.5);
  else
    xbin = PMT_NUMDIV/2 - (int)(fabs(x) + 0.5);

  if (y < 0.)
    ybin = PMT_NUMDIV/2 - (int)(fabs(y) + 0.5);
  else
    ybin = PMT_NUMDIV/2 + (int)(fabs(y) + 0.5);

  return tgm->pmt_map[xbin][ybin];
}

