#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "constants.h"
#include "control.h"
#include "event.h"
#include "fdconstants.h"
#include "toolbox.h"

#include "calibration.h"
#include "tacalibration.h"
#include "raytrace.h"

int hitTAMirror (geofd_dst_common *geo, TGeom *tgm, int cam, double u[]) {
  int i, j, seg;
  double dist[3], pdist[3], segdist[3];
  double dp, maxdot;

  /* Find closest mirror segment */
  seg = -1;
  maxdot = -1.;
  for (i=0; i<geo->nseg[cam]; i++) {
    dp = dotprod (u, geo->vseg3[cam][i]) / geo->rcurve3[cam];
    if (dp > maxdot) {
      maxdot = dp;
      seg = i;
    }
  }
//   printf("\n\n  (seg = %2d)\n",seg);
  for (i=0; i<3; i++)
    dist[i] = u[i] - geo->rcurve3[cam] * geo->vseg3[cam][seg][i];
//   printf("  (mag(dist) %f, tgm->maxseg %f)\n",magnitude(dist), tgm->maxseg[cam]);
//   printf("  (u %f %f %f\n)",u[0],u[1],u[2]);
  /* Check to see if it falls farther than maximum allowable distance */
  if ( magnitude(dist) > tgm->maxseg[cam] )
    return -1;

  /* Project dist onto plane perpendicular to geo->vseg vector */
  for (i=0; i<3; i++)
    pdist[i] = dist[i] - dotprod (dist, geo->vseg3[cam][seg]) * geo->vseg3[cam][seg][i];

  /*
   * Get vectors from each segment point to impact point
   * Dot hexmirpt and segdist - throw out if the angle
   * between them is greater than 60 degrees
   */
  for (i=0; i<6; i++) {
    for (j=0; j<3; j++)
      segdist[j] = (geo->seg_point2point/2) * tgm->hexmirpt[cam][seg][i][j] + 
      pdist[j];

    unitVector (segdist, segdist);

    if ( dotprod (tgm->hexmirpt[cam][seg][i], segdist) < 0.5 )
      return -1;
  }

  return seg;
}


int hitMDMirror (geofd_dst_common *geo, TGeom *tgm, int cam, double u[]) {
  int i, seg;
  double dist;
  double maxdot, dp;

  /* Check to see if it falls in cracks between mirrors */
  if (fabs(u[0]) < MD_CLOVER_GAP || fabs(u[1]) < MD_CLOVER_GAP)
    return -1;

//   /* Check to see if hits center button */
//   dist = u[0]*u[0] + u[1]*u[1];
//   if (dist < (HIRES_CLOVER_BUTTON*HIRES_CLOVER_BUTTON))
//     return FALSE;

  /* Find closest mirror segment */
  seg = -1;
  maxdot = -1.;
  for (i=0; i<geo->nseg[cam]; i++) {
    dp = dotprod (u, geo->vseg3[cam][i]) / geo->rcurve3[cam];
    if (dp > maxdot) {
      maxdot = dp;
      seg = i;
    }
  }

  if (maxdot > tgm->cosang[cam])
    return seg;
  else
    return -1;
}
