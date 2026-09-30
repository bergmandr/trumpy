/*
 * Routine to do raytracing for BR / LR FD
 * int num          number of photons to trace
 * int cam          camera id (0 - 11)
 * double vsite[]   vector from center of mirror to light source 
 *                  (in site coordinates)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "constants.h"
#include "control.h"
#include "event.h"
#include "fdconstants.h"
#include "toolbox.h"

#include "atmosphere-withdb.h"
#include "nkg.h"
#include "calibration.h"
#include "raytrace.h"

void doRayTrace (int num, int cam, double rpuv[],
		 double npln[], double vsite[], 
		 double age, double mol,
		 geofd_dst_common *geo, TGeom *tgeom, RayTrace *ray) {

  int i, hitsOnPMTPlane = 0;
  double v[3], disp[3];
  double alt, stheta;
  
// //   stheta = asin(vsite[2]/magvec(vsite)) * R2D;
  
// //   applyRotation(geo->site2clf,vsite,v);
// //   addvec(v,geo->local_vsite,v);
// //   alt=v[2] + CLF_ALTITUDE;
// //   printf("RAYTHROW mir %d alt %f stheta %f \n",cam, alt, stheta);
  ray->nthrow = num;
  ray->ngood  = 0;
  ray->cam = cam;

  /* v, vector from center of curvature to source */
  subvec(vsite, geo->rcurve[cam]*geo->vmir[cam], v);

  for (i=0; i<num; i++) {
    getNKGdr (age, mol, rpuv, npln, disp);

    /* Get final vector from center of curvature to source */
    addvec (v, disp, ray->vsite);

    /* Get vector to source in camera coordinates */
    applyRotation (geo->site2cam[cam], ray->vsite, ray->vsite);

     ray->pmt[i] = trace (geo, tgeom, ray);
/*    ray->pmt[i] = acptmap(tgeom, ray); */
    ray->tpe[i] = ray->time;

    if (ray->pmt[i] >= 0) {
      
// //       printf("RAYHIT mir %d tube %d row %d from alt %f stheta %f\n",cam,ray->pmt[i],ray->pmt[i] % 16, alt, stheta);
      ray->ngood++;
    }

    if ( (fabs(ray->xcam)<999.) && (fabs(ray->ycam)<999.) ) {
      ray->xave += ray->xcam;
      ray->yave += ray->ycam;
      hitsOnPMTPlane++;
    }
  }
  if (hitsOnPMTPlane > 0) {
    ray->xave /= (double)hitsOnPMTPlane;
    ray->yave /= (double)hitsOnPMTPlane;
  }
}
