/*
  Origin is center of curvature of spherical mirror
  (mirror axis aligned along z)
  when facing the camera, x is to the right, y is down.
  center of mirror at {0., 0., -rcurve}
  v is unit vector pointing back to source
  target plane is a "chord" on the sphere
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "constants.h"
#include "control.h"
#include "event.h"
#include "fdconstants.h"
#include "random.h"
#include "toolbox.h"

#include "calibration.h"
#include "tacalibration.h"
#include "raytrace.h"

int trace (geofd_dst_common *geo, TGeom *tgm, RayTrace *ray) {
  int i, tube, seg;

  double w[3];          // vector, to target plane
  double v[3];          // vector, target plane to light source
  double y[3];          // vector, to plane tangent to center of mirror
  double u[3];          // vector, to intersection point with mirror
  double n[3];          // vector, normal to mirror at intersection point
  double vp[3];         // vector, v with perpendicular component
                        //   to n subtracted
  double g[3];          // vector, intersection with mirror to camera box cover
  double ug[3];         // unit vector, g
  double f[3];          // vector, to intersection with PMT plane

  double gp[3];         // displacement vector inside camcover glass
  double gpp[3];  // displacement vector inside filter

  double cosvzen; // cosine of zenith angle of v[]
  double cosgzen; // cosine of zenith angle of g[]

  double dv[3];         // vector from intersection with mirror to light source

  double radius, the, phi;
  double B, C, t, factor, angle;
  double qeff;
  
  int num;

  /* Initialize some things */
  ray->xcam = 999.;
  ray->ycam = 999.;

  for (i=0; i<3; i++) {
    gp[i] = 0.;
    gpp[i] = 0.;
  }

  /* Random location on target plane */
  radius = (geo->diameters[ray->cam]/2.) * sqrt(RANDOM_NUMBER);
  phi = 2 * M_PI * RANDOM_NUMBER;
  w[0] = radius * cos(phi);
  w[1] = radius * sin(phi);
  w[2] = -(geo->rcurve3[ray->cam]) + tgm->h[ray->cam];

  /* Get vector to light source */
  for (i=0; i<3; i++)
    v[i] = ray->vsite[i] - w[i];
  unitVector (v, v);

  cosvzen = v[2];

  for (i=0; i<3; i++)
    y[i] = w[i] - (tgm->h[ray->cam]/cosvzen) * v[i];

  /* Find point where photon hits mirror.
   * For a sphere with center at 0, a ray starting at o
   * intersects the sphere at two points p.
   * p * p = R^2
   * if p = y + t * v, (y + t * v) * (y + t * v) - R^2 = 0
   * (v*v)t^2 + 2 * t * v * y + y*y - R^2 = 0
   * so v*v = 1., b = 2 * v * y, c = y^2 - R^2
   */

  B = 2 * dotProduct (v, y);
  double yy = magnitude(y);
  C = yy*yy - geo->rcurve3[ray->cam]*geo->rcurve3[ray->cam];

  /* The shorter t is the one we want */
  t = ( -B - sqrt(B*B - 4*C) ) / 2.;

  /* Normal vector is equivalent / opposite to u */
  for (i=0; i<3; i++) {
    u[i] = y[i] + t * v[i];
    n[i] = -u[i];
  }

  /* Getting mean travel time (to within a nanosecond) */           
  for (i=0; i<3; i++)
    dv[i] = ray->vsite[i] - u[i];
  ray->time = magnitude(dv) / SPEED_OF_LIGHT + tgm->meantime[ray->cam];

  /* Check camera box shadow */          
  if (USE_HIT_CAMERABOX && hitCameraBox (geo, tgm, ray, v, w))          
    return HIT_CAMERABOX;
  
  
  
  /* realistically speaking, an arriving photon will encounter
   *   the camera poles before striking the mirror.  However, 
   *   the rays must pass all of these tests regardless of 
   *   their order.  Since hitCameraPoles takes the longest,
   *   it make more sense computationally to put it after 
   *   the other tests to minimize the number of calls. */
  /* Check to see if photon hits camera poles */     
#ifdef RAYTRACE_VERBOSE
  switch (geo->camtype[ray->cam]) {
    case GEOFD_TA:
      num = NUMTAPOLES;
      break;
    case GEOFD_MD:
      num = NUMMDPOLES;
      break;
    case GEOFD_TALE:
      num = NUMTLPOLES;
      break;
  }
  if (USE_HIT_CAMERAPOLES && hitCameraPoles (num, tgm, ray, u, v))
    return HIT_CAMERAPOLES;
#endif
  /* Check to see if photon hits UVLED flasher at MD or TL */
  if (geo->siteid >= MIDDLE_DRUM_SITEID && USE_HIT_FLASHERBOX
      && hitFlasherBox(geo, tgm, ray, v, w))
    return HIT_FLASHERBOX;
  
  
  /* Check to see if photon hits mirror (or goes in crack) */
  switch (geo->camtype[ray->cam]) {
    case GEOFD_TA:
      if (USE_HIT_CRACK && (seg = hitTAMirror (geo, tgm, ray->cam, u))==-1)
        return HIT_CRACK;
      break;
    case GEOFD_MD:
    case GEOFD_TALE:
      if (USE_HIT_CRACK && (seg = hitMDMirror (geo, tgm, ray->cam, u))==-1)
        return HIT_CRACK;
      break;
  }
  
  for (i=0; i<3; i++)
    y[i] -= geo->seg_center[ray->cam][seg][i];
  B = 2 * dotProduct (v, y);
  yy = magnitude(y);
  C = yy*yy - geo->seg_rcurve[ray->cam][seg]*geo->seg_rcurve[ray->cam][seg];
  t = ( -B - sqrt(B*B - 4*C) ) / 2.;
  for (i=0; i<3; i++) {
    u[i] = y[i] + t * v[i];
    n[i] = -u[i];
  }
  
  for (i=0; i<3; i++) {
    u[i] += geo->seg_center[ray->cam][seg][i];
  }
  
  
  /* Check mirror absorption */
  if (USE_MIR_ABSORBED && RANDOM_NUMBER > tgm->mir_reflect)
    return MIR_ABSORBED;

#ifdef RAYTRACE_VERBOSE
    printf("uvec   %d %f %f\n", ray->cam, u[0], u[1]);
#endif

  /* 
   * Wiggle normal vector allow for mirror imperfections
   * Wiggle the vector {0, 0, R} and then rotate to it to the normal
   */
  if (USE_WIGGLE) {
    the = acos (n[2] / magnitude(n));
    phi = atan2 (n[1], n[0]);
//     factor = tgm->spot[ray->cam] * sqrt(-2. * log(RANDOM_NUMBER));
    factor = tgm->spot[ray->cam][seg] * sqrt(-2. * log(RANDOM_NUMBER));
    angle = 2.*M_PI * RANDOM_NUMBER;
    n[0] = factor * cos(angle);
    n[1] = factor * sin(angle);
    n[2] = geo->seg_rcurve[ray->cam][seg];
    roty (n, -the, n);
    rotz (n, -phi, n);
  }
  unitVector (n, n);

  for (i=0; i<3; i++) {
    vp[i] = v[i] - dotProduct (v, n) * n[i];
    g[i] = v[i] - 2 * vp[i];
  }
  unitVector (g, g);

  cosgzen = g[2];

  for (i=0; i<3; i++) {
    ug[i] = g[i];
    g[i] *= (-u[2] - (geo->rcurve3[ray->cam] - tgm->sep[ray->cam])) / cosgzen;
    f[i] = u[i] + g[i];
  }

  /* Check to make sure photon falls within camera box boundary */
  if (   USE_MISSED_CAMERA &&
       ( fabs(f[0]) > (geo->cam_width/2. ) ||
         fabs(f[1]) > (geo->cam_height/2.)
       )
     )
    return MISSED_CAMERA;

  if (USE_CAMCOVER_ABSORBED && RANDOM_NUMBER > tgm->camcover_trans)
    return CAMCOVER_ABSORBED;

  if (USE_CAMCOVER_SHIFT)
    snellShift (tgm, tgm->n_glass, tgm->thick_glass, ug, gp);

  if (USE_FILTER_ABSORBED && RANDOM_NUMBER > tgm->filter_trans)
    return FILTER_ABSORBED;

  if (geo->camtype[ray->cam] == GEOFD_TA && USE_FILTER_SHIFT)
    snellShift (tgm, tgm->n_filter, tgm->thick_filter, ug, gpp);

  /* Get g, f at PMT face */
  for (i=0; i<3; i++) {
    g[i] += ((tgm->thick_glass + tgm->thick_incam + 
            tgm->thick_filter)/cosgzen) * ug[i] + gp[i] + gpp[i];
    f[i] = u[i] + g[i];
  }

  /* camera coordinates of ray */
  ray->xcam =  f[0];
  ray->ycam = (geo->uniqID==0 || geo->uniqID==1320019200)?-f[1]:f[1];

#ifdef RAYTRACE_VERBOSE
    printf("oncam  %d %f %f\n", ray->cam, ray->xcam, ray->ycam);
#endif  
  /* Check PMT boundaries */
  tube = pmtDetect (geo, tgm, ray);

  if (tube < 0)
    return MISSED_TUBE;

  /* Get quantum efficiency at the hit location */
  qeff = 1.0;
  if (USE_PMT_UNIFORMITY)
    qeff *= pmtQE (geo, tgm, ray, tube);

  if (USE_PMT_QE)
    qeff *= tgm->pmtQE;

  if (RANDOM_NUMBER > qeff)
    return NO_PE_CONVERSION;
  
     
  
#ifndef RAYTRACE_VERBOSE

  /* see note near line 107 for explanation of why this is here. */  
  /* Check to see if photon hits camera poles */ 
    switch (geo->camtype[ray->cam]) {
      case GEOFD_TA:
        num = NUMTAPOLES;
        break;
      case GEOFD_MD:
        num = NUMMDPOLES;
        break;
      case GEOFD_TALE:
        num = NUMTLPOLES;
        break;
    }
    if (USE_HIT_CAMERAPOLES && hitCameraPoles (num, tgm, ray, u, v))
      return HIT_CAMERAPOLES;
#endif

#ifdef RAYTRACE_VERBOSE
  printf("intube %d %d %f %f\n", ray->cam, tube, ray->xcam, ray->ycam);
#endif

  /* Getting time of flight of ray */
  ray->time = magnitude(dv);
  ray->time += magnitude(g);
  ray->time /= SPEED_OF_LIGHT;

  return tube;
}
