#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "constants.h"
#include "control.h"
#include "event.h"
#include "toolbox.h"


#include "airshower.h"
#include "atmosphere-withdb.h"
#include "flash.h"
#include "nerling.h"
#include "scatter.h"
#include "showerlib.h"
#include "track.h"
#include "nkg.h"
#include "fdconstants.h"

//#define NOT_EXACTLY_ONEGRAM

/*
 *  Defines air shower segments
 *
 *  This function is slow, but not as slow as it used to be.  We are
 *  now using qsimp (TRACK_USE_QSIMP) to integrate the density between
 *  two points on a line and finding a distance of exactly 1 gram by
 *  successive approximations.
 */
int buildLaserTrack(double lambda, RuntimeParameters *par, AirShower *as,
                GaisserHillasParameters *gh, Track *track) {

  int nseg = buildUpwardTrackGeom(par, as, track);  
//   printf("done with bUTG; nseg=%d\n",nseg);
  fillLaserTrackLight(lambda, par, as, gh, track);
  return nseg;
}

int buildUpwardTrackGeom(RuntimeParameters *par, AirShower *as, Track *track) {
  int nseg, nsegest;
  double lat, lon, alt, vclf[3], J[3][3];
  double dl, x, dx, mz[3], r[3], uv[3], ri[3], rp[3];
  double dlm;
  double position;
  double end;
  static const double dxEps = DEPTH_RESOLUTION/10.;
  
  double CUTOFF_ALT, TOP_ALT;
  
  switch (par->siteid) {
    case TALE_SITEID:
      TOP_ALT = TOP_OF_ATMOSPHERE;
      CUTOFF_ALT = 0.9 * TOP_OF_ATMOSPHERE;
      break;
    default:
      TOP_ALT = 50000.0;
      CUTOFF_ALT = 40000.0;
      break;
  }
  
  TrackSegment *sgmt;

  /*
   * Determine the starting point in CLF and Earth Center coordinates
   */

  /* get the location of the CLF in Earth Center coords */
  lat = CLF_LATITUDE;
  lon = CLF_LONGITUDE;
  alt = CLF_ALTITUDE;
  lla2r(lat, lon, alt, vclf);

  /* get the rotation matrix that transforms the CLF to Earth Center coords */
  J[2][0] = cos(lon) * cos(lat);
  J[2][1] = sin(lon) * cos(lat);
  J[2][2] = sin(lat);
  zmatrix(J[2], J);
  matrixTranspose(J, J);

  /* approximate starting altitude (within 0.01%) */
  track->length = findTrackLength(as->zenith, alt=TOP_ALT);
//   track->length = findTrackLength(as->zenith, alt=50000.0);

//   printf("got track->length = %f\n",track->length);
  /* get time origin */
  track->t0 = 0.;

  /* Find starting point in ECC */
  applyRotation(J, as->impactv, r);
  addvec(r, vclf, r);

  /* Find ending point in ECC */
  applyRotation(J, as->trackuv, uv);
  addvec(r, track->length*uv, ri);

  /* get the atmosphere for this day */
//   fdatmos_param_dst_common *aparam = newInstanceOf(fdatmos_param_dst_common);
//   fdscat_dst_common *atrans = newInstanceOf(fdscat_dst_common);
  
  fdatmos_param_dst_common *aparam = &track->aparam;
  fdscat_dst_common *atrans = &track->atrans;
  
//   printf("about to getAtmosphere... ")     ;
  getAtmosphere(par, aparam, atrans);
//   printf("Got it!\n");
  /*
   * Find the real altitude at the top of the track and estimate the slant 
   * depth there by finding the vertical depth and comparing the track angle 
   * to the local vertical, mz.
   */
  r2lla(r, &lat, &lon, &alt);
  r2lla(ri, &lat, &lon, &end);
  mz[0] = J[0][2];
  mz[1] = J[1][2];
  mz[2] = J[2][2];
//   printf("about to getGrammageByAltitude\n");
  x = getGrammageByAltitude(aparam, end) / dotprod(mz, uv);
//   printf("and now, about to getGrammageBetweenPoints\n");
  x += getGrammageBetweenPoints(aparam, r, ri);
//   printf("done.\n");
  track->xoffset = x;

  /* 
   * Break shower segments down into steps of DEPTH_STEPSICE (1 g/cm2,
   * nominally) of slant depth.  Save the length of each segment in
   * km.  This is a numerical integration that takes into account the
   * ellipsoidal shape of the Earth and atmosphere. 
   */

  /*
   * Estimate the amount of memory needed for segment length array.
   * This expects a call to 'free()' once this information is no
   * longer needed.
   */
  
//   printf("calling it again...\n");
  nsegest = (int)(getGrammageBetweenPoints(aparam, r, ri)) + 1;
//   printf("okay, it worked that time too\n");
  track->segment = (TrackSegment*)malloc(nsegest*sizeof(TrackSegment));

  /*
   * Loop over segments until we reach the "ground" level at the
   * CLF
   */
  position = 0.;         // distance from top of first segment
  nseg = 0;
  while ( alt < end ) {
//     printf("alt %f end %f\n",alt,end);
    /* Allocate more segments if necessary */
    if ( nseg >= nsegest ) {
//       printf("realloc requested\n");
      nsegest += 25;
      track->segment = 
        (TrackSegment*)realloc(track->segment, nsegest*sizeof(TrackSegment));
    }
//     printf("past realloc\n");
    sgmt = &track->segment[nseg];

    /* Record position at top of segment */
    sgmt->position = position;

    /* Find altitiude at TOP of segment */
    sgmt->height = alt;

    /* Record slant depth at TOP of segment */
    sgmt->sdepth = x;

    /* Find how far, dl and dx, to the next segment */
#ifdef NOT_EXACTLY_ONEGRAM
    dxdl = 100.*getDensityByAltitude(aparam, alt);
    dl = DEPTH_STEPSIZE / dxdl;
    addvec(r, dl*uv, rp);
    dx = copysign(getGrammageBetweenPoints(aparam, r, rp), dl);
#else
    dx = DEPTH_STEPSIZE;
    dl = getDistanceByGrammage(aparam, r, uv, dx, dxEps); 
    
    while (par->siteid == TALE_SITEID && dl > 1000.0) {
//       printf("Halving dx for segment %d at altitude %f; currently %f (dl %f)\n",nseg,alt,dx,dl);
      dx /= 2.;
      dl = getDistanceByGrammage(aparam, r, uv, dx, dxEps);
    }
//     printf("segment altitude: %f\n",alt);
    
    addvec(r, dl*uv, rp);
#endif
    sgmt->dxseg = dx;
    sgmt->dlseg = dl;

    r2lla(rp, &lat, &lon, &alt);
//     if (alt > 40000.)
    if (alt > CUTOFF_ALT)
      break;

    /* Find midpoint in X of segment */
    dlm = getLRhoIntegralBetweenPoints(aparam, r, rp) / dx;
    sgmt->dlmid = dlm;

    //    printf("buildUpwardTrackGeom() : %d p h x dx dl %f %f %f %f %f\n", nseg, position, alt, x, dx, dl);
//        printf("alt_above_CLF %f gram %f\n", alt-CLF_ALTITUDE, x);

    /* Setup for next segment */
    nseg++;
    x -= dx;
    position += dl;
    memcpy(r, rp, 3*sizeof(double));
    r2lla(r, &lat, &lon, &alt);
//     printf("end of while-loop\n");
  }
  
//   abort();
  
//   printf("past while-loop\n");
  /* remove wasted space at end of arrays */
  if ( nseg < nsegest )
    track->segment = 
      (TrackSegment*)realloc(track->segment, nseg*sizeof(TrackSegment));

//   free(aparam);
//   free(atrans);
  
  track->nseg = nseg;
  return nseg;
}

void fillLaserTrackLight(double lambda, RuntimeParameters *par, AirShower *as,
                   GaisserHillasParameters *gh, Track *track) {
  int nseg, wlbin = 0;
  double lat, lon, alt, vclf[3], J[3][3];
  double dl, s, r[3], lr[3], uv[3], rp[3], rm[3];
  double ldl;
  double dlm;
  double am;
  double nch, att[NWAVELEN_BANDS];

  TrackSegment *sgmt;

  /* Get appropriate wavelength bin for laser */
/*   for (i=0; i<NWAVELEN_BANDS; i++) { */
/*     w0 = LAMBDA0 + (double)i*DLAMBDA; */
/*     w1 = w0 + DLAMBDA; */
/*     if (lambda >= w0 && lambda < w1) { */
/*       wlbin = i; */
/*       break; */
/*     } */
/*   } */
  wlbin = (int)( ( lambda - LAMBDA0 ) / DLAMBDA );

  /* get the current atmosphere */
//   fdatmos_param_dst_common *aparam = newInstanceOf(fdatmos_param_dst_common);
//   fdscat_dst_common *atrans = newInstanceOf(fdscat_dst_common);
  fdatmos_param_dst_common *aparam = &track->aparam;
  fdscat_dst_common *atrans = &track->atrans;
  getAtmosphere(par, aparam, atrans);

  /*
   * Determine the starting point in CLF and Earth Center coordinates
   */

  /* get the location of the CLF in Earth Center coords */
  lat = CLF_LATITUDE;
  lon = CLF_LONGITUDE;
  alt = CLF_ALTITUDE;
  lla2r(lat, lon, alt, vclf);

  /* get the rotation matrix that transforms the CLF to Earth Center coords */
  J[2][0] = cos(lon); J[2][0] *= cos(lat);
  J[2][1] = sin(lon); J[2][1] *= cos(lat);
  J[2][2] = sin(lat);
  zmatrix(J[2], J);
  matrixInverse(J, J);

  /*
   * Find vector, r, pointing from impact point to the top of the
   * first shower segment (CLF coord)
   */
  memcpy(r, as->impactv, 3*sizeof(double));

  /* Transform 'r' and 'uv' from CLF to Earth Center coordinates */
  applyRotation(J, r, r);
  addvec(r, vclf, r);
  applyRotation(J, as->trackuv, uv);
  
  r2lla(r, &lat, &lon, &alt);

  ldl = 0.0;            // length of previous segment
  nch = gh->nmax;
  for (nseg=0; nseg < track->nseg; nseg++) {
    sgmt = &track->segment[nseg];
    alt = sgmt->height; // Find altitiude at TOP of segment
    dl = sgmt->dlseg;   // Find the length of the segment
    addvec(r, dl*uv, rp);

    dlm = sgmt->dlmid;  // Find midpoint of the segment
    s = 0.;
    sgmt->age = s;
    sgmt->nch = nch;

    /* Find altitude at midpoint */
    addvec(r, dlm*uv, rm);
    r2lla(rm, &lat, &lon, &am);

    sgmt->molrad = 0.;

    /*
     * Find the optical attenuation between top previous segment and
     * here.
     */
    addvec(r, ldl*uv, lr);
    getAttenuation(aparam, lr, r, att);

    /* Find number of photons */
//     zeros(sgmt->pcv, NWAVELEN_BANDS);
//     zeros(sgmt->ncv, NWAVELEN_BANDS);
//     zeros(sgmt->nfl, NWAVELEN_BANDS);
    //    sgmt->ncv[wlbin] = nch;
  
    int i;
    for ( i=0; i<NWAVELEN_BANDS; i++) {
      sgmt->pcv[i] = 0.0;
      sgmt->ncv[i] = 0.0;
      sgmt->nfl[i] = 0.0;
    }
    //    printf("seg %4d N %e -> ", nseg, nch);

    nch *= att[wlbin];
    sgmt->ncv[wlbin] = nch;

    //  printf("%e (%e over %f meters)\n", sgmt->ncv[wlbin], att[wlbin], ldl);

    /* Save values for next segment */
    ldl = dl;
    cpyvec(rp, r);
  }
}
