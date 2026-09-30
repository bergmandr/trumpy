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

#include "histogram.h"

//#define NOT_EXACTLY_ONEGRAM

/*
 *  Defines air shower segments
 *
 *  This function is slow, but not as slow as it used to be.  We are
 *  now using qsimp (TRACK_USE_QSIMP) to integrate the density between
 *  two points on a line and finding a distance of exactly 1 gram by
 *  successive approximations.
 */
int buildTrack(RuntimeParameters *par, AirShower *as,
	       GaisserHillasParameters *gh, Track *track) {

  int nseg = buildTrackGeom(par, as, track);
  fillTrackLight(par, as, gh, track);
  return nseg;
}

int buildTrackGeom(RuntimeParameters *par, AirShower *as, Track *track) {
  int nseg, nsegest;
  double lat, lon, alt, vclf[3], J[3][3];
  double dl, x, dx, mz[3], r[3], uv[3], ri[3], rp[3];
  double dlm;
  double position;
  static const double dxEps = DEPTH_RESOLUTION/10.;
  TrackSegment *sgmt;

  fdatmos_param_dst_common *aparam = &track->aparam;
  fdscat_dst_common *atrans = &track->atrans;

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
  track->length = findTrackLength(as->zenith, alt=50000.0);

  /* get time origin */
  track->t0 = -track->length / SPEED_OF_LIGHT;

  /* Find vector, r, pointing from impact point to the top of the
   * first shower segment (CLF coord) */
  subvec(as->impactv, track->length*as->trackuv, r);

  /* Transform 'r' and 'uv' from CLF to Earth Center coordinates */
  applyRotation(J, r, r);
  addvec(r, vclf, r);
  applyRotation(J, as->trackuv, uv);

  /* Find impact point in ECC */
  applyRotation(J, as->impactv, ri);
  addvec(ri, vclf, ri);

  /* get the atmosphere for this day */
  getAtmosphere(par, aparam, atrans);

  /*
   * Find the real altitude at the top of the track and estimate the slant 
   * depth there by finding the vertical depth and comparing the track angle 
   * to the local vertical, mz.
   */
  r2lla(r, &lat, &lon, &alt);
  mz[0] = -J[0][2];
  mz[1] = -J[1][2];
  mz[2] = -J[2][2];
  x = getGrammageByAltitude(aparam, alt) / dotprod(mz, uv);
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
  nsegest = (int)(getGrammageBetweenPoints(aparam, r, ri)) + 1;
  track->segment = (TrackSegment*)malloc(nsegest*sizeof(TrackSegment));

  /*
   * Loop over segments until we reach the "ground" level at the
   * CLF
   */
  position = 0.;         // distance from top of first segment
  nseg = 0;
  while ( alt > CLF_ALTITUDE ) {
    /* Allocate more segments if necessary */
    if ( nseg >= nsegest ) {
      nsegest += 25;
      track->segment = 
	(TrackSegment*)realloc(track->segment, nsegest*sizeof(TrackSegment));
    }

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
    addvec(r, dl*uv, rp);
#endif
    sgmt->dxseg = dx;
    sgmt->dlseg = dl;

    /* Find midpoint in X of segment */
    dlm = getLRhoIntegralBetweenPoints(aparam, r, rp) / dx;
    sgmt->dlmid = dlm;

    /* Setup for next segment */
    nseg++;
    x += dx;
    position += dl;
    memcpy(r, rp, 3*sizeof(double));
    r2lla(r, &lat, &lon, &alt);
  }

  /* remove wasted space at end of arrays */
  if ( nseg < nsegest )
    track->segment =
      (TrackSegment*)realloc(track->segment, nseg*sizeof(TrackSegment));

  track->nseg = nseg;
  return nseg;
}

void fillTrackLight(RuntimeParameters *par, AirShower *as,
		    GaisserHillasParameters *gh, Track *track) {
  int i, nseg;
  double lat, lon, alt, vclf[3], J[3][3];
  double dl, x, dx, s, r[3], lr[3], uv[3], rp[3], rm[3];
  double ldl, ldx, ls;
  double xm, dlm;
  double am, lam;
  double nch, lnch, att[NWAVELEN_BANDS];
  double dnfl[NWAVELEN_BANDS], dncv[NWAVELEN_BANDS];
  double nfl[NWAVELEN_BANDS], ncv[NWAVELEN_BANDS];
  double pcv[NWAVELEN_BANDS];

  TrackSegment *sgmt, *lsgmt;

  /* get the current atmosphere */
  fdatmos_param_dst_common *aparam = &track->aparam;
  fdscat_dst_common *atrans = &track->atrans;
  getAtmosphere(par, aparam, atrans);  // may not be necessary now

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
  subvec(as->impactv, track->length*as->trackuv, r);

  /* Transform 'r' and 'uv' from CLF to Earth Center coordinates */
  applyRotation(J, r, r);
  addvec(r, vclf, r);
  applyRotation(J, as->trackuv, uv);
  
  r2lla(r, &lat, &lon, &alt);

  ldl = 0.0;		// length of previous segment
  ldx = 0.0;		// grammage of previous segment
  lam = alt;            // altitude at middle of prev segment
  ls = ageS(track->xoffset,gh->xmax); // age at middle of prev segment
  lnch = 0.0;           // num of charge particles at mid of prev segment
  lsgmt = NULL;         // pointer to previous segment
  for ( nseg=0; nseg<track->nseg; nseg++ ) {
    sgmt = &track->segment[nseg];
    alt = sgmt->height; // Find altitiude at TOP of segment
    x = sgmt->sdepth;   // Record slant depth at TOP of segment
    dl = sgmt->dlseg;   // Find the length of the segment
    addvec(r, dl*uv, rp);
    dx = sgmt->dxseg;   // Find the thickness of the segment

    dlm = sgmt->dlmid;  // Find midpoint of the segment
    xm = x + dx/2.;     // Find depth at midpoint
    s = ageS(xm,gh->xmax);       // Find shower age at midpoint
    sgmt->age = s;
    nch = (gh->lambda == -1)?gaussianInAgeFunction(gh,x):gaisserHillasFunction(gh, x); // Find shower size at midpoint
    sgmt->nch = nch;

    /* Find altitude at midpoint */
    addvec(r, dlm*uv, rm);
    r2lla(rm, &lat, &lon, &am);

    sgmt->molrad = getMoliereRadius (aparam, am);
    
    /*
     * Find the optical attenuation between top previous segment and
     * here.
     */
    subvec(r, ldl*uv, lr);
    getAttenuation(aparam, lr, r, att);

    /*
     * Get the rate, dN/dX, of production of Cerenkov photons, in 5 nm
     * bins, per shower particle, in the previous segment
     */
    getCvPhotonIncrRate(aparam, ls, lam, dncv);

    /*
     * Get the rate, dN/dX, of production of fluorescence photons, in
     * 5 nm bins (using FLASH spectrum), per eV deposited by shower
     */
    // getFluorescenceYieldByMeter(aparam, s, am, dnfl);
    getFluorescenceYieldByMeV(aparam, am, dnfl);

    /* Store energy deposited in the segement */
    sgmt->dedep = dx * nch * getAlphaEff(s);

#ifdef HISTOGRAM_MODE
    float yy = (float)( getAlphaEff(s) * 1.0e-6 );
    float sumfl = 0.00;
    if ( nch > 0.00 )
      for ( i=0; i<NWAVELEN_BANDS; i++ )
	sumfl += (float)( ( dnfl[i] * sgmt->dedep * 1.0e-6 ) / ( dl * nch ) );

    double rho, pres, temp;
    getPTDByAltitude(aparam, alt, &pres, &temp, &rho);
    pres /= 0.133322;

    hfill(402, (float)s, yy, 1.0);
    hfill(403, (float)(alt/1000.0), yy, 1.0);
    hfill(404, (float)pres, yy, 1.0);

    hfill(405, (float)s, sumfl, 1.0);
    hfill(406, (float)(alt/1000.0), sumfl, 1.0);
    hfill(407, (float)pres, sumfl, 1.0);
#endif

    /* Find number of fluorescence & Cherenkov photons produced */
    if (nseg==0) {
      // no need to be fancy (SS 10.4.11)
      for ( i=0; i<NWAVELEN_BANDS; i++ ) {
	pcv[i] = 0.0;
	ncv[i] = 0.0;
	nfl[i] = 0.0;
      }
    }
    else {
      for ( i=0; i<NWAVELEN_BANDS; i++ ) {
	pcv[i] = ldx * lnch * dncv[i];
	ncv[i] = att[i] * lsgmt->ncv[i] + pcv[i];
	// nfl[i] = nch * sgmt->dlseg * dnfl[i];
	nfl[i] = sgmt->dedep * dnfl[i] / 1.0e6;
      }
    }
    memcpy(sgmt->pcv, pcv, NWAVELEN_BANDS*sizeof(double));
    memcpy(sgmt->ncv, ncv, NWAVELEN_BANDS*sizeof(double));
    memcpy(sgmt->nfl, nfl, NWAVELEN_BANDS*sizeof(double));

    /* Save values for next segment */
    ldl = dl;
    ldx = dx;
    lam = am;
    ls = s;
    lnch = nch;
    lsgmt = sgmt;
    memcpy(r, rp, 3*sizeof(double));
  }
}

void clearTrack(Track *track) {
  track->nseg = 0;
  track->xoffset = 0.;
  track->length = 0.;
  free(track->segment);
}

/*
 * Uses the law of sines to compute the length of the shower track through the 
 *   atmosphere from the ground up to the given altitude.  i.e.
 *
 *                     sin(a)   sin(b)   sin(c)
 *                     ------ = ------ = ------
 *                        A        B        C
 *
 *   where A, B, & C are the lengths of the sides of a triangle, and a, b, & c 
 *   are the angles opposite to them ( a is angle between sides B & C, b is 
 *   angle between A & C, and c is angle between A & B ).
 */
double findTrackLength(double zenith, double altitude) {
  double A, B, C;
  double a, b, c;
  double r;

  A = EARTH_RADIUS - CLF_ALTITUDE + altitude;
  C = EARTH_RADIUS;

  if (zenith > (M_PI/2.))
    a = zenith;
  else
    a = M_PI - zenith;

  if ( a >= M_PI )
    return altitude - CLF_ALTITUDE;

  r = sin( a ) / A;
  c = asin( C * r );

  b = M_PI - a - c;
  B = sin( b ) / r;

  return B;
}

/*
 *  Input:
 *    AirShower *as    -- Extensive Air Shower object
 *    double *cvec     -- Observer coordinate wrt CLF (km)
 *    double *zvec     -- Observer viewing direction vector
 *    double tolerance -- Observer field of view (radians)
 *
 *  Output:
 *    int *segnum      -- Array of indices of shower segments in FOV
 *
 *  Returns:
 *    The number of shower segments in 
 */
int getSegmentsInView(AirShower *as, Track *track, double *cvec, 
		      double *zvec, double tolerance, int *segnum) {
  int i, nseg;
  double cq, ct, s[3], su[3];
  double dlseg;

  ct = cos( tolerance );

  /* start 's' to point to 'top' of air shower */
  s[0] = as->impactv[0] - cvec[0] - track->length*as->trackuv[0];
  s[1] = as->impactv[1] - cvec[1] - track->length*as->trackuv[1];
  s[2] = as->impactv[2] - cvec[2] - track->length*as->trackuv[2];

  nseg = 1;
  for ( i=1; i<track->nseg; i++ ) {
    unitVector(s, su);

    cq = dotprod(su, zvec);

    if ( cq >= ct )
      segnum[nseg++] = i;

    dlseg = track->segment[i].dlseg;
    s[0] += dlseg * as->trackuv[0];
    s[1] += dlseg * as->trackuv[1];
    s[2] += dlseg * as->trackuv[2];
  }

  if ( nseg == 1 )
    return 0;

  segnum[0] = segnum[1] - 1;

  return nseg;
}

/*
 *  Input:
 *    double r[3]      --- point in site's coordinate system.
 *    int nvec         --- number of 'cameras' in cvec anc zvec list
 *    double cvec[][3] --- location of observer in site coordinates
 *    double zvec[][3] --- observer viewing direction
 *    double tolerance --- opening angle of obersver's FOV (radians)
 *
 *  Output:
 *    int *camera      --- list of cameras that 'see' point 'r'
 *
 *  Returns:
 *    The number of cameras that 'see' the point 'r'.
 */
int getTelescopesInView(double r[3], int nvec, double cvec[][3], 
			double zvec[][3], double tolerance, int *camera) {
  int i, ncview;
  double ct, cq, s[3], su[3];

  ct = cos( tolerance );

  ncview = 0;
  for ( i=0; i<nvec; i++ ) {
    s[0] = r[0] - cvec[i][0];
    s[1] = r[1] - cvec[i][1];
    s[2] = r[2] - cvec[i][2];

    unitVector(s, su);

    cq = dotprod(su, zvec[i]);

    if ( cq >= ct )
      camera[ncview++] = i;
  }

  return ncview;
}
