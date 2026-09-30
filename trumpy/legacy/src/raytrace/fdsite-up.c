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

#include "airshower.h"
#include "atmosphere-withdb.h"
#include "nerling.h"
#include "scatter.h"
#include "showerlib.h"
#include "track.h"
#include "fdsite.h"

void getLaserObservedTrack(double lambda, RuntimeParameters *par,
			   AirShower *as, const Track *track, 
			   geofd_dst_common *site,
			   ObservedTrack *obtrack) {
  getLaserObservedTrackGeom   (par, as, track, site, obtrack);
  if (obtrack->nseg == 0)
    return;
  fillLaserObservedTrackLight (lambda, par, track, obtrack);
}

void getLaserObservedTrackGeom(RuntimeParameters *par, AirShower *as,
			       const Track *track, geofd_dst_common *site, 
			       ObservedTrack *obtrack) {
  int i, j, n;
  double t, cq, r;
  double R[3][3];
  double *uv;       // Unit vector along shower track direction in site coords
  double sv[3];     // Vector to top of shower segment in site coords
  double smv[3];    // Vector to middle of shower segment in site coords
  double snv[3];    // Vector to bottom of shower segment in site coords
  double msuv[3];   // Unit vector from mid of seg to detector in site coords
  double rpv[3];    // Vector to point-of-closest-approach in site coords
                    // (and vector from site to impact point before removing
                    // component along track)
  double rv[3];     // unit vector from segment to mirror in site coords
                    // (and vector from Rp to segment in site coords for TOF)
  double sv_ecc[3]; // Vector to top of shower seg in Earth center coords (ECC)
  double smv_ecc[3];// Vector to mid of shower seg in Earth center coords (ECC)
  double mv_ecc[3]; // Mirror location in ECC
  double uv_ecc[3]; // Shower track direction in ECC
  double snext[3];  // Vector to next shower segment in ECC
  double t0[3], t1[3], mt0, mt1; // Temp. vectors for finding angular extent
  double dx, lat, lon, ha, hb, hm;
  ObservedTrackSegment *sgmt;
  TrackSegment *tsgmt;

  double areamr[GEOFD_MAXMIR];
  for (i=0; i<site->nmir; i++)
    areamr[i] = M_PI * site->diameters[i] * site->diameters[i] / 4.0;

  /* Get Rp origin (from CLF in CLF coordinates) */
  getObservedTrackOrigin(site, as, obtrack);

  /* Get shower vector in site coordinates */
  matrixTranspose(site->site2clf, R);
  applyRotation(R, as->trackuv, obtrack->uv);
  uv = obtrack->uv;                         // uv is track vector in site coord

  /* Get vector from rp origin to impact point */
  subvec(as->impactv, obtrack->local_vsite, sv);      // still in OC
  applyRotation(R, sv, sv);                           // sv now in site coord
  memcpy(obtrack->rimp, sv, 3*sizeof(double));

  /* sv is from site to top of first segment */

  /* Get track geometry in terms of Earth-Center coordinates */
  applyRotation(par->oc2ecc, as->impactv, sv_ecc);// rotate OC impact to ECC...
  addvec(sv_ecc, par->origin, sv_ecc);            // ...and add OC origin to it
  applyRotation(par->oc2ecc, as->trackuv, uv_ecc);// rotate shower vec to ECC..
  /* sv_ecc is from site to top of first segment in ECC */

  /* Find Rp vector & Psi angle */
  obtrack->rp = getRpVector(site, as, rpv);
  obtrack->psi = getPsiAngle(site, as);

  rpv[0] = -rpv[0];
  rpv[1] = -rpv[1];
  rpv[2] = -rpv[2];

  /* Find rpuv and npln */
  unitVector(rpv, obtrack->rpuv);
  getSDPNVector(obtrack->rpuv, obtrack->uv, obtrack->npln)  ;

  /* Find vector along track from first segment to Rp, and set initial time */
  addvec(-1.*rpv, -1.*sv, rv);
  //  obtrack->t0 = -sqrt( dotprod(rv, rv) ) / SPEED_OF_LIGHT;
  //  obtrack->t0 += RANDOM_NUMBER * 12800.;	// randomize t0
  t = track->t0;

  /* Get the atmosphere for this date */
  fdatmos_param_dst_common *aparam = newInstanceOf(fdatmos_param_dst_common);
  fdscat_dst_common *atrans = newInstanceOf(fdscat_dst_common);
  getAtmosphere(par, aparam, atrans);

  /* Allocate memory for ObservedTrack and all ObservedTrackSegment's */
  obtrack->nseg = track->nseg;
  obtrack->segment = 
    (ObservedTrackSegment*)malloc(track->nseg*sizeof(ObservedTrackSegment));

  /* Main loop over segments */
  for ( i=0; i<track->nseg; i++ ) {
    sgmt  = &obtrack->segment[i];
    tsgmt = &track->segment[i];
    sgmt->dlseg = tsgmt->dlseg;

    /* Get vector and time to middle of segment */
    addvec(sv, tsgmt->dlmid*uv, smv);
    t += tsgmt->dlmid/SPEED_OF_LIGHT;

    /* Get vector to next segment ,and calculate angular extent */
    addvec(sv, tsgmt->dlseg*uv, snv);

    sgmt->tgen = t;                    // Time at middle of segment

    /* Find viewing angle at middle of segment */
    unitVector(smv, msuv);
    sgmt->qview = acos( dotprod(msuv, uv) );

    /* Find middle and next segment in ECC */
    addvec(sv_ecc, tsgmt->dlmid*uv_ecc, smv_ecc);
    addvec(sv_ecc, sgmt->dlseg*uv_ecc, snext);

    /* Find fraction of light scattered out of beam by aerosols */
    r2lla(sv_ecc, &lat, &lon, &ha);
    r2lla(snext, &lat, &lon, &hb);
    getMieScatterFraction(ha, hb, sgmt->dlseg, sgmt->fmie);

    /* Find fraction of light scattered out of beam by air molecules */
    dx = tsgmt->dxseg;
    getRayleighScatterFraction(dx, sgmt->fray);

    /* Set variables for finding Cerenkov angular distribution */
    r2lla(smv_ecc, &lat, &lon, &hm);
    sgmt->altmid = hm;

    n = 0;
    /* Loop over mirrors/cameras at site */
    for ( j=0; j<site->nmir; j++ ) {

      if (n == GEOFD_MAXMIR) {
        vperr("Maximum number of mirrors reached!\n");
        abort();
      }

      /*
       * Find unit vector (and distance) from segment to mirror in
       * site coordinates
       */
      subvec(smv, obtrack->local_vmir[j], rv);
      r = magnitude(rv);
      unitVector(rv,rv);
      
      /* Find angle between segment and mirror pointing direction */
      cq = dotprod(rv, site->vmir[j]);

      /* If segment is more than 15.0 degrees from the mirror skip it */
      if ( cq > COS15 ) {
	sgmt->costhe[n] = cq;

	/* Find angular extent */
	subvec(sv,  obtrack->local_vmir[j], t0);
	subvec(snv, obtrack->local_vmir[j], t1);
	mt0 = magnitude(t0);
	mt1 = magnitude(t1);
	sgmt->dtheta[n] = acos(dotprod(t0,t1)/mt0/mt1);

	/*
	 * get mirror location in terms of Earth-Center coordinates
	 * rotate vector to mirror into ecc
	 */
	applyRotation(site->site2earth, obtrack->local_vmir[j], mv_ecc);
	/* add vector to site in ecc */
	addvec (mv_ecc, obtrack->vsite, mv_ecc); // now it's vec to mir in ecc

	/* Find attenuation by wavelength between segment and mirror */
	getAttenuation(aparam, smv_ecc, mv_ecc, sgmt->att[n]);

	/* Set distance and time of light arrival */
	sgmt->dist[n] = r;
	sgmt->geoFactor[n] = areamr[j]/(4.0*M_PI*r*r);
	sgmt->mir[n] = j;
	n++;
      }
    }

    sgmt->nmir = n;

    /* Move to next segment */
    t += (tsgmt->dlseg-tsgmt->dlmid) / SPEED_OF_LIGHT;
    addvec(sv, tsgmt->dlseg*uv, sv);
    addvec(sv_ecc, tsgmt->dlseg*uv_ecc, sv_ecc);
  }

  free(aparam);
  free(atrans);
}

void fillLaserObservedTrackLight(double lambda, const RuntimeParameters *par,
				 const Track *track, 
				 ObservedTrack *obtrack) {
  int i, j, k, wlbin = 0;
  double geo;
  double *att;      // Attenuation between segment and mirror
  double *fmie;     // Fraction of Cerenkov beam scattered by aerosols
  double *fray;     // Fraction of Cerenkov bean scattered by molecules
  double angle;
  ObservedTrackSegment *sgmt;
  TrackSegment *tsgmt;

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

  fdatmos_param_dst_common *aparam = newInstanceOf(fdatmos_param_dst_common);
  fdscat_dst_common *atrans = newInstanceOf(fdscat_dst_common);
  getAtmosphere(par, aparam, atrans);

  /* Main loop over segments */
  for ( i=0; i<track->nseg; i++ ) {
    sgmt  = &obtrack->segment[i];
    tsgmt = &track->segment[i];

    /* Get age and scattering fractions */
    fmie = sgmt->fmie;
    fray = sgmt->fray;
    
    /* Loop over mirrors/cameras at site */
    for ( j=0; j<sgmt->nmir; j++ ) {

      for ( k=0; k<NWAVELEN_BANDS; k++ ) {
	sgmt->nfl[j][k] = 0.00;
	sgmt->ncvdir[j][k] = 0.00;
	sgmt->ncvmie[j][k] = 0.00;
	sgmt->ncvray[j][k] = 0.00;
      }

      /* Get attenuation between segment and mirror, and geometric factor */
      att = sgmt->att[j];
      geo = sgmt->geoFactor[j];

      angle = M_PI - sgmt->qview;

      sgmt->ncvray[j][wlbin] += tsgmt->ncv[wlbin] * 4.*M_PI*geo
	* att[wlbin] * fray[wlbin]*(1.-fmie[wlbin])
	* rayleighScatterFunction(angle);

      sgmt->ncvmie[j][wlbin] += tsgmt->ncv[wlbin] * 4.*M_PI*geo
	* att[wlbin] * fmie[wlbin]*(1.-fray[wlbin])
	* mieScatterFunction(angle);

      /*
      printf("seg %4d mir %2d q %.1f att %e RAY %e frac %e | AER %e frac %e\n",
	     i, j, angle*R2D, att[wlbin],
	     sgmt->ncvray[j][wlbin], fray[wlbin],
	     sgmt->ncvmie[j][wlbin], fmie[wlbin]);
      */
    }
  }

  free(aparam);
  free(atrans);
}

void getUpwardFluxInStepsOfX(const Track *trk, const ObservedTrack *in,
			     int nseg, double *xtop, double *xbot, 
			     double areamr[GEOFD_MAXMIR], ObservedTrack *out) {
  int i, mir, iseg, jseg, ic, iw;
  double x1, x2, y1, y2, d1, d2, f;
  ObservedTrackSegment *xseg, *yseg;
  double area;

  /*
   * initialize output Observed Track --- don't worry about time or site 
   *    geometry.
   */
  out->nseg = nseg;
  out->segment = 
    (ObservedTrackSegment*)malloc(nseg*sizeof(ObservedTrackSegment));
  for ( jseg=0; jseg<nseg; jseg++ ) {
    out->segment[jseg].nmir = 0;
    for ( ic=0; ic<GEOFD_MAXMIR; ic++ ) {
      out->segment[jseg].mir[ic] = 0;
      out->segment[jseg].costhe[ic] = 0.00;
      for ( iw=0; iw<NWAVELEN_BANDS; iw++ ) {
        out->segment[jseg].nfl[ic][iw] = 0.00;
        out->segment[jseg].ncvdir[ic][iw] = 0.00;
        out->segment[jseg].ncvmie[ic][iw] = 0.00;
        out->segment[jseg].ncvray[ic][iw] = 0.00;
      }
    }
  }

  for (jseg=out->nseg-1; jseg>=0; jseg--) {    // Looping over data grammages
    yseg = &out->segment[jseg];
    y2 = xtop[jseg];
    y1 = xbot[jseg];

    area = 0.;

    for (iseg=in->nseg-1; iseg>=1; iseg--) {   // Looping over MC grammages
      xseg = &in->segment[iseg];
      x2 = trk->segment[iseg-1].sdepth;
      x1 = trk->segment[iseg].sdepth;

      if (x2 < y1)    // MC has not yet gotten to front edge of data bin
	continue;

      if (x1 > y2)    // MC has gotten past back edge of data bin
	break;

      //      printf("%d x1 %f x2 %f | %d y1 %f y2 %f\n", iseg, x1, x2, jseg, y1, y2);

      d1 = 0.;
      d2 = 0.;
      f  = 0.;

      if (x1 <= y1 && x2 >= y1 && x2 <= y2) {
	d2 = x2;
	d1 = y1;
      }
      else if (x1 <= y1 && x2 >= y2) {
	d2 = y2;
	d1 = y1;
      }
      else if (x1 >= y1 && x2 <= y2) {
	d2 = x2;
	d1 = x1;
      }
      else if (x1 >= y1 && x1 <= y2 && x2 >= y2) {
	d2 = y2;
	d1 = x1;
      }

      f = (d2 - d1) / (x2 - x1);

      if (f > 0.) {
	for ( ic = 0; ic<xseg->nmir; ic++ ) {
	  mir = xseg->mir[ic];
	  area += f * areamr[mir] * xseg->costhe[ic];

	  i = addListItem(mir, &yseg->nmir, GEOFD_MAXMIR, yseg->mir);
	  if (i < 0) {
	    vperr("Maximum number of mirrors reached!\n");
	    abort();
	  }

	  for ( iw=0; iw<NWAVELEN_BANDS; iw++ ) {
	    yseg->nfl[i][iw]    += f * xseg->nfl[ic][iw]    / xseg->dtheta[ic];
	    yseg->ncvdir[i][iw] += f * xseg->ncvdir[ic][iw] / xseg->dtheta[ic];
	    yseg->ncvmie[i][iw] += f * xseg->ncvmie[ic][iw] / xseg->dtheta[ic];
	    yseg->ncvray[i][iw] += f * xseg->ncvray[ic][iw] / xseg->dtheta[ic];
	  }
	}

	//	printf("%d x1 %f x2 %f | %d y1 %f y2 %f | a %f : f %f\n", iseg, x1, x2, jseg, y1, y2, area, f);

      }  // f > 0.
    }    // Loop through MC segments

    /* Divide out active area */
    for ( ic=0; ic<yseg->nmir; ic++ ) {
      for ( iw=0; iw<NWAVELEN_BANDS; iw++ ) {
	if (area > 0.) {
	  yseg->nfl[ic][iw]    /= area;
	  yseg->ncvdir[ic][iw] /= area;
	  yseg->ncvmie[ic][iw] /= area;
	  yseg->ncvray[ic][iw] /= area;
	}
      }
    }
  }      // Loop through data segments
}
