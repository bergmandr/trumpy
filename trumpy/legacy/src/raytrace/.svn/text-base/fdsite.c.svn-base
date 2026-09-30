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

void getObservedTrack(RuntimeParameters *par, AirShower *as,
                  const Track *track, geofd_dst_common *site,
                  ObservedTrack *obtrack) {
  getObservedTrackGeom   (par, as, track, site, obtrack);
  if (obtrack->nseg == 0)
    return;
  fillObservedTrackLight (par, track, obtrack);
}

void getObservedTrackGeom(RuntimeParameters *par, AirShower *as,
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
  unitVector(obtrack->uv,obtrack->uv);  // shouldn't be necessary,
                                        // but may be in rare cases.
  uv = obtrack->uv;                         // uv is track vector in site coord

  /* Get vector from rp origin to impact point */
  subvec(as->impactv, obtrack->local_vsite, sv);      // still in OC
  applyRotation(R, sv, sv);                           // sv now in site coord
  memcpy(obtrack->rimp, sv, 3*sizeof(double));

  /* Find vector to top of track (top of first segment) */
  subvec(sv, track->length*uv, sv);   // sv is from rp origin to top seg

  /* Get track geometry in terms of Earth-Center coordinates */
  applyRotation(par->oc2ecc, as->impactv, sv_ecc);// rotate OC impact to ECC...
  addvec(sv_ecc, par->origin, sv_ecc);            // ...and add OC origin to it
  applyRotation(par->oc2ecc, as->trackuv, uv_ecc);// rotate shower vec to ECC..
  subvec(sv_ecc, track->length*uv_ecc, sv_ecc);   // and use it to translate
                                                  // impact pt to top of trk
  /* Find Rp vector & Psi angle */
  obtrack->rp = getRpVector(site, as, rpv);
  obtrack->psi = getPsiAngle(site, as);

  /* Find rpuv and npln */
  unitVector(rpv, obtrack->rpuv);
  getSDPNVector(obtrack->rpuv, obtrack->uv, obtrack->npln)  ;

  /* Find vector along track from first segment to Rp, and set initial time */
  subvec(rpv, sv, rv);

  /* Set the time reference (t=0) to be time of shower core hitting
     the ground */
  //  obtrack->t0 = -track->length / SPEED_OF_LIGHT;
  //  obtrack->t0 = -sqrt( dotprod(rv, rv) ) / SPEED_OF_LIGHT;
  //  obtrack->t0 += RANDOM_NUMBER * 12800.;    // randomize t0
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
    msuv[0] = -smv[0];
    msuv[1] = -smv[1];
    msuv[2] = -smv[2];
    unitVector(msuv, msuv);
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
        sgmt->geoFactor[n] = cq*areamr[j]/(4.0*M_PI*r*r);
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

void fillObservedTrackLight(const RuntimeParameters *par, const Track *track, 
                      ObservedTrack *obtrack) {
  int i, j, k;
  double s, geo;
  double *att;      // Attenuation between segment and mirror
  double *fmie;     // Fraction of Cerenkov beam scattered by aerosols
  double *fray;     // Fraction of Cerenkov bean scattered by molecules
  double pfmie;     // Mie phase function value for segment viewing angle
  double pfray;     // Rayleigh phase function value for segment viewing angle
  double pfcer;     // Ckv beam angular distribution for segment viewing angle
  double hm;
  ObservedTrackSegment *sgmt;
  TrackSegment *tsgmt;

  fdatmos_param_dst_common *aparam = newInstanceOf(fdatmos_param_dst_common);
  fdscat_dst_common *atrans = newInstanceOf(fdscat_dst_common);
  getAtmosphere(par, aparam, atrans);

  /* Main loop over segments */
  for ( i=0; i<track->nseg; i++ ) {
    sgmt  = &obtrack->segment[i];
    tsgmt = &track->segment[i];

    /* Get age and scattering fractions */
    s = tsgmt->age;
    hm = sgmt->altmid;
    fmie = sgmt->fmie;
    fray = sgmt->fray;

    /* Loop over mirrors/cameras at site */
    for ( j=0; j<sgmt->nmir; j++ ) {

      /* Get attenuation between segment and mirror, and geometric factor */
      att = sgmt->att[j];
      geo = sgmt->geoFactor[j];

      /* First calculate production of light */
      for ( k=0; k<NWAVELEN_BANDS; k++ ) {
      sgmt->nfl[j][k] = tsgmt->nfl[k];
      if ( i < track->nseg-1 ) {
        /*
         * Direct is what's created in the segment.
         * This already takes into account scattered Cv light.
         * Since pcv is what's generated in the previous segment,
         * we need to access the next segment's previous
         * production to find what was produced here
         */
        sgmt->ncvdir[j][k] = track->segment[i+1].pcv[k];

        /* Scattered is from average of beam in segment */
        sgmt->ncvmie[j][k] = (track->segment[i+1].ncv[k] +
                        tsgmt->ncv[k])/2.;
        sgmt->ncvray[j][k] = sgmt->ncvmie[j][k];
      }
      else {
        /* This may need to be computed eventually */
        sgmt->ncvdir[j][k] = 0.;
        sgmt->ncvmie[j][k] = 0.;
        sgmt->ncvray[j][k] = 0.;
      }

      /* Scattered cerenkov as fraction of beam */
      sgmt->ncvmie[j][k] *= fmie[k]*(1.-fray[k]);
      sgmt->ncvray[j][k] *= fray[k]*(1.-fmie[k]);
      }

      /*
       * Multiply by phase functions: 4Pi/N(?) dN/dOmega (theta) 
       * and calculate what reaches the mirror.
       */
      pfmie = 4.0*M_PI*geo * mieScatterFunction(sgmt->qview);
      pfray = 4.0*M_PI*geo * rayleighScatterFunction(sgmt->qview);
      pfcer = 4.0*M_PI*geo * cvPhaseFunction(aparam, s, hm, sgmt->qview);

      for ( k=0; k<NWAVELEN_BANDS; k++ ) {
      sgmt->nfl[j][k]    *= att[k] * geo;
      sgmt->ncvdir[j][k] *= att[k] * pfcer;
      sgmt->ncvmie[j][k] *= att[k] * pfmie;
      sgmt->ncvray[j][k] *= att[k] * pfray;
      }

      double rv[3];
      subvec(obtrack->rimp, track->length*obtrack->uv, rv);
      addvec(rv, tsgmt->position*obtrack->uv, rv);

/* #ifdef TESTPRINT */
/* #define DRB_PRINT */
/* //#define DRBFY_PRINT */
/*       double cvbeam = 0.; */
/*       double fy = 0.; */
/*       double fluxNA = 0.; */
/*       double attAve = 0.; */
/*       double at1Ave = 0.; */
/*       double at2Ave = 0.; */
/*       double fluxPh = 0.; */
/*       double fluxPE = 0.; */
/*       double mirPE = 0.; */
/*       TrackSegment *ntsgmt = (i<track->nseg-1)?track->segment[i+1]:track->segment[i]; */
/*       double cky = 0.; */
/*       double ckFluxNA = 0; */
/*       double ckFluxPh = 0; */
/*       double ckFluxPE = 0; */
/*       double fyAlt = 0.; */
/*       double fyAlt2 = 0.; */
/*       double dnfl[NWAVELEN_BANDS],dnfl2[NWAVELEN_BANDS]; */
/*       getFluorescenceYield(aparam, tsgmt->height, dnfl); */
/*       getFluorescenceYieldByMeter(aparam, tsgmt->age, tsgmt->height, dnfl2); */
/*       // Attenutation parts calculation */
/*       double J[3][3]; */
/*       J[2][0] = cos(CLF_LONGITUDE) * cos(CLF_LATITUDE); */
/*       J[2][1] = sin(CLF_LONGITUDE) * cos(CLF_LATITUDE); */
/*       J[2][2] = sin(CLF_LATITUDE); */
/*       zmatrix(J[2], J); */
/*       matrixInverse(J, J); */
/*       double vclf[3]; lla2r(CLF_LATITUDE, CLF_LONGITUDE, CLF_ALTITUDE, vclf); */
/*       double smv_ecc[3],smv[3]; */
/*       subvec(as->impactv,(track->length-(tsgmt->position+tsgmt->dlmid))*as->trackuv,smv); */
/*       applyRotation(J,smv,smv_ecc); */
/*       addvec(smv_ecc,vclf,smv_ecc); */
/*       double mv_ecc[3]; */
/*       applyRotation(site->site2earth,site->local_vmir[j],mv_ecc); */
/*       addvec(mv_ecc,site->vsite,mv_ecc); */
/*       // double att[NWAVELEN_BANDS]; */
/*       // getAttenuation(aparam, smv_ecc, mv_ecc, att); */
/*       double at1[NWAVELEN_BANDS]; */
/*       double at2[NWAVELEN_BANDS]; */
/*       // From getAttenuation */
/*       double dr[3],dl; */
/*       subvec(smv_ecc,mv_ecc,dr); */
/*       dl = sqrt(dotprod(dr,dr)); */
/*       double lat,lon,ha,hb; */
/*       r2lla(smv_ecc, &lat,&lon,&ha); */
/*       r2lla(mv_ecc, &lat,&lon,&hb); */
/*       double dx = getSlantDepth(aparam, smv_ecc, mv_ecc); */
/*       double atm[NWAVELEN_BANDS]; */
/*       double ato[NWAVELEN_BANDS]; */
/*       double atr[NWAVELEN_BANDS]; */
/*       getMieScatterFraction(ha, hb, dl, atm); */
/*       getOzoneAbsorbtionFraction(ha, hb, dl, ato); */
/*       getRayleighScatterFraction(dx, atr); */
/*       for ( k=0; k<NWAVELEN_BANDS; k++ ) { */
/*    at1[k] = (1.-atm[k]); */
/*    at2[k] = (1.-ato[k])*(1.-atr[k]); */
/*       } */
/*       // End of attenuation parts calculation */
/*       for ( k=0; k<NWAVELEN_BANDS; k++ ) { */
/*    cvbeam += tsgmt->ncv[k]; */
/*    fy +=     tsgmt->nfl[k]; */
/*    fluxNA += sgmt->nfl[j][k]/sgmt->att[j][k]; */
/*    // attAve += sgmt->att[j][k]; */
/*    // attAve += sgmt->att[j][k]*(sgmt->nfl[j][k]/sgmt->att[j][k]); */
/*    at1Ave += at1[k]*(sgmt->nfl[j][k]/sgmt->att[j][k]); */
/*    at2Ave += at2[k]*(sgmt->nfl[j][k]/sgmt->att[j][k]); */
/*    fluxPh += sgmt->nfl[j][k]; */
/*    fluxPE += sgmt->nfl[j][k]*0.81*0.306*pmt_qe[k]*bg3_avetrans[k]; */
/*    cky +=    ntsgmt->pcv[k]; */
/*    ckFluxNA += sgmt->ncvdir[j][k]; */
/*    ckFluxPh += sgmt->ncvdir[j][k]*sgmt->att[j][k]; */
/*    ckFluxPE += sgmt->ncvdir[j][k]*sgmt->att[j][k]*0.81*0.306*pmt_qe[k]*bg3_avetrans[k]; */
/*    fyAlt  += dnfl[k]; */
/*    fyAlt2 += dnfl2[k]; */
/*       } */
/*       GaisserHillasParameters ghp; */
/*       fy /= tsgmt->nch*tsgmt->dlseg; */
/*       double dedx = tsgmt->dedep/tsgmt->nch / 1.e6; */
/*       double areamr = M_PI * site->diameter * site->diameter / 4.0; */
/*       double geoFact = geo/areamr*(tsgmt->dlseg/sgmt->dtheta[j]); */
/*       double fluxGeoFact = geoFact/(180./M_PI); */
/*       // attAve /= NWAVELEN_BANDS; */
/*       // attAve /= fluxNA; */
/*       attAve = fluxPh/fluxNA; */
/*       at1Ave /= fluxNA; */
/*       at2Ave /= fluxNA; */
/*       fluxNA /= areamr*sgmt->dtheta[j]*(180./M_PI); */
/*       fluxPh /= areamr*sgmt->dtheta[j]*(180./M_PI); */
/*       fluxPE /= areamr*sgmt->dtheta[j]*(180./M_PI); */
/*       mirPE = fluxPE*areamr*sgmt->dtheta[j]*(180./M_PI); */
/*       double ckPhF = ckFluxNA*((r*r)/areamr)/cky; */
/*       double ckdNdQ = ckPhF*(2*3.14159*sin(sgmt->qview)); */
/*       cky /= tsgmt->nch*tsgmt->dlseg; */
/*       ckFluxNA /= areamr*sgmt->dtheta[j]*(180./M_PI); */
/*       ckFluxPh /= areamr*sgmt->dtheta[j]*(180./M_PI); */
/*       ckFluxPE /= areamr*sgmt->dtheta[j]*(180./M_PI); */
/*       fyAlt *= getAlphaEff(tsgmt->age) * 100.*getDensityByAltitude(aparam, tsgmt->height); */
/*       fyAlt2 *= getAlphaEff(tsgmt->age)/1.e6 / (5.06/2.08); // Compare to dE/dX of FLASH at 28 GeV */
/*       if (fluxPh>0.) { */
/* # ifdef DRB_PRINT */
/*    printf("drb " */
/*           "%6.1f " */
/*           "%4.1f %6.3f %6.3f " */
/*           "%6.2f %6.1f " */
/*           "%4.2f %4.2f " */
/*           "%5.1f %5.0f " */
/*           "%5.3f %5.3f %5.3f %5.0f %5.0f " */
/*           "%5.2f %5.0f " */
/*           "%4.2f %8.2e %8.2e %6.2f %6.2f %6.3f\n", */
/*           tsgmt->sdepth, */
/*           sgmt->qview*57.296, (sgmt->altmid-CLF_ALTITUDE)/1000., r/1000., */
/*           tsgmt->nch/1.e9, cvbeam/1.e12, */
/*           fy, dedx, */
/*           1./geoFact/1000., fluxNA,  */
/*           attAve, at1Ave, at2Ave, fluxPh, fluxPE, */
/*           areamr*sgmt->dtheta[j]*(180./M_PI),mirPE, */
/*           cky,ckPhF,ckdNdQ,ckFluxNA,ckFluxPh,ckFluxPE); */
/* # endif */
/* # ifdef DRBFY_PRINT */
/*    printf("  drbfy %6.1f %4.1f %6.3f %6.3f %6.2f: %6.3f?=%6.3f?=%6.3f\n", */
/*           tsgmt->sdepth, sgmt->qview*57.296, (sgmt->altmid-CLF_ALTITUDE)/1000., r/1000., */
/*           tsgmt->nch/1.e9, fy,fyAlt,fyAlt2); */
/* # endif */
/*       } */
/* #endif */
    }
  }

  free(aparam);
  free(atrans);
}

void clearObservedTrack(ObservedTrack *obt) {
  obt->nseg = 0;
  free(obt->segment);
}

void getFluxInStepsOfX(const Track *trk, const ObservedTrack *in,
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

  for (jseg=0; jseg<out->nseg; jseg++) {     // Looping over data grammages
    yseg = &out->segment[jseg];
    y1 = xtop[jseg];
    y2 = xbot[jseg];

    area = 0.;

    for (iseg=0; iseg<in->nseg-1; iseg++) {  // Looping over MC grammages
      xseg = &in->segment[iseg];
      x1 = trk->segment[iseg].sdepth;
      x2 = trk->segment[iseg+1].sdepth;

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

      //    printf("%d x1 %f x2 %f | %d y1 %f y2 %f | a %f : f %f\n", iseg, x1, x2, jseg, y1, y2, area, f);

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

/*
 *  getRpVector(geofd_dst_common, AirShowerGeometry, double rpv)
 *
 *  Author: Sean R. Stratton
 *          Rutgers University, Dept. of Physics & Astronomy
 *
 *  Computes the Rp vector, i.e. the vector pointing from the detector to the 
 *    closest point on the shower track, for a shower with the given geometry, 
 *    with respect to the given detector site.
 *
 *  Input:
 *    fdsg -- geofd dst bank.
 *    as -- EAS's impact point and track propogation direction
 *
 *  Output:
 *    rpv -- xyz vector of Rp in m in detector's coordinates
 *
 *  Returns:
 *    The magnitude of Rp.
 */
double getRpVector(geofd_dst_common *fdsg, const AirShower *as, double *rpv) {
  double rcosa, R[3][3], rp;
  ObservedTrack obt;

  getObservedTrackOrigin(fdsg, as, &obt);

  /* get shower impact point in terms of rp origin */
  subvec(as->impactv, obt.local_vsite, rpv);

  /* now find Rp vector */
  rcosa = dotprod(rpv, as->trackuv);
  subvec(rpv, rcosa*as->trackuv, rpv);

  /* rotate Rp vector into site's coordinate system */
  matrixTranspose(fdsg->site2clf, R);
  applyRotation(R, rpv, rpv);

  rp = magnitude(rpv);

  if (as->trackuv[2] > 0.)
    rp = -rp;

  return rp;
}

/*
 *  getPsiAngle(geofd_dst_common, AirShowerGeometry)
 *
 *  Author: Sean R. Stratton
 *          Rutgers University, Dept. of Physics & Astronomoy
 *
 *  Computes the Angle between the Shower track and the horizontal plane along 
 *    the SD plane (Psi angle) at the given detector site.
 *
 *  Input:
 *    fdsg -- Site's location and orientation wrt. the CLF
 *    as -- The EAS geometry (impact point and track direction)
 *
 *  Returns:
 *    Psi angle in radians.
 */
double getPsiAngle(geofd_dst_common *fdsg, AirShower *as) {
  double rx, ry, rr, cpsi, psi;
  double rcosa, rpv[3], uv[3];
  double R[3][3];
  ObservedTrack obt;

  getObservedTrackOrigin(fdsg, as, &obt);

  /* get shower impact point in terms of rp origin */
  subvec(as->impactv, obt.local_vsite, rpv);

  /* now find Rp vector */
  rcosa = dotprod(rpv, as->trackuv);
  subvec(rpv, rcosa*as->trackuv, rpv);

  /* rotate Rp vector into site's coordinate system */
  matrixTranspose(fdsg->site2clf, R);
  applyRotation(R, rpv, rpv);
  applyRotation(R, as->trackuv, uv);

  /* trace track from Rp to where it crosses z=0 plane */
  rx = rpv[0] - rpv[2]*uv[0]/uv[2];
  ry = rpv[1] - rpv[2]*uv[1]/uv[2];
  rr = sqrt( rx*rx + ry*ry );

  cpsi = ( rx*uv[0] + ry*uv[1] ) / rr;
  psi = acos(cpsi);

  if (uv[2] > 0.)
    psi = -psi;

  return psi;
}

/*
 *  getZenithAngle(geofd_dst_common, AirShowerGeometry)
 *
 *  Author: Sean R. Stratton
 *          Rutgers University, Dept. of Physics & Astronomy
 *
 *  Computes the Zenith angle of the shower track as seen by the given 
 *    detector site.
 *
 *  Input:
 *    fdsg -- geofd dst bank.
 *    as -- EAS's impact location and track propogation direction.
 *
 *  Returns:
 *    Shower zenith angle theta, as viewed by the given site.
 */
double getZenithAngle(geofd_dst_common *fdsg, AirShower *as) {
  double uv[3], R[3][3];

  matrixTranspose(fdsg->site2clf, R);
  applyRotation(R, as->trackuv, uv);

  return acos( -uv[2] );
}

void getSDPNVector(double rpuv[3], double uv[3], double n[3]) {
  crossProduct(rpuv, uv, n);
}

void getViewingAngles(geofd_dst_common *fdsg, const AirShower *as, 
                  const Track *track, ObservedTrack *obtrack,
                  double *viewang, double *chiang) {
  int i;
  double dl, x[3], ru[3], s[3], su[3], R[3][3];

  matrixTranspose(fdsg->site2clf, R);

  /* 'x' is impact coordinate wrt FD site */
  subvec(as->impactv, obtrack->local_vsite, x);

  /* find distance along track to ground level at site */
  dl = dotprod(x, R[2]) / dotprod(as->trackuv, R[2]);

  /* 'ru' is minus the unit vector pointing to the intersection of the shower 
   *   track with the ground level plane. */
  subvec(dl*as->trackuv, x, ru);
  unitVector(ru, ru);

  /* 's' is minus the vector pointing to the shower segment.  'su' is the 
   *   corresponding unit vector. */
  subvec(track->length*as->trackuv, x, s);

  for ( i=0; i<track->nseg; i++ ) {
    unitVector(s, su);

    viewang[i] = acos( dotprod(su, as->trackuv) ); /* = acos( (s*uv)/|s| ) */
    chiang[i] = acos( dotprod(su, ru) );           /* = acos( (s*ru)/|s| ) */

    subvec(s, track->segment[i].dlseg*as->trackuv, s);
  }
}

#if 0
int getAngularSegments(const geofd_dst_common *fdsg, const AirShower *as, 
                   const Track *track, ObservedTrack *obtrack,
                   double *dlseg) {
  int i, nseg;
  double s[3], su[3], x[3], xu[3];
  double dq;

  /* 'x' is impact coordinate wrt rp origin */
  subvec(as->impactv, obtrack->local_vsite, x);
  unitVector(x, xu);

  /* start 's' at the "top" of the air shower */
  subvec(x, track->length*as->trackuv, s);
  unitVector(s, su);

  dq = acos( dotprod(su, xu) );
  nseg = (int)( dq / HALF_DEGREE ) + 1;

  for ( i=0; i<nseg; i++ ) {
    dq = acos(dotprod(su, as->trackuv)) - HALF_DEGREE;
    dlseg[i] = sqrt(dotprod(s, s)) * sin(HALF_DEGREE) / sin(dq);

    addvec(s, dlseg[i]*as->trackuv, s);
    unitVector(s, su);
  }

  return nseg;
}
#endif

int applyEvsRpCut(double logE, double Rp) {
  double rpkm = Rp/1000.;
  /* This version has linear segments 
  if (logE < 18.5)
    // return (logE > 16.+rpkm*(18.5-16.)/40.)? 1 : 0;
    return (logE > 16.+rpkm*0.0625)? 1 : 0;
  else if (logE < 19.0) 
    // return (logE > 18.5+(rpkm-40.)*(19.0-18.5)/(60.-40.))? 1 : 0;
    return (logE > 18.5+(rpkm-40.)*0.025)? 1 : 0;
  */

  /* This version uses a sqrt
  if (rpkm <= 5.)
    return 1;
  else
    return (logE > 16.+sqrt(rpkm-5)/2.3)? 1 : 0;
  */
    
  if (rpkm <= 7.5)
    return TRUE;
  else
    return ( logE > (16.+sqrt(rpkm - 7.5) / 2.5) ) ? TRUE : FALSE;

}

void getObservedTrackOrigin(geofd_dst_common *fdsg, const AirShower *as,
                      ObservedTrack *obtrack) {

  double rcosa, R[3][3], rpv[3], rpuv[3], uv[3], npln[3];

  cpyvec(as->trackuv, uv);

  /* get shower impact point in terms of site's origin */
  subvec(as->impactv, fdsg->local_vsite, rpv);

  /* now find Rp vector */
  rcosa = dotprod(rpv, uv);
  subvec(rpv, rcosa*uv, rpv);

  /* get rp unit vector and find SDP normal */
  unitVector(rpv, rpuv);
  getSDPNVector(rpuv, uv, npln);

  /* rotate npln into site's coordinate system */
  matrixTranspose(fdsg->site2clf, R);
  applyRotation(R, npln, npln);

  /*
   * Find mirror which "sees" the event best (the mirror whose pointing
   * direction is closest to perpendicular to the SDPn).
   * This will be the origin of Rp.
   */
  getRpOrigin(fdsg, npln, obtrack);
}

void getRpOrigin(geofd_dst_common *fdsg, double *npln,
             ObservedTrack *obtrack) {
  int i, mir, bestmir;
  double cosang, mincosang;
  double siteOrigin[3];

  mincosang = M_PI/2.;
  bestmir = -1;
  for (mir=0; mir<fdsg->nmir; mir++) {
    cosang = fabs( dotProduct(fdsg->vmir[mir], npln) );
    if (cosang < mincosang) {
      mincosang = cosang;
      bestmir = mir;
    }
  }

  /* get rp origin from CLF in CLF coordinate system */
  if (bestmir >= 0) {
    for (i=0; i<3; i++)
      siteOrigin[i] = fdsg->local_vmir[bestmir][i];
    applyRotation(fdsg->site2clf, siteOrigin, obtrack->local_vsite);
    addvec(obtrack->local_vsite, fdsg->local_vsite, obtrack->local_vsite);
  }
  else {
    for (i=0; i<3; i++) {
      siteOrigin[i] = 0.;
      obtrack->local_vsite[i] = fdsg->local_vsite[i];
    }
  }

  /* Get Rp origin in ECC */
  applyRotation(fdsg->site2earth, siteOrigin, obtrack->vsite);
  addvec(obtrack->vsite, fdsg->vsite, obtrack->vsite);

  /* Get mirror locations wrt Rp origin */
  for (i=0; i<fdsg->nmir; i++)
    subvec(fdsg->local_vmir[i], siteOrigin, obtrack->local_vmir[i]);
}
