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
#include "calibration.h"
#include "tacalibration.h"
#include "raytrace.h"
#include "showerlib.h"
#include "track.h"
#include "fdsite.h"
#include "acpttrack.h"
#include "taelectronics.h"

#define BLOCK_SIZE 1024
#define SATURATION_RATE 8192.0  /* Accepted PE / 100ns */
#define MAX_SATURATED_TUBES 16

#ifndef COS15
#define COS15 0.965925826
#endif

#ifndef COS13
#define COS13 0.974370065
#endif

/*
 *  getAcceptedTrack
 *
 *  Inputs:
 *    const AirShower *as    -- Used for vector to shower impact coordinate
 *    const Track *tk        -- EAS track structure (CLF coordinates)
 *    geofd_dst_common *site -- FD Site geometry
 *    ObservedTrack *obtk    -- EAS as observed by *site
 *    Calibration *calib     -- FD Calibration (used for tube live flags)
 *    TGeom *tgeom           -- Passed directly to 'doRayTrace()'
 *
 *  Output:
 *    PETimes *pet           -- Arrival times of PE indexed by camera & tube
 *
 *
 *  This function computes the EAS observed at an FD site's cameras.  It is the
 *    connection between the shower track as an array of segments in CLF space 
 *    and the electronics simulation.
 *
 */
int getAcceptedTrack(const AirShower *as, const Track *tk, 
                  geofd_dst_common *site, ObservedTrack *obtk, 
                  TACalibration *calib, TGeom *tgeom, PETimes *pet) {

  int i, j, k, m;        /* general use iterators                            */
  int icam;              /* iterator for loop over cameras (also camera ID#) */
  int jcam;              /* iterator for participating cameras               */
  int iseg;              /* iterator for loop over shower segments           */
  int nrtp;              /* number of angular pieces in shower segment       */
  int nphot;             /* number of photons in angular piece               */
  int tube;              /* tube ID number                                   */
  int len;               /* alias for pet->len[][] array                     */
  double *npe;           /* alias for *pet->n[][] array                      */
  double *tpe;           /* alias for *pet->t[][] array                      */
  double C2S[3][3];      /* CLF to site rotation matrix                      */
  double tkuv[3];        /* Shower direction vector in site coordinates      */
  double rtop[3];        /* vector from site to "top" of air shower          */
  double rmid[3];        /* vector from site to middle of shower segment     */
  double rbeg[3];        /* vector from site to beginning of shower segment  */
  double rend[3];        /* vector from site to end of shower segment        */
  double rgen[ACPTTRACK_MAX_NRTP][3];  /* vectors to angular pcs of AS seg   */
  double tgen[ACPTTRACK_MAX_NRTP];  /* gen time of light from angular piece  */
  double tflux;          /* total photon flux at mirror                      */
  TrackSegment *tksgmt;  /* alias for track segment                          */
  ObservedTrackSegment *obsgmt;  /* alias for observed track segment         */
  RayTrace ray;          /* output from doRayTrace() func                    */

  pet->t1 = 0x7FFFFFFF;
  pet->t2 = -pet->t1;

  /* Get CLF-to-site rotation matrix */
  matrixTranspose(site->site2clf, C2S);

  /* Find top of shower (in site coordinates) */
  addvec(as->impactv, -tk->length*as->trackuv, rtop);
  subvec(rtop, obtk->local_vsite, rtop);
  applyRotation(C2S, rtop, rtop);

  /* get shower direction vector in site coordinates */
  cpyvec(as->trackuv, tkuv);   /* preserves "const" on as */
  applyRotation (C2S, tkuv, tkuv);

  /*
   *  Loop over cameras
   */
  int nstubes = 0;
  pet->nmir = 0;
  for ( icam=0; icam<site->nmir; icam++ ) {
    double t1[NTUBES_CAMERA];
    double t2[NTUBES_CAMERA];
    for ( j=0; j<NTUBES_CAMERA; j++ ) {
      pet->n[icam][j] = NULL;
      pet->t[icam][j] = NULL;
      pet->len[icam][j] = 0;

      t1[j] = 1.0e+137;
      t2[j] = -1.0e+137;
    }

    /*
     *  Loop over shower segments.
     */
    for ( iseg=0; iseg<tk->nseg; iseg++ ) {
      tksgmt = &tk->segment[iseg];
      obsgmt = &obtk->segment[iseg];

      /* Find vector to middle of segment (site) */
      addvec(rtop, (tksgmt->position+tksgmt->dlmid)*tkuv, rmid);

      double dp = dotprod(rmid, site->vmir[icam]) / magvec(rmid);
      if ( dp < COS13 )  /* skip segment if not in view */
      continue;

      /* 
       *  Find angular extent of segment, and split into <0.1 deg steps 
       */
      addvec(rtop, tksgmt->position*tkuv, rbeg);
      addvec(rbeg, tksgmt->dlseg   *tkuv, rend);
      double den = magvec(rbeg);
      den *= magvec(rend);
      double dTheta = acos(dotprod(rbeg,rend)/den);
      nrtp = (int)ceil( dTheta / ACPTTRACK_DTH );

      if ( nrtp > 1 ) {
        /* Add protection on too large nrtp */
        nrtp = min(nrtp, ACPTTRACK_MAX_NRTP);

        /* For multiple raytrace split up evenly in X (?) */
        double r[3];
        double dl = tksgmt->dlseg / (double)nrtp;
        double t = obsgmt->tgen - (tksgmt->dlmid-dl/2.)/SPEED_OF_LIGHT;
        addvec(rbeg, (dl/2.)*tkuv, r);
        for ( k=0; k<nrtp; k++ ) {
          tgen[k] = t;
          subvec(r, obtk->local_vmir[icam], rgen[k]);
          t += dl/SPEED_OF_LIGHT;
          addvec(r, dl*tkuv, r);
        }
      }
      else if ( nrtp == 1 ) {
        /* For only one raytrace in a segment use the "mid" (in X) point */
        tgen[0] = obsgmt->tgen;
        addvec(rbeg, tksgmt->dlmid*tkuv, rgen[0]);
      }
      else {
        vperr("Bad nrtp value: %d\n", nrtp);
        abort();
      }

      /*
       * Find the total flux seen by this telescope from this segment.
       */
      tflux = 0.00;
      for ( m=0; m<obsgmt->nmir; m++ ) {
        if ( obsgmt->mir[m] == icam ) {
          for ( k=0; k<NWAVELEN_BANDS; k++ ) {
            double lflux = obsgmt->nfl[m][k] + obsgmt->ncvdir[m][k] +
              obsgmt->ncvmie[m][k] + obsgmt->ncvray[m][k];
            if (site->siteid == BLACK_ROCK_SITEID || site->siteid == LONG_RIDGE_SITEID)  
              lflux *= calib->reductFactor[icam][k];
            tflux += lflux;
          }
          tflux /= (double)nrtp;
          break;
        }
      }

      for ( k=0; k<nrtp; k++ ) {
        if ( tflux < 20.0 ) {
          nphot = prand(tflux);
        }
        else {
          double u = grand();
          nphot = (int)( u*sqrt(tflux) + tflux );
        }

        // pout("cam=%d  seg=%d/%d  vang=%6.2lf  nphot=%d\n", icam, iseg+1, 
        //      tk->nseg, 180.0*acos(dp)/M_PI, nphot);

        if ( nphot <= 0 )  /* continue to next angular segment */
          continue;

        int ii, nn = 1;
        while ( nphot > MAXNPE ) {
          nphot >>= 1;
          nn <<= 1;
        }

        for ( ii=0; ii<nn; ii++ ) {
          doRayTrace(nphot, icam, obtk->rpuv, obtk->npln, rgen[k], 
              tksgmt->age, tksgmt->molrad, site, tgeom, &ray);

          for ( m=0; m<ray.nthrow; m++ ) {
            tube = ray.pmt[m];
            if ( (tube>=0) && (
            (site->siteid==BLACK_ROCK_SITEID || site->siteid==LONG_RIDGE_SITEID)?
            calib->liveflag[icam][tube]:TRUE ) ) {
              double time = tgen[k] + ray.tpe[m];

              t1[tube] = min(t1[tube], time);
              t2[tube] = max(t2[tube], time);

              jcam = addListItem(icam, &pet->nmir, NCAMERAS_SITE, pet->mir);
              if ( pet->len[jcam][tube] == 0 ) {
          pet->t[jcam][tube]=(double *)malloc(BLOCK_SIZE*sizeof(double));
          if ( jcam < 0 ) {
            vperr("Maximum number of mirrors reached!\n");
            abort();
          }
              }
              else if ( (pet->len[jcam][tube]%BLOCK_SIZE) == 0 ) {
                pet->t[jcam][tube] = (double *)realloc(pet->t[jcam][tube], 
              (pet->len[jcam][tube]+BLOCK_SIZE)*sizeof(double));
              }

              len = pet->len[jcam][tube];
              tpe = pet->t[jcam][tube];
              /* this type of insertion sort is inefficient for large data
              * sets (O(^2) with N>1000). implemented a heap sort
              * algorithm below with runs ~500x faster */
      //          for ( j=len; (time<tpe[j-1])&&(j>0); j-- )
      //           tpe[j] = tpe[j-1];
      //          tpe[j] = time;
              tpe[len] = time;

              double perate;
              if ( len > MAXNPE )
                perate = 100.0 * (double)len / ( t2[tube] - t1[tube] );
              else
                perate = 0.00;

              if ( perate < SATURATION_RATE ) {
                pet->len[jcam][tube]++;
              }
              else {
                /* saturation... shut the tube down */
                pout("cam %2d  tube %3d  saturated!\n", icam, tube);
                if (site->siteid == BLACK_ROCK_SITEID || site->siteid == LONG_RIDGE_SITEID)
                  calib->liveflag[icam][tube] = FALSE;
                pet->len[jcam][tube] = 0;
                free(pet->t[jcam][tube]);
                pet->t[jcam][tube] = NULL;
                nstubes++;
                if ( nstubes >= MAX_SATURATED_TUBES ) {
                  pout("ERROR!  light saturation!\n");
                  return 1;
                }
              }

            }  //  if ( (tube>=0) && calib->liveflag[icam][tube] )
          }  //  for ( m=0; m<ray.nthrow; m++ )
        }  //  for ( ii=0; ii<nn; ii++ )
      }  //  for ( k=0; k<nrtp; k++ )
    }  //  for ( iseg=0; iseg<tk->nseg; iseg++ )
  }  //  for ( icam=0; icam<site->nmir; icam++ )

  /*
   *  collapse PE arrays
   */
  for ( jcam=0; jcam<pet->nmir; jcam++ ) {
    for ( tube=0; tube<NTUBES_CAMERA; tube++ ) {
      len = pet->len[jcam][tube];
      if ( len <= 0 ) continue;

      pet->n[jcam][tube] = (double *)malloc(BLOCK_SIZE*sizeof(double));
      pet->n[jcam][tube][0] = 1.0;

      npe = pet->n[jcam][tube];
      tpe = pet->t[jcam][tube];
      
      heapsort(len, tpe); // sort the Tpe array
      
      i = 1;
      while ( i < len ) {
        if ( (i%BLOCK_SIZE) == 0 )
          pet->n[jcam][tube] = (double *)realloc(pet->n[jcam][tube], 
              (i+BLOCK_SIZE)*sizeof(double));

//         npe = pet->n[jcam][tube];
//         tpe = pet->t[jcam][tube];
        if ( fabs(tpe[i]-tpe[i-1]) < DOUBLE_FP_PRECISION ) {
          len--;
          npe[i-1] += 1.0;
          for ( j=i; j<len; j++ )
            tpe[j] = tpe[j+1];
        }
        else {
          npe[i] = 1.0;
          i++;
        }
      }

      pet->len[jcam][tube] = len;
//       pet->t1 = min(pet->t1, pet->t[jcam][tube][0]);
//       pet->t2 = max(pet->t2, pet->t[jcam][tube][len-1]);
      pet->t1 = min(pet->t1, tpe[0]);
      pet->t2 = max(pet->t2, tpe[len-1]);
    }  // for ( tube=0; tube<NTUBES_CAMERA; tube++ )
  }  //  for ( jcam=0; jcam<pet->nmir; jcam++ )

  return 0;
}

void clearPETimes(PETimes* pe) {
  int j,k;
  for (j=0;j<GEOFD_MAXMIR;j++)
    for (k=0;k<GEOFD_MIRTUBE;k++) {
      /*
       * pe->n and pe-> are allocated in fillPETimes from a PEList.
       * They are only allocated when a given camera/tube appear, so
       * need to check before free'ing
       */
      if (pe->n[j][k] != NULL) free(pe->n[j][k]);
      if (pe->t[j][k] != NULL) free(pe->t[j][k]);
      pe->n[j][k] = NULL;
      pe->t[j][k] = NULL;
      pe->len[j][k] = 0;
    }
  pe->t1 = 0.;
  pe->t2 = 0.;
}

/*
void clearAcceptedTrack(AcceptedTrack *at) {
  int i, j;
  for (i=0; i<GEOFD_MAXMIR; i++) {
    for (j=0; j<GEOFD_MIRTUBE; j++) {   
      if (at->len[i][j] > 0) {           
        free(at->npe[i][j]);
        free(at->tpe[i][j]);
      }
    }
  } 
}
*/
