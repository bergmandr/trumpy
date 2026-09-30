#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "event.h"
#include "random.h"
#include "constants.h"
#include "control.h"
#include "toolbox.h"

#include "atmosphere-withdb.h"
#include "airshower.h"
#include "showerlib.h"
#include "eventlist.h"
#include "track.h"

#include "fdconstants.h"
#include "fdsite.h"
#include "calibration.h"
#include "tacalibration.h"
#include "raytrace.h"
#include "acpttrack.h"
#include "taelectronics.h"
//#include "pelist.h"

/* #include "histogram.h" */

#include "bank-manager.h"

/*
 *  Removes extra cameras & tubes from a fdraw bank.
 *
 *  Author: Sean R. Stratton
 *          Rutgers University, Dept. of Physics & Astronomy
 *
 *  No arguments or return value.
 */
void compactFDRawBank(fdraw_dst_common* fdraw) {
  short ntube;
  int i, j, k;
  // This caused problems in cygwin
  // fdraw_dst_common fdtemp;
  fdraw_dst_common* fdtemp = newInstanceOf(fdraw_dst_common);

  *fdtemp = *fdraw;

  for ( i=0; i<fdraw->num_mir; i++ ) {
    fdraw->num_chan[i] = 0;
    for ( j=0; j<fdraw_nchan_mir; j++ )
      fdraw->hit_pt[i][j] = fdtemp->hit_pt[i][j];

    ntube = 0;
    for ( j=0; j<fdraw_nchan_mir; j++ ) {
      int wf[fdraw_nt_chan_max];
      for ( k=0; k<fdraw_nt_chan_max; k++ )
        wf[k] = (int)fdtemp->m_fadc[i][j][k];

      double cl;
      if ( fdtemp->channel[i][j] >= 0 )
        cl = fadcSig(fdtemp->mean[i][j][0], fdtemp->disp[i][j][0], 
              fdraw_nt_chan_max, wf);
      else
        cl = 0.00;

/* #ifdef HISTOGRAM_MODE */
/*    double dmean = 0.00; */
/*    double dvari = 0.00; */
/*    for ( k=0; k<512; k++ ) { */
/*      dmean += (double)fdtemp->m_fadc[i][j][k]; */
/*      dvari += (double)fdtemp->m_fadc[i][j][k] *  */
/*        (double)fdtemp->m_fadc[i][j][k]; */
/*    } */
/*    dmean = dmean/512.0; */
/*    dvari = dvari/512.0 - dmean*dmean; */

/*    dmean -= (double)fdtemp->mean[i][j][0]/16.0; */
/*    dvari -= (double)fdtemp->disp[i][j][0]/16.0; */

/*    hf1(211, (float)cl, 1.0); */
/*    hf1(212, (float)dmean, 1.0); */
/*    hf1(213, (float)dvari, 1.0); */
/* #endif */
      if ( cl >= MIN_READOUT_CL ) {
        fdraw->channel[i][ntube] = fdtemp->channel[i][j];
        for ( k=0; k<4; k++ ) {
          fdraw->mean[i][ntube][k] = fdtemp->mean[i][j][k];
          fdraw->disp[i][ntube][k] = fdtemp->disp[i][j][k];
        }

        for ( k=0; k<fdraw_nt_chan_max; k++ )
          fdraw->m_fadc[i][ntube][k] = fdtemp->m_fadc[i][j][k];

        ntube++;
/* #ifdef HISTOGRAM_MODE */
/*      hf1(221, (float)cl, 1.0); */
/*      hf1(222, (float)dmean, 1.0); */
/*      hf1(223, (float)dvari, 1.0); */
/* #endif */
      }
/* #ifdef HISTOGRAM_MODE */
/*    else { */
/*      hf1(231, (float)cl, 1.0); */
/*      hf1(232, (float)dmean, 1.0); */
/*      hf1(233, (float)dvari, 1.0); */
/*    } */
/* #endif */
    }

    fdraw->num_chan[i] = ntube;
  }

  free(fdtemp);
}
void clearTrumpMCBank(trumpmc_dst_common *tmc) {
  int i,j,k,m;
  
  tmc->nDepths = 0;
  for (m=0; m<TRUMPMC_MAXDEPTHS; m++)
    tmc->depth[m] = 0;
    
  tmc->nSites = 0;
  for (i=0; i<TRUMPMC_MAXSITES; i++) {
    tmc->siteid[i] = 0;
    tmc->nMirrors[i] = 0;

    for (j=0; j<TRUMPMC_MAXMIRRORS; j++) {
      tmc->mirror[i][j] = 0;
      tmc->totalNPEMirror[i][j] = 0;
      tmc->nTubes[i][j] = 0;
      for (m=0; m<TRUMPMC_MAXDEPTHS; m++) {
        tmc->fluoFlux[i][j][m] = 0;
        tmc->aeroFlux[i][j][m] = 0;
        tmc->raylFlux[i][j][m] = 0;
        tmc->dirCFlux[i][j][m] = 0;
      }
      
      for (k=0; k<TRUMPMC_MAXTUBES; k++) {
        tmc->tube[i][j][k] = 0;
        tmc->aveTime[i][j][k] = 0;
        tmc->totalNPE[i][j][k] = 0;
      }
    }
  }
}

void buildTrumpMCBank(int siteid, const AirShower* as, 
                  const geofd_dst_common* geo,
                  const GaisserHillasParameters* gh, const Track* t, 
                  const ObservedTrack* ot, const PETimes* pet,  
                  const TACalibration* calib, const RuntimeParameters *par, trumpmc_dst_common* tmc) {
  int i, j, k, m, n, icam; 
  FDSiteGeometry *fdsg = newInstanceOf(FDSiteGeometry);
  *fdsg = *geo;

  double areamr[GEOFD_MAXMIR];
  for (i=0; i<fdsg->nmir; i++)
    areamr[i] = M_PI * fdsg->diameters[i] * fdsg->diameters[i] / 4.0;

  tmc->julian = par->jday;
  tmc->jsec = par->jsec;
  tmc->nano = par->nsec;
  for (i=0; i<3; i++) {
    tmc->impactPoint[i] = (real4)as->impactv[i];
    tmc->showerVector[i] = (real4)as->trackuv[i];
  }
  tmc->energy = (real4)pow(10.,as->loge);
#ifdef LASER_OPTION
  tmc->primary = 0;
#else
  tmc->primary = as->species;
#endif
  tmc->ghParm[0] = (real4)gh->x0;
  tmc->ghParm[1] = (real4)gh->xmax;
  tmc->ghParm[2] = (real4)gh->nmax;
  tmc->ghParm[3] = (real4)gh->lambda;

  /* Only one site at moment (no loop yet) */
  tmc->nSites = 1;
  tmc->siteid[0] = siteid;
  for (i=0; i<3; i++)
    tmc->siteLocation[0][i] = (real4)fdsg->local_vsite[i];
  tmc->psi[0] = (real4)ot->psi;
  double rp[3];
  getRpVector(fdsg, as, &rp[0]);
  for (i=0; i<3; i++)
    tmc->rp[0][i] = (real4)rp[i];

  /* Find depths to use */
  int idMod = (t->nseg-1)/TRUMPMC_MAXDEPTHS + 1;
  tmc->nDepths = t->nseg/idMod + 1;
  for (i=0, j=0; i<t->nseg; i+=idMod, j++)
#ifdef LASER_OPTION
    tmc->depth[j] = (real4)( t->segment[i].height);
#else
    tmc->depth[j] = (real4)( t->segment[i].sdepth + t->segment[i].dxseg/2. );
#endif
//     tmc->depth[j] = t->segment[i].height;
  /* Find mirrors which viewed track */
  int nMirror = 0;
  int mirror[TRUMPMC_MAXMIRRORS] = {0};
  for (i=0; i<ot->nseg; i++) {
    for (j=0; j<ot->segment[i].nmir; j++) {
      icam = addListItem(ot->segment[i].mir[j], &nMirror,
                   TRUMPMC_MAXMIRRORS, mirror);
      if (icam < 0) {
        vperr("Maximum number of mirrors reached!\n");
        abort();
      }
      //      pout("mirror %2d sees segment %4d\n", icam, j);
    }
  }

  tmc->nMirrors[0] = nMirror;
  for (i=0; i<nMirror; i++) {
    tmc->mirror[0][i] = mirror[i];
    for (j=0; j<TRUMPMC_MAXDEPTHS; j++) {
      tmc->fluoFlux[0][i][j] = 0.;
      tmc->aeroFlux[0][i][j] = 0.;
      tmc->raylFlux[0][i][j] = 0.;
      tmc->dirCFlux[0][i][j] = 0.;
    }
  }

  /* Find fluxes */
// //   for (i=0; i<12; i++)
// //     fprintf(stderr,"mirror %d calibration reduction (photons / pe) at 355nm %f\n",i,calib->reductFactor[i][21]);
  for (i=0, j=0; i<t->nseg; i+=idMod, j++) {
    ObservedTrackSegment *osg = &ot->segment[i];
    for (k=0; k<osg->nmir; k++) {
      int kk = osg->mir[k];
      for (n=0; n<tmc->nMirrors[0]; n++) {
        if (kk == tmc->mirror[0][n]) {
          for (m=0; m<NWAVELEN_BANDS; m++) {
            if (geo->siteid == BLACK_ROCK_SITEID || geo->siteid == LONG_RIDGE_SITEID) {          
              tmc->fluoFlux[0][n][j] += (real4)( osg->nfl[k][m] * 
                calib->reductFactor[kk][m] );
              tmc->aeroFlux[0][n][j] += (real4)( osg->ncvmie[k][m] *
                calib->reductFactor[kk][m] );
              tmc->raylFlux[0][n][j] += (real4)( osg->ncvray[k][m] *
                calib->reductFactor[kk][m] );
              tmc->dirCFlux[0][n][j] += (real4)( osg->ncvdir[k][m] *
                calib->reductFactor[kk][m] );
            }
            else {
              tmc->fluoFlux[0][n][j] += (real4)( osg->nfl[k][m]);
              tmc->aeroFlux[0][n][j] += (real4)( osg->ncvmie[k][m]);
              tmc->raylFlux[0][n][j] += (real4)( osg->ncvray[k][m]);
              tmc->dirCFlux[0][n][j] += (real4)( osg->ncvdir[k][m]);
            }
          }
          tmc->fluoFlux[0][n][j] /= (real4)(osg->dtheta[k] * areamr[kk]);
          tmc->aeroFlux[0][n][j] /= (real4)(osg->dtheta[k] * areamr[kk]);
          tmc->raylFlux[0][n][j] /= (real4)(osg->dtheta[k] * areamr[kk]);
          tmc->dirCFlux[0][n][j] /= (real4)(osg->dtheta[k] * areamr[kk]);
          break;
        }
      }
    }
  }

  /* Get mirror info */
//   printf("%s: nMirror = %d\n",__FUNCTION__,nMirror);
  for (i=0; i<pet->nmir; i++) {
    int ij = -1;  // TAS: j-index of mirror i
    int tNPE = 0;
    int nTube = 0;
    for (j=0; j<nMirror; j++) {
      // ELB, bug below that causes only one mirror NPE data to be filled
      //      if (tmc->mirror[0][i] == pet->mir[i]) {
      if (tmc->mirror[0][j] == pet->mir[i]) {
        ij = j;
//           printf("%s: filling NPE for mirror %d\n",__FUNCTION__,pet->mir[i]);
        for (k=0; k<GEOFD_MIRTUBE; k++) {
//           printf("%s: tube %d len = %d\n",__FUNCTION__,k,pet->len[i][k]);
          if (pet->len[i][k] > 0) {
            double npe = 0.0;
            double time = 0.0;
            for (m=0; m<pet->len[i][k]; m++) {
              npe  += pet->n[i][k][m];
              time += pet->t[i][k][m]*pet->n[i][k][m];
            }
            time /= npe;
            tmc->tube[0][j][nTube]     = k;
            tmc->aveTime[0][j][nTube]  = (real4)time;
            tmc->totalNPE[0][j][nTube] = (real4)npe;
            tNPE += (int)npe;
            nTube++;
          }
        }
        break;
      }
    }
    // TAS, noticed ELB's indexing bug appears here too.
//     tmc->nTubes[0][i] = nTube;
//     tmc->totalNPEMirror[0][i] = tNPE;
    tmc->nTubes[0][ij] = nTube;
    tmc->totalNPEMirror[0][ij] = tNPE;
  }

  free(fdsg);
}
