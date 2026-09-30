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

#include "atmosphere-withdb.h"
#include "airshower.h"
#include "showerlib.h"
#include "eventlist.h"
#include "track.h"

#include "fdconstants.h"
#include "adjcamera.h"
#include "fdsite.h"
#include "tacalibration.h"
#include "calibration.h"
#include "raytrace.h"
#include "acpttrack.h"


#include "tlelectronics.h"

// #include "parameters.h"

// #include "petime.h"



// typedef struct {
static double ut1[101], ut2[101], ut3[101], ut4[101];
static double vt1[101], vt2[101];
static double uutn1[101], uvtn1[101], uutn2[101];
static double fu1, fu2, fu3, fu4;
static double fv1, fv2;
static double fnu1, fnu2, fnv1;
static double fv2u1, fv2v1;
static double a1, a2, b1, b2;
static double unorm, vnorm;
// } tltf;

// tltf tltf_;
// static integer4 ifirst = 0;

// taufadc taufadc_ = { .tau_neg = 100000.};

void initTLElectronics(geofd_dst_common *geo) {
  
//   double gain_adc, gain_16, ped_16, ped_adc;
  int i,j,k;
  int mir;
  int idpmt;
//   int 

  i = setTLPmtID();
  printf("Assigned idpmt for %d PMTs\n",i);
//   gain_adc = 1.0;
//   ped_adc = 10.50;
//   ped_16 = 20.50;
//   gain_16 = 1.0;
  
  for (mir=0; mir<GEOTL_MAXMIR; mir++) {
//     gain_16 = (geo->ring[mir]==2)?1.50:1.89;
//     for (j=0; j<GEOFD_MIRTUBE; j++) {
// //       pedfadc_.ped_fadc[mir][j] = ped_adc;
//       tes->ped_fadc[mir][j] = ped_adc;
// //       gainfadc_.gain_fadc[mir][j] = gain_adc;
//       tes->gain_fadc[mir][j] = gain_adc;
//     }
//     for (j=0; j<16; j++) {
// //       pedfadc_.ped_fadc[mir][256+j] = ped_16;
//       tes->ped_fadc[mir][256+j] = ped_16;
// //       pedfadc_.ped_fadc[mir][288+j] = ped_16;
//       tes->ped_fadc[mir][288+j] = ped_16;
// //       gainfadc_.gain_fadc[mir][256+j] = gain_16;
//       tes->gain_fadc[mir][256+j] = gain_16;
// //       gainfadc_.gain_fadc[mir][288+j] = gain_16;
//       tes->gain_fadc[mir][288+j] = gain_16;
// //       pedfadc_.ped_fadc[mir][272+j] = ped_adc;
//       tes->ped_fadc[mir][272+j] = ped_adc;
// //       pedfadc_.ped_fadc[mir][304+j] = ped_adc;
//       tes->ped_fadc[mir][304+j] = ped_adc;
// //       gainfadc_.gain_fadc[mir][272+j] = gain_adc/9.8;
//       tes->gain_fadc[mir][272+j] = gain_adc/9.8;
// //       gainfadc_.gain_fadc[mir][304+j] = gain_adc/9.8;
//       tes->gain_fadc[mir][304+j] = gain_adc/9.8;
//     }
    
    for (j=0; j<256; j++) {

      idpmt = mirrs_.mirtub[mir][j] - 1;
      for (k=0; k<3; k++)
        centre_.tbv[idpmt][k] = (int)(1.e9*geo->vtube[mir][j][k] + 0.5);
//       printf("assigning vectors... tbv[%d] = %d %d %d\n",idpmt,
//              centre_.tbv[idpmt][0],centre_.tbv[idpmt][1],centre_.tbv[idpmt][2]);
//              abort();
    }    
  }

#ifndef NO_NOISE
  init_noise_();
#endif
  
//   taufadc_.tau_fadc = 35.;
//   taufadc_.tau_16 = 500.;
//   taufadc_.tau_neg = 500000.;
  
  initTLWaveformGenerator();
  
//   float v16[max_tim] = {0.}, v[max_tim]={0.}, tpe[10000]={0.}, tstart;
//   int ntpe, incr_npe[10000]={0}, nt;
//   
//   tstart=1000.;
//   tpe[0] = 1901.;
//   incr_npe[0] = 1;
//   nt = 400;
//   ntpe = 1;
//   glb_prm_.xmu_noise = 6.6;
//   printf("Test waveform generation for 1 pe\n");
//   generateTLWaveform(&ntpe, tpe, incr_npe, &tstart, &nt, v, v16);
// 
//   for (i=0; i<nt; i++) printf("%f ",v[i]);
//   printf("\n");
//   abort();
  
  
  
  return;
}





numcon numcon_ = { 
  .pi = PI,
  .radian = R2D,
  .root2 = 1.4142135,
  .root3 = 1.7320508,
  .twopi = 6.283185307,
  .root_2pi = 2.5066282,
  .halfpi = 1.5707963268,
  .cinv = .333564e4
};




pemir pemir_;


// whicheye whicheye_ = { .neye = 1};
mir_trig mir_trig_;
rawcheck rawcheck_;
// eye_hit eye_hit_;
// db_info db_info_ /*= { .db_trig_version = 3, .db_adj_opt = 2}*/;
// prtlev prtlev_;



rawcheck rawcheck_;
jim jim_;
prescale prescale_ = { .nsingle = 0} ;
// glb_prm glb_prm_ /*= { .thresh_1 = 32.0, .thresh_2 = 22.0}*/;
mir_hit mir_hit_;
idtemp idtemp_;
mirrs mirrs_;
petime petime_;
// control control_ = { .inoise_mir = 1};
// pedfadc pedfadc_;
filter_m filter_m_;
// trig_def trig_def_;
centre centre_;
// calib calib_ = { .calib_corr = 1.0} ;
petemp petemp_;
// gainfadc gainfadc_;

int setTLPmtID(void) {
//   printf("running setTLPmtID\n");
  int nt = 1;
  int i,j;
  for (i=0; i<GEOTL_MAXMIR; i++)
    for (j=0; j<GEOFD_MIRTUBE; j++)
      mirrs_.mirtub[i][j] = nt++;
//       printf("assigned mirtub[%d][%d][%d]=%d\n"
    
  return --nt;
}



double convertPETimesToTLPE(PETimes *pt) {
//   printf("%s::%s (%d): inside\n",__FILE__,__FUNCTION__,__LINE__);
  int mir, pmt, imir;
  int len,tlen,it,tnpe;
  int it0mir,it1mir,itav;
  double t0, t1;
  double mir_cur_max, t_cur_max;
  
//   double omni_t0 = 1.e9;    // earliest arrival time of first PE anywhere
  
  double offset;
  int us;
  
  
  pemir_.nmir_view = pt->nmir;
  
  idtemp_.ntba = 0;
  
 
  us = (int)(pt->t1-15000.)/1000;
  offset = (double)us * 1000.;
  
  for (mir=0; mir<pt->nmir; mir++) {
    

    imir = pt->mir[mir];
    pemir_.imir_view[mir] = imir + 1;
    
    len = 0;
    t0 = 1.e9;
    t1 = -t0;
    mir_cur_max=0.;
    t_cur_max=0.;
    for (pmt=0; pmt<GEOFD_MIRTUBE; pmt++) {
      tlen = pt->len[mir][pmt];
      len += tlen;
      if (tlen == 0)
        continue;
//       printf("len = %d; t0 = %f\n",len,pt->t[mir][pmt][0]);
      t0 = min(t0, pt->t[mir][pmt][0] - offset);
//       printf("t0 %f\n",t0);
      t1 = max(t1, pt->t[mir][pmt][tlen-1] - offset);
      tnpe = 0;
      for (it=0; it<tlen; it++) {
        tnpe += pt->n[mir][pmt][it];
        petemp_.tpe_temp[idtemp_.ntba][it] = pt->t[mir][pmt][it] - offset;
        petemp_.incr_npe_temp[idtemp_.ntba][it] = pt->n[mir][pmt][it];
        if (pt->n[mir][pmt][it] > mir_cur_max) {
          mir_cur_max = pt->n[mir][pmt][it];
          t_cur_max = pt->t[mir][pmt][it] - offset;
        }
      }
      
      idtemp_.itba[idtemp_.ntba] = mirrs_.mirtub[imir][pmt];
      petime_.npe_pmt[idtemp_.ntba] = tnpe;
      petemp_.ntpe_temp[idtemp_.ntba] = tlen;

//       if (imir == 6 && pmt == 248) {
//         printf("convPE: mr.m %d id.ntba %d id.itba %d npe %d tlen %d\n",
//                mirrs_.mirtub[imir][pmt],idtemp_.ntba,idtemp_.itba[idtemp_.ntba],
//                petime_.npe_pmt[idtemp_.ntba],petemp_.ntpe_temp[idtemp_.ntba]);
//       }
      
      idtemp_.ntba++;
    }
    
    // t0 is the time since event-t0; let's offset it to time within the microsecond
//     omni_t0 = min(omni_t0, t0);
    
    pemir_.pemir[imir] = (double)len;
    pemir_.t0_mir[imir] = t0 - 150. * DT_ADC;
    pemir_.t1_mir[imir] = t1 + 150. * DT_ADC;
    it0mir = (int)(pemir_.t0_mir[imir]/DT_ADC);
    it1mir = (int)(pemir_.t1_mir[imir]/DT_ADC) + 1;
    pemir_.nt_mir[imir] = it1mir - it0mir + 1;
    if (pemir_.nt_mir[imir] > max_tim) {
      printf("truncating nt_mir for mirror %d\n",imir);
      itav = (int)(t_cur_max / DT_ADC);
      itav = max(itav, it0mir + max_tim/2);
      itav = min(itav, it1mir - max_tim/2);
      it0mir = itav - max_tim/2 + 1;
      it1mir = itav + max_tim/2 - 1;
      pemir_.nt_mir[imir] = it1mir - it0mir + 1;
    }
    
//     pemir_.nt_mir[imir] = 
    printf("%d PE for mirror %d\n",len,pt->mir[mir]);
  }  
  
  return offset;
    
}





int simTLTriggerResponse(/*TLCalibration *cb, */PETimes *pt, 
                         fraw1_dst_common *fraw1,TLCalibration *tes) {
  double offset;
//   int i;
//     printf("%s::%s (%d): inside\n",__FILE__,__FUNCTION__,__LINE__);
//   printf("mirrs_.mirtub[6][248] = %d\n",mirrs_.mirtub[6][248]);
  offset = convertPETimesToTLPE(pt);
//   printf("%s::%s (%d): done with conversion\n",__FILE__,__FUNCTION__,__LINE__);
//   printf("pemir_.
//   abort();
  electronics_fadc_(tes);

  
  printf("value of mir_trig_.nmir_trig_2: %d\n",mir_trig_.nmir_trig_2);
//   printf("Trigger, by fiat, DOES NOT RESPOND\n");


  if (mir_trig_.nmir_trig_2) { // if the event did trigger
    fraw1->event_code = 1;
    fraw1->site = TALE_SITEID;
//     fraw1->part = 
    fraw1->jclkcnt += (int)offset;

    while (fraw1->jclkcnt >= 1e9 ) {
      fraw1->jclkcnt -= 1e9;
      fraw1->jsecond++;
    }
    
    while (fraw1->jsecond >= 86400) {
      fraw1->jsecond -= 86400;
      fraw1->julian++;
    }
    
    
//     for (i=0; i<fraw1->nmir; i++) {
// //       fraw1->second[i] = fraw1->jsecond;
// //       fraw1->clkcnt[i] = 
//     }
  }

//   abort();
  return mir_trig_.nmir_trig_2;
}

void initTLWaveformGenerator() {
  double tau_neg_16 = 100000.;
  double utau, vtau, ltau, ktau;
  double utau2, vtau2, duvtau, duvtau2, duvtau3;
//   double a1, a2, b1, b2;
  double un = 0., vn = 0., uacn = 0., vacn = 0.;

  double gt, gtu, gtu2, gtu3, gtl, gtk, gtv;
  double eu, ev, el1, el2;
  double utn, vtn;
  double dt, t;
  int idt, it;
  double xu, xv, xl1, xl2;
  double fu, fv;
  double tu, tv, tl1, tl2;
  double ut, vt;
//   float fu1, fv1, fnu1, fnv1, fv2u1, fv2v1;
//   float fu2, fv2, fu3, fu4;
  
  utau = TAU_FADC/*taufadc_.tau_fadc*/;
  vtau = TAU_16/*taufadc_.tau_16*/;
  ltau = TAU_NEG/*taufadc_.tau_neg*/;
  ktau = tau_neg_16;
  
  utau2 = utau * utau;
  vtau2 = vtau * vtau;
  duvtau = vtau - utau;
  duvtau2 = duvtau * duvtau;
  duvtau3 = duvtau * duvtau2;

  // Constants used for v16 term (4 poles, 2*utau + 2*vtau)
  
  a1 =      DT_ADC * 2 * utau2 * vtau  / duvtau3;
  a2 =      DT_ADC *     utau2         / duvtau2;
  b1 = -2 * DT_ADC *     utau  * vtau2 / duvtau3;
  b2 =      DT_ADC *             vtau2 / duvtau2;
  
  for (idt=0; idt<1000; idt++) {
    for (it=0; it<101; it++) {
      gt = DT_ADC * (double)it + 0.1* ((double)idt + 0.5);
      //f77: idt ranges from 1..1000, it 1..101, so "it-1" and "idt-0.5" appear
      gtu = gt / utau;
      gtu2 = 0.5 * gtu * gtu;
      gtu3 = gtu * gtu2 / 3.0;
      
      if (gtu < 25.) {
        eu = exp(-gtu)/utau;
      }
      else {
        eu = 0.;
      }
      
      gtl = gt / ltau;
      
      if (gtl < 25.) {
        el1 = exp(-gtl)/ltau;
        vtn = exp(-gtl);
      }
      else {
        vtn = 0.;
        el1 = 0.;
      }
      
      gtk = gt / ktau;
      if (gtk < 25.) {
        vtn -= exp(-gtk);
        el2 = exp(-gtk)/ktau;
      }
      else
        el2 = 0.;
      
      vtn /= (ktau - ltau);
      gtv = gt/vtau;
      
      if (gtv < 25.) {
        ev = exp(-gtv)/vtau;
      }
      else
        ev = 0;
      
      utn = DT_ADC * (2.0 - gtl) * el1;
      vtn += el1 + el2;
      ut = DT_ADC * gtu3 * eu - utn;
      vt = (a1 + a2 * gtu) * eu + (b1 + b2*gtv)*ev - vtn;
      if (ut > 0) {
        uacn += utn;
        un += ut;
      }
      if (vt > 0) {
        vn += vt;
        vacn += vtn;
      }
    }
  }
  unorm = un / 1000.;
  vnorm = vn / 1000.;      
  
  
  dt = 1.;
  for (it=0; it<101; it++) {
  // t is time from pe to measurement
    t = dt * (double)it;
    tv = t / vtau;
    tu = t / utau;
    tl1 = t / ltau;
    tl2 = t / ktau;
    fv = exp(-tv)/vtau;
    fu = exp(-tu)/utau;
    vt1[it] = fv;
    vt2[it] = fv*tv;
    uutn1[it] = exp(-tl1)/ltau;
    uutn2[it] = (2. - tl1)*exp(-tl1)/ltau;
    uvtn1[it] = exp(-tl2)/ktau;
    ut1[it] = fu;
    ut2[it] = fu*tu;
    ut3[it] = ut2[it]*tu/2.;
    ut4[it] = ut3[it]*tu/3.;
  }
  
  
  // Constants used to get voltages u1,u2,u3,u4,v1,v2
  // based on same voltages at previous clock cycle
  
  xu = DT_ADC/TAU_FADC/*taufadc_.tau_fadc*/;
  xv = DT_ADC/TAU_16/*taufadc_.tau_16*/;
  xl1 = DT_ADC/TAU_NEG/*taufadc_.tau_neg*/;
  xl2 = DT_ADC/tau_neg_16;
  fu1 = exp(-xu);
  fv1 = exp(-xv);
  fnu1 = exp(-xl1);
  fnv1 = exp(-xl2);
  fnu2 = exp(-xl1)*xl1;
  fv2u1 = tau_neg_16/(tau_neg_16 - TAU_NEG/*taufadc_.tau_neg*/);
  fv2v1 = TAU_NEG/*taufadc_.tau_neg*/ / (TAU_NEG/*taufadc_.tau_neg*/ - tau_neg_16);
  fu2 = fu1*xu;
  fv2 = fv1*xv;
  fu3 = fu2*xu/2.;
  fu4 = fu3*xu/3.;
  
  return;
}

void generateTLWaveform(int *ntpe, double tpe[100000], int incr_npe[100000], 
            double *tstart, int *nt, double v[max_tim], double v16[max_tim],
                        TLCalibration *tes) {
  
  // formerly new_pe.f in MCRU
  // note: tpe must be time ordered
  
  double u1[3200], u2[3200], u3[3200], u4[3200]; // sum1 voltages
  double v1[3200], v2[3200];                     // sum16 voltages
  double uun1[3200], uun2[3200];                 // Neg (AC) voltages
  double uvn1[3200], uvn2[3200];                 // Neg (AC) voltages
      

  // arguments for the random-number generator
//   integer4 lrvec=1;
//   real4 rvec[lrvec];

  int it, jt;
  
  double dt,t;
  double xmu_noise = tes->xmu_noise;
  int nt_more, ipe;
  int npe_n;
  int ipen,ix;
  double xt;
  double wpe;

//   if (*ntpe == 1433) printf("Got a tube with ntpe = 1433; is it 248?\n");
  
  nt_more = (int)(2.5 * xmu_noise * TAU_16/DT_ADC) + 1;

  v1[0] = xmu_noise/DT_ADC;
  v2[0] = xmu_noise/DT_ADC;
  u1[0] = xmu_noise/DT_ADC;
  u2[0] = xmu_noise/DT_ADC;
  u3[0] = xmu_noise/DT_ADC;
  u4[0] = xmu_noise/DT_ADC;
  uun1[0] = xmu_noise/DT_ADC;
  uvn1[0] = xmu_noise/DT_ADC;
  uun2[0] = xmu_noise/DT_ADC;
  uvn2[0] = xmu_noise/DT_ADC;

  ipe = 0;    // we're going to walk through all the 
              // photo-electron arrival times in this PMT
  if (*nt > max_tim) {
    printf(" hr_new_pe called for nt,ntpe = %d %d\n",*nt,*ntpe);
    *nt = max_tim;
  }

//   if (*ntpe == 1433) printf("nt %d nt_more %d\n",*nt,nt_more);

// The value of t is different in each iteration of the loop over "it" ...
// first, nt_more is subtracted from "it." The difference is the number of
// 100ns time samples away from tstart where the boundary "t" will be. Note that
// for early steps, t may be earlier than tstart. t is the time right now!



  for (it=1; it<(*nt + nt_more); it++) {
    jt = it - nt_more;
    t = *tstart + DT_ADC * (double)jt;
    // first update voltages:
    u1[it] = fu1*u1[it-1];
    u2[it] = fu1*u2[it-1] + fu2*u1[it-1];
    u3[it] = fu1*u3[it-1] + fu2*u2[it-1] + fu3*u1[it-1];
    u4[it] = fu1*u4[it-1] + fu2*u3[it-1] + fu3*u2[it-1] + fu4*u1[it-1];

    v1[it] = fv1*v1[it-1];
    v2[it] = fv1*v2[it-1] + fv2*v1[it-1];
    uun1[it] = fnu1*uun1[it-1];
    uun2[it] = fnu1*uun2[it-1] - fnu2*uun1[it-1];
    uvn1[it] = fnv1*uvn1[it-1];
    
   // now add in noise pe for this time bin
   
    if (xmu_noise > 0) {
//       rnpssn_(&glb_prm_.xmu_noise,&npe_n,&ierr);
      npe_n = prand(xmu_noise);
    }
    else
      npe_n = 0;
#ifdef NO_NOISE    
    npe_n = 0; // remove this line!
#endif    
    
    
    // correct for calibration correction:
//     npe_n = (int)(npe_n/calib_.calib_corr + 0.5);

    if (dt <= 0.)
      dt = 1.;

    for (ipen=0; ipen<npe_n; ipen++) {

//       ranlux_(&rvec,&lrvec);

//       xt = DT_ADC*rvec[0];
      xt = DT_ADC * RANDOM_NUMBER;
      ix = (int)(xt/dt + 0.5);

      u1[it] += ut1[ix];
      u2[it] += ut2[ix];
      u3[it] += ut3[ix];
      u4[it] += ut4[ix];
      v1[it] += vt1[ix];
      v2[it] += vt2[ix];
      uun1[it] += uutn1[ix];
      uun2[it] += uutn2[ix];
      uvn1[it] += uvtn1[ix];
    }

//      for debugging, check for NaN:
//     if (u1[it]*0 != 0)
//       abort();

// now do signal pe

    while (ipe < *ntpe) { // we haven't run out of PEs yet
      if (tpe[ipe] < t) { // does this PE land before "t"?
        wpe = (real4)(incr_npe[ipe]); // how many PEs together?
        ix = (integer4)( (t - tpe[ipe])/dt + 0.5); // how many slices earlier?
        if (ix > 100) {
          break;  // more than 100 slices... that one's too long ago. Exit loop?
        }
        if (ix < 0 || ix > 100) {
          ix = max(ix,0);
          ix = min(ix,100);;
        }
        u1[it] += wpe*ut1[ix];
        u2[it] += wpe*ut2[ix];
        u3[it] += wpe*ut3[ix];
        u4[it] += wpe*ut4[ix];
        v1[it] += wpe*vt1[ix];
        v2[it] += wpe*vt2[ix];
        uun1[it] += uutn1[ix];
        uun2[it] += uutn2[ix];
        uvn1[it] += uvtn1[ix];
        ipe++;    
        continue;
      }     
      break;
    }

    uvn2[it] = fv2u1*uun1[it] + fv2v1*uvn1[it];

    // now save answer
    
    if (jt >= 0 && jt < max_tim) { 
      v[jt] = DT_ADC * u4[it];
      v[jt] -= DT_ADC * uun2[it];
      v16[jt] = a1*u1[it] + b1*v1[it] + a2*u2[it] + b2*v2[it];
      v16[jt] -= DT_ADC*uvn2[it];
      v[jt] /= unorm;
      v16[jt] /= vnorm;
      // Add a factor to correct the calibration constant between
      // npe and FADC counts:
      v[jt] *= tes->calib_corr;
      v16[jt] *= tes->calib_corr;
    }
  }   // end loop over "it" bins
//   if (*ntpe == 1433) printf("it = %d\n",it);
  return;
}

void filterIsolatedHits(void) {
  // identify PMTs that pass the DSP scan but have no passing neighbors, and exclude them
  // from contributing to... something. TBD. This was previously filter_1_m.f


  // external structures: mir_hit_, filter_m_, glb_prm_



  int nhits;                  // number of hit channels for this mirror
  nhits = mir_hit_.nhit_mir;
  
  int i,j,k;                  // loop indices
  int i1, j1, k1;             // i,j,k + 1 for Fortran-based functions
  
  double tmin_av = 0.;         // for calculating the mean and RMS minimum interval between
  double tminsq_av = 0.;       // neighboring pairs of hit PMTs
  double ntmin = 0;     

  double tmin_sig;
  
  int itmin, itmax = 100;
  int itdiff, distance;
  int nnt;
  
  double cut_tim, cut_tim_2;
  double sig_cut_tim = 5.;
  double cut_tim_max = 30.;
  int mhit, khit;
  int itc, itc2;
  
  int nother;
  
  // First, give all PMTs a clean slate:
  for (i=0; i<nhits; i++)
    filter_m_.lcut_m[i] = 0;
  
   
  
  // Investigate each hit tube's relationship with its neighbors
  for (i=0; i<nhits; i++) {
//     mirror_xy_(&(mir_hit_.id_hit_mir[i]),&iy,&ix); // determine position on camera
    
    itmin = itmax; // least time difference observed between neighboring tubes
    
    for (j=0; j<nhits; j++) {
      if (filter_m_.lcut_m[j] > 0)        // has this tube already been cut?
        continue;
      
      if (i==j)                           // no autocomparison
        continue;
      
      itdiff = (int)(fabs(mir_hit_.t_hit_mir[i] - mir_hit_.t_hit_mir[j]) / DT_ADC);
      
      if (itdiff < itmin) { // is the time difference smaller than the previous minimum?
        i1 = i + 1;
        j1 = j + 1;
        distance = iang_m_(&i1, & j1);
        
        if (distance <= 1) // are the tubes neighbors?
          itmin = itdiff;
      }
    }
      
    if (itmin >= itmax) {   // no neighbors found for tube i
      filter_m_.lcut_m[i] = 11; // this tube could get another chance
      continue;                 // but will not contribute to calculation below
    }
      
    tmin_av += (double)itmin;
    tminsq_av += (double)itmin * (double)itmin;
    ntmin++;
  }
  
  if (ntmin < 2) {            // not enough surviving pairs
    for (i=0; i<nhits; i++) {
      if (filter_m_.lcut_m[i] == 0) // tube has not yet been cut
        filter_m_.lcut_m[i] = 12;   // these tubes don't get another chance
    }
    return;                   // nothing left to do with this mirror
  }


  nnt = ntmin;
  tmin_av /= (double)nnt;
  tminsq_av /= (double)nnt;
  
  tmin_sig = sqrt(tminsq_av - tmin_av*tmin_av);
  
  tmin_sig = max(tmin_sig, 2.);
  
  cut_tim = tmin_av + tmin_sig*sig_cut_tim; // set an interval cut based on distribution
                                            // of observed intervals: if a pair exceeds the
                                            // average by more than sig_cut_tim sigmas,
                                            // it will be discarded.
                                            
  cut_tim = min(cut_tim, cut_tim_max);      // let's not be too permissive, though.

  mhit = 0;
  
  for (i=0; i<nhits; i++) {   // looping over all hit tubes again
    if (filter_m_.lcut_m[i] > 0)     // tube has been cut already
      continue;
    
    khit = 0;
    for (j=0; j<nhits; j++) {
      if (filter_m_.lcut_m[i] > 0)   // tube's potential neighbor has been cut already
        continue;
      
      itdiff = (int)(fabs(mir_hit_.t_hit_mir[i] - mir_hit_.t_hit_mir[j]) / DT_ADC);
      if ((double)itdiff < cut_tim) { // interval is small enough
        i1 = i + 1;
        j1 = j + 1;
        distance = iang_m_(&i1, & j1);
        
        if (distance <= 2) // are the tubes somewhat neighborly?
          khit = 1;
      }
    }
    if (khit == 0)            // no general-vicinity tubes nearby in time
      filter_m_.lcut_m[i] = 13;
  }
  
  // Now, give good tubes' filtered neighbors the benefit of the doubt.
  
  cut_tim_2 = cut_tim + tmin_av; // relax the timing criterion
  
  for (i=0; i<nhits; i++) {
    if (filter_m_.lcut_m[i] != 11) // tube was filtered for other reason
      continue;

    // look for good neighbors. If they're like State Farm, they're there.
    itmin = itmax;
    for (j=0; j<nhits; j++) {
      if (filter_m_.lcut_m[j] > 0) // not a good neighbor. Fie!
        continue;
      
      if (i==j)
        continue;
      
      itdiff = (int)(fabs(mir_hit_.t_hit_mir[i] - mir_hit_.t_hit_mir[j]) / DT_ADC);
      if (itdiff < itmin) {
        i1 = i + 1;
        j1 = j + 1;
        distance = iang_m_(&i1, & j1);
        
        if (distance <= 2) // are the tubes somewhat neighborly?
          itmin = itdiff;
      }
    }
    
    if (itmin < (int)cut_tim_2)     // tube does have a good neighbor somewhere
      filter_m_.lcut_m[i] = 0;      // so it gets a free pass. Maybe.
      
  }
  
  // Pairs of tubes that are adjacent only to each other, like that couple at a party that
  // doesn't talk to anyone else, will be cut.
  
  itc = (int)cut_tim;
  itc2 = (int)cut_tim_2;
  
  for (i=0; i<nhits; i++) {
    if (filter_m_.lcut_m[i] > 0)
      continue;
    
    for (j=0; j<nhits; j++) {
      if (filter_m_.lcut_m[j] > 0 || i==j)
        continue;
        
      i1 = i + 1;
      j1 = j + 1;
      distance = iang_m_(&i1, & j1);   
      itdiff = (int)(fabs(mir_hit_.t_hit_mir[i] - mir_hit_.t_hit_mir[j]) / DT_ADC);
      
      if (distance <= 1 && itdiff < itc) {      // i & j are a pair

        nother = 0; // counter for other neighboring pairs
        for (k=0; k<nhits; k++) {
          if (filter_m_.lcut_m[k] > 0)
            continue;
          
          if (i==k || j==k)
            continue;
          
          i1 = i + 1;
          k1 = k + 1;
          distance = iang_m_(&i1, & k1);   
          itdiff = (int)(fabs(mir_hit_.t_hit_mir[i] - mir_hit_.t_hit_mir[k]) / DT_ADC); 
          
          if (distance <= 2 && itdiff < itc2)
            nother++;
          
          j1 = j + 1;
          k1 = k + 1;
          distance = iang_m_(&j1, & k1);   
          itdiff = (int)(fabs(mir_hit_.t_hit_mir[j] - mir_hit_.t_hit_mir[k]) / DT_ADC); 

          if (distance <= 2 && itdiff < itc2)
            nother++;
        }
        
        if (nother==0) {      // pair is isolated
          filter_m_.lcut_m[i] = 41;
          filter_m_.lcut_m[j] = 41;
        }
      }
    }
  }
  return;
}       

void sortHitsByTime(void) {
  // from time_order_m.f, this function does what its name would imply.
  // external structures: mir_hit_, filter_m_
  
  
  int inew[nhit_mir_max],iold[nhit_mir_max];
  int id_new[nhit_mir_max],it0_new[nhit_mir_max];
  int nt_new[nhit_mir_max],npe_new[nhit_mir_max];
  double t_new[nhit_mir_max],dt_new[nhit_mir_max];
  double sig_new[nhit_mir_max],sig_t_new[nhit_mir_max];
  int ich_new[nhit_mir_max],ich_t_new[nhit_mir_max];
  int it0_t_new[nhit_mir_max],nfadc_new[nhit_mir_max];
  int max_fadc_new[nhit_mir_max];
  int lcut_m_new[nhit_mir_max];

  int i,j,jmin = -1;
  double t,tmin;

  for (i=0; i<mir_hit_.nhit_mir; i++)
    inew[i] = 0;
  
  for (i=0; i<mir_hit_.nhit_mir; i++) {
    tmin = 999999.;
    for (j=0; j<mir_hit_.nhit_mir; j++) {
      if (inew[j] == 0) {
        t = mir_hit_.t_hit_mir[j];
        if (t < tmin) {
          tmin = t;
          jmin = j;
        }
      }
    }
    iold[i] = jmin;
    inew[jmin] = i;
    id_new[i] = mir_hit_.id_hit_mir[jmin];
    it0_new[i] = mir_hit_.it0_hit_mir[jmin];
    nt_new[i] = mir_hit_.it0_hit_mir[jmin];
    npe_new[i] = mir_hit_.npe_hit_mir[jmin];
    nfadc_new[i] = mir_hit_.nfadc_hit_mir[jmin];
    max_fadc_new[i] = mir_hit_.max_fadc_hit_mir[jmin];
    sig_new[i] = mir_hit_.sig_hit_mir[jmin];
    it0_t_new[i] = mir_hit_.it0_trig_mir[jmin];
    t_new[i] = mir_hit_.t_hit_mir[jmin];
    dt_new[i] = mir_hit_.dt_hit_mir[jmin];
    sig_t_new[i] = mir_hit_.sig_trig_mir[jmin];
    ich_new[i] = mir_hit_.ichan_hit_mir[jmin];
    ich_t_new[i] = mir_hit_.ichan_trig_mir[jmin];
    lcut_m_new[i] = filter_m_.lcut_m[jmin];
  }
  
  for (i=0; i<mir_hit_.nhit_mir; i++) {
    mir_hit_.id_hit_mir[i] = id_new[i];
    mir_hit_.nt_hit_mir[i] = nt_new[i];
    mir_hit_.npe_hit_mir[i] = npe_new[i];
    mir_hit_.sig_trig_mir[i] = sig_t_new[i];
    mir_hit_.it0_hit_mir[i] = it0_new[i];
    mir_hit_.nfadc_hit_mir[i] = nfadc_new[i];
    mir_hit_.max_fadc_hit_mir[i] = max_fadc_new[i];
    mir_hit_.it0_trig_mir[i] = it0_t_new[i];
    mir_hit_.t_hit_mir[i] = t_new[i];
    mir_hit_.dt_hit_mir[i] = dt_new[i];
    mir_hit_.sig_hit_mir[i] = sig_new[i];
    mir_hit_.ichan_trig_mir[i] = ich_t_new[i];
    mir_hit_.ichan_hit_mir[i] = ich_new[i];
    filter_m_.lcut_m[i] = lcut_m_new[i];
  }
  
  return;
}

int tentativeRayleighFilter(void) {
  // This code looks at hits to whether there is indeed a track. I think.
  
  int nhits;
  nhits = mir_hit_.nhit_mir;
  
  int i,ifi,j,il;
  int nhit_filt, nbad, n1;
  int ibad, jbad;
  double uhat[3],duhat[nhit_mir_max][3],dangle[nhit_mir_max];
  double vec[nhit_mir_max][3];
  int ihit_filt[nhit_mir_max];
  
  int nhit_unorm;
  int idpmt, ipmt;
  
  double dd;
  
  
  // Iterative process: run this loop once; repeat as necessary until nbad = 0
  
  do {

    nhit_filt = 0;
    
    // first, get the indices of surviving tubes
    for (i=0; i<nhits; i++) {
      if (filter_m_.lcut_m[i] == 0) {     // tube survived filters to this point
        ihit_filt[nhit_filt++] = i;
      }
    }
    
    if (nhit_filt < 4)        // not enough hits
      return 3;
    
    for (j=0; j<3; j++)
      uhat[j] = 0.;
    
    nhit_unorm = 0;

    // loop over time-sorted hits, get tube pointing vector and difference from previous hit
    for (ifi=0; ifi<nhit_filt; ifi++) {
      i = ihit_filt[ifi];
      idpmt = mir_hit_.id_hit_mir[i] - 1;
      for (j=0; j<3; j++) {
        vec[ifi][j] = 1.e-9 * (double)centre_.tbv[idpmt][j];
//         printf("ifi %d vec = %f %f %f\n",ifi,vec[ifi][0],vec[ifi][1],vec[ifi][2]);
        if (ifi > 0)
          duhat[ifi][j] = vec[ifi][j] - vec[ifi-1][j]; // difference vector
      }
      
      if (ifi > 0) {
        dd = dotprod(vec[ifi],vec[ifi-1]);
        if (fabs(dd) < .9999) {     // this condition will always be satisfied for angles
                                    // larger than 1 degree
          unitVector(duhat[ifi],duhat[ifi]);
          dangle[ifi] = R2D*acos(dd);     // what is the angle between these two hit tubes?
          nhit_unorm++;             // how many tubes are participating in the
                                    // cumulative difference vector?
          addvec(uhat,duhat[ifi],uhat);
        }
      }
    }
//     abort();
    if (nhit_unorm < 2)       // not enough contributions to the difference vector
      return 4;
    
    // now, uhat is the sum of many difference vectors (each normalized prior to addition)
    // normalize uhat:
    
    if (magvec(uhat) < 1.)    // no clear direction: not a track?
      return 5;
    
    unitVector(uhat,uhat);
    
    // now we loop over all hits again, seeing which ones are consistent with this track
    
    n1 = 0;
//     ssum = 0.;          // will be a measure of seconds per projected degree along track
    jbad = 1;
    nbad = 0;
    
    for (ifi=1; ifi<nhit_filt; ifi++) {
      i = ihit_filt[ifi];
      il = ihit_filt[ifi-1];
      ipmt = mir_hit_.id_hit_mir[i];
      if (dangle[ifi] < 1.85) {     // hit tube points within 1.85 degrees of previous hit
        dd = dotprod(uhat,duhat[ifi]);
        if (dd > 0.3) {             // difference vector is within 72.5 degrees of track
          n1++;
//           s = (mir_hit_.t_hit_mir[i] - mir_hit_.t_hit_mir[il])/(dd * dangle[ifi]);
//           ssum += s;
          ibad = 0;     // this is a good tube so far
        }
        else if (dd < -0.5) {       // difference vector is at least 120 degrees off track
//           s = 0.;
          ibad = 1;     // this is a bad tube
        }
        else {          // sort of off to the side, neither forward nor backward
          ibad = 0;
//           s = 0.;
        }
      }
      else {      // hit tube does not point within 1.85 degrees of previous hit
        dd = 0.;
//         s = 0.;
        ibad = 1;
      }
      
      // Consecutive bad tubes? Lose a hit!
      if (ibad == 1 && jbad == 1) {
        nbad++;
        filter_m_.lcut_m[il] = 61;
      }
      if (ibad == 1 && (ifi+1) == nhit_filt) { // last tube is bad
        nbad++;
        filter_m_.lcut_m[i] = 62;
      }
      
      jbad = ibad;
    }
    if (n1 < 3)        // too few would-be participants in the track
      return 6;
    
//     sav = ssum / (float)n1;   // seconds per degree
    
  } while (nbad > 0);
  
  return 0;
}