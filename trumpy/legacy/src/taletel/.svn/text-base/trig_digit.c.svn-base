#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "control.h"
#include "event.h"
#include "fdconstants.h"
#include "tlcalibration.h"
#include "tlelectronics.h"

#define TRUE 1
#define FALSE 0

typedef struct {
  integer4 jloop;
} to_trig_digit;

extern to_trig_digit to_trig_digit_;



void trig_digit_(TLCalibration *tes, int *mir_id, int *istat) {
//   printf("trig_digit active\n");

  double v_sum16[64][max_tim];
  double tpe[100000], v1[max_tim], v16_new[max_tim];
  int incr_npe[100000];
  
  static int ifirst = 0; // static may be unnecessary; pending...
  
  int idpmt,ipmt,mirid,jtba;
  int ixcol,iyrow;
  int ntpe, /*mpe, */itpe;
  int it, nt;
  
  int j,jch,jmraw;
  int madci;
  
  double t0/*, t1*/;
  double hit_tubes;
  int debug_check;
  int inoise,icnt,iall,id_noise[256];

  double rn_prob,prob;
  int nntub,idntub,incr_noise;

  double tpe_noise;
  int coinc_tubes,iloop;

  int extranoise = TRUE; // f77: logical extranoise/.true./
                              // set to .false. to turn off high-amplitude sky noise

  int fiarg; // integer argument to be passed by reference to Fortran function
//   real4 ffarg1,ffarg2; // floating-point arguments to be passed by reference to F function


// ************************************************************
// @     Initialization
// ************************************************************

// This function operates on one mirror at a time. Steps taken during
// initialization include setting noise parameters: incr_noise, the contribution
// to the PMT voltage from noise PEs, and glb_prm_.xmu_noise, based on
// database ambient noise values,
//
// The mirror index jmraw refers to the position in the FRAW1 bank of this mirror.
// Currently the values stored in fraw1 are Fortran-indexed, hence the +1 when
// assigning fraw1_.mir_num[jmraw].
//
// Important: pemir_.t0_mir[mirid] is not the time of the first PE to arrive at the
// mirror. It is 15 microseconds earlier. Likewise, t1 is 15 microseconds later.
// Right now, fraw1_.it0_chan is the same for all mirrors, and fraw1_.nt_chan
// contains potentially >>100 samples, but later this will be reduced to the
// 100 most relevant samples for those channels that are kept, with it0_chan
// indicating where these begin.

  incr_noise = (int)(THRESH_2) + 1;
  incr_noise = (int)(incr_noise/tes->calib_corr + 1.);
#ifdef NO_NOISE  
  incr_noise = 0; // remove this line
#endif  
  
  mirid = *mir_id - 1;
  printf("*mir_id = %d, mirid = %d\n",*mir_id,mirid);
  if (FALSE /*param_src_.use_database == 0*/)
    extranoise = FALSE;
  
  if (TRUE /*param_src_.use_database > 0*/) {
    tes->xmu_noise = bg_noise_.amb_noise[mirid]/(tes->calib_corr * tes->calib_corr);
  }
  
  if (ifirst == 0) {
    if (FALSE/*param_src_.use_database <= 0*/)
      tes->xmu_noise /= (tes->calib_corr * tes->calib_corr);
  }

#ifdef NO_NOISE
  tes->xmu_noise = 0.; // remove this line
#endif  
  
  for (j=0; j<64; j++) {
    for (it=0; it<max_tim; it++) {
      v_sum16[j][it] = 0.;
    }
  }
  
  t0 = pemir_.t0_mir[mirid];
//   t1 = pemir_.t1_mir[mirid];
  nt = pemir_.nt_mir[mirid];
  // define this mirror
  
//   fraw1_.num_mir++;
//   jmraw = fraw1_.num_mir - 1; // f77: -0
  jmraw = fraw1_.num_mir++;
  fraw1_.mir_num[jmraw] = mirid + 1; // f77: mirid + 0
  fraw1_.num_chan[jmraw] = 320;
  
  for (jch=0; jch<320; jch++) {
    fraw1_.channel[jmraw][jch] = (integer2)jch + 1; // f77: + 0
    fraw1_.it0_chan[jmraw][jch] = (integer2)(t0 / DT_ADC);
    fraw1_.nt_chan[jmraw][jch] = nt;
  }

  
// ******************************************************************
// @       Defining the noise tubes
// ******************************************************************
//      
//      Additionally to the Poisson-distributed noise in 'new_pe', 
//      extra noise tubes are added to simulate high amplitude sky
//      noise.

  nntub = 0;
  if (extranoise && /*db_info_.db_trig_version*/DB_TRIG_VERSION < 3) {
    
    // noise distribution taken from DEC '99 - MAY '00 HiRes2 data
    
    coinc_tubes = 0;
    rn_prob = 1.0;
    prob = 0.0;
    
    while (rn_prob > prob) {
//       ranlux_(&rvec2,&lrvec2);
//       rn_nntub = rvec2[0];
//       rn_nntub = RANDOM_NUMBER;
      nntub = (integer4)(nnbins1*/*rn_nntub*/RANDOM_NUMBER);
      prob = skynoise_.skynoise1[nntub]/skynoise_.norm1; //f77: nntub+1
//       rn_prob = rvec2[1];
      rn_prob = RANDOM_NUMBER;
    }
  }
  else if (extranoise && /*db_info_.db_adj_opt*/DB_ADJ_OPT < 2) {
    // noise distribution taken from SEP '00 - Mar '01 HiRes2 data
    
    coinc_tubes = 0;
    rn_prob = 1.0;
    prob = 0.0;
    
    while (rn_prob > prob) {
//       ranlux_(&rvec2,&lrvec2);
//       rn_nntub = rvec2[0];
//       rn_nntub = RANDOM_NUMBER;
      nntub = (integer4)(nnbins2 * /*rn_nntub*/RANDOM_NUMBER);
      if (mirid <=3 || (mirid >=6 && mirid <= 9) ) // f77: 4, 7, 10 respectively
        prob = skynoise_.skynoise2d[nntub]/skynoise_.norm2d;  //f77: nntub+1
      else {
        if (to_trig_digit_.jloop == 1)
          prob = skynoise_.skynoise2t[nntub]/skynoise_.norm2t; //f77: nntub+1
        else
          prob = skynoise_.skynoise2a[nntub]/skynoise_.norm2a; //f77: nntub+1
      }
//       rn_prob = rvec2[1];
      rn_prob = RANDOM_NUMBER;
    }
  }
  else if (extranoise) {
    // noise distribution taken from March '01 - Sep '01 HiRes-2 data
    
    coinc_tubes = 0;
    rn_prob = 1.0;
    prob = 0.0;
    while (rn_prob > prob) {
//       ranlux_(&rvec2,&lrvec2);
//       rn_nntub = rvec2[0];
//       rn_nntub = RANDOM_NUMBER;
      nntub = (integer4)(nnbins3*/*rn_nntub*/RANDOM_NUMBER);
      if (mirid <= 3 || (mirid >= 6 && mirid <= 9) ) {
        prob = skynoise_.skynoise3d[nntub]/skynoise_.norm3d;
      }
      else {
        if (to_trig_digit_.jloop == 1)
          prob = skynoise_.skynoise3t[nntub]/skynoise_.norm3t;
        else
          prob = skynoise_.skynoise3a[nntub]/skynoise_.norm3a;
      }
//       rn_prob = rvec2[1];
      rn_prob = RANDOM_NUMBER;
    }
  }
  
// *****************************************************************
// @    select the noise tubes and mark them in array id_noise(256)
// *****************************************************************

  for (icnt=0; icnt<256; icnt++)
    id_noise[icnt] = 0;
  
  icnt = 1; //f77: = 1
  iloop = 0;

  while (icnt <= nntub && iloop < 10000) {
    iloop++;
//     ranlux_(&rvec,&lrvec);
//     rn_idntub = rvec[0];
//     rn_idntub = RANDOM_NUMBER;
    idntub = (integer4)(/*rn_idntub*/RANDOM_NUMBER*256 + 0.5)-0; // f77: -0

    if (id_noise[idntub] == 0) {
      icnt++;
      id_noise[idntub] = 1;
      idpmt = mirrs_.mirtub[mirid][idntub];

      // a noise-tube can coincide with a signal tube
      // in that case, do not count the tube as noise tube
      
      for (jtba=0; jtba<idtemp_.ntba; jtba++) {
        if (idpmt == idtemp_.itba[jtba]) {
          icnt--;
          /* uncomment next line if you do not want to allow noise in signal tube */
          // id_noise[idntub] = 0;
        }
      }
    }
  }
  
// c***********************************************************************    
// c@   Loop over tubes, add signal plus noise tubes to row/col sums
// c***********************************************************************
  hit_tubes = 0;
  
  for (ipmt=0; ipmt<256; ipmt++) {
    idpmt = mirrs_.mirtub[mirid][ipmt];
    ntpe = 0;
//     mpe = 0;

//      add noise hits  (new arrays tpe_noise and incr_noise)
//      time for noise hits is distributed randomly in a time
//      window starting 50 time slices after t0_mir and ending
//      50 time slices before t1_mir. The extra 2 * 50 time slices
//      will contain electronic noise only.
  
    inoise = 0;
    tpe_noise = 0.;
    if (id_noise[ipmt] == 1) {
//       ranlux_(&rvec,&lrvec);
//       rn_time = rvec[0];
//       rn_time = RANDOM_NUMBER;
      tpe_noise = t0 + 50*DT_ADC + /*rn_time*/RANDOM_NUMBER*DT_ADC * (nt - 100);
      inoise = 1;
    }

    // put all hits (signal + noise) into the right time order

    debug_check = 0;
    
    iall = -1;
    for (jtba=0; jtba<idtemp_.ntba; jtba++) {
      if (idpmt == idtemp_.itba[jtba]) {

        if (debug_check == 0) {
          hit_tubes++;
          debug_check = 1;
        }
        if (inoise == 1)
          coinc_tubes++;
//         mpe = petime_.npe_pmt[jtba];
        ntpe = petemp_.ntpe_temp[jtba];

  // c  signal and noise-pe's both present for this tube - get them in time order:
        for (itpe=0; itpe<ntpe; itpe++) {
          if (petemp_.tpe_temp[jtba][itpe] > tpe_noise && inoise == 1) {
            iall++;
            tpe[iall] = tpe_noise;            // first add in noise
            incr_npe[iall] = incr_noise;
            inoise = 0;
            iall++;                           // then add in signal
            tpe[iall] = petemp_.tpe_temp[jtba][itpe];
            incr_npe[iall] = petemp_.incr_npe_temp[jtba][itpe];
          }
          else if (petemp_.tpe_temp[jtba][itpe] == tpe_noise && inoise == 1) {
            iall++;
            tpe[iall] = petemp_.tpe_temp[jtba][itpe];// add noise+signal increments together
            incr_npe[iall] = petemp_.incr_npe_temp[jtba][itpe] + incr_noise;
            inoise = 0;
          }
          else {
            iall++;
            tpe[iall] = petemp_.tpe_temp[jtba][itpe];   // add in signal
            incr_npe[iall] = petemp_.incr_npe_temp[jtba][itpe];
          }
        }
        
    // add in remaining noise-PEs
    
        if (inoise == 1) {
          iall++;
          tpe[iall] = tpe_noise;
          incr_npe[iall] = incr_noise;
          inoise = 0;
        }
      }
    }
    
    // only noise-PEs for this tube:
    if (inoise == 1 && iall == -1) {
      iall++;
      tpe[iall] = tpe_noise;
      incr_npe[iall] = incr_noise;
      inoise = 0;
    }
    
    ntpe = iall + 1; // = iall in F77

    fiarg = ipmt + 1;
    mirror_xy_(&fiarg,&iyrow,&ixcol);
//     printf("%s (%d): calling waveform generator\n",__FILE__,__LINE__);
    generateTLWaveform(&ntpe, tpe, incr_npe, &t0, &nt, v1, v16_new, tes);

// ***********************************************************            
// @        Add signal to analog sums:
// ***********************************************************

    for (it=0; it<nt; it++) {
      // trigger sums:
      v_sum16[ixcol-1][it] += v16_new[it];
      v_sum16[iyrow+31][it] += v16_new[it];
      
      // low gain sums:
      
      v_sum16[ixcol+15][it] += v1[it];
      v_sum16[iyrow+47][it] += v1[it];
      
// ***************************************************************
// @    Save raw data for individual pmt
// ***************************************************************

      madci = (integer4)(v1[it]*tes->gain_fadc[mirid][ipmt] +
                          tes->ped_fadc[mirid][ipmt]); // f77: neye-0
      madci = min(madci,255);
      madci = max(madci,0);
//       printf("%3d ",madci);
      fraw1_.m_fadc[jmraw][ipmt][it] = (integer1)madci;
    }
  }
    
// ****************************************************************
// @    Now digitize sums
// ****************************************************************

  for (j=0; j<64; j++) {
    for (it=0; it<nt; it++) {
      madci = (integer4)(v_sum16[j][it]*tes->gain_fadc[mirid][256+j] +
                          tes->ped_fadc[mirid][256+j]); // f77: neye-0
      madci = min(madci,255);
      madci = max(madci,0);
      fraw1_.m_fadc[jmraw][256+j][it] = (integer1)madci;
    }
  }

  
  *istat = 0;
  
  return;
}