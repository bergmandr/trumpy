#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "tlelectronics.h"



typedef struct {
  integer4 nhit_trig;
  integer4 flag[256];
  integer4 slave[256];
  integer4 bdr_chn[256];
  integer4 npe_trig[256];
  integer4 tav_trig[256];
} to_mirtrig;

extern to_mirtrig to_mirtrig_;
void dumpFraw(fraw1_dst_common *fraw1, int line);

integer4 trigger_mir_(void);
;

int sec_trigger_(TLCalibration *tes, integer4 *mirid1, 
                 integer4 *itrig1/*, integer4 *iret*/) {
  int iret;
  integer4 trig_version;
  integer4 im,ih,j,k,jj,nj,it,imt;
  
  integer4 f77im,f77imt,cmirid;
  
  
  double t0,t1,ped,tt,tt0;
  integer4 nt,mped,ipmt,mhit,npe_min;
  integer4 madc[max_tim];
  double tau_2 = 1000.;
  integer4 scan_pulse,scan_itime;
  integer4 scan_it0, scan_nt;
  double scan_area,scan_time;
  

  integer4 idead;
  
  integer4 trig_info;
  
  // variables needed for trigger version 2
  integer4 cluster_found, hit_found;
  integer4 ichan, irow, icol, jrow, jcol;
  integer4 drow,dcol,itube;
  
  // variables needed for trigger version 3
  // the arguments passing into trigger_mir()
  // -- see the struct typedef prior to this function definition
  
  integer4 jt0[320], jadc[320][100], jt;
  integer4 it0_store,nt_store,madc_short[max_tim];
  integer4 it_min,it_max,first_slice;
  
  integer4 fiarg;
  dumpFraw(&fraw1_,__LINE__);
//   c********************************************************************
//   c@         Initialization
//   c********************************************************************

  /***/iret = -1;
  trig_info = 0;
  
  f77im = fraw1_.num_mir;
  
  im = f77im - 1;
  
  // first check if there is an entry in the ftrg1-block for this mirror
  
  f77imt = ftrg1_.num_mir;
  
  imt = f77imt - 1;
  
//   printf("ftrg1_.mir_num %d\n",ftrg1_.mir_num[imt]); // f77: (imt)
//   printf("fraw1_.mir_num %d\n",fraw1_.mir_num[im]); // f77: (im)
  
  if (f77imt != 0 && f77im != 0) {
    if (ftrg1_.mir_num[imt] == fraw1_.mir_num[im])
      trig_info = 1;
  }
  
  cmirid = *mirid1 - 1;
// c*************************************************************************     
// c@  rewriting the ftrg1- and fraw1- blocks to include only triggered hits:
// c*************************************************************************  

  // make a copy of the fraw1 variables for this mirror
  
  fraw1_copy_second = fraw1_.second[im];
  fraw1_copy_clkcnt = fraw1_.clkcnt[im];
  fraw1_copy_mir_num = fraw1_.mir_num[im];
  fraw1_copy_num_chan = fraw1_.num_chan[im];
  for (j=0; j<fraw1_.num_chan[im]; j++) {
    fraw1_copy_channel[j] = fraw1_.channel[im][j];
    fraw1_copy_it0_chan[j] = fraw1_.it0_chan[im][j];
    fraw1_copy_nt_chan[j] = fraw1_.nt_chan[im][j];
    if (fraw1_copy_nt_chan[j] > nt_chan_max) {
      printf("fraw1_copy_nt_chan[%d] > nt_chan_max: %d!\n",j,fraw1_copy_nt_chan[j]);
      abort();
    }
    for (k=0; k<fraw1_.nt_chan[im][j]; k++) {
      fraw1_copy_m_fadc[j][k] = fraw1_.m_fadc[im][j][k];
    }
  }
  
  // remove this mirror from the fraw1-block
  
  fraw1_.num_mir--;
  fraw1_.second[im] = 0;
  fraw1_.clkcnt[im] = 0;
  fraw1_.mir_num[im] = 0;
  fraw1_.num_chan[im] = 0;
  for (j=0; j<fraw1_copy_num_chan; j++) {
    fraw1_.channel[im][j] = 0;
    fraw1_.it0_chan[im][j] = 0;
    fraw1_.nt_chan[im][j] = 0;
    for (k=0; k<fraw1_copy_nt_chan[j]; k++)
      fraw1_.m_fadc[im][j][k] = 0;
  }
  // delete the fraw1-block and leave the program if there is no trigger
  // info for this irror
  
  if (trig_info == 0 && *itrig1 == 2)
    printf("trigger_az: ftrg1_.num_mir maxed out!!!!!!!!!!\n");
  
  if (trig_info == 0)
    return iret;
  
  // we have an ftrg1-entry: start fixing the ftrg1-block
  // first mke a copy of the ftrg1-variables for this mirror
  ftrg1_copy_trig_code = ftrg1_.trig_code[imt];
  ftrg1_copy_mir_num = ftrg1_.mir_num[imt];
  ftrg1_copy_nhit_dsp = ftrg1_.nhit_dsp[imt];
//   printf("nhit_dsp[%d] in sec_trigger: %d\n",imt,ftrg1_copy_nhit_dsp);
  for (ih=0; ih<ftrg1_.nhit_dsp[imt]; ih++) {
    ftrg1_copy_ichan_dsp[ih] = ftrg1_.ichan_dsp[imt][ih];
    ftrg1_copy_it0_dsp[ih] = ftrg1_.it0_dsp[imt][ih];
    ftrg1_copy_nt_dsp[ih] = ftrg1_.nt_dsp[imt][ih];
    ftrg1_copy_tav_dsp[ih] = ftrg1_.tav_dsp[imt][ih];
    ftrg1_copy_sigt_dsp[ih] = ftrg1_.sigt_dsp[imt][ih];
    ftrg1_copy_nadc_dsp[ih] = ftrg1_.nadc_dsp[imt][ih];
  }
  
  ftrg1_copy_t_pld_start = ftrg1_.t_pld_start[imt];
  ftrg1_copy_t_pld_end = ftrg1_.t_pld_end[imt];
  ftrg1_copy_row_pattern = ftrg1_.row_pattern[imt];
  ftrg1_copy_col_pattern = ftrg1_.col_pattern[imt];
  for (ih=0; ih<16; ih++) {
    ftrg1_copy_wp[ih] = ftrg1_.wp[imt][ih];
  }
  ftrg1_copy_delay = ftrg1_.delay[imt];
  
  // remove this mirror from the ftrg1-block
  
  ftrg1_.num_mir--;
  ftrg1_.trig_code[imt] = 0;
  ftrg1_.mir_num[imt] = 0;
  ftrg1_.nhit_dsp[imt] = 0;
  for (ih=0; ih<ftrg1_copy_nhit_dsp; ih++) {
    ftrg1_.ichan_dsp[imt][ih] = 0;
    ftrg1_.it0_dsp[imt][ih] = 0;
    ftrg1_.nt_dsp[imt][ih] = 0;
    ftrg1_.tav_dsp[imt][ih] = 0;
    ftrg1_.sigt_dsp[imt][ih] = 0;
    ftrg1_.nadc_dsp[imt][ih] = 0;
  }
  ftrg1_.t_pld_start[imt] = 0;
  ftrg1_.t_pld_end[imt] = 0;
  ftrg1_.row_pattern[imt] = 0;
  ftrg1_.col_pattern[imt] = 0;
  for (ih=0; ih<15; ih++) {
    ftrg1_.wp[imt][ih] = 0;
  }
  ftrg1_.delay[imt] = 0;
  
// c****************************************************************   
// c@   DATABASE OPTION: read dead mirrors from the database files
// c                    and skip them when writing the ftrg1 and 
// c                    fraw1 banks; also read the trigger version
// c                    from the database file
// c****************************************************************

//   if ( TRUE/*(param_src_.use_database > 0)*/ || FALSE/*(param_src_.use_data_set == 3)*/) {
//     
//     for (idead=0; idead<db_info_.db_deadmir; idead++) {
//       if (db_info_.db_dead[idead] == *mirid1) {
//         printf("dead mirror skipped: # %d\n",*mirid1);
//         return;
//       }
//     }
//     trig_version = DB_TRIG_VERSION/*db_info_.db_trig_version*/;
//   }
//   else
//     trig_version = trig_def_.def_sec_trig;

  for (idead=0; idead<tes->db_deadmir; idead++) {
    if (tes->db_dead[idead] == *mirid1) {
      printf("dead mirror skipped: # %d\n",*mirid1);
      return iret;
    }
  }
  trig_version = DB_TRIG_VERSION/*db_info_.db_trig_version*/;
  if (*itrig1 != 2) {  // let adjacent triggers pass
/*c*********************************************************************
c@         check PRIMARY TRIGGER 
c*********************************************************************

c   the Primary Trigger is applied in electronics_fadc.f; this subroutine
c    only contains a check of the trigger conditions...
c
c  trigger-conditions for each mirror (valid for all trigger versions):
c
c      at least two 3-fold coincidences of neighbor or next-to-neighbor
c      trigger channels
c      (satisfied by two 3-folds in one view (horizontal or vertical) or 
c       one 3-fold in each view)


c  determine the number of 3-fold coincidences for each view
c  from ftrg1_copy_row_pattern and ftrg1_copy_col_pattern

c         vcnt=0
c         hcnt=0
c         row_pattern=0
c         col_pattern=0
c         row_pattern=ftrg1_copy_row_pattern
c         col_pattern=ftrg1_copy_col_pattern
c         do icount=1,14
c            if(mod(col_pattern,2) .eq. 1) vcnt = vcnt + 1
c            if(mod(row_pattern,2) .eq. 1) hcnt = hcnt + 1  
c            row_pattern=row_pattern/2
c            col_pattern=col_pattern/2
c         enddo
c         cntsum = hcnt + vcnt
c         if (itrig1 .ne. 1 .and. cntsum .ge. 2) then
c            print*,'error in trigger_az (primary trigger) !'
c            print*,'itrig1,cntsum',itrig1,cntsum
c            stop
c         endif   */ 

    if (*itrig1 != 1) {
      /***/iret = -2;
      printf("iret %d\n",/***/iret);
      return iret;
    }
//     c  now state the trigger conditions for the different trigger versions
// 
// c*********************************************************************
// c@  ###################  SECONDARY TRIGGER  #########################
// c*********************************************************************
//     printf("trig_version %d\n",trig_version);
    if (trig_version == 1) {
// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
// @
//   TRIGGER VERSION 1  (valid until 01/05/2000) :
// 
// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//          
//    DSP scan (confirming scan) shows at least 6 hits above threshold
      
      if ( ftrg1_copy_nhit_dsp < 6 ) {
        /***/iret = -3;
        printf("iret %d\n",/***/iret);
        return iret;
      }
      /***/iret = 1;
      
    }
    else if (trig_version == 2) {
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c@
c  TRIGGER VERSION 2  (valid from 01/05/2000 until 06/01/2000) :
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c          New version of Host/Trigger code, looks for cluster of
c           at least 4 hits in a 6x6 sub-cluster -- Checks all 121
c           possible 6x6 subclusters. Scheme is as follows: Less than
c           4 hits, reject. Less than 16 hits, look for cluster.
c          [ Less than 100 hits, accept. More than 100 hits, use
c           prescale. Accept every 64th trigger regardless of
c           results of confirming scan. ] <- no prescale in the MC */       
      if (ftrg1_copy_nhit_dsp < 4) {
        /***/iret = -3;
        printf("iret %d\n",/***/iret);
        return iret;
      }
      if (ftrg1_copy_nhit_dsp < 16) {
        
        cluster_found = 0;
        jcol = 0;
        while (jcol < 11 && cluster_found == 0) {
          jcol++;
          jrow = 0;
          while (jrow < 11 && cluster_found == 0) {
            jrow++;
            
            hit_found = 0;
            for (itube=0; itube<ftrg1_copy_nhit_dsp; itube++) {
              ichan = ftrg1_copy_ichan_dsp[itube];
              printf("calling mirror_xy\n");
              mirror_xy_(&ichan,&irow,&icol);
              drow = irow - jrow;
              dcol = icol - jcol;
              if (drow >= 0 && drow < 6 && dcol >= 0 && dcol < 6) {
                hit_found++;
              }
            }
            
            if (hit_found >= 4)
              cluster_found = 1;
          }
        }
        
        if (cluster_found == 0) {
          /***/iret = -3;
          printf("iret %d (cluster_found = 0)\n",/***/iret);
          return iret;
        }
      }               // more than 16 hits -> accept
      
      /***/iret = 2;
      
    }
    
    else if (trig_version == 3) {
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c@
c  TRIGGER VERSION 3  (valid from 06/01/2000 until ?) :
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
c  the latest trigger version uses an algorithm as described in 
c  S. Westerhoff's routine 'trigger_mir.c'     */ 

      /***/iret = trigger_mir_();
      

      if (/***/iret > 0) {
        /***/iret = 3;
      }
      else {
        /***/iret = -3;
        printf("iret %d (trigger_mir_ <= 0))\n",/***/iret);
        return iret;
      }
    }               // END TRIGGER VERSIONS
    
  }             // end "if(itrig1.ne.2)" loop
  // This miror has been triggered! (or is adjacent trigger)

// c**********************************************************************
// c@         Set up for DSP scan (READOUT SCAN)
// c**********************************************************************

  t0 = pemir_.t0_mir[cmirid];
  t1 = pemir_.t1_mir[cmirid];
  
              
// c calculate store window for this mirror (store windows for adjacent 
// c triggers have already been calculated in electronics_fadc; here I'm 
// c only subtracting the 2*50 time slices to get the original store windows)

  if (*itrig1 == 2) {
    it0_store = (integer4)(t0/DT_ADC) + 50;
    nt_store = pemir_.nt_mir[cmirid] - 100;
  }
  else {
    it0_store = ftrg1_copy_t_pld_start;
    nt_store = ftrg1_copy_t_pld_end - it0_store + 1;
    if (/*db_info_.db_adj_opt*/DB_ADJ_OPT < 2) {
      it0_store -= (nt_store + 128);
      nt_store = 3*nt_store + 192;
    }
    else {
      it0_store += nt_store/2;
      nt_store = 4*nt_store + 200;
      it0_store -= nt_store/2;
    }
  }

  // Loop over all pmt in mirror:
//         printf("it0_store: %d, nt_store %d\n",it0_store,nt_store);
  for (ipmt=0; ipmt<fraw1_copy_num_chan; ipmt++) {
    nt = fraw1_copy_nt_chan[ipmt];
    first_slice = 0;
    for (it=0; it<nt; it++) {
      jt = fraw1_copy_it0_chan[ipmt] + (it);// - 1; // f77: it - 1
      
      // only allow data inside store window.

      if (jt >= it0_store && jt < (it0_store + nt_store)) {
        if (first_slice == 0)
          it_min = it;
        first_slice = 1;
        it_max = it;
      }

      madc[it] = (integer4)fraw1_copy_m_fadc[ipmt][it];
    }

    // Now scan for hits
    ped = /*pedfadc_.*/tes->ped_fadc[cmirid][ipmt];
    mped = (integer4)(16. * (ped - 0.5) + 0.5);
    
    scan_pulse = 0;
    
    
// c  Scan only from time slice 50 to time slice nt-50 (the first 50 and
// c   last 50 slices contain only electronic noise) and only if inside
// c   store window.
    it_min = max(it_min,50);        // f77: 51
    it_max = min(it_max,nt-50-1); // f77: nt - 50

    for (it=it_min; it<=it_max; it++) {
      madc_short[it-it_min] = madc[it];
    }

    fiarg = it_max - it_min + 1;
    TLdsp_scan_(&tau_2,madc_short,&fiarg,&mped,&scan_pulse,&scan_itime,
              &scan_it0,&scan_nt,&scan_area,&scan_time);
              
    scan_itime += it_min;// - 1;
    scan_it0 += it_min;// - 1;
      
    if (scan_pulse < THRESH_2) 
      continue;
    tt = t0 + scan_time;
    tt0 = DT_ADC * (real4)(scan_it0 - 1) + t0;
    
    if (ftrg1_.nhit_dsp[imt] < nhit_dsp_max) {
      ftrg1_.nhit_dsp[imt]++;
      mhit = ftrg1_.nhit_dsp[imt] - 1; // f77: - 0
    }
    else {
      npe_min = 1.e6;
      mhit = 0;
      for (j=0; j<ftrg1_.nhit_dsp[imt]; j++) {
        if (ftrg1_.nadc_dsp[imt][j] < npe_min) {
          npe_min = ftrg1_.nadc_dsp[imt][j];
          mhit = j;
        }
      }
    }
    
    jt0[mhit] = fraw1_copy_it0_chan[ipmt] + scan_itime - 50;
    
    for (jt=0; jt<100; jt++) {
      it = jt0[mhit] - fraw1_copy_it0_chan[ipmt] + jt;
      if (it >= 0 && it < nt) {
        jadc[mhit][jt] = madc[it];
      }
      else {
        printf("DSP scan error!\n");
        abort();
      }
    }
    
    // Load into trig_dsp_info
    
    ftrg1_.ichan_dsp[imt][mhit] = ipmt + 1;
    ftrg1_.it0_dsp[imt][mhit] = (integer2)(tt0/DT_ADC);
    ftrg1_.tav_dsp[imt][mhit] = (real4)scan_itime;
    ftrg1_.nt_dsp[imt][mhit] = scan_nt;
    ftrg1_.sigt_dsp[imt][mhit] = (real4)scan_pulse;
    ftrg1_.nadc_dsp[imt][mhit] = (integer2)scan_area;
  }
    
  ftrg1_.trig_code[imt] = ftrg1_copy_trig_code;
  ftrg1_.mir_num[imt] = ftrg1_copy_mir_num;
  ftrg1_.t_pld_start[imt] = ftrg1_copy_t_pld_start;
  ftrg1_.t_pld_end[imt] = ftrg1_copy_t_pld_end;
  ftrg1_.row_pattern[imt] = ftrg1_copy_row_pattern;
  ftrg1_.col_pattern[imt] = ftrg1_copy_col_pattern;
  for (ih=0; ih<16; ih++) {
    ftrg1_.wp[imt][ih] = ftrg1_copy_wp[ih];
  }
  ftrg1_.delay[imt] = ftrg1_copy_delay;
  
  ftrg1_.num_mir++;
  
//   printf("trigger: ftrg1_.num_mir 2 %d\n",ftrg1_.num_mir);
// c*************************************************************** 
// c@        now fix the fraw1-block
// c***************************************************************
// 
// c  Only record channel information if it is seen in ftrg1 too   
  fraw1_.num_mir++;
  fraw1_.second[im] = fraw1_copy_second;
  fraw1_.clkcnt[im] = fraw1_copy_clkcnt;
  fraw1_.mir_num[im] = fraw1_copy_mir_num;

  nj = -1;
//   printf("fraw1_copy_num_chan %d, ftrg1_.nhit_dsp[%d] = %d\n",fraw1_copy_num_chan,imt,ftrg1_.nhit_dsp[imt]);
  for (j=0; j<fraw1_copy_num_chan; j++) {
    jj = 0;
    while (jj < ftrg1_.nhit_dsp[imt]) {

      if (ftrg1_.ichan_dsp[imt][jj] == fraw1_copy_channel[j]) {

        nj++;
        break;
      }
      jj++;
    }
    if (jj < ftrg1_.nhit_dsp[imt]) {
      fraw1_.channel[im][nj] = fraw1_copy_channel[j];
      fraw1_.it0_chan[im][nj] = jt0[jj];
      fraw1_.nt_chan[im][nj] = 100;
      for (k=0; k<100; k++) {
        fraw1_.m_fadc[im][nj][k] = jadc[jj][k];
      }
    }
    else if (fraw1_copy_channel[j] > 256) {
/*c  trigger & LG channels: if there's no signal, write them out
c  anyway to preserve the structure of the fraw1-bank. Start after
c  the first 150 time slices (only el. + sky noise!)   */   
      nj++;
      fraw1_.channel[im][nj] = fraw1_copy_channel[j];
      fraw1_.it0_chan[im][nj] = fraw1_copy_it0_chan[j];
      fraw1_.nt_chan[im][nj] = 100;
      for (k=0; k<100; k++) {
        fraw1_.m_fadc[im][nj][k] = fraw1_copy_m_fadc[j][k+150];
      }
    }
  }
  fraw1_.num_chan[im] = nj + 1; // f77: = nj
  return iret;
}
    
    