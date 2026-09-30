#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#include "tlelectronics.h"


#include "geotl_dst.h"
#include "adjinit.h"



void dumpFraw(fraw1_dst_common *fraw1, int line);


// caozh      the arguments passing into trigger_mir()
typedef struct {
  integer4 nhit_trig;
  integer4 flag[256];
  integer4 slave[256];
  integer4 bdr_chn[256];
  integer4 npe_trig[256];
  integer4 tav_trig[256];
} to_mirtrig;

extern to_mirtrig to_mirtrig_;



typedef struct {
  integer4 mpe;
} cztmp;

cztmp cztmp_;


// extern prescale prescale_;

typedef struct {
  integer4 jloop;
} to_trig_digit;

to_trig_digit to_trig_digit_;


// void sec_trigger_(TLElectronicsState *pesinteger4 *mirid1, integer4 *itrig1, integer4 *iret);
void electronics_fadc_(TLCalibration *tes) { 
  double ped,pemir_test[GEOTL_MAXMIR];
  integer4 madc[max_tim];
  integer4 mped,scan_pulse,scan_itime;
  integer4 scan_it0,scan_nt;
  double scan_area, scan_time;

  integer4 num_mir_max,nhit_filt_max,istatr;
  
  integer4 i,j,mr=-1,mirid1,jm,km,imir_max,idpmt,ipmt;
  integer4 irow,icol,istat_16,istat_old,ihit,nhit_filt;
  integer4 iret_td;
  integer4 itrig_16,ntw,jtba;
  
  integer4 ititle,mhit,npe_min;
  integer4 jmt=-1,col,row;
  integer4 it,miridj;
  
  double pemir_max,pemir_tot,time_mir,time_mir_max;
  double pemir_sc_max,t0,t1;
  double tt,tt0;

  double tau_1=1000.;


  double pemir_max_min = 100.;
  static integer4 ifirst = 0;
  
  integer2 v_pattern,h_pattern,start,stop,code,mask32 = 32;

  
  integer4 iret_trig,trig1;
  
  integer4 adj_opt,adjid,adj_view,nt_prim;
  double t0_prim;
  integer4 adj_mir[GEOTL_MAXMIR],adj_loop;
  integer4 nmir_loop,nmir_trig_a;
  integer4 idmir_loop[GEOTL_MAXMIR],idmir_trig_a[GEOTL_MAXMIR];
  integer4 ii,jj,kk,ipr,jpr;
  
  integer4 fiarg;
 
// **************************************************************************
// @                   Initialization
// **************************************************************************

  if (ifirst == 0) {
    num_mir_max = 0;
    ifirst = 1;
  }
  
  pemir_max = 0.;

  pemir_tot = 0.;
  rawcheck_.nmir_try = 0;

  mir_trig_.nmir_trig_1 = 0;
  mir_trig_.nmir_trig_2 = 0;
  nhit_filt_max = 0;

  ftrg1_.num_mir = 0;
  fraw1_.num_mir = 0;
  
  for (i=0; i<GEOTL_MAXMIR; i++) {
    adj_mir[i] = -1;
    idmir_trig_a[i] = 0;
    idmir_loop[i] = 0;
  }
  
//   if ( TRUE/*(param_src_.use_database > 0)*/ || FALSE/*(param_src_.use_data_set == 3)*/ )
//     adj_opt = DB_ADJ_OPT/*db_info_.db_adj_opt*/;
//   else
//     adj_opt = trig_def_.def_adj_opt;
//   
//   if (FALSE /*param_src_.use_database <= 0*/)
//     db_info_.db_prescale = 0;

  adj_opt = DB_ADJ_OPT;
  
  // Initialize trigger status
  istat_16 = 1;
  istat_old = 1;

// ***************************************************************************
// @  We check trigger for mirrors which view more than a minimum number of pe
//    Find mirid_max and all other mirrors above minimum
// ***************************************************************************

  for (jm=0; jm<pemir_.nmir_view; jm++) {
    mirid1 = pemir_.imir_view[jm];
    i = mirid1 - 1; // f77: -0

    pemir_tot += pemir_.pemir[i];
    time_mir = pemir_.t1_mir[i] - pemir_.t0_mir[i] - 6000.; // why?
    if (pemir_.pemir[i] > pemir_max) {
      pemir_max = pemir_.pemir[i];
      time_mir_max = time_mir;
      pemir_sc_max = pemir_.pemir_sc[i];
      imir_max = jm;
    }
   
    if (time_mir < 8000.)
      pemir_test[i] = pemir_max_min*sqrt(time_mir/8000.);
    else
      pemir_test[i] = pemir_max_min;
    
    pemir_test[i] = max(100.,pemir_test[i]);
   
    if (pemir_.pemir[i] >= pemir_test[i])
      rawcheck_.nmir_try++;   
  }
  printf("%s (%d): pemir_tot = %f\n",__FILE__,__LINE__,pemir_tot);

  if (rawcheck_.nmir_try == 0) {
    return;
  }

//         if adj_opt > 0, perform two loops:
//      first over all viewing mirrors to determine primary triggers,          
//      then over adjacent triggers to add them to the triggered
//      mirrors...

  adj_loop = 1;
  if (adj_opt > 0)
    adj_loop = 2;  // adjacent mirror trigger ON


// ********************************************************************
//      @   start loop over loop over mirrors
// ********************************************************************
  for (to_trig_digit_.jloop = 1; to_trig_digit_.jloop <= adj_loop; to_trig_digit_.jloop++) {
    if (to_trig_digit_.jloop == 2) {
//       printf("now looking for adjacent mirrors?\n");
      nmir_trig_a = 0;
      for (ii=0; ii<GEOTL_MAXMIR; ii++) {
//         printf("adj_mir[%d] = %d\n",ii,adj_mir[ii]);
        if (adj_mir[ii] >= 0) {
          ipr = adj_mir[ii];
          jpr = 0;
          for (jj=0; jj<ftrg1_.num_mir; jj++) {
            if (ftrg1_.mir_num[jj] == ipr) {
              jpr = jj;
            }
          }
          
//      create store window for this mirror; overwrite the original window, 
//      if this is a viewing mirror.

          t0_prim = ftrg1_.t_pld_start[jpr];
          nt_prim = ftrg1_.t_pld_end[jpr] - t0_prim + 1;
//           printf("jpr %d, adj_opt = %d, t0_prim %f, nt_prim %d tps %d\n",
//                  jpr,adj_opt,t0_prim,nt_prim,ftrg1_.t_pld_start[jpr]);
//           abort();
          if (adj_opt == 1) {
            pemir_.nt_mir[ii] = 3*nt_prim + 192;
            pemir_.t0_mir[ii] = 100*(t0_prim - nt_prim - 128);
            pemir_.t1_mir[ii] = pemir_.t0_mir[ii] + 100*(pemir_.nt_mir[ii] - 1);
          }
          else if (adj_opt == 2) {
            pemir_.nt_mir[ii] = 4*nt_prim + 200;
            pemir_.t0_mir[ii] = 100*(t0_prim + nt_prim/2 - pemir_.nt_mir[ii]/2);
            pemir_.t1_mir[ii] = pemir_.t0_mir[ii] + 100*(pemir_.nt_mir[ii]-1);
          }
//           printf("before adjustment: t0 %f, t1 %f\n",pemir_.t0_mir[ii],pemir_.t1_mir[ii]);
          // add 2 * 50 extra time slices for ambient noise
          
          pemir_.t0_mir[ii] -= 50*DT_ADC;
          pemir_.t1_mir[ii] += 50*DT_ADC;
          pemir_.nt_mir[ii] += 100;
          
          
          // prevent overflow
          pemir_.nt_mir[ii] = min(pemir_.nt_mir[ii],max_tim);
          
          
//           nmir_trig_a++;
          idmir_trig_a[nmir_trig_a++] = ii + 1; // f77: +0
          adj_view = 0;
          for (kk=0; kk<pemir_.nmir_view; kk++) {
            if (ii == pemir_.imir_view[kk]-1) { // f77: -0
              adj_view = 1;
            }
          }
          
          if (adj_view == 0) {
            pemir_.pemir[ii] = 0.;
            pemir_.pemir_sc[ii] = 0.;
          }
        }
      }
      nmir_loop = nmir_trig_a;
      for (i=0; i</*nmir_eye*/GEOTL_MAXMIR; i++)
        idmir_loop[i] = idmir_trig_a[i];
    }
    else {
      nmir_loop = pemir_.nmir_view;
      for (i=0; i<pemir_.nmir_view; i++) {
        idmir_loop[i] = pemir_.imir_view[i];
      }
    }
    
// ***************************************************************
// @      BEGIN LOOP OVER MIRRORS
// ***************************************************************

    for (jm=0; jm<nmir_loop; jm++) {
      mirid1 = idmir_loop[jm];
      
      
      
      if (ftrg1_.num_mir == ftrg1_nmir_max)
        return;
      
      // set up for digitization of rows and cols
      
      
      if (pemir_.pemir[mirid1-1] >= pemir_test[mirid1-1] || to_trig_digit_.jloop == 2) {
        
                    
// **************************************************************
// @   Generate raw data and digitize channels and sums.
//     Also, add ambient noise and extra hits.
// **************************************************************
    

        trig_digit_(tes,&mirid1,&iret_td);
        dumpFraw(&fraw1_,__LINE__);
        if (iret_td != 0) {
          printf(" IRET_TD = %d: SKIP THIS EVENT!\n",iret_td);
          return;
        }
        

// **************************************************************
// @  Check the PRIMARY TRIGGER for this mirror
//    ( row / column sums )
// **************************************************************
        itrig_16 = fadc_trig_pld_(&mirid1,&v_pattern,&h_pattern,&start,&stop,&code);
//         printf("Got past fadc_trig_pld; code=%d, start %d, stop %d\n",code,start,stop);

          if (ftrg1_.num_mir < ftrg1_nmir_max) {
            ftrg1_.num_mir++;
            mr = ftrg1_.num_mir - 1; // f77: -0
            ftrg1_.mir_num[mr] = mirid1;
            ftrg1_.t_pld_start[mr] = start;
            ftrg1_.t_pld_end[mr] = stop;
            ftrg1_.row_pattern[mr] = h_pattern;
            ftrg1_.col_pattern[mr] = v_pattern;
            if (to_trig_digit_.jloop == 2)
              ftrg1_.trig_code[mr] = 20;
            else
              ftrg1_.trig_code[mr] = code;
            
            ftrg1_.nhit_dsp[mr] = 0;
          }
        
  
// **********************************************************************
// @       # # # #          PRIMARY TRIGGER         # # # #
// **********************************************************************

//  PRESCALE 1: one out of 8 triggers with a single 3-fold coincidence
//  ----------   is allowed to pass the primary trigger.       

        if (code == 512 || code == 32) {
          prescale_.nsingle++;
          if ( (prescale_.nsingle % 8) == 0) {
            mir_trig_.idmir_trig_1[mir_trig_.nmir_trig_1++] = mirid1;
            printf("PRESCALE 1 !\n");
            ftrg1_.trig_code[mr] = 21;
          }
        }

//    At least 2 three-fold coincidences are required (code .gt. 0).
//     mark primary triggered mirrors in 'adj_mir' with '-1'    
                    
        code %= mask32;
        if (code > 0) {
          mir_trig_.idmir_trig_1[mir_trig_.nmir_trig_1++] = mirid1;
          
        }
        
        // begin LOOP over primary triggered and adjacent mirrors
        
        if ( ((mir_trig_.nmir_trig_1 != 0) && (mir_trig_.idmir_trig_1[mir_trig_.nmir_trig_1-1] == mirid1))
                        || to_trig_digit_.jloop == 2) {
          

          
// ***********************************************************************
// @    Set up for DSP scan  ( CONFIRMING SCAN )
// ***********************************************************************

          t0 = pemir_.t0_mir[mirid1-1];
          t1 = pemir_.t1_mir[mirid1-1];
          ntw = (integer4)( (t1 - t0)/DT_ADC) + 1;
          
          if (ntw > max_tim) 
            ntw = max_tim;
          
//           if (prtlev_.iprtlev >= 2)
//           if (mirid1 == 7)  printf(" t0,t1_window: %f %f; ntw %d \n",t0,t1,ntw);
          
          mir_hit_.nhit_mir = 0;
          
          // loop over all pmt in mirror:
          ititle = 0;
          for (ipmt = 0; ipmt<256; ipmt++) {
            idpmt = mirrs_.mirtub[mirid1-1][ipmt];
            fiarg = ipmt + 1;
            mirror_xy_(&fiarg,&irow,&icol);
            
            cztmp_.mpe = 0;
            for (jtba=0; jtba<idtemp_.ntba; jtba++) {
              if (idpmt == idtemp_.itba[jtba])
                cztmp_.mpe = petime_.npe_pmt[jtba];
            }
            
            // scan if inoise_mir > 0 .or. signal > 10 pe
            
            if (cztmp_.mpe < 10) {
              if (/*control_.inoise_mir == 0*/FALSE || tes->xmu_noise <= 0)
                continue;
            }
            
            // Scan only from time slice 50 (f77: 51) to time slice nt-51 (f77: 50)
            // (the first 50 and last 50 slices contain only ambient noise in the MC)
//             if (mirid1 == 7 && ipmt == 248) printf("waveform to be scanned:\n");
            for (it=50; it<ntw-50; it++) {
//               if (whicheye_.neye == 1)
//                 madc[it-50] = jstar1_(&(fraw1_.m_fadc[fraw1_.num_mir-1][ipmt][it]));
              madc[it-50] = (integer4)fraw1_.m_fadc[fraw1_.num_mir-1][ipmt][it];
//               if (mirid1 == 7 && ipmt == 248) printf("%3d ",madc[it-50]);
            }
            
//             if (mirid1 == 7 && ipmt == 248) printf("\n");
            
            // Now scan for hits
            
            ped = /*pedfadc_.*/tes->ped_fadc[mirid1-1][ipmt];
            mped = (integer4)(16. * (ped - 0.5) + 0.5);
//             if (prtlev_.iprtlev >= 3) {
//               if (whicheye_.neye == 1) {
//                 printf(" Scan imir, mirid, ipmt, nt %d %d %d %d %d\n",fraw1_.num_mir,
//                        mirid1,ipmt,cztmp_.mpe,ntw);
//               }
//             }
            
            fiarg = ntw-100;
            TLdsp_scan_(&tau_1,madc,&fiarg,&mped,&scan_pulse,
                      &scan_itime,&scan_it0,&scan_nt,&scan_area,&scan_time);
                      
            scan_itime +=50;
            scan_it0 +=50;
            
//             printf("scan_pulse %3d, THRESH_1 %f\n",scan_pulse,THRESH_1);
            if (scan_pulse < THRESH_1)
              continue;
//             printf("did not 'continue'; nhit_dsp[%d] right now is %d\n",jmt,ftrg1_.nhit_dsp[jmt]);
//             if (prtlev_.iprtlev >= 3)
//               printf(" found hit, npe = %f\n",scan_area);
            
            // put this info into /mir_hit/ common
            //nhit_mir++; // no increment yet; move to end
            tt = t0 + scan_time;
            mir_hit_.t_hit_mir[mir_hit_.nhit_mir] = tt;
            tt0 = DT_ADC*(real4)(scan_it0) + t0; // f77: scan_it0 - 1
            mir_hit_.it0_hit_mir[mir_hit_.nhit_mir] = (integer4)(tt0/DT_ADC) + 1; // f77: +0
            mir_hit_.id_hit_mir[mir_hit_.nhit_mir] = ipmt + 1; // f77: +0
            mir_hit_.sig_trig_mir[mir_hit_.nhit_mir] = (real4)scan_pulse;
            mir_hit_.npe_hit_mir[mir_hit_.nhit_mir] = (integer4)scan_area;
            mir_hit_.nt_hit_mir[mir_hit_.nhit_mir] = scan_nt;
            mir_hit_.nhit_mir++; 
        
            // save hits in ftrg1. First check if mirror defined yet
            

            jmt = -1; // f77: = 0
            for (km=0; km<ftrg1_.num_mir; km++) {
              miridj = ftrg1_.mir_num[km];
              if (mirid1 == miridj)
                jmt = km;
            }
            if (jmt == -1) { // f77: if jmt.eq.0
              printf("this mirror was not defined before scan!\n");
              if (ftrg1_.num_mir < ftrg1_nmir_max) {
                ftrg1_.num_mir++;
                jmt = ftrg1_.num_mir - 1; // f77: -0
                ftrg1_.mir_num[jmt] = mirid1;
                ftrg1_.nhit_dsp[jmt] = 0;
                if (to_trig_digit_.jloop == 2) {
                  ftrg1_.trig_code[jmt] = 20;
                  printf("adj. mirror saved in ftrg1\n");
                }
                else
                  ftrg1_.trig_code[jmt] = -1;
                ftrg1_.row_pattern[jmt] = 0;
                ftrg1_.col_pattern[jmt] = 0;
                ftrg1_.t_pld_start[jmt] = 0;
                ftrg1_.t_pld_end[jmt] = 0;
              }
            }
            
            if (ftrg1_.nhit_dsp[jmt] < nhit_dsp_max) {
              ftrg1_.nhit_dsp[jmt]++;
              mhit = ftrg1_.nhit_dsp[jmt] - 1; // f77: -0
            }
            else { // replace smallest hit
              npe_min = 1.e6;
              mhit = 0; // f77: =1
              for (j=0; j<ftrg1_.nhit_dsp[jmt]; j++) {
                if (ftrg1_.nadc_dsp[jmt][j] < npe_min) {
                  npe_min = ftrg1_.nadc_dsp[jmt][j];
                  mhit = j;
                }
              }
            }
    
            // load into /trig_dsp_info/
            
            ftrg1_.ichan_dsp[jmt][mhit] = ipmt + 1; // f77: +0
            ftrg1_.it0_dsp[jmt][mhit] = (integer4)(tt0/DT_ADC);
            // for confirming scan, this time is supposed to be the time
            // of peak of digital scan in clkcnts (100 nsec = 1 unit):
            
            ftrg1_.tav_dsp[jmt][mhit] = (real4)scan_itime;
            ftrg1_.nt_dsp[jmt][mhit] = scan_nt;
            ftrg1_.sigt_dsp[jmt][mhit] = (real4)scan_pulse;
            ftrg1_.nadc_dsp[jmt][mhit] = (integer4)scan_area;
          }
  
// *******************************************************************
// @     Filter hits for secondary trigger version 3
// *******************************************************************
          

          filterIsolatedHits();

          sortHitsByTime();

          
          if (tes->xmu_noise >= 0 && /*control_.inoise_mir > 0*/TRUE) {
            istatr = tentativeRayleighFilter();
          }
          
  
          // Count filtered hits:

          nhit_filt = 0;
//           printf("%s (%d): istatr = %d\n",__FILE__,__LINE__,istatr);
          if (istatr == 0) {
            printf("mir_hit_.nhit_mir = %d,ftrg1_.nhit_dsp[%d] = %d\n",
                   mir_hit_.nhit_mir,jmt,ftrg1_.nhit_dsp[jmt]);
            for (ihit=0; ihit</*mir_hit_.nhit_mir*/ftrg1_.nhit_dsp[jmt]; ihit++) {
              if (filter_m_.lcut_m[ihit] == 0) {
                
                
                ipmt = ftrg1_.ichan_dsp[jmt][ihit];
                if (ipmt < 1) printf("%s (%d): warning! bad ipmt (%d) at jmt %d ihit %d\n",
                  __FILE__,__LINE__,ipmt,jmt,ihit);
                mirror_xy_(&ipmt,&row,&col);
                to_mirtrig_.slave[nhit_filt] = col - 1;
                to_mirtrig_.bdr_chn[nhit_filt] = row - 1;
                to_mirtrig_.npe_trig[nhit_filt] = ftrg1_.nadc_dsp[jmt][ihit];
                to_mirtrig_.tav_trig[nhit_filt] = ftrg1_.tav_dsp[jmt][ihit];
                to_mirtrig_.flag[nhit_filt] = 0;
                nhit_filt++;  //increment at end               
              }
            }
          }
//           printf("%s (%d): nhit_filt = %d\n",__FILE__,__LINE__,nhit_filt);
          to_mirtrig_.nhit_trig = nhit_filt;
          for (ihit=nhit_filt; ihit<256; ihit++) {
            to_mirtrig_.slave[ihit] = 0;
            to_mirtrig_.bdr_chn[ihit] = 0;
            to_mirtrig_.npe_trig[ihit] = 0;
            to_mirtrig_.npe_trig[ihit] = 0;
            to_mirtrig_.tav_trig[ihit] = 0;
            to_mirtrig_.flag[ihit] = 0;
          }
          

          if (nhit_filt > nhit_filt_max)
            nhit_filt_max = nhit_filt;
          
        }         // End loop over prim. triggers and adjacent
   
// **********************************************************************
// @        # # # #        SECONDARY TRIGGER           # # # #
// **********************************************************************

//   If this mirror survived the primary trigger, apply the secondary
//   trigger. Determine adjacent triggers and let them pass through the
//   secondary trigger. Perform a readout scan on each mirror. 
//   Take the prescale constant from the database and apply PRESCALE 2
//   in sec_trigger.      
        
        
        iret_trig = -1;
        trig1 = 0;
                
        if (to_trig_digit_.jloop == 1) {
          if ( (mir_trig_.nmir_trig_1 != 0 && 
                  mir_trig_.idmir_trig_1[mir_trig_.nmir_trig_1-1] == mirid1) )
            trig1 = 1;
          
//           printf("%s (%d): trig1 = %d\n",__FILE__,__LINE__, trig1);
          iret_trig = sec_trigger_(tes, &mirid1,&trig1/*,&iret_trig*/);
          
          
          if (iret_trig > 0) {
            printf("MIRROR TRIGGERED!\n");
            
            sortHitsByTime();
            mir_trig_.idmir_trig_2[mir_trig_.nmir_trig_2] = mirid1;
            mir_trig_.nmir_trig_2++;
            
            adj_mir[mirid1 - 1] = -2;
            

/*c    check for adjacent triggers:   
c    * first adj. mir. trig. version: if there are at least 3
c    three-fold coincidences in the primary trigger (code .gt. 2), 
c    mark the adjacent triggers in 'adj_mir' with the ID of the
c    primarily triggered mirror.   */          
            
            if (adj_opt > 0 && iret_trig != 99 && code > 2) {
              for (ii=0; ii<5; ii++) {
                adjid = adjinit[mirid1-1][ii+1];
                if (adj_mir[adjid-1]==-1)
                  adj_mir[adjid-1] = mirid1;
              }
            }
/*c    * second adj. mir. trig. version: activates trigger if there
c      are 2 three-folds in either view with nt > 10.     */       
            
            if (adj_opt == 2 && iret_trig != 99 && code == 1) {
              if (ftrg1_.t_pld_end[ftrg1_.num_mir-1] - 
                  ftrg1_.t_pld_start[ftrg1_.num_mir-1] > 10) {
                for (ii=0; ii<5; ii++) {
                  adjid = adjinit[mirid1-1][ii+1];
                  if (adj_mir[adjid-1] == -1)
                    adj_mir[adjid-1] = mirid1;
                }
              }
            }
          }
        }
        else if (to_trig_digit_.jloop == 2) {
          // for adjacent triggers, perform only READOUT SCAN:
          
          trig1 = 2;
//           printf("calling sec_trigger at end\n");
          iret_trig = sec_trigger_(tes, &mirid1,&trig1/*,&iret_trig*/);
//           printf("finished sec_trigger at end\n");
          mir_trig_.idmir_trig_2[mir_trig_.nmir_trig_2] = mirid1;
          mir_trig_.nmir_trig_2++;
        }
      }
    }
    
  }
  return;
}
            
void dumpFraw(fraw1_dst_common *fraw1, int line) {        
//   int i,j,it;
//   printf("dump, called on line %d\n",line);
//   for (i=6; i<9; i++) {
//     printf("f.m[0][%d][20:29] ",i);
//     for(it=20;it<30;it++)printf("%3d ",fraw1->m_fadc[0][i][it]);
//     printf("\n");
//   }
  return;
}