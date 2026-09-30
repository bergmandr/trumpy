/*  Fills trigger patterns, start/stop times, and code for mirror ID = mirid
    Modify to remove references to fraw2 etc. 3/21/02 (jhb)
    Also remove fadc_trig_init. Only constants left are set in data statements. */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "tlelectronics.h"
// #include "maxmin.h"

// #ifndef max
// #define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
// #endif

// #ifndef min
// #define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
// #endif

// #include "dst_std_types.h"
// #include "fraw1_dst.h"

integer2 fadc_trig_pld_(integer4 *mirid, integer2 *v_pattern, integer2 *h_pattern,
                       integer2 *start, integer2 *stop, integer2 *code) {
  
  integer2 return_value;
  
  integer4 trigger_ver = 1;
  
  integer2 v_patterni,h_patterni,starti=-1,stopi/*,codei*/;
  
  integer4 mir,/*status,*/cd,vl,hl;
  
  static integer4 first = 0;
  
  // constants  
  
  integer2 trig_pld_thresh = 32;
  integer2 mask2 = 2;  
  integer4 fadc_trig_quiet = 32;
  
  integer2 h_ctr_1[16], h_ctr_2[16], v_ctr_1[16], v_ctr_2[16], vi,hi;
  integer2 i,j,ich,ch,ih,iv,nt,vj,hj,v_coin,h_coin,quiet_time=-1;
  integer4 hcnt,vcnt,madc,trigi,cnt,cnti,nm,nch,ntch;
  integer4 it0,it1,jt0,it,jt,vmax,hmax,vor,hor;
  integer2 vadc[16],hadc[16],vadcl[16],hadcl[16];
  integer2 vadcn[16],hadcn[16];
  integer4 ht0/*,k,l*/,nti;
  
/*  countdown controls logic for decrementing counter measuring confidence
    countdown = 0   is original logic: count down if coincidence ends
    countdown = 1   is alternative logic: count down only when both channels end. */
  integer4 countdown = 0;
  
//   real4 y;
//   printf("Done with declarations! un-initialized: it,it0,it1 %d %d %d\n",it,it0,it1);
//   fraw1_common_to_dumpf_(stdout,0);
//   printf("ecode %d site %d part %d nm %d\n",fraw1_.event_code,fraw1_.site,
//          fraw1_.part,fraw1_.num_mir);
//   printf("julian %d jsecond %d jclkcnt %d\n",fraw1_.julian,fraw1_.jsecond,fraw1_.jclkcnt);
//   for (i=0; i<fraw1_.num_mir; i++) {
//     printf("i %d sec %d clk %d mn %d nch %d\n",i,fraw1_.second[i],fraw1_.clkcnt[i],
//            fraw1_.mir_num[i],fraw1_.num_chan[i]);
//   }
  if (trigger_ver > 0)
    trig_pld_thresh = 32;
  else
    trig_pld_thresh = 28;
  
  first++;
  
  /*if (fraw1_.fraw1_julian >= 20030618) fadc_trig_quiet = 64;*/ 
  
  return_value = 0;           // default no trigger
  
  *v_pattern = 0;
  *h_pattern = 0;
  *code = -1;
//   printf("successfully wrote v_pattern, h_pattern, code: %d %d %d\n",*v_pattern,*h_pattern,*code);
  cnt = 0;
  vmax = 0;
  hmax = 0;
  
//   it=-1; // not in F77
  mir = -1; // starts at 0 in F77
  it1 = -1; // starts at 0 in F77
  it0 = 99999;

  nm = fraw1_.num_mir;
//   printf("nm: %d; mirid=%d,fmn[0]=%d\n",nm,*mirid,fraw1_.mir_num[0]);
  for (i=0; i<nm; i++) {
//     printf("i: %d, nm: %d, i<nm: %d\n",i,nm,(i<nm));
    if (fraw1_.mir_num[i] == *mirid) {
//       printf("found a match when i=%d (but mir=%d)\n",i,mir);
      mir = i; 
//          printf("found a match when i=%d (but mir=%d)\n",i,mir);    
    }
  }
//   printf("%s (%d): mir: %d\n",__FILE__,__LINE__,mir);    
  if (mir < 0)
    return return_value;
  
  
  nch = fraw1_.num_chan[mir];
//   printf("nch: %d\n",nch);
  for (ich=0; ich<nch; ich++) {
    ch = fraw1_.channel[mir][ich];
    iv = ch - 256; // compared with stored F77 index "ch", no decrement
    ih = ch - 288; // compared with stored F77 index, " " 
    if ((iv >= 1 && iv <= 16) || (ih >= 1 && ih <= 16)) {
      jt0 = fraw1_.it0_chan[mir][ich] - 1; // index of first time; no -1 in F77
      it0 = min(it0,jt0);
      nt = fraw1_.nt_chan[mir][ich]; // number of times, no shift
      it1 = max(it1,jt0 + nt); // F77: jt0 + nt - 1, but now jt0 can = 0
//       printf(" ch, it0, it1, nt, jt0: %d %d %d %d %d\n",ch,it0,it1,nt,jt0);
    }
  }
  
  
  
  

  // Find time limits for mirror
  
  
  
  
  
  nt = it1 - it0 + 1; // probably no + 1, but this quantity is not used again!
  ht0 = it0; // this quantity is also not used.
//   printf("ht0, it1: %d %d\n",ht0,it1);
  
  
//   abort();
 
  
  
  
/*   the next part was formerly structured like this:
<FORTRAN 77>
50    continue
      ...
      
      if (stopi .lt. it1 - fadc_trig_quiet) then
        it0 = stopi + fadc_trig_quiet
        go to 50
      end if
</FORTRAN 77>      
 I'll try replacing it with:
  
  while (stopi < it1 - fadc_trig_quiet) {
    ...
    
    it0 = stopi + fadc_trig_quiet;
  } 
  */




  
  stopi = 0;
//   printf("about to begin -while- loop\n");
  while (stopi < it1 - fadc_trig_quiet) {
//     printf("in the loop\n");
    trigi = 0;
    v_patterni = 0;
    h_patterni = 0;
    
    for (i=0; i<16; i++) {
      h_ctr_1[i] = 0;
      v_ctr_1[i] = 0;
      h_ctr_2[i] = 0;
      v_ctr_2[i] = 0;
      vadcl[i] = 0;
      hadcl[i] = 0;
      vadc[i] = 0;
      hadc[i] = 0;
      vadcn[i] = 0;
      hadcn[i] = 0;
    }
//     printf("about to loop over -it- until %d (it0 %d)\n",it1,it0);
//     printf("it %d, it0 %d, it1 %d\n",it,it0,it1);
//     abort();

    
    for (it=it0; it<it1; it++) {
//       printf("it first assigned to it0 (%d) currently = %d...\n",it0,it);
      for (ich=0; ich<nch; ich++) {
//         printf("inner loop iteration: ich=%d\n",ich);
        ch = fraw1_.channel[mir][ich]; // stored F77 index
        ntch = fraw1_.nt_chan[mir][ich]; 
//         printf("about to call jstar2 with arg based on mir %d ich %d \n",mir,ich);
//         printf("just to reiterate, mir has the value %d\n",mir);
//         printf("the arg is fraw1_.it0_chan[%d][%d] = %d\n",mir,ich,fraw1_.it0_chan[mir][ich]); 
//         printf("now, what was mir? Oh yeah, %d. Time to run jstar2?\n",mir);
//         jt0 = jstar2_(&(fraw1_.it0_chan[mir][ich])) - 1; // -0 in F77
        jt0 = (integer4)fraw1_.it0_chan[mir][ich] - 1; // -0 in F77
//         printf("complete; jt0 = %d\n",jt0);
        iv = ch - 256 - 1; // -0 in F77
        ih = ch - 288 - 1; // -0 in F77
        jt = it - jt0; // +1 in F77
//         printf("jt %d, ntch %d ...",jt,ntch);
        if (jt>=0 && jt<ntch)
//           madc = jstar1_(&(fraw1_.m_fadc[mir][ich][jt]));
          madc = (integer4)fraw1_.m_fadc[mir][ich][jt];
        else
          madc = 0;
//         printf("assigned madc = %d\n",madc);
        if (iv >= 0 && iv <16) {
//               printf("fraw1_.m_fadc[mir %d][ich %d][jt %d] (ch %d) = %d\n",
//                      mir,ich,jt,ch,madc);
//               abort();
          vmax = max(vmax,madc);
          if (madc >= trig_pld_thresh)
            vadcn[iv] = 1;
          else
            vadcn[iv] = 0;
          
          if (vadcl[iv] > 0) {
            if (v_ctr_1[iv] < 63)
              v_ctr_1[iv]++;
          }
          else {
            if (v_ctr_1[iv] > 0)
              v_ctr_1[iv]--;
          }
        }
        
        if (ih >= 0 && ih < 16) {
          hmax = max(hmax,madc);
          if (madc >= trig_pld_thresh)
            hadcn[ih] = 1;
          else
            hadcn[ih] = 0;
          
          if (hadcl[ih] > 0) {
            if (h_ctr_1[ih] < 63)
              h_ctr_1[ih]++;
          }
          else {
            if (h_ctr_1[ih] > 0)
              h_ctr_1[ih]--;
          }
        }
      }
//       printf("ctr values set\n");
      for (i=0; i<15; i++) {
        j = i + 1;
        vj = 0; // 0 in F77
        hj = 0; // 0 in F77
        if (j==15) {
          if (vadc[i] > 0 || vadc[j] > 0)
            vj = vadc[j] + v_ctr_1[j];
          if (hadc[i] > 0 || hadc[j] > 0)
            hj = hadc[j] + h_ctr_1[j];
        }
        else {
          if (vadc[i] > 0 || vadc[j] > 0 || vadc[j+1] > 0)
            vj = vadc[j] + v_ctr_1[j] + vadc[j+1] + v_ctr_1[j+1];
          if (hadc[i] > 0 || hadc[j] > 0 || hadc[j+1] > 0)
            hj = hadc[j] + h_ctr_1[j] + hadc[j+1] + h_ctr_1[j+1];
        }

        vi = vadc[i] + v_ctr_1[i];
        hi = hadc[i] + h_ctr_1[i];

        if (v_ctr_2[i] == 0) {
          if (vi > 0 && vj > 0)
            v_ctr_2[i]++;
        }
        else {
          if (v_ctr_2[i] < 63) {
            if (vi > 0 && vj > 0)
              v_ctr_2[i]++;
            else if ( countdown == 1 && (vi > 0 || vj > 0))
              v_ctr_2[i]++;
          }
          if (countdown == 1) {
            if (vi == 0 && vj == 0)
              v_ctr_2[i]--;
          }
          else {
            if (vi == 0 || vj == 0)
              v_ctr_2[i]--;
          }
        }
          
        if (h_ctr_2[i] == 0) {
          if (hi > 0 && hj > 0)
            h_ctr_2[i]++;
        }
        else {
          if (h_ctr_2[i] < 63) {
            if (hi > 0 && hj > 0)
              h_ctr_2[i]++;
            else if ( countdown == 1 && (hi > 0 || hj > 0))
              h_ctr_2[i]++;
          }
          if (countdown == 1) {
            if (hi == 0 && hj == 0)
              h_ctr_2[i]--;
          }
          else {
            if (hi == 0 || hj == 0)
              h_ctr_2[i]--;
          }
        }
      }
    
    
      v_coin = 0;
      h_coin = 0;
      for (i=0; i<14; i++) {
        vor = 0;
        hor = 0;
        for (j=0; j<5; j++) { // subtracting 1 from j here compared to F77
          if (j != 1 && i+j < 16) { // (j .ne. 2 .and. i+j-1 .le 16) in F77
            if (vadcn[i+j] > 0) // i+j-1 in F77
              vor = 1;
            if (hadcn[i+j] > 0) // i+j-1 in F77
              hor = 1;
          }
        }
        j = i + 1;
        vj = v_ctr_2[j];
        hj = h_ctr_2[j];
        if (j < 14 ) {
          vj += v_ctr_2[j+1];
          hj += h_ctr_2[j+1];
        }
        
        if (vor > 0 && v_ctr_2[i] > 0 && vj > 0) {
          v_coin += pow(2,i); // (i-1) in F77
//           printf("i %d, j %d: increasing v_coin to %d\n",i,j,v_coin);
//           abort();
        }
        if (hor > 0 && h_ctr_2[i] > 0 && hj > 0)
          h_coin += pow(2,i); // (i-1) in F77
      }
      
      vl = v_coin;
      hl = h_coin;
      
      if (v_coin > 0 || h_coin > 0) {
        if (trigi == 0) {
          starti = it; // index alert
          trigi = 1;
        }
//         printf("v_patterni %d, v_coin %d\n",v_patterni,v_coin);
        v_patterni |= v_coin;
        h_patterni |= h_coin;
        quiet_time = 0;
      }
      else {
        if (trigi == 1) {
          quiet_time++;
          if (quiet_time == fadc_trig_quiet) {
  
            stopi = it + 1 - quiet_time; // index alert
            // go to 100 -- how to implement? will use "break" to exit the for-loop,
            // arriving a point (just before 100)
            // where stopi will still be set the same way;
            // and comment out the assignment of stopi right here
//             printf("qt=ftq; stopi %d it %d\n",stopi,it);
            
     
            break;
          }
        }
      }
      for (j=0; j<16; j++) {
        vadcl[j] = vadc[j];
        hadcl[j] = hadc[j];
        vadc[j] = vadcn[j];
        hadc[j] = hadcn[j];
      }
    }
    if (it == it1)
      stopi = it1 - quiet_time; // +1 in F77
//     printf("past -for- loop; index it=%d, it1=%d\n",it,it1);
    
//     printf("stopi = %d\n",stopi);
//     abort();
    // 100 continue -- was here
//     printf("moving beyond 100, stopi=%d\n",stopi);
    
    if (trigi == 0)
      return return_value;
    
    return_value = 1;
    
    // Now count bits in v_pattern and h_pattern
    
    hcnt = 0;
    vcnt = 0;
    v_coin = v_patterni;
    h_coin = h_patterni;
    
    for (i=0; i<14; i++) {
      if (v_coin % mask2 == 1)
        vcnt++;
      if (h_coin % mask2 == 1)
        hcnt++;
      
      v_coin /= 2;
      h_coin /= 2;
    }
    cnti = vcnt + hcnt;
    if (cnti > cnt || (cnti == cnt && stopi - starti > *stop - *start)) {
      cnt = cnti;
      *v_pattern = v_patterni;
      *h_pattern = h_patterni;
      *start = starti + 1; // +0 in F77
      *stop = stopi + 1; // +0 in F77
//       cd = vcnt * (2**9) + hcnt * (2**5);
      cd = 512 * vcnt + 32 * hcnt; // used 512 and 32 instead of 2**9 and 2**5
      
      if (hcnt >= 2 || vcnt >= 2)
        cd++;
      if (hcnt >= 1 && vcnt >= 1)
        cd += 2;
      if (hcnt >= 3 || vcnt >= 3)
        cd += 4;
      if (hcnt >= 2 && vcnt >= 2)
        cd += 8;
      if (hcnt >= 3 && vcnt >= 3)
        cd += 16;
      
      if (fraw1_.event_code >= 2)
        cd -= 32768;
      
      *code = cd;
      nti = (stopi-starti)*8 + 200;
    }
    it0 = stopi + fadc_trig_quiet;
//     printf("loop end: stopi = %d\n",stopi);
  }
//   printf("past the -while- loop. stopi = %d\n",stopi);
  return return_value;
}
      
      
    
    
    
      