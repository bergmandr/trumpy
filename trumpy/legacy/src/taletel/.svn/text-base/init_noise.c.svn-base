#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "tlelectronics.h"
// #include "maxmin.h"
// #include "dst_std_types.h"

// #include "numcon.h"
// #include "glb_prm.h"
// #include "skynoise.h"
// #include "db_mcru.h"


// param_src param_src_ = {.use_database = 1};
// which_set which_set_;

skynoise skynoise_;
bg_noise bg_noise_;


/* not declared here in MCRU, but yes for TRUMP */




/* end of TRUMP-only declaration */


void init_noise_(void) {
  
  
  
  


//   real4 ds1_bgnoise[42] = {
//     5.3,5.8,4.7,4.8,5.1,4.5,4.6,4.8,5.3,
//     4.8,4.3,4.7,5.3,5.2,5.2,5.3,5.5,5.2,
//     4.6,4.5,4.5,4.4,5.5,4.5,5.5,5.3,4.3,
//     4.1,5.3,4.9,3.9,4.2,5.1,5.1,4.6,4.2,
//     5.3,5.4,4.4,4.2,5.8,5.2};
//                       
//                       
//   real4 ds2_bgnoise[42] = {
//     7.0,8.5,5.8,6.7,6.4,5.5,5.8,5.6,6.9,
//     7.1,5.9,6.7,7.3,6.9,7.8,6.9,8.0,6.2,
//     5.9,5.7,6.2,5.3,7.6,6.2,6.5,7.0,5.0,
//     5.4,6.1,5.8,4.5,4.6,6.2,5.9,4.5,4.6,
//     6.9,6.5,5.8,5.4,7.9,7.0};
    
  real4 ds3_bgnoise[42] = {
    7.2,7.6,6.2,6.8,6.8,5.9,6.6,6.8,6.7,
    6.2,5.7,6.0,6.3,6.2,6.4,6.0,6.4,6.5,
    5.6,5.4,5.3,5.6,6.1,5.4,6.4,6.6,5.8,
    5.8,7.1,7.2,4.5,5.6,6.7,5.6,6.0,5.7,
    5.6,6.9,6.0,5.5,6.4,6.1};
    
  int i/*,j*/;
  real4 temp13;
  real4 tempt, tempa, tempd;
  real4 temp3t, temp3a, temp3d;
  
  if (TRUE/*param_src_.use_database > 0*/) {
//     switch (which_set_.setnumber) {
//       case 1:
//         for (i=0; i<42; i++)
//           bg_noise_.amb_noise[i] = ds1_bgnoise[i];
//         break;
//       case 2:
//         for (i=0; i<42; i++)
//           bg_noise_.amb_noise[i] = ds2_bgnoise[i];
//         break;
//       default:
//         for (i=0; i<42; i++)
//           bg_noise_.amb_noise[i] = ds3_bgnoise[i];
//         break;
//     }
    for (i=0; i<42; i++)
      bg_noise_.amb_noise[i] = ds3_bgnoise[i];
    
  }
  

  
  
  skynoise_.norm1 = 0.;
  skynoise_.norm2t = 0.;
  skynoise_.norm2a = 0.;
  skynoise_.norm2d = 0.;
  skynoise_.norm3a = 0.;
  skynoise_.norm3t = 0.;
  skynoise_.norm3d = 0.;
  
  for (i=0; i<nnbins1; i++) {
    skynoise_.skynoise1[i] = exp(-0.5*pow( ( (double)i + 50.)/20.,2) ) / 
      (numcon_.root_2pi * 20.);
    
    temp13 = skynoise_.skynoise1[i];
    skynoise_.norm1 = max(skynoise_.norm1,temp13);
    
  }
  
  for (i=0; i<nnbins2; i++) {
    skynoise_.skynoise2t[i] = exp(-0.5*pow( ((double)i + 4.44)/(20. + 10. + 7.),2))/
      (numcon_.root_2pi * (20. + 10. + 7.) );
      
    if (i==0)
      skynoise_.skynoise2t[i] += 0.04;
      
    skynoise_.skynoise2a[i] = exp(-0.5*pow( ( (double)i + 21.184+5.)/(21.775+15.+7.),2)) /
      (numcon_.root_2pi * (21.775+15.+7.));
    
    if (i==0)
      skynoise_.skynoise2a[i] += 0.1;
    
    if (i==1) {
      skynoise_.skynoise2a[i] = skynoise_.skynoise2a[i] - 0.12 * skynoise_.skynoise2a[i];
      skynoise_.skynoise2a[i-1] = skynoise_.skynoise2a[i-1] + 0.12*skynoise_.skynoise2a[i];
    }
    
    skynoise_.skynoise2d[i] = exp(-0.5*pow( ( (double)i + 15.)/(30.+7.),2)) / 
      (numcon_.root_2pi * (30.+7.));
      
    if (i==0)
      skynoise_.skynoise2d[i] += 0.04;
  }
  
  for (i=0; i<nnbins2; i++) {
    tempt = skynoise_.skynoise2t[i];
    skynoise_.norm2t = max(skynoise_.norm2t, tempt);
    tempa = skynoise_.skynoise2a[i];
    skynoise_.norm2a = max(skynoise_.norm2a, tempa);
    tempd = skynoise_.skynoise2d[i];
    skynoise_.norm2d = max(skynoise_.norm2d, tempd);
  }
  
  for (i=0; i<nnbins3; i++) {
    skynoise_.skynoise3t[i] = exp(-0.5*pow( ((double)i + 4.44) / (37.+37.),2)) /
      (numcon_.root_2pi * (37.+37.));
    
    if (i==0)
      skynoise_.skynoise3t[i] += 0.04;
    
    skynoise_.skynoise3a[i] = exp(-0.5*pow( ((double)i + 26.) / (44.+44.),2)) /
      (numcon_.root_2pi * (44.+44.));
    if (i==0)
      skynoise_.skynoise3a[i] += 0.1;
    if (i==1) {
      skynoise_.skynoise3a[i] = skynoise_.skynoise3a[i] - 0.12*skynoise_.skynoise3a[i];
      skynoise_.skynoise3a[i-1] += 0.12*skynoise_.skynoise3a[i];
    }
    
    skynoise_.skynoise3d[i] = exp(-0.5*pow( ((double)i + 15.) / (37.+37.),2)) /
      (numcon_.root_2pi * (37.+37.));
    
    if (i==0)
      skynoise_.skynoise3d[i] += 0.04;
    
    temp3a = skynoise_.skynoise3a[i];
    skynoise_.norm3a = max(skynoise_.norm3a,temp3a);
    temp3t = skynoise_.skynoise3t[i];
    skynoise_.norm3t = max(skynoise_.norm3t,temp3t);
    temp3d = skynoise_.skynoise3d[i];
    skynoise_.norm3d = max(skynoise_.norm3d,temp3d);
  }
  return;
}
      
