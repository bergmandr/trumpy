#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "tlelectronics.h"
// #include "maxmin.h"

// #include "dst_std_types.h"

// #include "glb_prm.h"
// 

void TLdsp_scan_(double *tau, int madc[1024], int *nbin, int *mped,
               int *pulse, int *itime, int *i1, int *nt,
               double *area, double *time) {

  // for compatibility with parent F77 code
  // the itime argument needs to be incremented
  // likewise, when reporting back i1, needs to be incremented
  
  
  integer4 i,effped, fsum, fmax;
  double ped, tsum, adc;
  integer4 ifilter, i2, j;
  
  double filter;
  
  
  
  filter = exp(-DT_ADC/ *tau);
  ifilter = (integer4)(filter*32768 + 0.5);
  filter = (double)ifilter / 32768.;
  
  effped = (integer4)( (double)(*mped)/(1. - filter) + 0.5);
  
  fsum = effped;
  fmax = 0;
  *pulse = 0;
  *itime = 0;
//   printf("Detailed peak-finding begins here (%d entries):\n",*nbin);
  for (i=0; i<*nbin; i++) {
    fsum = (integer4)( (double)fsum * filter + 0.5) + 16* (madc[i]);
//     printf("(%d,%d,%d), ",i,fsum,madc[i]);
    if (fsum > fmax) {
//       printf("(%d,%d,%d)[<--],",i,fsum,madc[i]);
      fmax = fsum;
      *itime = i;
    }
  }
//   printf("\n");
//   printf("fmax %d, effped %d\n",fmax,effped);
  *pulse = (fmax - effped)/16;
  *pulse = min(*pulse,255);
  
  *i1 = *itime;
  i2 = *itime;
  
  ped = (double)(*mped)/16.;
  for (i=1; i<=*nbin; i++) { // preserving F77 loop since it is added or subtracted
    j = *itime - i;
    if (j >= 0 && j < *nbin) { // in F77: j .ge. 1 .and. j .le. nbin
      adc = madc[j] - ped;
      if (adc > 0.)
        *i1 = j;
      else
        break; // F77: go to 10
    }
  }
  // "10 continue" was here  
  
  for (i=0; i<*nbin; i++) {
    j = *itime + i;
    if (j >= 0 && j < *nbin) {
      adc = madc[j] - ped;
      if (adc > 0.)
        i2 = j;
      else
        break; // F77: go to 20
    }
  }
  // "20 continue" was here
//   printf("DSPSCAN: fm %d p %d,fp %d, i1 %d, i2 %d\n",fmax,*pulse,effped,*i1,i2);
  tsum = 0.;
  *area = 0.;
  
  for (i=*i1; i<=i2; i++) {
    adc = madc[i] - ped;
    *area += adc;
    tsum += adc*DT_ADC*(double)(i); // F77: (i-1)
  }
  
  *nt = i2 - *i1 + 1;
  
  if (*nt > 1)
    *time = tsum/ *area;
  else
    *time = 0.;
  *itime +=1;
  *i1 += 1;
  return;
}