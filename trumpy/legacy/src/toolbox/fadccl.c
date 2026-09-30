#include <math.h>
#include "toolbox.h"

double fadccl(int mean, int vari, int nbin, short *wf) {
  int j, k, hsum, sum, sig, maxsig, psrs;
  double cl, maxcl;

  /* check the 8-bin set */
  hsum = 0;
  for ( j=0; j<4; j++ )
    hsum += wf[j];

  sum = hsum;
  maxsig = 0;
  for ( j=4; j<nbin; j+=4 ) {
    for ( k=0; k<4; k++ )
      sum += wf[j+k];

    psrs = ( sum << 1 ) - mean;
    sig = psrs * (psrs<<1);
    if ( sig > maxsig )
      maxsig = sig;

    sum -= hsum;
    hsum = sum;
  }
  maxcl = sqrt( (double)maxsig / (double)vari ) / 2.0;

  /* check the 16-bin set */
  hsum = 0;
  for ( j=0; j<8; j++ )
    hsum += wf[j];

  sum = hsum;
  maxsig = 0;
  for ( j=8; j<nbin; j+=8 ) {
    for ( k=0; k<8; k++ )
      sum += wf[j+k];

    psrs = sum - mean;
    sig = (psrs<<1) * (psrs<<1);
    if ( sig > maxsig )
      maxsig = sig;

    sum -= hsum;
    hsum = sum;
  }
  cl = sqrt( (double)maxsig / (double)vari ) / 2.0;
  if ( cl > maxcl )
    maxcl = cl;

  mean *= 8;

  /* check the 32-bin set */
  hsum = 0;
  for ( j=0; j<16; j++ )
    hsum += wf[j];

  sum = hsum;
  maxsig = 0;
  for ( j=16; j<nbin; j+=16 ) {
    for ( k=0; k<16; k++ )
      sum += wf[j+k];

    psrs = (sum<<2) - mean;
    sig = (psrs>>2) * (psrs>>1);
    if ( sig > maxsig )
      maxsig = sig;

    sum -= hsum;
    hsum = sum;
  }
  cl = sqrt( (double)maxsig / (double)vari ) / 2.0;
  if ( cl > maxcl )
    maxcl = cl;

  /* check the 64-bin set */
  hsum = 0;
  for ( j=0; j<32; j++ )
    hsum += wf[j];

  sum = hsum;
  maxsig = 0;
  for ( j=32; j<nbin; j+=32 ) {
    for ( k=0; k<32; k++ )
      sum += wf[j+k];

    psrs = ( sum<<1 ) - mean;
    sig = (psrs>>1) * (psrs>>1);
    if ( sig > maxsig )
      maxsig = sig;

    sum -= hsum;
    hsum = sum;
  }
  cl = sqrt( (double)maxsig / (double)vari ) / 2.0;
  if ( cl > maxcl )
    maxcl = cl;

  /* check the 128-bin set */
  hsum = 0;
  for ( j=0; j<64; j++ )
    hsum += wf[j];

  sum = hsum;
  maxsig = 0;
  for ( j=64; j<nbin; j+=64 ) {
    for ( k=0; k<64; k++ )
      sum += wf[j+k];

    psrs = sum - mean;
    sig = (psrs>>1) * psrs;
    if ( sig > maxsig )
      maxsig = sig;

    sum -= hsum;
    hsum = sum;
  }
  cl = sqrt( (double)maxsig / (double)vari ) / 2.0;
  if ( cl > maxcl )
    maxcl = cl;

  return maxcl;
}
