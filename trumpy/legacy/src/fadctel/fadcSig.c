#include <assert.h>

#include <stdio.h>
#include <math.h>

#include "constants.h"
#include "fdconstants.h"
#include "event.h"
#include "control.h"
#include "airshower.h"
#include "showerlib.h"
#include "track.h"
#include "fdsite.h"
#include "calibration.h"
#include "raytrace.h"
#include "acpttrack.h"
#include "taelectronics.h"

/*
 *  Finds the peak significance of a tube signal.
 *
 *  Author:  Sean R. Stratton
 *           Rutgers University, Dept. of Physics & Astronomy
 *
 *
 *  input:
 *    sdf   -  an FD_Waveform structure containing raw SDF data
 *
 *  returns:
 *    The largest significance value using 16 & 32-bin scan windows.
 *
 */
double fadcSig(int mean, int var, int nwf, int *wf) {
  int i, sum16, sum32, smax;

  assert( nwf >= 32 );

  /* initialize the 16-bin sum */
  sum16 = 0;
  for ( i=0; i<16; i++ )
    sum16 += wf[i];

  /* initialize the 32-bin sum, while scanning the first of the 16-bin sums */
  smax = sum16;
  sum32 = sum16;
  for ( i=16; i<32; i++ ) {
    sum32 += wf[i];

    sum16 += wf[i] - wf[i-16];
    if ( sum16 > smax )
      smax = sum16;
  }

  /* scan the rest of the fadc trace for the largest sum */
  for ( i=32; i<nwf; i++ ) {
    sum16 += wf[i] - wf[i-16];
    sum32 += wf[i] - wf[i-32];

    if ( sum16 > smax )
      smax = sum16;

    if ( sum32 > 2*smax )
      smax = sum32 / 2;
  }

  return (double)( smax - mean ) / sqrt( (double)var );
}
