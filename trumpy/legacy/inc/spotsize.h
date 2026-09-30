#ifndef _SPOTSIZE_H_
#define _SPOTSIZE_H_

// These are the values for H. Tokuno's "spot size parameter" to direct
// the 2-D Gaussian smearing of the reflection normal vector on the mirror.

// These are in units of degrees. Translation to TRUMP values requires multiplication
// by conversion factor: 0.01 m per 0.195 degrees (see initRayTrace.c).

static double ssp[2][12] = {
  {0.187,  // br 00
   0.112,  // br 01 
   0.189,  // br 02
   0.098,  // br 03
   0.077,  // br 04
   0.139,  // br 05
   0.189,  // br 06
   0.087,  // br 07
   0.196,  // br 08
   0.080,  // br 09
   0.215,  // br 10
   0.073}, // br 11
  {0.175,  // lr 00
   0.055,  // lr 01
   0.195,  // lr 02
   0.072,  // lr 03
   0.183,  // lr 04
   0.077,  // lr 05
   0.187,  // lr 06
   0.038,  // lr 07
   0.211,  // lr 08
   0.069,  // lr 09
   0.213,  // lr 10
   0.077}  // lr 11
};

#endif
