#ifndef _FDMC_ELECTRONICS_H_
#define _FDMC_ELECTRONICS_H_

#include "fdsite.h"

#define MIN_CL_FOR_HITPT 6.00
#define MIN_CL_FOR_WRITE 3.00
#define PE_MIN 100
#define VOFFSET 0.00

typedef struct {
  int pedestal[GEOFD_MAXMIR][GEOFD_MIRTUBE];
  int trigth[GEOFD_MAXMIR][GEOFD_MIRTUBE];
  int mean[GEOFD_MAXMIR][GEOFD_MIRTUBE];
  int disp[GEOFD_MAXMIR][GEOFD_MIRTUBE];
  double gain[GEOFD_MAXMIR][GEOFD_MIRTUBE];
} PMTState;

/* will be moving to 'raytrace' component soon */
typedef struct {
  double t1;  /* first pe arrival time (in ns) */
  double t2;  /* final pe arrival time (in ns) */
  double *n[GEOFD_MAXMIR][GEOFD_MIRTUBE]; /* weight of pe for each entry */
  double *t[GEOFD_MAXMIR][GEOFD_MIRTUBE]; /* array of pe arrival times (ns) */
  int len[GEOFD_MAXMIR][GEOFD_MIRTUBE];   /* number of entries in arrays */
} PETimes;

int initTriggerSimulator(void);
int simElectronics(PMTState *pmtstate, PETimes *pe);
unsigned int getHitBitPattern(int threshold, int *fadc);

#endif
