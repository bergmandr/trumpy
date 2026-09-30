#ifndef _ACPTTRACK_H_
#define _ACPTTRACK_H_

// 0.00017453 = 0.1 * M_PI/180
// #define ACPTTRACK_DTH 1.7453e-3
#define ACPTTRACK_DTH 1.7453e-4
// #define ACPTTRACK_DTH 1.7453e-5
#define ACPTTRACK_MAX_NRTP 50
#define ACPTTRACK_DT 100  // Number of nanoseconds per bin

// #include "track.h"

typedef struct {
  double t1;    // Time of first pe
  double t2;    // Time of last pe
  int nmir;                      
  int mir[GEOFD_MAXMIR];        
  double *n[GEOFD_MAXMIR][GEOFD_MIRTUBE];
  double *t[GEOFD_MAXMIR][GEOFD_MIRTUBE];
  int len[GEOFD_MAXMIR][GEOFD_MIRTUBE];
} PETimes;

/* in "acpttrack.c" */
#ifdef __cplusplus
extern "C"{
#endif
int getAcceptedTrack(const AirShower *as, 
                 const Track *tk,
                 geofd_dst_common *site,
                 ObservedTrack *obtk,
                 TACalibration *calib,
                 TGeom *tgeom, 
                 PETimes *pet);
int getUpwardAcceptedTrack(const AirShower *as,
                     const Track *tk,
                     geofd_dst_common *site,
                     ObservedTrack *obtk,
                     TACalibration *calib,
                     TGeom *tgeom, 
                     PETimes *pet);
//void sortAcceptedTrack (AcceptedTrack *at);
//void clearAcceptedTrack(AcceptedTrack *at);

void clearPETimes(PETimes* pe);
#ifdef __cplusplus
}
#endif
#endif
