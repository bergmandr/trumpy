#ifndef _FDSITE_H_
#define _FDSITE_H_

typedef geofd_dst_common FDSiteGeometry;

typedef struct {
  int nmir;
  int mir[GEOFD_MAXMIR];

  double dlseg;
  double qview;
  double tgen;
  double altmid;     /* altitude at middle of segment */
  double dist[GEOFD_MAXMIR];
  double dtheta[GEOFD_MAXMIR]; /* angle fraction (radians) */
  double geoFactor[GEOFD_MAXMIR];
  double fmie[NWAVELEN_BANDS];
  double fray[NWAVELEN_BANDS];
  double att[GEOFD_MAXMIR][NWAVELEN_BANDS];
  double costhe[GEOFD_MAXMIR];
  double nfl[GEOFD_MAXMIR][NWAVELEN_BANDS];
  double ncvdir[GEOFD_MAXMIR][NWAVELEN_BANDS];
  double ncvmie[GEOFD_MAXMIR][NWAVELEN_BANDS];
  double ncvray[GEOFD_MAXMIR][NWAVELEN_BANDS];
} ObservedTrackSegment;

typedef struct {
  int nseg;
  double rp;
  double psi;
  double uv[3];   // Shower direction unit vector in site coordinates
  double rimp[3]; // Shower impact position in site coordinates
  double rpuv[3]; // unit vector to point of closest approach in site coords
  double npln[3]; // shower-detector plane normal in site coordinates

  double local_vsite[3];  // CLF to Rp origin (CLF coordinates)
  double vsite[3];        // Location of Rp origin in ECC
  double local_vmir[GEOFD_MAXMIR][3];   // Location of mirrors 
                                        // wrt Rp origin (site coords)
  
  ObservedTrackSegment *segment;
  int siteid;
} ObservedTrack;

/* in "fdsite.c" */
#ifdef __cplusplus
extern "C"{
#endif
void getObservedTrack(RuntimeParameters *par, AirShower *as,
                      const Track *track, geofd_dst_common *site,
                      ObservedTrack *obtrack);
void getObservedTrackGeom(RuntimeParameters *par, AirShower *as,
                          const Track *track, geofd_dst_common *site, 
                          ObservedTrack *obtrack);
void fillObservedTrackLight(const RuntimeParameters *par, const Track *track, 
                            ObservedTrack *obtrack);
void getObservedTrack2(RuntimeParameters *par, AirShower *as,
                       const Track *track, geofd_dst_common *site,
                       ObservedTrack *obtrack);
void getObservedTrackGeom2(RuntimeParameters *par, AirShower *as,
                           const Track *track, geofd_dst_common *site, 
                           ObservedTrack *obtrack);
void fillObservedTrackLight2(const RuntimeParameters *par, const Track *track, 
                             ObservedTrack *obtrack);
void clearObservedTrack(ObservedTrack *obt);

void getFluxInStepsOfX(const Track *trk, const ObservedTrack *in,
                       int nseg, double *xtop, double *xbot, 
                       double areamr[GEOFD_MAXMIR], ObservedTrack *out);
double getRpVector(geofd_dst_common *fdsg, const AirShower *as, double *rpv);
double getPsiAngle(geofd_dst_common *fdsg, AirShower *as);
double getZenithAngle(geofd_dst_common *fdsg, AirShower *as);
void getSDPNVector(double rpuv[3], double uv[3], double n[3]);
void getViewingAngles(geofd_dst_common *fdsg, const AirShower *as,
                      const Track *track, ObservedTrack *obtrack,
                      double *viewang, double *chiang);
int getAngularSegments(const geofd_dst_common *fdsg, const AirShower *as, 
                       const Track *track, ObservedTrack *obtrack,
                       double *dlseg);
int applyEvsRpCut(double logE, double Rp);
void getObservedTrackOrigin(geofd_dst_common *fdsg, const AirShower *as,
                            ObservedTrack *obtrack);
void getRpOrigin(geofd_dst_common *fdsg, double *npln,
                 ObservedTrack *obtrack);

/* in "fdsite_up.c" */
void getLaserObservedTrack(double lambda, RuntimeParameters *par,
                          AirShower *as, const Track *track, 
                          geofd_dst_common *site,
                          ObservedTrack *obtrack);
void getLaserObservedTrackGeom(RuntimeParameters *par, AirShower *as,
                              const Track *track, geofd_dst_common *site, 
                              ObservedTrack *obtrack);
void fillLaserObservedTrackLight(double lambda, const RuntimeParameters *par, 
                                const Track *track, 
                                ObservedTrack *obtrack);
void getUpwardFluxInStepsOfX(const Track *trk, const ObservedTrack *in,
                              int nseg, double *xtop, double *xbot, 
                              double areamr[GEOFD_MAXMIR], ObservedTrack *out);

/* in "fdsitegeometry.c" */
int loadFDSiteGeometry(char *path, geofd_dst_common *fdsg);
void setFDSiteGeometry(const RuntimeParameters *par, geofd_dst_common *fdsg);
void getDefaultFDGeometry(int siteid);
#ifdef __cplusplus
}
#endif

#endif
