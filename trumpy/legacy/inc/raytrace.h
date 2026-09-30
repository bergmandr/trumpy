#ifndef _RAYTRACE_H_
#define _RAYTRACE_H_

// don't use this one anymore!
// #define VERBOSE                   0

// use this one instead! (see trace.c for effects)
// #define RAYTRACE_VERBOSE


#define NUMTAPOLES              56           // Added 2 poles - SBT 20150430
#define NUMMDPOLES              11
#define NUMTLPOLES              NUMMDPOLES

#define MAXNPE                4096

#define USE_WIGGLE             1
#define USE_HIT_CAMERABOX      1
#define USE_HIT_CAMERAPOLES      1
#define USE_HIT_FLASHERBOX      1
#define USE_MIR_ABSORBED       0
#define USE_HIT_CRACK          1
#define USE_MISSED_CAMERA      1
#define USE_CAMCOVER_ABSORBED  0
#define USE_CAMCOVER_SHIFT     1
#define USE_FILTER_ABSORBED      0
#define USE_FILTER_SHIFT       1
#define USE_MISSED_TUBE        1
#define USE_PMT_UNIFORMITY       1
#define USE_PMT_QE             0

#define NUMERRS                10
#define HIT_CAMERABOX         -10
#define HIT_CAMERAPOLES         -9
#define HIT_FLASHERBOX          -8
#define MIR_ABSORBED          -7
#define HIT_CRACK       -6
#define MISSED_CAMERA         -5
#define CAMCOVER_ABSORBED     -4
#define FILTER_ABSORBED         -3
#define MISSED_TUBE           -2
#define NO_PE_CONVERSION      -1

typedef struct {
  double spot[GEOFD_MAXMIR][GEOFD_SEGMENT];                          // spot size on camera (meters)
  double mir_reflect;                   // mirror reflectivity
  double camcover_trans;                // camera cover transmission
  double filter_trans;                  // BG3 filter transmission
  double n_air, n_glass, n_filter;      // indices of refraction
  double pmtQE;                         // average PMT quantum efficiency
  double thick_glass;                   // thickness of cam box cover (meters)
  double thick_incam;                   // distance between cam box cover 
                                        //   glass and PMT filter (meters)
  double thick_filter;                  // distance between front of filter 
                                        //   and PMT face (meters)

  double sep[GEOFD_MAXMIR];             // distance between mirror and camera 
                                        //   box cover (meters)
  double h[GEOFD_MAXMIR];               // height of target plane above 
                                        //   center of mirror (meters)
  double maxseg[GEOFD_MAXMIR];          // maximum dist. allowed for a ray to
                                        //   land from center of mir seg (m)
  double cosang[GEOFD_MAXMIR];          // maximum angle that a ray u vector
                                        //   is allowed to have with respect
                                        //   to geofd_.vseg
  double meantime[GEOFD_MAXMIR];        // mean travel time from mirror
                                        //   to camera (accurate to within a
                                        //   nanosecond)


  double hexpt[6][3];                   // unit vectors from each point on a 
                                        //   hexagon to the hexagon's 
                              //   center (with hexagon in xy plane)
  double seg_the[GEOFD_MAXMIR][GEOFD_SEGMENT];        // angles of centers of mirror segments
                                        //   (spherical polar theta, phi)
  double seg_phi[GEOFD_MAXMIR][GEOFD_SEGMENT];
  double hexmirpt[GEOFD_MAXMIR][GEOFD_SEGMENT][6][3]; // unit vectors from mir. seg. point 
                                        //   to mir. seg. center

  double pmt_map[PMT_NUMDIV][PMT_NUMDIV]; // PMT uniformity map
  double max_sensitivity;             // value of maximum PMT sensitivity

  /* These variables are for the camera mounting structure */
  double base[GEOFD_MAXMIR][NUMTAPOLES][3];   // vector to base of each pole
  double pole[GEOFD_MAXMIR][NUMTAPOLES][3];   // vector in direction of each pole
  double radius[NUMTAPOLES];                  // radius of each pole (meters)
  double length[NUMTAPOLES];                  // length of each pole (meters)
} TGeom;

typedef struct {
  int nthrow;                           // number of photons / npe to throw
  int ngood;                            // number of hit tubes
  int cam;                              // camera to shoot at
  double vsite[3];                      // vector from cam center of
                                        // curvature to light source (meters)
                                        // (camera coordinate system)

  double xcam, ycam;                    // location of ray on camera (camera
                                        // coordinates [see geofd_dst.h])
  double time;                          // time of flight of photon from
                                        // source to PMT (ns)

  int pmt[MAXNPE];                      // array of tube IDs or error codes
  double tpe[MAXNPE];                   // array of travel times (ns)

  double xave, yave;                    // average location for rays on
                                        // PMT plane
} RayTrace;
      
/* in "initRayTrace.c" */
#ifdef __cplusplus
extern "C"{
#endif
void initRayTrace (geofd_dst_common *g, TGeom *tgm, UCalibration *calib);
void clearRayTrace (RayTrace *ray);

/* in "doRayTrace.c" */
void doRayTrace (int num, int cam, double rpuv[],
             double npln[], double vsite[],
             double age, double mol,
             geofd_dst_common *geo, TGeom *tgm, RayTrace *ray);

/* in "trace.c" */
int trace (geofd_dst_common *geo, TGeom *tgm, RayTrace *ray);

/* in "acptmap.c" */
int acptmap(TGeom *tgm, RayTrace *ray);

/* in "shadow.c" */
void buildTACameraPoles (geofd_dst_common *geo, TGeom *tgm);
void buildMDCameraPoles (geofd_dst_common *geo, TGeom *tgm);
int hitCameraBox (geofd_dst_common *geo, TGeom *tgm, RayTrace *ray,
                  double v[], double w[]);
int hitFlasherBox (geofd_dst_common *geo, TGeom *tgm, RayTrace *ray,
                   double v[], double w[]);
                   
int hitCameraPoles (int num, TGeom *tgm, RayTrace *ray, double u[], double v[]);

/* in "hitMirror.c" */
int hitTAMirror (geofd_dst_common *geo, TGeom *tgm, int cam, double u[]);
int hitMDMirror (geofd_dst_common *geo, TGeom *tgm, int cam, double u[]);

/* in "camera.c" */
void snellShift (TGeom *tgm, double n2, double thick, double g[], double s[]);
int pmtDetect (geofd_dst_common *geo, TGeom *tgm, RayTrace *r);
double pmtQE (geofd_dst_common *geo, TGeom *tgm, RayTrace *r, int tube);

#ifdef __cplusplus
}
#endif

#endif
