#ifndef _CONTROL_H_
#define _CONTROL_H_

#define newInstanceOf(type) malloc(sizeof(type))

#ifndef QUIET_PRINT
#define  pout(format, args...) fprintf(stdout, "%20s(): "format, __FUNCTION__, ##args)
#else
#define  pout(format, args...) fprintf(stderr, "%20s(): "format, __FUNCTION__, ##args)
#endif
#define  perr(format, args...) fprintf(stderr, "%20s(): "format, __FUNCTION__, ##args)
#define vperr(format, args...) fprintf(stderr, "%s:%d: %s(): "format, __FILE__, __LINE__, __FUNCTION__, ##args)

#define NCONFIG_MAX 3

#define MAX_ARGS   5
#define MAX_STRLEN 256
#define MAX_NBREAK 5
#define MAX_PART   100
#define DEFAULT_DT 1.0  /* seconds */

#define MAX_LIBRARIES   3   /* will typically be '2' (p, Fe) */

#define DEFAULT_ATMOSDB    "atmosphere/raobslc_por20141106.dst.gz"
#define DEFAULT_SCATTERDB  "atmosphere/fdscat.dst.gz"

#ifdef TRACK_MEMORY
#include "memtrack.h"
#endif

/* runtime parameters */
typedef struct {
  /*
   *  RUNTIME-SPECIFIC PARAMETERS
   *
   *  This includes output DST file names, number of trials/successes, random 
   *    number seeds, generally parameters that are independent of detector 
   *    type/location and the EAS in general.
   */
  char outfn[MAX_STRLEN];   /* the output file name */
  char pdirn[MAX_STRLEN];   /* the parent directory of the output file(s) */
  int trigid;               /* ID number for next trigger */
  int ntry;                 /* number of trials to take before quitting */
  int nevt;                 /* number of events to collect before quitting */
  int seed;                 /* random number seed */
  double origin[3];
  char infn[MAX_STRLEN];    /* user-input event list */
  char ontime[MAX_STRLEN];  /* fdped file for on-times by part number */
  double dt;

  /*
   *  SHOWER-SPECIFIC PARAMETERS
   *
   *  This includes parameters that are independent of detector type/location. 
   *    These will define the EAS as far as primary particle species, and 
   *    provides TRUMP with limits on particle energy & geometry.  Note that 
   *    'rplo' & 'rphi' are stored here.  This variable exists in both SD and 
   *    FD analysis, but with different meanings.  Therefore, it is 
   *    appropriately placed here.
   */
  int species[MAX_LIBRARIES];  /* primary particle species (uses Corsika convention) */
  char showlibFile[MAX_LIBRARIES][MAX_STRLEN];  /* the name of the shower library directory */
  char atmosdbFile[MAX_STRLEN];
  char scatterFile[MAX_STRLEN];
  double logelo, logehi;
  double ebreak[MAX_NBREAK+1];
  double eslope[MAX_NBREAK+1];

#ifdef LASER_OPTION
  double laser_energy;
  double laser_lambda;
  double laser_latitude;
  double laser_longitude;
  double laser_altitude;
#endif

  double rplo, rphi;
  double philo, phihi;
  double thetalo, thetahi;
  double lat;                     /* angular distance from CLF to throw evts */

  /*
   *  SD-SPECIFIC PARAMETERS
   *
   *  None as of yet...
   */

  /*
   *  FD-SPECIFIC PARAMETERS
   *
   *  This includes site identification & geometry, calibration file names, 
   *    stricter geometry boundaries like limits on psi & impact locations 
   *    (as an angle as seen from a FD site).
   */
  int siteid;  /* 0=BR, 1=LR, 2=MD, 3=TL */

  char geometryFile[MAX_STRLEN];  /* the name of the geometry file to use */

  /* calibration file names */
  char mirrefFile[MAX_STRLEN];    /* name of the fdmir_ref DST file */
  char paraglasFile[MAX_STRLEN];  /* name of the fdparaglas_trans DST file */
  char bg3transFile[MAX_STRLEN];  /* name of the fdbg3_trans DST file */
  char pmtQEFile[MAX_STRLEN];     /* name of the fdpmt_qece DST file */
  char pmtGainFile[MAX_STRLEN];   /* name of the fdpmt_gain DST file */
  char pmtUnifFile[MAX_STRLEN];   /* name of the fdpmt_uniformity DST file */
  char pmtCalFile[MAX_STRLEN];    /* name of the fdped DST file */

  double psilo, psihi;
  double phiimplo, phiimphi;

  double nsbackground;  /* night sky background level in 1/100 ns */

  /*
   *  SPECIAL PARAMETERS
   *
   *  These are variables that need to be kept in memory, but are not set in 
   *    the configuration file per se.  These are not organized in any 
   *    particular order.
   */

  int nshowlib;  /* number of attached shower library files */

  /* Flags for DST banks */
  int flag_mirref;
  int flag_paraglas;
  int flag_bg3;
  int flag_pmtQE;
  int flag_pmtGain;
  int flag_pmtUnif;
  int flag_pmtCal;

  /* Flags for date/event format */
  int flag_eventlist;
  int flag_ontime;
  int flag_today;

  int jday;     // integer julian day
  int jsec;     // seconds after jday
  int year, month, day, hour, min, sec, nsec;

  int tpart;                /* target part number (-1 for any) */
  int part[MAX_PART], ipart, numparts;
  int t0day[MAX_PART], t0sec[MAX_PART], t0nsec[MAX_PART];
  int t1day[MAX_PART], t1sec[MAX_PART], t1nsec[MAX_PART];

  double oc2ecc[3][3];  /* transformation matrix to go from "Origin Coordinates" to "Earth-Center Coordinates" */

  /*
   *  Variables used to throw particles on a broken energy spectrum
   */
  int nbreak;

  /* shower track geometry parameters */
} RuntimeParameters;

/* from "control.c" */
#ifdef __cplusplus
extern "C"{
#endif
int loadConfigurationFile(const char *cfn, RuntimeParameters *par);


int initializeParameters(RuntimeParameters *par);

int resolveConfiguration(RuntimeParameters *par);

void dumpConfTags(void);

void writeHeader (char *conf, RuntimeParameters *par);
void dumpConfiguration(const RuntimeParameters *par);

int setupParameters(RuntimeParameters *par);
#ifdef __cplusplus
}
#endif

#endif
