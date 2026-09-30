#ifndef _TACALIBRATION_H_
#define _TACALIBRATION_H_

typedef struct {
  /*
   * pmtgain requires bg3 and QE information.
   * pmt pedestal values require pmtgain information.
   * These are flags to indicate that bg3, pmtqe and pmtgain 
   * information exists.
   */
  int mask;

  double mirref[GEOFD_MAXMIR][NWAVELEN_BANDS];  /* mirror reflectivity */
  double paraglas[NWAVELEN_BANDS];  /* Paraglas transmission */
  double bg3[NWAVELEN_BANDS];  /* BG3 filter transmission */
  double pmtmaxqe;  /* Maximum PMT QE */
  double pmtrelqe[NWAVELEN_BANDS];  /* Relative PMT QE */
  double pmtgain[GEOFD_MAXMIR][GEOFD_MIRTUBE];  /* PMT Gain (FADC/NPE) */
  double pmtmaxunif;  /* Maximum PMT uniformity */
  double pmtrelunif[PMT_NUMDIV][PMT_NUMDIV];  /* Relative PMT uniformity */
  int mean[GEOFD_MAXMIR][GEOFD_MIRTUBE];  /* PMT mean */
  int vari[GEOFD_MAXMIR][GEOFD_MIRTUBE];  /* PMT variance */
  int pedestal[GEOFD_MAXMIR][GEOFD_MIRTUBE];  /* PMT pedestal */
  int liveflag[GEOFD_MAXMIR][GEOFD_MIRTUBE];  /* PMT live flag */

  double reductFactor[GEOFD_MAXMIR][NWAVELEN_BANDS];
  
  int nmir;         // number of mirrors in site geometry file
} TACalibration;

typedef struct {
  int unit;
  int wbl;
  int hbl;
  char path[MAX_STRLEN];
  boolean isOpen;
} TACalibrationFile;

#define MAX_CALIBRATION_UNITS 10
#define MAX_WARN 10
#define MAX_REL_GAINERROR 0.25

#define PARAGLAS_FLAG 0x1
#define BG3_FLAG      0x2
#define PMTQE_FLAG    0x4
#define PMTUNIF_FLAG  0x8
#define MIRREF_FLAG   0x10
#define PMTGAIN_FLAG  0x20
#define PEDESTAL_FLAG 0x40
#define REDUCT_FLAG   0x80

/* in "calibration.c" */
#ifdef __cplusplus
extern "C"{
#endif

/* Global loaders */
int getTADefaultCalibration(RuntimeParameters *par, TACalibration *cal);
int getTATimeIndependentCalibration(RuntimeParameters *par, TACalibration *cal);
int getTATimeDependentCalibration(RuntimeParameters *par, TACalibration *cal);
int getTACalibrationReduction(TACalibration *cal);
void detachTACalibration(void);

/* lower level routines */
TACalibrationFile *attachTACalibrationFile(int fd, char *path);
int detachTACalibrationFile(int fd);

int getTAMirrorReflectivity(TACalibrationFile *cf, RuntimeParameters *par, TACalibration *cal);
int getTAPmtGains(TACalibrationFile *cf, RuntimeParameters *par, TACalibration *cal);
int getTAPedestals(TACalibrationFile *cf, RuntimeParameters *par, TACalibration *cal);

int getTADefaultParaglasTrans(TACalibration *cal);
int getTADefaultBG3Trans(TACalibration *cal);
int getTADefaultPmtQE(TACalibration *cal);
int getTADefaultPmtUniformity(TACalibration *cal);

int getTADefaultMirrorReflectivity(RuntimeParameters *par, TACalibration *cal);
int getTADefaultPmtGains(RuntimeParameters *par, TACalibration *cal);
int getTADefaultPedestals(RuntimeParameters *par, TACalibration *cal);

#ifdef __cplusplus
}
#endif
#endif
