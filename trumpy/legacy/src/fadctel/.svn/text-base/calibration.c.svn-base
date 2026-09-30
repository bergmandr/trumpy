#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "control.h"
#include "constants.h"
#include "fdconstants.h"
#include "event.h"
#include "convtime.h"
#include "calibration.h"

#include "paraglas.h"
#include "bg3.h"
#include "pmt_qe.h"
#include "pmt_uniformity.h"
#include "mirref.h"
#include "pmt_gain.h"
#include "pmt_calibration.h"

#ifdef MIRREF_CORRECTION
#warning "USING MIRROR REFLECTIVITY CORRECTION!!!"
#endif

static int ncalFile = 0;
static CalibrationFile calFile[MAX_CALIBRATION_UNITS] = {
  { 0, 0, 0, "", FALSE }
};
static boolean globalLiveFlag[2][GEOFD_MAXMIR][GEOFD_MIRTUBE];

/*
 *  if this works, then I'll replace it in stools.c
 */
static int _endswith(char *str, const char *sfx) {
  int n = strlen(sfx);
  char *p = strstr(str, sfx);
  if ( p == NULL )
    return 0;

  if ( strlen(p) == n )
    return 1;
  else
    return 0;
}

/*
 *  These should be the most commonly called functions
 */

/*
 *  Fills Calibration structure pointed to by 'cal' with the default 
 *    calibration
 */
int getDefaultCalibration(RuntimeParameters *par, Calibration *cal) {

  getDefaultParaglasTrans(cal);
  getDefaultBG3Trans(cal);
  getDefaultPmtQE(cal);
  getDefaultPmtUniformity(cal);

  getDefaultMirrorReflectivity(par, cal);
  getDefaultPmtGains(par, cal);
  getDefaultPedestals(par, cal);

  getCalibrationReduction(cal);

  return cal->mask;
}

/*
 *  Does one-time loading of time-independent calibration data
 */
int getTimeIndependentCalibration(RuntimeParameters *par, Calibration *cal) {
  int i, j, k, buf;
  char cdir[MAX_STRLEN], path[MAX_STRLEN];
  CalibrationFile *cf;

  sprintf(cdir, "%s/calibration/%s", RTDATA, CALIBRATION_VERSION);

  /*
   *  LOAD PARAGLAS TRANSMITTANCE
   */
  if ( strcmp(par->paraglasFile, "(none)") == 0 )
    sprintf(path, "%s/paraglasDSTBank.dst.gz", cdir);
  else
    strcpy(path, par->paraglasFile);

  for ( i=0; i<NWAVELEN_BANDS; i++ )
    cal->paraglas[i] = 0.0;

  if ( (cf=attachCalibrationFile(PARAGLAS_FILEDES, path)) != NULL ) {
    eventRead(PARAGLAS_FILEDES, cf->wbl, cf->hbl, &buf);
    double nent[NWAVELEN_BANDS] = { 0.0 };

    for ( i=0; i<fdparaglas_trans_.nLambda; i++ ) {
      double lam = fdparaglas_trans_.minLambda + 
	(double)i*fdparaglas_trans_.deltaLambda;

      k = (int)((1.0e-9*lam-LAMBDA0)/DLAMBDA);
      if ( (k>=0) && (k<NWAVELEN_BANDS) ) {
	cal->paraglas[k] += fdparaglas_trans_.transmittance[i];
	nent[k] += 1.0;
      }
    }

    for ( i=0; i<NWAVELEN_BANDS; i++ ) {
      if ( nent[i] > 0.5 )
	cal->paraglas[i] /= nent[i];
      else
	cal->paraglas[i] = 0.0;
    }

    detachCalibrationFile(PARAGLAS_FILEDES);
    cal->mask |= PARAGLAS_FLAG;
  }

  /*
   *  LOAD BG3 TRANSMITTANCE
   */
  if ( strcmp(par->bg3transFile, "(none)") == 0 )
    sprintf(path, "%s/bg3DSTBank.dst.gz", cdir);
  else
    strcpy(path, par->bg3transFile);

  for ( i=0; i<NWAVELEN_BANDS; i++ )
    cal->bg3[i] = 0.0;

  if ( (cf=attachCalibrationFile(BG3_FILEDES, path)) != NULL ) {
    eventRead(BG3_FILEDES, cf->wbl, cf->hbl, &buf);
    double nent[NWAVELEN_BANDS] = { 0.0 };

    for ( i=0; i<fdbg3_trans_.nLambda; i++ ) {
      double lam = fdbg3_trans_.minLambda + 
	(double)i*fdbg3_trans_.deltaLambda;

      k = (int)((1.0e-9*lam-LAMBDA0)/DLAMBDA);
      if ( (k>=0) && (k<NWAVELEN_BANDS) ) {
	cal->bg3[k] += fdbg3_trans_.transmittance[i];
	nent[k] += 1.0;
      }
    }

    for ( i=0; i<NWAVELEN_BANDS; i++ ) {
      if ( nent[i] > 0.5 )
	cal->bg3[i] /= nent[i];
      else
	cal->bg3[i] = 0.0;
    }

    detachCalibrationFile(BG3_FILEDES);
    cal->mask |= BG3_FLAG;
  }

  /*
   *  LOAD PMT QE AND CE
   */
  if ( strcmp(par->pmtQEFile, "(none)") == 0 )
    sprintf(path, "%s/pmtQECEDSTBank.dst.gz", cdir);
  else
    strcpy(path, par->pmtQEFile);

  cal->pmtmaxqe = 0.0;
  for ( i=0; i<NWAVELEN_BANDS; i++ )
    cal->pmtrelqe[i] = 0.0;

  if ( (cf=attachCalibrationFile(PMTQE_FILEDES, path)) != NULL ) {
    eventRead(PMTQE_FILEDES, cf->wbl, cf->hbl, &buf);
    double nent[NWAVELEN_BANDS] = { 0.0 };

    for ( i=0; i<fdpmt_qece_.qeNLambda; i++ ) {
      double lam = fdpmt_qece_.qeMinLambda + 
	(double)i*fdpmt_qece_.qeDeltaLambda;

      k = (int)((1.0e-9*lam-LAMBDA0)/DLAMBDA);
      if ( (k>=0) && (k<NWAVELEN_BANDS) ) {
	cal->pmtrelqe[k] += fdpmt_qece_.ce * fdpmt_qece_.qe[i];
	nent[k] += 1.0;
      }
    }

    for ( i=0; i<NWAVELEN_BANDS; i++ ) {
      if ( nent[i] > 0.5 )
	cal->pmtrelqe[i] /= nent[i];
      else
	cal->pmtrelqe[i] = 0.0;

      if ( cal->pmtrelqe[i] > cal->pmtmaxqe )
	cal->pmtmaxqe = cal->pmtrelqe[i];
    }

    for ( i=0; i<NWAVELEN_BANDS; i++ )
      cal->pmtrelqe[i] /= cal->pmtmaxqe;

    detachCalibrationFile(PMTQE_FILEDES);
    cal->mask |= PMTQE_FLAG;
  }

  /*
   *  LOAD PMT UNIFORMITY MAP
   */
  if ( strcmp(par->pmtUnifFile, "(none)") == 0 )
    sprintf(path, "%s/pmtUniformityDSTBank.dst.gz", cdir);
  else
    strcpy(path, par->pmtUnifFile);

  if ( (cf=attachCalibrationFile(PMTUNIF_FILEDES, path)) != NULL ) {
    eventRead(PMTUNIF_FILEDES, cf->wbl, cf->hbl, &buf);

    cal->pmtmaxunif = 0.0;
    for ( i=0; i<fdpmt_uniformity_.xNDivision; i++ ) {
      for ( j=0; j<fdpmt_uniformity_.yNDivision; j++ ) {
	cal->pmtrelunif[i][j] = fdpmt_uniformity_.uniformity[i][j];
	if ( cal->pmtrelunif[i][j] > cal->pmtmaxunif )
	  cal->pmtmaxunif = cal->pmtrelunif[i][j];
      }
    }

    for ( i=0; i<fdpmt_uniformity_.xNDivision; i++ )
      for ( j=0; j<fdpmt_uniformity_.yNDivision; j++ )
	cal->pmtrelunif[i][j] /= cal->pmtmaxunif;

    detachCalibrationFile(PMTUNIF_FILEDES);
    cal->mask |= PMTUNIF_FLAG;
  }

  return cal->mask;
}

/*
 *  This routine needs to be a little smarter...
 *    needs to check if new file will need to be opened
 *    needs to check if input file is a directory, then decide which new file
 *      to open.
 */
int getTimeDependentCalibration(RuntimeParameters *par, Calibration *cal) {
  int i;
  char cdir[MAX_STRLEN], path[MAX_STRLEN];
  CalibrationFile *cf;

  // get the mirror reflectivity first
  sprintf(cdir, "%s/calibration/cal_%s", RTDATA, CALIBRATION_VERSION);

  /*
   *  LOAD MIRROR REFLECTIVITY
   */
  // count on the year, month, day and time to be defined in *par already
  if ( _endswith(par->mirrefFile, ".dst.gz") || 
       _endswith(par->mirrefFile, ".dst") ) {
    strcpy(path, par->mirrefFile);
  }
  else if ( strcmp(par->mirrefFile, "(none)") != 0 ) {
    // mirrefFile points to a containing folder
    sprintf(path, "%s/mirrorDSTBank.%4d.dst.gz", par->mirrefFile, par->year);
  }
  else if ( strcmp(par->mirrefFile, "DISABLED") != 0 ) {
    // sprintf(path, "%s/fdmirror_ref/mirrorDSTBank.%4d.dst.gz", cdir, par->year);  // this has problems, because events may span from one year into the next.
    sprintf(path, "%s/calibration/%s/mirrorDSTBank.dst.gz", RTDATA, CALIBRATION_VERSION);
  }

  // check if the file channel already exists
  for ( i=0; i<ncalFile; i++ )
    if ( calFile[i].unit == MIRROR_FILEDES )
      break;

  // file channel is open and points to the correct file
  if ( strcmp(par->mirrefFile, "DISABLED") != 0 ) {
    if ( calFile[i].isOpen && (strcmp(calFile[i].path,path)==0) ) {
      cal->mask |= getMirrorReflectivity(&calFile[i], par, cal);
    }
    else {
      // file channel is open but points to the wrong file
      if ( calFile[i].isOpen )
	detachCalibrationFile(MIRROR_FILEDES);

      if ( (cf=attachCalibrationFile(MIRROR_FILEDES, path)) != NULL )
	cal->mask |= getMirrorReflectivity(cf, par, cal);
      else
	cal->mask ^= MIRREF_FLAG;   // if MIRREF_FLAG bit is on, turn it off
    }

    if ( (cal->mask&MIRREF_FLAG) == 0 ) {
      perr("error getting mirror reflectivity.\n");
      // revert to default?
    }
  }

  /*
   *  LOAD PMT GAINS
   */
  if ( _endswith(par->pmtGainFile, ".dst.gz") || 
       _endswith(par->pmtGainFile, ".dst") ) {
    strcpy(path, par->pmtGainFile);
  }
  else if ( strcmp(par->pmtGainFile, "(none)") != 0 ) {
    // pmtGainFile points to a containing folder
    sprintf(path, "%s/pmtGainDSTBank.%4d%02d.dst.gz", 
	    par->pmtGainFile, par->year, par->month);
  }
  else {
    sprintf(path, "%s/fdpmt_gain/pmtGainDSTBank.%4d%02d.dst.gz", 
	    cdir, par->year, par->month);
  }

  // check if the file channel already exists
  for ( i=0; i<ncalFile; i++ )
    if ( calFile[i].unit == PMTGAIN_FILEDES )
      break;

  // file channel is open and points to the correct file
  if ( calFile[i].isOpen && (strcmp(calFile[i].path,path)==0) ) {
    cal->mask |= getPmtGains(&calFile[i], par, cal);
  }
  else {
    // file channel is open but points to the wrong file
    if ( calFile[i].isOpen )
      detachCalibrationFile(PMTGAIN_FILEDES);

    if ( (cf=attachCalibrationFile(PMTGAIN_FILEDES, path)) != NULL )
      cal->mask |= getPmtGains(&calFile[i], par, cal);
    else
      cal->mask ^= PMTGAIN_FLAG;  // if PMTGAIN_FLAG bit is on, turn it off
  }

  if ( (cal->mask&PMTGAIN_FLAG) == 0 ) {
    perr("error getting PMT gains.\n");
    // revert to default?
  }

  /*
   *  LOAD PEDESTALS
   */
  int fd;
  char site[MAX_STRLEN];
  if ( par->siteid == 0 ) {
    fd = PMTCAL0_FILEDES;
    strcpy(site, "black-rock");
  }
  else if ( par->siteid == 1 ) {
    fd = PMTCAL1_FILEDES;
    strcpy(site, "long-ridge");
  }
/*   else if ( par->siteid == 2 ) { */
/*     fd = PMTCAL2_FILEDES; */
/*     perr("error: pedestals not defined for MD yet.\n"); */
/*     return cal->mask; */
/*   } */
  else {
    vperr("invalid site ID (%d).  Pedestals not read.\n", par->siteid);
    return cal->mask;
  }

  if ( _endswith(par->pmtCalFile, ".dst.gz") || 
       _endswith(par->pmtCalFile, ".dst") ) {
    strcpy(path, par->pmtCalFile);
  }
  else if ( strcmp(par->pmtCalFile, "(none)") != 0 ) {
    // pmtCalFile points to a containing folder
    sprintf(path, "%s/%s/y%4dm%02dd%02d.ped.dst.gz", 
	    par->pmtCalFile, site, par->year, par->month, par->day);
  }
  else if ( strcmp(par->pmtCalFile, "DISABLED") != 0 ) {
    sprintf(path, "%s/calibration/fdped/%s/y%4dm%02dd%02d.ped.dst.gz", 
	    RTDATA, site, par->year, par->month, par->day);
  }

  if ( strcmp(par->pmtCalFile, "DISABLED") != 0 ) {
    // check if the file channel already exists
    for ( i=0; i<ncalFile; i++ )
      if ( calFile[i].unit == fd )
	break;

    // file channel is open and points to the correct file
    if ( calFile[i].isOpen && (strcmp(calFile[i].path,path)==0) ) {
      cal->mask |= getPedestals(&calFile[i], par, cal);
    }
    else {
      // file channel is open but points to the wrong file
      if ( calFile[i].isOpen )
	detachCalibrationFile(fd);

      if ( (cf=attachCalibrationFile(fd, path)) != NULL )
	cal->mask |= getPedestals(cf, par, cal);
      else
	cal->mask ^= PEDESTAL_FLAG;  // if PEDESTAL_FLAG bit is on, turn it off
    }

    if ( (cal->mask&PEDESTAL_FLAG) == 0 ) {
      perr("error getting pedestal data.\n");
      // get defaults?
    }
  }

  return cal->mask;
}

int getCalibrationReduction(Calibration *cal) {
  int i, j;

  for (i=0; i<GEOFD_MAXMIR; i++)
    for (j=0; j<NWAVELEN_BANDS; j++)
      cal->reductFactor[i][j] = cal->mirref[i][j] * 
                                cal->paraglas[j]  * 
                                cal->bg3[j]       *
                                cal->pmtmaxunif   *
                                cal->pmtrelqe[j]  *
                                cal->pmtmaxqe;

  cal->mask |= REDUCT_FLAG;
  return cal->mask;
}

void detachCalibration(void) {
  int i;
  for ( i=0; i<ncalFile; i++ ) {
    if ( calFile[i].isOpen ) {
      dstCloseUnit(calFile[i].unit);
      delBankList(calFile[i].wbl);
      delBankList(calFile[i].hbl);
    }
    calFile[i].unit = 0;
    calFile[i].wbl = 0;
    calFile[i].hbl = 0;
    strcpy(calFile[i].path, "(none)");
    calFile[i].isOpen = FALSE;
  }
}

/*
 *  These functions should not be called outside of this source file, but they 
 *    will be open so that users have more control of lower level behavior.
 */
CalibrationFile *attachCalibrationFile(int fd, char *path) {
  int i;
  CalibrationFile *cf = NULL;

  // check if this file channel is already set up.
  for ( i=0; i<ncalFile; i++ )
    if ( calFile[i].unit == fd )
      break;

  // if the file channel already exists, check if it's open
  if ( (i<ncalFile) && calFile[i].isOpen ) {
    perr("warning: this file channel is already open!\n");
    // should it be automatically closed and re-opened?
    return NULL;
  }

  switch ( fd ) {
  case PARAGLAS_FILEDES:
  case BG3_FILEDES:
  case PMTQE_FILEDES:
  case PMTUNIF_FILEDES:
    // these must be direct paths to the DST file
    if ( dstOpenUnit(fd, path, MODE_READ_DST) == 0 ) {
      calFile[i].unit = fd;
      calFile[i].wbl = newBankList(100);
      calFile[i].hbl = newBankList(100);
      strcpy(calFile[i].path, path);
      calFile[i].isOpen = TRUE;
      cf = &calFile[i];
      ncalFile++;
    }
    break;

  case MIRROR_FILEDES:
  case PMTGAIN_FILEDES:
  case PMTCAL0_FILEDES:
  case PMTCAL1_FILEDES:
  case PMTCAL2_FILEDES:
    // these may be the direct paths to a DST file or containing folder
    if ( _endswith(path, ".dst.gz") || _endswith(path, ".dst") ) {
      if ( dstOpenUnit(fd, path, MODE_READ_DST) == 0 ) {
	calFile[i].wbl = newBankList(100);
	calFile[i].hbl = newBankList(100);
	calFile[i].isOpen = TRUE;
	cf = &calFile[i];
	if ( i == ncalFile ) {
	  calFile[i].unit = fd;
	  strcpy(calFile[i].path, path);
	  ncalFile++;
	}
	// and make sure we've read in the first event
	int buf;
	eventRead(fd, calFile[i].wbl, calFile[i].hbl, &buf);
      }
    }
    else {
      // if it's only a containing folder, then save a place for it and wait 
      //   for the next call to this function to actually open it.
      calFile[i].unit = fd;
      strcpy(calFile[i].path, path);
      ncalFile++;
      cf = NULL;
    }
    break;

  default:
    perr("file desc. %d does not represent a valid calibration type.\n", fd);
    cf = NULL;
    ;;
  }

  return cf;
}

int detachCalibrationFile(int fd) {
  int i;

  for ( i=0; i<ncalFile; i++ )
    if ( calFile[i].unit == fd )
      break;

  dstCloseUnit(fd);
  delBankList(calFile[i].wbl);
  delBankList(calFile[i].hbl);
  calFile[i].unit = 0;
  sprintf(calFile[i].path, "(none)");
  calFile[i].isOpen = FALSE;

  ncalFile--;
  for ( ; i<ncalFile; i++ )
    calFile[i] = calFile[i+1];

  return 0;
}

int getDefaultParaglasTrans(Calibration *cal) {
  pout("setting default Paraglas transmittance data.\n");
  int i;
  for ( i=0; i<NWAVELEN_BANDS; i++ )
    cal->paraglas[i] = paraglas_default[i];
  return cal->mask |= PARAGLAS_FLAG;
}

int getDefaultBG3Trans(Calibration *cal) {
  pout("setting default BG3 filter tranmittance data.\n");
  int i;
  for ( i=0; i<NWAVELEN_BANDS; i++ )
    cal->bg3[i] = bg3_default[i];
  return cal->mask |= BG3_FLAG;
}

int getDefaultPmtQE(Calibration *cal) {
  pout("setting default PMT QE data.\n");
  int i;
  cal->pmtmaxqe = pmtmaxqe_default;
  for ( i=0; i<NWAVELEN_BANDS; i++ )
    cal->pmtrelqe[i] = pmtrelqe_default[i];
  return cal->mask |= PMTQE_FLAG;
}

int getDefaultPmtUniformity(Calibration *cal) {
  pout("setting default PMT uniformity data.\n");
  int i, j;

  cal->pmtmaxunif = 0.00;
  for ( i=0; i<PMT_NUMDIV; i++ ) {
    for ( j=0; j<PMT_NUMDIV; j++ ) {
      cal->pmtrelunif[i][j] = pmtunif_default[i][j];
      if ( cal->pmtrelunif[i][j] > cal->pmtmaxunif )
	cal->pmtmaxunif = cal->pmtrelunif[i][j];
    }
  }

  for ( i=0; i<PMT_NUMDIV; i++ )
    for ( j=0; j<PMT_NUMDIV; j++ )
      cal->pmtrelunif[i][j] /= cal->pmtmaxunif;

  return cal->mask |= PMTUNIF_FLAG;
}

int getDefaultMirrorReflectivity(RuntimeParameters *par, Calibration *cal) {
  pout("setting default mirror reflectivity data.\n");
  int i, j;
  for ( i=0; i<GEOFD_MAXMIR; i++ )
    for ( j=0; j<NWAVELEN_BANDS; j++ )
      cal->mirref[i][j] = 
#ifdef MIRREF_CORRECTION
	MIRREF_CORRECTION * 
#endif
	mirref_default[par->siteid][i][j];
  return cal->mask |= MIRREF_FLAG;
}

int getDefaultPmtGains(RuntimeParameters *par, Calibration *cal) {
  pout("setting default PMT gain data.\n");
  int i, j, k;
  double bg3qe = 0.223393;
  if ( (cal->mask&(BG3_FLAG|PMTQE_FLAG)) == (BG3_FLAG|PMTQE_FLAG) ) {
    k = (int)((1e-9*CALIB_LAMBDA-LAMBDA0)/DLAMBDA);
    if ( (k>=0) && (k<NWAVELEN_BANDS) )
      bg3qe = cal->bg3[k] * cal->pmtmaxqe * cal->pmtrelqe[k];
  }

  for ( i=0; i<GEOFD_MAXMIR; i++ )
    for ( j=0; j<GEOFD_MIRTUBE; j++ )
      cal->pmtgain[i][j] = pmtgain_default[par->siteid][i][j] * 
	RXF_FACTOR / bg3qe;

  return cal->mask |= PMTGAIN_FLAG;
}

int getDefaultPedestals(RuntimeParameters *par, Calibration *cal) {
  int i, j;

  if ( (cal->mask&PMTGAIN_FLAG) == 0 ) {
    perr("ERROR: PMT gain data must be defined before calling this function.\n");
    return 0x0;
  }

  pout("setting default pedestals.\n");

  for ( i=0; i<GEOFD_MAXMIR; i++ ) {
    for ( j=0; j<GEOFD_MIRTUBE; j++ ) {
      cal->mean[i][j] = pedestal_default[par->siteid][i][j];
      cal->vari[i][j] = pedrms_default[par->siteid][i][j];
      cal->pedestal[i][j] = pedestal_default[par->siteid][i][j] - 
	(int)(pedrms_default[par->siteid][i][j]*RHO_ALPHA/cal->pmtgain[i][j]);
      cal->pedestal[i][j] = cal->pedestal[i][j]/16 + 1;
    }
  }

  return cal->mask |= PEDESTAL_FLAG;
}

int getMirrorReflectivity(CalibrationFile *cf, RuntimeParameters *par, Calibration *cal) {
  int i, j, k, buf;
  int eposec, istat, iter = 0;
  double mjlday;
  fdmirror_ref_dst_common *bank = &fdmirror_ref_;
  static boolean needsUpdate = TRUE;

  if ( cf == NULL ) {
    perr("ERROR: Mirror reflectivity file is not attached\n");
    return 0x0;
  }

  if ( ! cf->isOpen ) {
    perr("ERROR: Mirror reflectivity file is not open (path=%s)\n", cf->path);
    return 0x0;
  }

  mjlday = (double)par->jsec/SECPDAY + (double)par->jday - MJLDOFF;
  eposec = (int)mjlday2eposec(mjlday);

  while ( iter < 2 ) {

    if ( eposec >= bank->dateTo ) {
      // advance forward
      istat = eventRead(cf->unit, cf->wbl, cf->hbl, &buf);
      needsUpdate = TRUE;
    }
    else if ( eposec >= bank->dateFrom ) {
      // we're in the right time frame now
      if ( ! needsUpdate )
	return MIRREF_FLAG;

      double ref[GEOFD_MAXMIR][300] = {{ 0.00 }};

      /* find av. reflectivity from mirror segments */
      for ( i=0; i<bank->nTelescope; i++ ) {
	for ( j=0; j<bank->nMirror; j++ ) {
	  fdmirror_ref_data *data = &(bank->mirrorData[par->siteid][i][j]);
	  for ( k=0; k<bank->nLambda; k++ )
	    ref[i][k] += data->reflection[k] / (double)bank->nMirror;
	}
      }

      for ( i=0; i<bank->nTelescope; i++ )
	for ( j=0; j<NWAVELEN_BANDS; j++ )
	  cal->mirref[i][j] = 0.00;

      /*
       * Find average reflectivity for an entire telescope
       *   in NWAVELEN_BANDS bins.
       */
      for ( i=0; i<bank->nTelescope; i++ ) {
	double lam;
	double nent[NWAVELEN_BANDS] = { 0.00 };
	for ( j=0; j<bank->nLambda; j++ ) {
	  lam = 1e-9 * ( bank->minLambda + (double)j*bank->deltaLambda );
	  k = (int)((lam-LAMBDA0)/DLAMBDA);
	  if ( (k>=0) && (k<NWAVELEN_BANDS) ) {
	    cal->mirref[i][k] += ref[i][j];
	    nent[k]++;
	  }
	}
	for ( k=0; k<NWAVELEN_BANDS; k++ ) {
	  if ( nent[k] > 0.5 )
	    cal->mirref[i][k] = cal->mirref[i][k]
#ifdef MIRREF_CORRECTION
	      * MIRREF_CORRECTION 
#endif
	      / nent[k];
	  else
	    cal->mirref[i][k] = 0.0;
	}
      }

      needsUpdate = FALSE;
      return MIRREF_FLAG;
    }
    else {
      // DST file has advanced too far, force it to re-open and read from the
      //   beginning.
      istat = -1;
      needsUpdate = TRUE;
    }

    if ( istat < 0 ) {
      // end of DST file reached.  try re-reading it from the beginning
      dstCloseUnit(cf->unit);
      dstOpenUnit(cf->unit, cf->path, MODE_READ_DST);
      iter++;
    }
  }

  // only gets here when data is not found
  perr("WARNING: Mirror reflectivity data for current time not found in attached DST file.\n");
  return 0x0;
}

int getPmtGains(CalibrationFile *cf, RuntimeParameters *par, Calibration *cal) {
  int i, j, k;
  int eposec, istat, buf, iter = 0;
  double mjlday;
  fdpmt_gain_dst_common *bank = &fdpmt_gain_;
  static boolean needsUpdate = TRUE;
  static double bg3qe = 0.223393;  // default BG3*QE*CE @ 335nm

  if ( cf == NULL ) {
    perr("ERROR: PMT gains file is not attached\n");
    return 0x0;
  }

  if ( ! cf->isOpen ) {
    perr("ERROR: PMT gains file is not open (path=%s)\n", cf->path);
    return 0x0;
  }

  if ( (cal->mask&(BG3_FLAG|PMTQE_FLAG)) == (BG3_FLAG|PMTQE_FLAG) ) {
    // i.e. at least the BG3 trans and QE calibrations have been read
    k = (int)((1e-9*CALIB_LAMBDA-LAMBDA0)/DLAMBDA);
    if ( (k>=0) && (k<NWAVELEN_BANDS) )
      bg3qe = cal->bg3[k] * cal->pmtmaxqe * cal->pmtrelqe[k];
  }

  mjlday = (double)par->jsec/SECPDAY + (double)par->jday - MJLDOFF;
  eposec = (int)mjlday2eposec(mjlday);

  while ( iter < 2 ) {

    if ( eposec >= bank->dateTo ) {
      // advance forward
      istat = eventRead(cf->unit, cf->wbl, cf->hbl, &buf);
      needsUpdate = TRUE;
    }
    else if ( eposec >= bank->dateFrom ) {
      // we're in the right time frame now
      if ( ! needsUpdate )
	return PMTGAIN_FLAG;

      // reset the live flags
      for ( i=0; i<GEOFD_MAXMIR; i++ )
	for ( j=0; j<GEOFD_MIRTUBE; j++ )
	  cal->liveflag[i][j] = TRUE;

      int nwarn = 0;
      for ( i=0; i<bank->nTelescope; i++ ) {
	for ( j=0; j<bank->nPmt; j++ ) {
	  fdpmt_gain_data *data = &(bank->gainData[par->siteid][i][j]);
	  if ( data->badFlag ) {
	    cal->pmtgain[i][j] = 0.0;
	    cal->liveflag[i][j] = FALSE;
	    if ( nwarn < MAX_WARN )
	      perr("WARNING: badFlag!=0 : cam %02d, tube %3d.\n", i, j);
	    nwarn++;
	  }
	  else if ( (data->gainError/data->gain) > MAX_REL_GAINERROR ) {
	    cal->pmtgain[i][j] = 0.0;
	    cal->liveflag[i][j] = FALSE;
	    if ( nwarn < MAX_WARN )
	      perr("WARNING: relative gain error too high (=%.6g): "
		   "cam %02d, tube %3d.\n", data->gainError/data->gain, i, j);
	    nwarn++;
	  }
	  else
	    cal->pmtgain[i][j] = RXF_FACTOR * data->gain / bg3qe;
	}
      }

      if ( nwarn >= MAX_WARN )
	perr("%d warnings follow.\n", nwarn-MAX_WARN);

      /* save the live flag array */
      for ( i=0; i<GEOFD_MAXMIR; i++ )
	for ( j=0; j<GEOFD_MIRTUBE; j++ )
	  globalLiveFlag[par->siteid][i][j] = cal->liveflag[i][j];

      needsUpdate = FALSE;
      return PMTGAIN_FLAG;
    }
    else {
      // DST file has advanced too far, force it to re-open and read from the
      //   beginning.
      istat = -1;
      needsUpdate = TRUE;
    }

    if ( istat < 0 ) {
      // end of DST file reached.  try re-reading it from the beginning
      dstCloseUnit(cf->unit);
      dstOpenUnit(cf->unit, cf->path, MODE_READ_DST);
      iter++;
    }
  }

  // only gets here when data is not found
  perr("WARNING: PMT gains data for current time not found in attached DST file.\n");
  return 0x0;
}

int getPedestals(CalibrationFile *cf, RuntimeParameters *par, Calibration *cal) {
  int i, j, buf;
  int eposec, istat, iter = 0;
  double mjlday, mjd1, mjd2;
  fdped_dst_common *bank;

  if ( cf == NULL ) {
    perr("ERROR: pedestals file is not attached\n");
    return 0x0;
  }

  if ( ! cf->isOpen ) {
    perr("ERROR: pedestals file is not open (path=%s)\n", cf->path);
    return 0x0;
  }

  if ( (cal->mask&PMTGAIN_FLAG) == 0 ) {
    perr("ERROR: PMT gains need to be defined before calling this function.\n");
    return 0x0;
  }

  if ( par->siteid == 0 )
    bank = &brped_;
  else if ( par->siteid == 1 )
    bank = &lrped_;
  else {
    perr("invalid site ID for pedestals (=%d)\n", par->siteid);
    return 0x0;
  }
  
  mjlday = (double)par->jsec/SECPDAY + (double)par->jday - MJLDOFF;
  eposec = (int)mjlday2eposec(mjlday);

  while ( iter < 2 ) {
    // convert julian time to seconds after midnight, 1/1/1970
    mjd1 = bank->jsecfrac_start/NSECPDAY + 
      bank->jsecond_start/SECPDAY + 
      (double)bank->julian_start - MJLDOFF;

    mjd2 = bank->jsecfrac_end/NSECPDAY + 
      bank->jsecond_end/SECPDAY + 
      (double)bank->julian_end - MJLDOFF;

    int dateFrom = (int)mjlday2eposec(mjd1);
    int dateTo = (int)mjlday2eposec(mjd2);

    if ( eposec >= dateTo ) {
      // advance forward
      istat = eventRead(cf->unit, cf->wbl, cf->hbl, &buf);
    }
    else if ( eposec >= dateFrom ) {
      // we're in the right time frame now
      int k = ( eposec - dateFrom ) / 60;
      int nwarn = 0;
      for ( i=0; i<GEOFD_MAXMIR; i++ ) {
	for ( j=0; j<GEOFD_MIRTUBE; j++ ) {
	  cal->liveflag[i][j] = globalLiveFlag[par->siteid][i][j];  // retrieve saved live flag
	  if ( ! cal->liveflag[i][j] )
	    continue;

	  cal->mean[i][j] = bank->pedestal[i][j][k];
	  cal->vari[i][j] = bank->pedrms[i][j][k];
	  cal->pedestal[i][j] = bank->pedestal[i][j][k] - 
	    (int)(bank->pedrms[i][j][k]*RHO_ALPHA/cal->pmtgain[i][j]);
	  cal->pedestal[i][j] = cal->pedestal[i][j]/16 + 1;

	  if ( cal->pedestal[i][j] < 0 ) {
	    if ( nwarn < MAX_WARN )
	      perr("cam %02d, tube %3d: pedestal < 0 (=%d)\n", i, j, cal->pedestal[i][j]);
	    cal->liveflag[i][j] = FALSE;
	    nwarn++;
	  }
	  else if ( cal->vari[i][j] <= 0 ) {
	    if ( nwarn < MAX_WARN )
	      perr("cam %02d, tube %3d; vari <= 0 (=%d)\n", i, j, cal->vari[i][j]);
	    cal->liveflag[i][j] = FALSE;
	    nwarn++;
	  }
	  else if ( (double)cal->vari[i][j] > 1.75*(double)cal->mean[i][j] ) {
	    if ( nwarn < MAX_WARN )
	      perr("cam %02d, tube %3d, vari/mean > 7/4 (=%f)\n", i, j, 
		   (double)cal->vari[i][j]/(double)cal->mean[i][j]);
	    cal->liveflag[i][j] = FALSE;
	    nwarn++;
	  }
	}
      }

      if ( nwarn >= MAX_WARN )
	perr("%d warnings follow.\n", nwarn);

      return PEDESTAL_FLAG;
    }
    else {
      // DST file has advanced too far, force it to re-open and read from the
      //   beginning.
      istat = -1;
    }

    if ( istat < 0 ) {
      // end of DST file reached.  try re-reading it from the beginning
      dstCloseUnit(cf->unit);
      dstOpenUnit(cf->unit, cf->path, MODE_READ_DST);
      iter++;
    }
  }

  // only gets here when data is not found
  perr("WARNING: pedestals data for current date and time not found in attached DST file.\n");
  return 0x0;
}
