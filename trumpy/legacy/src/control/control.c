#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "constants.h"
#include "control.h"
#include "event.h"
#include "fdconstants.h"
#include "toolbox.h"

#include "histogram.h"

#include "airshower.h"
#include "showerlib.h"
#include "eventlist.h"

#define REQUIRED_OPTS 0x000008C7

double _MOLIERE_SCALE_FACTOR_ = 0.;

enum {
  PARAM_FILE, 
  PARAM_PARENT, 
  PARAM_STARTID, 
  PARAM_NTRIALS, 
  PARAM_NEVENTS, 
  PARAM_SEED, 
  PARAM_ORIGIN, 
  PARAM_EVENTLIST, 
  PARAM_ONTIME, 
  PARAM_DTIME,
  PARAM_TARGETPT,
#ifdef HISTOGRAM_MODE
  PARAM_HISTOGRAM, 
#endif

  PARAM_SPECIES, 
  PARAM_SHOWLIB,
  PARAM_ATMOSDB, 
  PARAM_SCATTER, 
  PARAM_ENERGY, 
  PARAM_NBREAK, 
  PARAM_EBREAK, 
  PARAM_ESLOPE, 
#ifdef LASER_OPTION
  PARAM_LENERGY,
  PARAM_LLAMBDA,
  PARAM_LLOCATION,
#endif
  PARAM_RP, 
  PARAM_PHI, 
  PARAM_THETA, 
  PARAM_LAT,

  PARAM_SITEID, 
  PARAM_GEOFILE, 
  PARAM_MIRREF, 
  PARAM_PARAGLAS, 
  PARAM_BG3, 
  PARAM_PMTQE, 
  PARAM_PMTGAIN,
  PARAM_PMTUNIF, 
  PARAM_PMTCAL, 
  PARAM_PSI, 
  PARAM_PHIIMP, 
  PARAM_SKYBG,  
  PARAM_MRSCALE
};

static const char *PARAM_NAME[] = {
  "FILE",       /* output file name */
  "PARENT",     /* parent directory of output file */
  "STARTID",    /* ID number of first successful event */
  "NTRIALS",    /* number of trigger attempts before quitting */
  "NEVENTS",    /* number of successful triggers before quitting */
  "SEED",       /* initial random number seed */
  "ORIGIN",     /* global origin latitude longitude and altitude (deg,deg,m) */
  "EVENTLIST",  /* list of user-input events */
  "ONTIME",     /* fdped file for on-times by part number */
  "DTIME",      /* Average dt for consecutive events (seconds) */
  "TARGETPT",   /* Target part number for run (part=run ID mod 100) */
#ifdef HISTOGRAM_MODE
  "HISTOGRAM",  /* path to output histogram file */
#endif

  "SPECIES",    /* particle species ID, using Corsika convention */
  "SHOWLIB",    /* path to shower library DST file */
  "ATMOSDB",    /* path to atmospheric parameter database file */
  "SCATTER",    /* path to molecular scattering parameter database file */
  "ENERGY",     /* primary particle energy ( log(eV) ) */
  "NBREAK",     /* number of breaks in energy spectrum */
  "EBREAK",     /* break point energies in log(eV) */
  "ESLOPE",     /* slope(s) in energy spectrum */
#ifdef LASER_OPTION
  "LENERGY",    /* Laser energy in millijoules (mJ) */
  "LLAMBDA",    /* Laser wavelength (nm) */
  "LLOCATION",  /* Laser coordinates in lat lon alt (deg, deg, m) */
#endif
  "RP",         /* impact parameter of shower track (km) */
  "PHI",        /* azm direction of shower development (deg) */
  "THETA",      /* zenith angle of shower (deg) */
  "LAT",        /* angular distance from CLF to throw events */

  "SITEID",     /* site ID number: 0=BR, 1=LR, 2=MD */
  "GEOFILE",    /* path to geometry DST file */
  "MIRREF",     /* path to fdmir_ref DST file */
  "PARAGLAS",   /* path to fdparaglas_trans DST file */
  "BG3",        /* path to fdbg3_trans DST file */
  "PMTQE",      /* path to fdpmt_qece DST file */
  "PMTGAIN",    /* path to fdpmt_gain DST file */
  "PMTUNIF",    /* path to fdpmt_uniformity DST file */
  "PMTCAL",     /* path to fdped calibration file */
  "PSI",        /* angle between shower track & ground in SD plane (deg) */
  "PHIIMP",     /* azm angle to impact point of shower track w/ ground (deg) */
  "SKYBG",      /* night sky background level in 1/(100ns) */
  "MRSCALE",    /* scale factor to Moliere radius to make data match MC */
  NULL
};

static const char *PARAM_DESC[] = {
  "[%s] output file name", 
  "[%s] parent directory of output file", 
  "[%d] ID number of first successful event", 
  "[%d] number of trigger attempts before quitting", 
  "[%d] number of successful triggers before quitting", 
  "[%d] initial random number seed", 
  "[%f %f %f] global origin latitude longitude and altitude (deg,deg,m)", 
  "[%s] list of user-input events", 
  "[%s] fdped file for on-times by part number", 
  "[%f] Average dt for consecutive events (seconds)", 
  "[%d] Target part number for run (part=runID mod 100)",
#ifdef HISTOGRAM_MODE
  "[%s] path to histogram output file", 
#endif

  "[%d] particle species ID, using Corsika convention", 
  "[%s] path to shower library DST file", 
  "[%s] path to atmospheric parameter database file", 
  "[%s] path to molecular scattering parameter database file", 
  "[%f %f] primary particle energy bounds (lo,hi in log(eV))", 
  "[%d] number of breaks in energy spectrum", 
  "[%f ... %f] break point energies in log(eV)", 
  "[%f ... %f] slope(s) in energy spectrum", 
#ifdef LASER_OPTION
  "[%f] Laser energy in millijoules (mJ)",
  "[%f] Laser wavelength (nm)", 
  "[%f %f %f] Laser location in lat(N) lon(W) alt (deg deg m)", 
#endif
  "[%f %f] bounds of impact parameter of shower track (lo,hi in km)", 
  "[%f %f] bounds of azm direction of shower development (lo,hi in deg)", 
  "[%f %f] bounds of zenith angle of shower (lo,hi in deg)", 
  "[%f] angular distance from CLF to throw events (deg)",

  "[%d] site ID number: 0=BR, 1=LR, 2=MD", 
  "[%s] path to geometry DST file", 
  "[%s] path to fdmir_ref DST file", 
  "[%s] path to fdparaglas_trans DST file", 
  "[%s] path to fdbg3_trans DST file", 
  "[%s] path to fdpmt_qece DST file", 
  "[%s] path to fdpmt_gain DST file", 
  "[%s] path to fdpmt_uniformity DST file", 
  "[%s] path to fdped calibration file", 
  "[%f %f] bounds in angle between shower track & ground in SDP (deg)", 
  "[%f %f] bounds in azm angle to impact point of shower track w/ ground (deg)", 
  "[%f] night sky background level in 1/(100ns)", 
  "[%f] scale factor to Moliere radius to make data match MC",
  NULL
};

int loadConfigurationFile(const char *cfn, RuntimeParameters *par) {
  int i, idx, narg, missingOpts;
  double lat, lon, alt, z[3];
  char *param, line[MAX_STRLEN], args[MAX_ARGS][MAX_STRLEN];
  FILE *cfp = fopen(cfn, "r");

  missingOpts = REQUIRED_OPTS;

  while ( (narg=readline(cfp, line, args)) >= 0 ) {
    if ( narg == 0 )
      continue;

    param = args[0];

    idx = 0;
    while ( PARAM_NAME[idx] != NULL && strcmp(PARAM_NAME[idx], param) != 0 )
      idx++;

    switch ( idx ) {
      /*
       *  TESTS FOR RUNTIME-SPECIFIC OPTIONS 
       */
    case PARAM_FILE:
      strcpy(par->outfn, args[1]);
      break;

    case PARAM_PARENT:
      strcpy(par->pdirn, args[1]);
      break;

    case PARAM_STARTID:
      par->trigid = atoi(args[1]);
      break;

    case PARAM_NTRIALS:
      par->ntry = atoi(args[1]);
      break;

    case PARAM_NEVENTS:
      par->nevt = atoi(args[1]);
      break;

    case PARAM_SEED:
      par->seed = atoi(args[1]);
      break;

    case PARAM_ORIGIN:
      lat = M_PI*atof(args[1])/180.0;
      lon = M_PI*atof(args[2])/180.0;
      alt = atof(args[3]);

      lla2r(lat, lon, alt, par->origin);

      z[0] = cos(lon); z[0] *= cos(lat);
      z[1] = sin(lon); z[1] *= cos(lat);
      z[2] = sin(lat);

      zmatrix(z, par->oc2ecc);
      matrixInverse(par->oc2ecc, par->oc2ecc);
      break;

    case PARAM_EVENTLIST:
      strcpy(par->infn, args[1]);
      par->flag_eventlist = TRUE;
      break;

    case PARAM_ONTIME:
      strcpy(par->ontime, args[1]);
      par->flag_ontime = TRUE;
      break;

    case PARAM_DTIME:
      par->dt = atof(args[1]);
      break;

    case PARAM_TARGETPT:
      par->tpart = atoi(args[1]);
      break;

#ifdef HISTOGRAM_MODE
    case PARAM_HISTOGRAM:
      initHistograms(args[1]);
      break;
#endif

      /*
       *  TESTS FOR SHOWER-SPECIFIC OPTIONS
       */
    case PARAM_SPECIES:
      par->species[par->nshowlib] = atoi(args[1]);
      break;

    case PARAM_SHOWLIB:
      if ( par->nshowlib >= MAX_LIBRARIES )
	perr("Too many shower library files.");
      else
	strcpy(par->showlibFile[par->nshowlib++], args[1]);
      break;

    case PARAM_ATMOSDB:
      strcpy(par->atmosdbFile, args[1]);
      break;

    case PARAM_SCATTER:
      strcpy(par->scatterFile, args[1]);
      break;

    case PARAM_ENERGY:
      par->logelo = atof(args[1]);
      par->logehi = atof(args[2]);
      break;

    case PARAM_NBREAK:
      par->nbreak = atoi(args[1]);
      if ( par->nbreak > MAX_NBREAK )
	par->nbreak = MAX_NBREAK;
      break;

    case PARAM_EBREAK:
      /* check that nbreak is already defined */
      if ( (missingOpts&(0x1<<PARAM_NBREAK)) != 0x0 ) {
	perr("NBREAK not defined before EBREAK.");
	return -1;
      }

      for ( i=0; i<par->nbreak; i++ )
	par->ebreak[i] = atof(args[i+1]);

      break;

    case PARAM_ESLOPE:
      /* check that nbreak is already defined */
      if ( (missingOpts&(0x1<<PARAM_NBREAK)) != 0x0 ) {
	perr("NBREAK not defined before ESLOPE.");
	return -1;
      }

      for ( i=0; i<(par->nbreak+1); i++ )
	par->eslope[i] = atof(args[i+1]);

      break;

#ifdef LASER_OPTION
    case PARAM_LENERGY:
      par->laser_energy = atof(args[1]) * 1.0e-3 / JOULE_PER_EV;
      break;

    case PARAM_LLAMBDA:
      par->laser_lambda = atof(args[1]) * 1.0e-9;
      break;

    case PARAM_LLOCATION:
      par->laser_latitude = atof(args[1]) * M_PI / 180.0;
      /* NOTE: longitude assumes W */
      par->laser_longitude = -atof(args[2]) * M_PI / 180.0;
      par->laser_altitude = atof(args[3]);
      break;
#endif

    case PARAM_RP:
      par->rplo = atof(args[1]);
      par->rphi = atof(args[2]);
      break;

    case PARAM_PHI:
      par->philo = atof(args[1]) * M_PI / 180.0;
      par->phihi = atof(args[2]) * M_PI / 180.0;
      break;

    case PARAM_THETA:
      par->thetalo = atof(args[1]) * M_PI / 180.0;
      par->thetahi = atof(args[2]) * M_PI / 180.0;
      break;

    case PARAM_LAT:
      par->lat = atof(args[1]) * M_PI / 180.0;
      break;

      /*
       *  TESTS FOR SD-SPECIFIC OPTIONS  (none have been defined yet)
       */

      /*
       *  TESTS FOR FD-SPECIFIC OPTIONS
       */
    case PARAM_SITEID:
      par->siteid = atoi(args[1]);
      break;

    case PARAM_GEOFILE:
      strcpy(par->geometryFile, args[1]);
      break;

    case PARAM_MIRREF:
      strcpy(par->mirrefFile, args[1]);
      par->flag_mirref = TRUE;
      break;

    case PARAM_PARAGLAS:
      strcpy(par->paraglasFile, args[1]);
      par->flag_paraglas = TRUE;
      break;

    case PARAM_BG3:
      strcpy(par->bg3transFile, args[1]);
      par->flag_bg3 = TRUE;
      break;

    case PARAM_PMTQE:
      strcpy(par->pmtQEFile, args[1]);
      par->flag_pmtQE = TRUE;
      break;

    case PARAM_PMTGAIN:
      strcpy(par->pmtGainFile, args[1]);
      par->flag_pmtGain = TRUE;
      break;

    case PARAM_PMTUNIF:
      strcpy(par->pmtUnifFile, args[1]);
      par->flag_pmtUnif = TRUE;
      break;
      
    case PARAM_PMTCAL:
      strcpy(par->pmtCalFile, args[1]);
      par->flag_pmtCal = TRUE;
      break;

    case PARAM_PSI:
      par->psilo = atof(args[1]) * M_PI / 180.0;
      par->psihi = atof(args[2]) * M_PI / 180.0;
      break;

    case PARAM_PHIIMP:
      par->phiimplo = atof(args[1]) * M_PI / 180.0;
      par->phiimphi = atof(args[2]) * M_PI / 180.0;
      break;

    case PARAM_SKYBG:
      par->nsbackground = atof(args[1]);
      break;

    case PARAM_MRSCALE:
      _MOLIERE_SCALE_FACTOR_ = atof(args[1]);
      break;

    default:
      perr("unrecognized parameter (%s)\n", param);
      ;;
    }

    missingOpts ^= ( 0x1 << idx ) & REQUIRED_OPTS;
  }

  fclose(cfp);

  return missingOpts;
}

void dumpConfTags(void) {
  int i=-1;

  while ( PARAM_NAME[++i] != NULL ) {
    printf("%10s : %s\n", PARAM_NAME[i], PARAM_DESC[i]);
  }
}

void dumpConfiguration(const RuntimeParameters *par) {
  printf("output file: %s/%s\n", par->pdirn, par->outfn);
  printf("energy bounds: %lf  %lf\n", par->logelo, par->logehi);
  printf("Rp bounds: %lf  %lf\n", par->rplo, par->rphi);
  printf("psi bounds:  %lf %lf\n", 180.0/M_PI*par->psilo, 
  	 180.0/M_PI*par->psihi);
  printf("distance from CLF (deg): %lf\n", par->lat);
}

int initializeParameters(RuntimeParameters *par) {
  /*
   *  Initialize runtime-specific parameters.
   */
  strcpy(par->outfn, "");
  strcpy(par->pdirn, "");
  par->trigid = 1;
  par->ntry = 100;
  par->nevt =  10;
  par->seed = 1;
  lla2r(CLF_LATITUDE, CLF_LONGITUDE, CLF_ALTITUDE, par->origin);
  strcpy(par->infn, "");
  strcpy(par->ontime, "");
  par->dt = 0.25;
  par->tpart = -1;

  /*
   *  Initialize shower-specific parameters.
   */
  par->species[0] = PSPEC_PROTON;
  strcpy(par->showlibFile[0], "");
  sprintf(par->atmosdbFile, "%s/%s", RTDATA, DEFAULT_ATMOSDB);
  sprintf(par->scatterFile, "%s/%s", RTDATA, DEFAULT_SCATTERDB);

  par->logelo = 17.0;
  par->logehi = 20.0;
  par->ebreak[0] = par->logehi;
  par->eslope[0] = 3.00;
#ifdef LASER_OPTION
  par->laser_energy = 2.0 * 1.0e-3 / JOULE_PER_EV;
  par->laser_lambda = 3.550011e-07;
  par->laser_latitude = CLF_LATITUDE;
  par->laser_longitude = CLF_LONGITUDE;
  par->laser_altitude = CLF_ALTITUDE;
#endif

  par->rplo =   100.00;
  par->rphi = 35000.00;
  par->philo = -M_PI;
  par->phihi = M_PI;
  par->thetalo = 0.00;
  par->thetahi = 70.0 * M_PI / 180.0;

  /*
   *  Initialize SD-specific parameters (none yet)
   */

  /*
   *  Initialize FD-specific parameters
   */
  par->siteid = -1;
  strcpy(par->geometryFile, "(none)");
  strcpy(par->mirrefFile, "(none)");
  strcpy(par->paraglasFile, "(none)");
  strcpy(par->bg3transFile, "(none)");
  strcpy(par->pmtQEFile, "(none)");
  strcpy(par->pmtGainFile, "(none)");
  strcpy(par->pmtUnifFile, "(none)");
  strcpy(par->pmtCalFile, "(none)");

  par->psilo = M_PI / 6.00;
  par->psihi = 5.0 * M_PI / 6.00;
  par->phiimplo = -M_PI;
  par->phiimphi = M_PI;
  par->lat = 1.0 * M_PI / 180.0;

  par->nsbackground = 9.00;

  /*
   *  Initialize special parameters
   */
  par->nshowlib = 0;

  /* Setting flags to OFF for optional parameters */
  par->flag_mirref    = FALSE;
  par->flag_paraglas  = FALSE;
  par->flag_bg3       = FALSE;
  par->flag_pmtQE     = FALSE;
  par->flag_pmtGain   = FALSE;
  par->flag_pmtUnif   = FALSE;
  par->flag_pmtCal    = FALSE;

  /* Setting flags to OFF for date/event formats */
  par->flag_eventlist = FALSE;
  par->flag_ontime    = FALSE;
  par->flag_today     = FALSE;

  /* Setting defaults */
  par->nsbackground = NIGHT_SKY_BG;

  int ymd, isec;
  double sec;
  time_t tt = time(&tt);

  eposec2ymdsec((double)tt, &ymd, &sec);
  par->year = ymd / 10000;
  par->month = ( ymd / 100 ) % 100;
  par->day = ymd % 100;
  isec = tt % 86400;
  par->hour = isec / 3600;
  par->min = ( isec / 60 ) % 60;
  par->sec = isec % 60;
  par->nsec = 0;

  /* get jday, jsec */
  cal2jul(par);

  par->numparts  = 1;
  par->ipart     = 1;
  par->part[0]   = 1;
  par->t0day[0]  = 0;
  par->t0sec[0]  = 0;
  par->t0nsec[0] = 0;
  par->t1day[0]  = 9999999;
  par->t1sec[0]  = (int)SECPDAY;
  par->t1nsec[0] = 999999999;

  par->oc2ecc[2][0] = cos(CLF_LONGITUDE);
  par->oc2ecc[2][0] *= cos(CLF_LATITUDE);
  par->oc2ecc[2][1] = sin(CLF_LONGITUDE);
  par->oc2ecc[2][1] *= cos(CLF_LATITUDE);
  par->oc2ecc[2][2] = sin(CLF_LATITUDE);
  zmatrix(par->oc2ecc[2], par->oc2ecc);
  matrixInverse(par->oc2ecc, par->oc2ecc);
  
  par->nbreak = 0;

  return 0;
}

/*
 * Eventually this function will check for inconsistencies between dependent 
 *   parameters.  For example, limits in zenith angle will affect psi 
 *   distribution.
 */
int setupParameters(RuntimeParameters *par) {
  int i, j, rc = 0;

  /* 
   *  check for consistency of energy parameters 
   */

  /* make sure break points fall in overall energy bounds */
  for ( i=0; i<par->nbreak; i++ ) {
    if ( par->ebreak[i] < par->logelo ) {
      fprintf(stderr, "setupParameters(): spectrum breakpoint below minimum "
	      "energy bound! (%.2lf)\n", par->logelo);
      par->nbreak -= 1;
      for ( j=i; j<par->nbreak; j++ ) {
	par->ebreak[j] = par->ebreak[j+1];
	i--;
      }
    }

    if ( par->ebreak[i] > par->logehi ) {
      fprintf(stderr, "setupParameters(): spectrum breakpoint above maximum "
	      "energy bound! (%.2lf)\n", par->logelo);
      par->nbreak -= 1;
      for ( j=i; j<par->nbreak; j++ ) {
	par->ebreak[j] = par->ebreak[j+1];
	i--;
      }
    }
  }

  /* check that breakpoint energies occur in ascending order */
  for ( i=0; i<par->nbreak-1; i++ ) {
    if ( par->ebreak[i+1] < par->ebreak[i] ) {
      fprintf(stderr, "setupParameters(): breakpoint energies not in "
	      "ascending order!\n");
      return -1;
    }
  }
  par->ebreak[par->nbreak] = par->logehi;

  return rc;
}
 
void writeHeader (char *conf, RuntimeParameters *par) {
  int i, len, pos;
  char sconf[80], sout[80], str[80];

  printf(  "\n");
  printf(  "  *********************************************\n");
  printf(  "  *                    TRUMP                  *\n");
  printf(  "  *********************************************\n");

  pos = -1;
  len = strlen(conf);
  for (i=0; i<len; i++)
    if (conf[i] == '/')
      pos = i;
  strcpy(sconf, conf+pos+1);

  pos = -1;
  len = strlen(par->outfn);
  for (i=0; i<len; i++)
    if (par->outfn[i] == '/')
      pos = i;
  strcpy(sout, par->outfn+pos+1);

  if (par->flag_today) {
    printf("  *     Running in TODAY mode                 *\n");
    printf("  *     Config file : %-24s*\n", sconf);
    printf("  *     Output DST  : %-24s*\n", sout);
    printf("  *     Throwing    : %-6d events           *\n", par->ntry);
    printf("  *     Collecting  : %-6d events           *\n", par->nevt);
  }

  else if (par->flag_ontime) {

    pos = -1;
    len = strlen(par->ontime);
    for (i=0; i<len; i++)
      if (par->ontime[i] == '/')
	pos = i;
    strcpy(str, par->ontime+pos+1);

    printf("  *    Running in ON-TIME mode                *\n");
    printf("  *    Config file  : %-24s*\n", sconf);
    printf("  *    Output DST   : %-24s*\n", sout);
    printf("  *    Average Rate : %.3e Hz            *\n", 1./par->dt);
    printf("  *    On-time file : %-24s*\n", str);
  }

  else if (par->flag_eventlist) {

    pos = -1;
    len = strlen(par->infn);
    for (i=0; i<len; i++)
      if (par->infn[i] == '/')
	pos = i;
    strcpy(str, par->infn+pos+1);

    printf("  *  Running in EVENTLIST mode                *\n");
    printf("  *  Config file    : %-24s*\n", sconf);
    printf("  *  Output DST     : %-24s*\n", sout);
    printf("  *  Throwing       : %-6d events           *\n", par->nevt);
    printf("  *  Eventlist file : %-24s*\n", str);
  }

  printf(  "  *********************************************\n\n");
}

int resolveConfiguration(RuntimeParameters *par) {
  par->ebreak[par->nbreak] = par->logehi;
  return 0;
}
