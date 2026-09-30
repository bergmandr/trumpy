#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "constants.h"
#include "control.h"
#include "event.h"
#include "fdconstants.h"
#include "random.h"
#include "toolbox.h"

#include "airshower.h"
#include "showerlib.h"
#include "eventlist.h"

/*
 * All event info must be on one line.
 * "#" comments out the line.
 * getNextEvent() Reads next line of event list file.
 * returns 1 if all NUMPARAMS values are accounted for
 * returns 0 if line begins with "#", if not all NUMPARAMS
 *   values are present or if it reaches the end of the file.
 * Any amount of whitespace on a line is allowed.
 */
int getNextEvent (int nc, FILE *infile, RuntimeParameters *par, AirShower *as,
		  GaisserHillasParameters *gh, fdraw_dst_common *f[]) {
  int i, argc;
  char argv[18][MAX_STRLEN];
  char line[MAX_STRLEN];

  while ( (argc=readline(infile, line, argv)) == 0 ) {
    if ( argc < 0 )
      return -1;
  }

  as->species = par->species[0];

  par->year      = atoi(argv[0]);
  par->month     = atoi(argv[1]);
  par->day       = atoi(argv[2]);
  par->hour      = atoi(argv[3]);
  par->min       = atoi(argv[4]);
  par->sec       = atoi(argv[5]);
  par->nsec      = atoi(argv[6]);
  as->loge       = atof(argv[7]);
  as->impactv[0] = atof(argv[8]);
  as->impactv[1] = atof(argv[9]);
  as->impactv[2] = atof(argv[10]);
  as->trackuv[0] = atof(argv[11]);
  as->trackuv[1] = atof(argv[12]);
  as->trackuv[2] = atof(argv[13]);
  gh->x0         = atof(argv[14]);
  gh->lambda     = atof(argv[15]);
  gh->xmax       = atof(argv[16]);
  gh->nmax       = atof(argv[17]);

  unitVector (as->trackuv, as->trackuv);

  if (par->year < 1900 || par->year > 2100) {
    perr("  Bad value for year\n");
    return 1;
  }
  if (par->month < 1 || par->month > 12) {
    perr("  Bad value for month\n");
    return 1;
  }
  if (par->day < 1 || par->day > 31) {
    perr("  Bad value for day\n");
    return 1;
  }
  if (par->hour < 0 || par->hour > 24) {
    perr("  Bad value for hour\n");
    return 1;
  }
  if (par->min < 0 || par->min > 60) {
    perr("  Bad value for minute\n");
    return 1;
  }
  if (par->sec < 0 || par->sec > 60) {
    perr("  Bad value for second\n");
    return 1;
  }
  if (par->nsec < 0 || par->nsec > 1000000000) {
    perr("  Bad value for nanosecond\n");
    return 1;
  }
#ifndef LASER_OPTION
  if (as->loge < 15. || as->loge > 22.) {
    perr("  Bad value for logE\n");
    return 1;
  }
#endif
  if (fabs(gh->x0) > 1000.) {
    perr("  Bad value for x0\n");
    return 1;
  }
  if (gh->xmax < 0. || gh->xmax > 2000.) {
    perr("  Bad value for Xmax\n");
    return 1;
  }
  if (gh->nmax < 1.0e6 || 
#ifdef LASER_OPTION
    gh->nmax > 1.0e17)
#else
    gh->nmax > 1.0e13)
#endif
  {
    perr("  Bad value for Nmax : %e\n", gh->nmax);
    return 1;
  }

  /* Get jday and jsec from calendar date */
  cal2jul(par);

  for (i=0; i<nc; i++) {
    f[i]->julian = par->jday;
    f[i]->jsecond = par->jsec;
    f[i]->ctdclock = par->nsec / 25;
    f[i]->gps1pps_tick = 0;
  }  

  as->zenith = acos (-as->trackuv[2]);

  return 0;
}

int getEventListNum (char *infilename) {
  int num=0, narg;
  char line[MAX_STRLEN], args[18][MAX_STRLEN];
  FILE *infile;

  infile = fopen(infilename, "r");
  if ( infile == NULL ) {
    perr("unable to open '%s' for reading.\n", infilename);
    return -1;
  }

  while ( (narg=readline(infile, line, args)) >= 0 ) {
    if ( narg == 18 )
      num++;
  }

  return num;
}

/*
 * ***ELB Added this to check that fdped file exists for ontime
 * checking.  Something simple like this could be really useful
 * for dst too, dst open doesn't check if the files open through
 * bz2 or gz libs
 */
static int file_exists(const char *filename){
    FILE *fp;
    fp = fopen(filename, "r");
    if (fp != NULL){
        fclose(fp);
        return TRUE;
  }
  return FALSE;
}

int isOnTimeEvent(RuntimeParameters *par) {
  int rc, u, m=MODE_READ_DST;
  int w, g, s=100, e;
  char str[MAX_STRLEN];

  double j0, j1;

  fdped_dst_common *p = NULL;

  if (par->siteid == BLACK_ROCK_SITEID) {
    u = PMTCAL0_FILEDES;
    p = &brped_;
    sprintf(str, "dat/calibration/fdped/black-rock/y%4dm%02dd%02d.ped.dst.gz",
	    par->year, par->month, par->day);
  }
  else if (par->siteid == LONG_RIDGE_SITEID) {
    u = PMTCAL1_FILEDES;
    p = &lrped_;
    sprintf(str, "dat/calibration/fdped/long-ridge/y%4dm%02dd%02d.ped.dst.gz",
	    par->year, par->month, par->day);
  }
  else{
    return FALSE;
  }
  // this check added by ***ELB
  if(file_exists(str) == FALSE){
    perr("File doesn't exist.\n");
    return FALSE;
  }
  rc = dstOpenUnit (u, str, m);
  if ( rc != 0 ) {
    perr("unable to open '%s' for reading (%d).\n", str, rc);
    dstCloseUnit(u);
    return FALSE;
  }

  w = newBankList(s);
  g = newBankList(s);

  while ( eventRead(u, w, g, &e) >= 0 ) {
    j0 = (double)p->julian_start +
      (double)p->jsecond_start / SECPDAY +
      (double)p->jsecfrac_start / NSECPDAY;
    j1 = (double)p->julian_end +
      (double)p->jsecond_end / SECPDAY +
      (double)p->jsecfrac_end / NSECPDAY;

    if ((par->jday+par->jsec/SECPDAY) >= j0 &&
          (par->jday+par->jsec/SECPDAY) <= j1){

      dstCloseUnit(u);
      delBankList(w);
      delBankList(g);
      return TRUE;
    }
  }
  
  delBankList(w);
  delBankList(g);
  rc = dstCloseUnit(u);

  return FALSE;
}

