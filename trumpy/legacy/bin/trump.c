#include <ctype.h>
#include <sys/types.h>
#include <sys/times.h>

#include "trump.h"
#include "histogram.h"  /* for histogramming */
// #include "iomonitor.h" /* comment this line out to deactivate I/O monitoring */

#define UT2EASTERN 18000

#ifdef QUIET_PRINT
#include "rts.h"  /* for runtime specifications file */
#endif

static void printUsage(void) {
  printf(
  "\nTRUMP - Telescope-array Reversible+Updateable Monte carlo Program.\n"
  "\n"
  "usage:         trump config1 [ config2 ... configN ]\n"
  "\n"
  "  where 'configX' is the path of a configuration file (max of %d).  One \n"
  "  configuration file should be used for each site and detector type.\n"
  "  Overall runtime parameters (random number seed, atmospheric parameter\n"
  "  database file, etc.) will be initialized according to the first valid\n"
  "  configuration file (master) and ignored in all others (slaves).\n"
  "  Configuration files are tested for consistency.  If inconsistencies\n"
  "  cannot be automatically resolved, then the configuration will be\n"
  "  ignored.  TRUMP will run if there is at least one valid configuration.\n"
  "\n"
  "  Configuration files are read line-by-line.  A line ends when either a\n"
  "  newline character ('\\n') or a '#' is reached (for commenting).  Each\n"
  "  line consists of a 'tag' (see list of tags below), followed by a list\n"
  "  of arguments, depending on what the tag does.  The arguments may be\n"
  "  strings (%%s), integers (%%d), or floating point numbers (%%f).\n"
  "  The table below describes the tag, types of arguments, and a brief\n"
  "  description of the tag's function.\n"
  "\n"
  "configuration file tags:\n"
  "\n", 
  NCONFIG_MAX);

  dumpConfTags();
  printf("\n");
}

boolean updateOutfilePath(RuntimeParameters *par) {
  char *p, *q, sfx[MAX_STRLEN], newpath[MAX_STRLEN];
  strcpy(newpath, par->outfn);

  p = strstr(newpath, ".");
  strcpy(sfx, p+1);

  // first, look for pattern p??, where ? is a digit
  if ( (q=strstr(newpath, "p")) != NULL ) {
    if ( isdigit(q[1]) && isdigit(q[2]) )
      sprintf(q, "p%02d.%s", par->part[par->ipart], sfx);
  }
  else {
    sprintf(p, "p%02d.%s", par->part[par->ipart], sfx);
  }

  if ( strcmp(newpath, par->outfn) != 0 ) {
    strcpy(par->outfn, newpath);
    return TRUE;
  }
  else
    return FALSE;
}

#ifdef IDENTIFY
#include <sys/utsname.h>
static void dumpIdentity(void) {
  struct utsname name;
  uname(&name);
#if _UTSNAME_DOMAIN_LENGTH - 0
  printf("%s %s %s %s %s %s\n", 
       name.sysname, name.nodename, 
#  ifdef __USE_GNU
       name.domainname, 
#  else
       name.__domainname, 
#  endif
       name.release, name.version, name.machine
       );
#else
  printf("%s %s %s %s %s\n", 
       name.sysname, name.nodename, name.release, name.version, name.machine
       );
#endif
  fflush(stdout);
  return;
}
#endif

#ifdef LASER_OPTION
static void createLaser(RuntimeParameters *par, AirShower *as) {
  double rl[3], rclf[3], rotm[3][3];

  as->zenith     = M_PI;
  as->loge       = log10(par->laser_energy);
  as->trackuv[0] = 0.;
  as->trackuv[1] = 0.;
  as->trackuv[2] = 1.;

  lla2r(CLF_LATITUDE, CLF_LONGITUDE, CLF_ALTITUDE, rclf);

  rotm[2][0] = cos(CLF_LONGITUDE);
  rotm[2][0] *= cos(CLF_LATITUDE);
  rotm[2][1] = sin(CLF_LONGITUDE);
  rotm[2][1] *= cos(CLF_LATITUDE);
  rotm[2][2] = sin(CLF_LATITUDE);
  zmatrix(rotm[2], rotm);
  //  matrixTranspose(rotm, rotm);

  lla2r(par->laser_latitude, par->laser_longitude, par->laser_altitude, rl);
  subvec(rl, rclf, rl);

  applyRotation(rotm, rl, as->impactv);
  pout("Laser at (%lf, %lf, %lf)\n", as->impactv[0], as->impactv[1], as->impactv[2]);
}
#endif

/*
 *  MAIN routine for TRUMP
 *
 *  SS --  1.07.09
 *  DRB    1.15.09 Added getAcceptedTrack, makePEList & sortPEList calls
 */
int main(int argc, char *argv[] ) {

  fprintf(stdout, "This version built using DEDX_MODEL %d\n", DEDX_MODEL);

#ifdef IDENTIFY
  dumpIdentity();
#endif

  int i, nconfig=0, n[NCONFIG_MAX], v[NCONFIG_MAX], rc;
  int verp[NCONFIG_MAX];
  int unit[NCONFIG_MAX];
  int banklist[NCONFIG_MAX], size=10;
  int timerc=0;
  //  double second;
  boolean triggered;
  int numGoodConfigs;

  double rp;
  double rpv[3]; 
  double psi;

  int t0, t1, t_beginEl, t_endEl;
  struct tms tbuff;
  double total_time;
  int total_hour, total_min, total_sec, total_dec;

  RuntimeParameters *par[NCONFIG_MAX], runtime;
  UCalibration *calib[NCONFIG_MAX];
  TACalibration *tacalib = NULL;
  TLCalibration *tlcalib = NULL;
  
  FDSiteGeometry *fdsg[NCONFIG_MAX];
  TGeom *tgeom[NCONFIG_MAX];

  AirShower *as = newInstanceOf(AirShower);
  GaisserHillasParameters *gh = newInstanceOf(GaisserHillasParameters);
  Track *track = newInstanceOf(Track);
  ObservedTrack *otrack = newInstanceOf(ObservedTrack);
  // AcceptedTrack *atrack = newInstanceOf(AcceptedTrack);
  // PEList *pelist = newInstanceOf(PEList);
  PETimes *petimes = newInstanceOf(PETimes);

  brraw_dst_common *brraw = &brraw_;
  lrraw_dst_common *lrraw = &lrraw_;
  fdraw_dst_common *fdraw = &fdraw_;
  
  fraw1_dst_common *fraw1 = &fraw1_;
//   ftrg1_dst_common *ftrg1 = &ftrg1_;
  
  
  fdraw_dst_common *fdary[NCONFIG_MAX];
#ifdef _OPENMP
  int rs;
#endif
  FILE *eventlistfile = NULL;
//   int j,k,l;
#ifdef HISTOGRAM_MODE
  float dv[NCONFIG_MAX];
  float xtup[NCONFIG_MAX][NTUPLE_NENT];
#endif

  if ( argc == 1 ) {
    printUsage();
    return -1;
  }

#ifndef QUIET_PRINT
  printf("\n"
  "*************************************************************************\n"
  "*                                                                       *\n"
  "*           ********  ******    **    **  **    **  ******              *\n"
  "*              **     **    **  **    **  ***  ***  **    **            *\n"
  "*              **     **    **  **    **  ********  **    **            *\n"
  "*              **     ******    **    **  ** ** **  ******              *\n"
  "*              **     **  **    **    **  **    **  **                  *\n"
  "*              **     **   **   **    **  **    **  **                  *\n"
  "*              **     **    **   ******   **    **  **                  *\n"
  "*                                                                       *\n"
  "*       Telescope-array Reversible+Updateable Monte carlo Program       *\n"
  "*                                                                       *\n"
  "*************************************************************************\n"
  "\n");
#endif

  t0 = times(&tbuff);

  /***************************************************************************
   *
   *  INITIALIZATION 
   *
   ***************************************************************************/

  pout("  * INITIALIZING\n");

  /*
   *  Loop over the command line arguments and attach configuration files.  Let 
   *    the first valid configuration file be the "master", meaning the runtime 
   *    and shower specific parameters will be set by that file, regardless of 
   *    what the others have set.  Inconsistent configuration files will be 
   *    disregarded.
   */
  pout("      - Scanning configuration files.\n");
  for ( i=1; i<argc; i++ ) {
    pout("      - Reading configuration file '%s'... ", argv[i]);

    par[nconfig] = newInstanceOf(RuntimeParameters);
    initializeParameters(par[nconfig]);

    loadConfigurationFile(argv[i], par[nconfig]);
    rc = resolveConfiguration(par[nconfig]);
    if ( rc != 0 ) {
      perr("inconsistencies found (rc=%d)... configuration discarded.\n", rc);
      free(par[nconfig]);
    }
    else {
      pout("done. (set to %s)\n", (nconfig==0)?"master":"slave");
      nconfig++;
    }

#ifdef HISTOGRAM_MODE
    /* histogram file is loaded from 'control.h' */
    addConfigHists(nconfig);
    dv[nconfig-1] = 0.0f;
#endif

    if ( nconfig >= NCONFIG_MAX )
      break;
  }

#ifdef QUIET_PRINT
  char rtspath[256];
  sprintf(rtspath, "%s/%s", par[0]->pdirn, par[0]->outfn);
  char *p = strstr(rtspath, ".dst.gz");
  sprintf(p, ".rts");
  FILE *rtsfp = fopen(rtspath, "w");
  RuntimeSpecsBitField rtsbf[NCONFIG_MAX];
#endif

  /*
   *  Initialize the site/detector independent parameters.
   */
  pout("      - Initializing site/detector independent parameters...\n");
  rc = initTRUMP(nconfig, par);
//   abort();
  if ( rc != 0 ) {
    perr("  fatal error... aborting.\n");
    return rc;
  }
  pout("        ...done.\n");

  /* 
   *  Check eventlist, if it exists.
   */
  if ( par[0]->flag_eventlist ) {
    /* Checking validity / usefulness of event list */
    if ( par[0]->nevt < 1) {
      perr("No good events found in %s !\n", par[0]->infn);
      return 0;
    }

    eventlistfile = fopen(par[0]->infn, "r");
  }

  /*
   *  Initialize the geometry & trigger routines for each detector type & site.
   */
  for ( i=0; i<nconfig; i++ ) {
    if ( par[0]->flag_eventlist ) {
      par[i]->nevt = getEventListNum(par[0]->infn);
      par[i]->ntry = par[i]->nevt + 1;
    }

    pout("      - Initializing geometry + trigger routines for site %d...\n", par[i]->siteid);

    banklist[i] = newBankList(size);
    addBankList(banklist[i], TRUMPMC_BANKID);

    calib[i] = newInstanceOf(UCalibration);
    fdsg[i] = newInstanceOf(FDSiteGeometry);
    tgeom[i] = newInstanceOf(TGeom);

    switch ( par[i]->siteid ) {
    case BLACK_ROCK_SITEID:
      unit[i] = FDDSTOUT0_FILEDES;
      
      calib[i]->tc = newInstanceOf(TACalibration);
      tacalib = calib[i]->tc;
      fdary[i] = 
      initFDSite(par[i], calib[i], fdsg[i], tgeom[i], banklist[i], brraw);
      break;

    case LONG_RIDGE_SITEID:
      unit[i] = FDDSTOUT1_FILEDES;
      
      calib[i]->tc = newInstanceOf(TACalibration);
      tacalib = calib[i]->tc;
      fdary[i] = 
      initFDSite(par[i], calib[i], fdsg[i], tgeom[i], banklist[i], lrraw);
      break;

    case MIDDLE_DRUM_SITEID:
      unit[i] = FDDSTOUT2_FILEDES;
      
      break;
      
    case TALE_SITEID:
      unit[i] = FDDSTOUT3_FILEDES;
      calib[i]->tl = newInstanceOf(TLCalibration);
      tlcalib = calib[i]->tl;
      fdary[i] = 
      initFDSite(par[i], calib[i], fdsg[i], tgeom[i], banklist[i], fdraw);
      break;
      
      
    default:
      unit[i] = FDDSTOUT4_FILEDES;
      
      calib[i]->tc = newInstanceOf(TACalibration);
      tacalib = calib[i]->tc;
      fdary[i] = 
      initFDSite(par[i], calib[i], fdsg[i], tgeom[i], banklist[i], fdraw);
      
    }
    if ( rc != 0 ) {
      perr("fatal error setting up FD Site... aborting.\n");
      return rc;
    }
    pout("        ...done.\n");

    n[i] = 0;
    v[i] = 0;
    verp[i] = 0;
  }

  /* 
   *  Ready to go!
   */

  fprintf(stdout, "* SIMULATION IS NOW RUNNING *\n");
  fflush(stdout);

#ifdef _OPENMP
  int ntmp = omp_get_max_threads();
  unsigned long long ranq2V[ntmp], ranq2W[ntmp];
  for ( rs=0; rs<omp_get_max_threads(); rs++ )
    ranq2getState(rs, &ranq2V[rs], &ranq2W[rs]);
#else
  unsigned long long ranq2V, ranq2W;
  ranq2getState(&ranq2V, &ranq2W);
#endif  
  
  /***************************************************************************
   * 
   *  MAJOR LOOP OVER SHOWERS 
   *
   ***************************************************************************/
#ifdef CHECK_TRIAL
  n[0] = CHECK_TRIAL - 1;
#endif
  boolean go;
  do { /* loop until satisfied */
#ifdef TRACK_MEMORY
    dumpLiveTrackers();
#endif

#ifdef CHECK_TRIAL
//     if ( v[0] == 1 )
//       break;
#endif

    go = TRUE;

    if ( n[0] > 0 ) {
#ifdef _OPENMP

#ifndef CHECK_TRIAL

      for ( rs=0; rs<omp_get_max_threads(); rs++ ) {
        ranq2setState(rs, ranq2V[rs], ranq2W[rs]);
        ranq2doub(rs);
        ranq2getState(rs, &ranq2V[rs], &ranq2W[rs]);
        ranq2setup(rs, ranq2V[rs]);
      }
#else
    int ii;
    for (rs=0; rs<omp_get_max_threads(); rs++) {
      ranq2setup(rs, par[0]->seed+rs);
    
      
      for (ii=0; ii<n[0]; ii++)
        ranq2doub(rs);
      ranq2getState(rs, &ranq2V[rs], &ranq2W[rs]);
      ranq2setup(rs, ranq2V[rs]);
    }
#endif
#else
#ifndef CHECK_TRIAL
      ranq2setState(ranq2V, ranq2W);
      ranq2doub();
      ranq2getState(&ranq2V, &ranq2W);
      ranq2setup(ranq2V);
#else
    ranq2setup(par[0]->seed);
    int ii;
    for (ii=0; ii<n[0]; ii++)
      ranq2doub();
    ranq2getState(&ranq2V, &ranq2W);
    ranq2setup(ranq2V);
#endif
#endif
    }
    time_t now = time(&now);
    eposec2jul((int)now-UT2EASTERN, &runtime);
    jul2cal(&runtime);

    t1 = times(&tbuff);
    total_time = (double)(t1 - t0) / (double)(CLOCKS_PER_SEC/10000);
    total_hour = (int)total_time / 3600;
    total_min  = ((int)total_time % 3600) / 60;
    total_sec  = ((int)total_time % 3600) % 60;
    total_dec  = (int)(100.0*(total_time - floor(total_time)));

#ifndef QUIET_PRINT
    printf("\n- TRIAL %6d : %4d/%02d/%02d %02d:%02d:%02d "
         "[ %02d:%02d:%02d.%02d ]\n",
         n[0]+1, runtime.year, runtime.month, runtime.day,
         runtime.hour, runtime.min, runtime.sec,
         total_hour, total_min, total_sec, total_dec);
    printf("  v/n so far (iconfig:V/V_ERp/N): ");
    for ( i=0; i<nconfig; i++ )
      printf("  %1d:%d/%d/%d", i, v[i], verp[i], n[i]);
    printf("\n");
#else
    for (i=0; i<nconfig; i++ ) {
      rtsbf[i].itry = n[0] + 1;
      rtsbf[i].triglev = 0;
      rtsbf[i].rtday = runtime.day;
      rtsbf[i].rtmonth = runtime.month;
      rtsbf[i].rtyear = runtime.year;
      rtsbf[i].rtsec100 = 100*runtime.sec;
      rtsbf[i].rtminute = runtime.min;
      rtsbf[i].rthour = runtime.hour;
    }
#endif

    /* DRB: Print out random number generator state */
#ifndef QUIET_PRINT
    printf("  ranq2 state: ");
#ifdef _OPENMP
    for (rs=0; rs<omp_get_num_threads(); rs++)
      ranq2printState(rs);
#else
    ranq2printState();
#endif

    printf("  jday before getEventTime: %f\n",
         (double)par[0]->jday +
         (double)par[0]->jsec / SECPDAY +
         (double)par[0]->nsec / NSECPDAY);
#endif

    /*
     *  GENERATE AN EAS 
     *  Read from list or create random event
     */
    //    printf("  creating air shower... ");
    //    fflush(stdout);
    if ( par[0]->flag_eventlist ) {
      getNextEvent(nconfig, eventlistfile, par[0], as, gh, fdary);
    }
    else {
#ifdef LASER_OPTION
      createLaser(par[0], as);
      gh->x0     = 0.;
      gh->lambda = 0.;
      gh->xmax   = 0.;
      gh->nmax   = (par[0]->laser_energy)*par[0]->laser_lambda/HC;
#else
      createAirShower(par[0], as);
      getGaisserHillasParameters(as, gh);
#endif

      // timerc = getEventTime(par[0], fdraw);
      if ( par[0]->tpart < 0 ) {
#ifdef LASER_OPTION        
        timerc = nextEventTime(nconfig, par, fdary);
#else        
        if ( nconfig > 1 )
          timerc = nextTandemEventTime(nconfig, par, fdary);
        else
          timerc = nextEventTime(nconfig, par, fdary);
#endif      
      }
      else {
        do {
          timerc = nextEventTime(nconfig, par, fdary);
          if ( par[0]->part[par[0]->ipart] > par[0]->tpart )
            timerc = -1;

          if ( timerc != 0 )
            break;

        } while ( par[0]->tpart != par[0]->part[par[0]->ipart] );
      }
    }

    //    printf("done.\n");

#ifndef QUIET_PRINT
    printf("  date %4d/%02d/%02d %02d:%02d:%02d.%09d\n",
         par[0]->year, par[0]->month, par[0]->day, 
         par[0]->hour, par[0]->min, par[0]->sec, par[0]->nsec);
#ifdef LASER_OPTION
    printf("  primary: LASER, %5.2f mJ, %5.1f nm\n",
           par[0]->laser_energy * 1.0e3 * JOULE_PER_EV,
           par[0]->laser_lambda * 1.0e9);
#else
    printf("  primary: %s\n",
         (as->species==PSPEC_IRON)?"IRON":"PROTON");
#endif         
    printf("  {E,Nmax,X0,Xmax,Lam} %8.2e %8.2e %5.1f %5.1f %5.1f\n",
         pow(10.,as->loge),gh->nmax,gh->x0,gh->xmax,gh->lambda);
    fflush(stdout);
#else
    for (i=0; i<nconfig; i++) {
      rtsbf[i].pspec = (as->species==PSPEC_IRON)?2:1;

      rtsbf[i].stday = par[0]->day;
      rtsbf[i].stmonth = par[0]->month;
      rtsbf[i].styear = par[0]->year;
      rtsbf[i].stsec100 = 100*par[0]->sec;
      rtsbf[i].stminute = par[0]->min;
      rtsbf[i].sthour = par[0]->hour;
      rtsbf[i].stnano = par[0]->nsec;

      rtsbf[i].loge1E4 = (int)( 1.0e+4 * as->loge );
      double xtmp = as->impactv[0]-as->impactv[2]*as->trackuv[0]/as->trackuv[2];
      double ytmp = as->impactv[1]-as->impactv[2]*as->trackuv[1]/as->trackuv[2];
      rtsbf[i].xcore = (int)(100.*xtmp);
      rtsbf[i].ycore = (int)(100.*ytmp);

      rtsbf[i].logn1E4 = (int)( 1.0e+4 * log10(gh->nmax) );
      rtsbf[i].xmax1E4 = (int)( 1.0e+4 * gh->xmax );
      rtsbf[i].lambda10 = (int)(10.*gh->lambda);
      rtsbf[i].x010 = (int)(10.*gh->x0);
    }
#endif

    numGoodConfigs = 0;
    for ( i=0; i<nconfig; i++ ) {
      if ( n[i] == par[i]->ntry || v[i] == par[i]->nevt )
      continue;

      /* Apply E vs Rp Cut */
      rp  = getRpVector(fdsg[i], as, rpv);
      psi = getPsiAngle(fdsg[i], as);
      int goodEvsRp = applyEvsRpCut(as->loge ,rp);

#ifndef QUIET_PRINT
      printf("  {Conf,Rp,Psi,Zen,Phi} %d %5.1f %5.1f %5.1f %5.1f\n",
           i, rp/1000., psi*R2D, as->zenith*R2D,
           atan2(as->trackuv[2],as->trackuv[1])*R2D);
  printf("EVLISTITEM %04d %02d %02d %02d %02d %02d %d %f %f %f %f %f %f %f %f %f %f %e\n",
         par[i]->year,par[i]->month,par[i]->day,par[i]->hour,par[i]->min,par[i]->sec,
         par[i]->nsec,as->loge,as->impactv[0],as->impactv[1],as->impactv[2],
         as->trackuv[0],as->trackuv[1],as->trackuv[2],
         gh->x0,gh->lambda,gh->xmax,gh->nmax);
#else
      rtsbf[i].configID = i;
      rtsbf[i].rpcm = (int)( 100.0 * rp );
      rtsbf[i].psirad1E4 = (int)( 1.0e+4 * psi );
      rtsbf[i].zenrad1E4 = (int)( 1.0e+4 * as->zenith );
      rtsbf[i].azmrad1E4 = (int)( 1.0e+4 * atan2(as->trackuv[1], as->trackuv[0]) );
#endif

      n[i]++;

      if ( goodEvsRp )
      if ( (nconfig==1) || (timerc>0) )
        numGoodConfigs++;

#ifdef HISTOGRAM_MODE
      xtup[i][0] = (float)as->loge;
      xtup[i][1] = (float)as->species;
      xtup[i][2] = (float)gh->xmax;
      xtup[i][3] = (float)gh->nmax;
      xtup[i][4] = (float)( as->zenith * 180.0/M_PI );
      xtup[i][5] = (float)as->impactv[0];
      xtup[i][6] = (float)as->impactv[1];
      xtup[i][7] = (float)rp;
      xtup[i][8] = (float)( psi * 180.0/M_PI );
      xtup[i][9] = 0.0f;

      dv[i] += 1.00f;
#endif
    }

    if ( numGoodConfigs == 0 ) {
#ifndef QUIET_PRINT
      if ( (nconfig>1) && (timerc==0) )
      pout("  FAILED TIME CUT.\n");
      else
      printf("    FAILED E vs. Rp cut\n");
#else
      fwrite(&rtsbf[0], sizeof(RuntimeSpecsBitField), 1, rtsfp);
#endif

      fflush(stdout);
      if ( timerc < 0 ) { /* no more parts to read */
        go = FALSE;
      }

      for ( i=0; i<nconfig; i++ ) {
        if ( n[i] >= par[i]->ntry || v[i] >= par[i]->nevt )
          go = FALSE;                                     

#ifdef HISTOGRAM_MODE
        hfn(10*(i+1), xtup[i]);
      if ( par[i]->ntry > 1000000 ) {
        if ( (n[i]+1)%1000000 == 0 ) {
          perr("exceeded 1M events, opening new histogram file.\n");
          openNewHistogramFile();
        }
      }  
#endif
      }      

      continue;              
    }               

    //    printf("          building EAS track... ");
    //    fflush(stdout);
#ifdef LASER_OPTION
    buildLaserTrack(par[0]->laser_lambda, par[0], as, gh, track);
#else
    buildTrack(par[0], as, gh, track);
#endif
    //    printf("done.\n");
//   printf("track->nseg = %d\n",track->nseg);
    /*
     *  LOOP OVER DETECTOR CONFIGURATIONS
     */
    //    printf("          looping over detector configurations:\n");
    for ( i=0; i<nconfig; i++ ) {
#ifndef QUIET_PRINT
      printf("  Configuration %d:\n", i);
#endif
      triggered = FALSE;

      if ( n[i] > par[i]->ntry || v[i] > par[i]->nevt )
      continue;

      /* Apply E vs Rp Cut */
      rp = getRpVector(fdsg[i], as, rpv);
      //      double psi = getPsiAngle(fdsg[i], as);
      int goodEvsRp = applyEvsRpCut(as->loge ,rp);

      //      printf("            {E,Rp,Psi,Zen,Phi} %8.2e %5.1f %5.1f %5.1f %5.1f\n", pow(10.,as->loge), rp/1000., psi*R2D, as->zenith*R2D, atan2(as->trackuv[2],as->trackuv[1])*R2D);

      if ( ! goodEvsRp ) {
#ifndef QUIET_PRINT
	printf("    FAILED E vs. Rp Cut\n");
#endif
	continue;
      }
#ifndef QUIET_PRINT
      printf("    PASSED E vs. Rp Cut\n");
#endif

      if ( nconfig > 1 ) {
	if ( ((0x1<<i)&timerc) == 0x0 ) {
#ifndef QUIET_PRINT
	  printf("    FAILED time window cut.\n");
#endif
	  continue;
	}
#ifndef QUIET_PRINT
	printf("    PASSED time window cut.\n");
#endif
      }

      verp[i]++;
#ifdef HISTOGRAM_MODE
      xtup[i][9] = 1.0f;
#endif

      /** uncomment this if you want to revert to traditional 'dtime' mode *
      if ( par[0]->flag_ontime ) {
      par[i]->ntry++;
      par[i]->nevt++;
      }
      **/

      /* copy the event times across the configurations to keep the 
       calibration routines from failing. */
      if ( i > 0 ) {
        par[i]->jday = par[0]->jday;
        par[i]->jsec = par[0]->jsec;
        par[i]->year = par[0]->year;
        par[i]->month = par[0]->month;
        par[i]->day = par[0]->day;
        par[i]->hour = par[0]->hour;
        par[i]->sec = par[0]->sec;
        par[i]->nsec = par[0]->nsec;
      }

      /* Load time-dependent DST banks */
      switch (par[i]->siteid) {
        case BLACK_ROCK_SITEID:
        case LONG_RIDGE_SITEID:
//           printf("Entering getTATimeDep...\n");
          getTATimeDependentCalibration(par[i], tacalib);
//           printf("completed!\n");
          getTACalibrationReduction(tacalib);
          break;
          
        case MIDDLE_DRUM_SITEID:
          
        case TALE_SITEID:
          break;
      }
      
      pout("    computing observed track... \n");

#ifdef LASER_OPTION
      getLaserObservedTrack(par[i]->laser_lambda, par[i], as, track, fdsg[i], otrack);
#else
      getObservedTrack(par[i], as, track, fdsg[i], otrack);
#endif
      // pout("done.\n");
//     printf("otrack->nseg = %d\n",otrack->nseg);
      /* 
       *  DO RAYTRACING 
       */
      pout("    doing ray tracing... \n");
      fflush(stdout);
#ifdef LASER_OPTION
      rc = getUpwardAcceptedTrack(as, track, fdsg[i], otrack, tacalib, 
                          tgeom[i], petimes);
#else
      rc = getAcceptedTrack(as, track, fdsg[i], otrack, tacalib, tgeom[i], 
             petimes);
#endif
//     printf("petimes->nmir: %d\n",petimes->nmir );
//     printf("t1, t2: %f %f\n",petimes->t1,petimes->t2);
//     for (j=0; j<petimes->nmir; j++) {
//       printf("mirror at index j=%2d: mirror %2d\n",j,petimes->mir[j]);
//       for (k=0; k<GEOFD_MIRTUBE; k++)
//         if (petimes->len[j][k]) {
// //           printf("  cam %2d pmt %3d: len = %d\n",j,k,petimes->len[j][k]);
//           for (l=0; l<petimes->len[j][k]; l++)
//             printf("DUMPPET %2d %3d %f %f\n",petimes->mir[j],k,
//                    petimes->t[j][k][l],petimes->n[j][k][l]);
//         }
//     }
//     abort();
      if ( rc != 0 ) {
        perr("trial %6d: config %d: getAcceptedTrack() returned icode %d.\n", 
            n[0], i, rc);
        clearObservedTrack(otrack);
        clearPETimes(petimes);
        continue;
      }

      // pout("done.\n");

      /*
       *  SIMULATE ELECTRONICS RESPONSE 
       */
      pout("    doing electronics simulation... \n");//\n");
      t_beginEl = times(&tbuff);
      switch (par[i]->siteid) {
        case BLACK_ROCK_SITEID:
        case LONG_RIDGE_SITEID:
      
          fdraw->julian = fdary[i]->julian;
          fdraw->jsecond = fdary[i]->jsecond;
          fdraw->ctdclock = fdary[i]->ctdclock;
          fdraw->gps1pps_tick = fdary[i]->gps1pps_tick;

          triggered = simTATriggerResponse(tacalib, petimes, fdraw);
          break;
        case MIDDLE_DRUM_SITEID:
          break;
        case TALE_SITEID:
          fraw1->julian = fdary[i]->julian;
          fraw1->jsecond = fdary[i]->jsecond;
          fraw1->jclkcnt = fdary[i]->ctdclock * 25;
          
          triggered = simTLTriggerResponse(petimes, fraw1, tlcalib);
          printf("triggered = %d\n",triggered);
          break;
      }
      // pout("    done.\n");
      t_endEl = times(&tbuff);
#ifndef QUIET_PRINT      
      printf("Electronics elapsed time: %F\n",(double)(t_endEl - t_beginEl)/((double)(CLOCKS_PER_SEC/10000)));
#endif
      if ( triggered ) {
      /*
       *  WRITE EVENT TO DST FILES
       */
/*
      switch (par[i]->siteid) {
        case BLACK_ROCK_SITEID:
          for (j=0; j<512; j++)
            printf("BRWAVE 7 113 %d\n",fdraw->m_fadc[1][113][j]);
          for (j=0; j<512; j++)
            printf("BRWAVE 6 113 %d\n",fdraw->m_fadc[0][113][j]);
          break;
        case TALE_SITEID:
          for (j=0; j<100; j++)
            printf("TLWAVE 6 248 %d\n",fraw1->m_fadc[0][31][j]);
          for (j=0; j<100; j++)
            printf("TLWAVE 7 248 %d\n",fraw1->m_fadc[1][16][j]);  
          break;
      }*/
  
#ifdef HISTOGRAM_MODE
      xtup[i][9] = 2.0f;
      hf1(TRIAL_DIST_BASE+i+1, dv[i], 1.0);
      dv[i] = 0.0f;
#endif
//   printf("about to buildTrumpMCBank...\n");
      clearTrumpMCBank(&trumpmc_);
      buildTrumpMCBank(par[i]->siteid, as, fdsg[i], gh, track, otrack, petimes, tacalib, par[i], &trumpmc_);
//   printf("done.\n");
  switch (par[i]->siteid) {
    case BLACK_ROCK_SITEID:
    case LONG_RIDGE_SITEID:
      compactFDRawBank(fdraw);
      fdraw->part = fdary[i]->part;
      fdraw->event_num = v[i] + par[i]->trigid;

      while ( (fdraw->ctdclock - fdraw->gps1pps_tick) >= 40000000 ) {
        fdraw->ctdclock -= 40000000;
        fdraw->jsecond++;
      }
      while ( fdraw->jsecond >= 86400 ) {
        fdraw->jsecond -= 86400;
        fdraw->julian++;
      }
      break;
    case MIDDLE_DRUM_SITEID:
      break;
    case TALE_SITEID:
      fraw1->event_num = v[i] + par[i]->trigid;
      break;
  }

#ifndef QUIET_PRINT
      printf("    SUCCESS writing event %d : site %d part %d event %d\n", v[i], i, fdraw->part, fdraw->event_num);
#else
      rtsbf[i].triglev = 2;
      fwrite(&rtsbf[i], sizeof(RuntimeSpecsBitField), 1, rtsfp);
#endif
      // this will only apply to tandem mode
      if ( nconfig > 1 ) {
        if ( updateOutfilePath(par[i]) ) {
          char newpath[MAX_STRLEN];
          dstCloseUnit(unit[i]);
          sprintf(newpath, "%s/%s", par[i]->pdirn, par[i]->outfn);
          dstOpenUnit(unit[i], newpath, MODE_WRITE_DST);
        }
      }
#ifndef QUIET_PRINT
      pout("  OUTFILE = %s\n", par[i]->outfn);
#endif

//    if ( par[i]->siteid == BLACK_ROCK_SITEID ) {
//      *brraw = *fdraw;
//    }
//    else if ( par[i]->siteid == LONG_RIDGE_SITEID ) {
//      *lrraw = *fdraw;
//    }
      switch (par[i]->siteid) {
      case BLACK_ROCK_SITEID:
	*brraw = *fdraw;
	break;
      case LONG_RIDGE_SITEID:
	*lrraw = *fdraw;
	break;
      }
      eventWrite(unit[i], banklist[i], TRUE);
      v[i]++;
      }
#ifdef QUIET_PRINT
      else {
	rtsbf[i].triglev = 1;
	fwrite(&rtsbf[i], sizeof(RuntimeSpecsBitField), 1, rtsfp);
      }
#endif

      /*
       *  CLEANUP AND PREPARE FOR NEXT EAS
       */
      pout("    collecting garbage... \n");

      clearObservedTrack(otrack);
      // clearAcceptedTrack(atrack);
      // clearPEList(pelist);
      clearPETimes(petimes);
      pout("done.\n");
    }

    for ( i=0; i<nconfig; i++ ) {
      if ( timerc < 0 ) { /* no more parts to read */
      go = FALSE;
      }

      if ( par[0]->flag_eventlist && n[i] == (par[i]->ntry-1) )
      go = FALSE;

      if ( !par[0]->flag_eventlist && 
         (n[i] == par[i]->ntry || v[i] == par[i]->nevt) )
      go = FALSE;

#ifdef HISTOGRAM_MODE
      hfn(10*(i+1), xtup[i]);
      if ( par[0]->ntry > 1000000 ) {
      if ( (n[i]+1)%1000000 == 0 ) {
        perr("exceeded 1M trials.  opening new histogram file.\n");
        openNewHistogramFile();
      }
      }
#endif
    }

    clearTrack(track);
    //    printf("          jday after getEventTime:  %f\n\n",par[0]->jday);

    fflush(stdout);
  } while ( go );

  /***************************************************************************
   *
   *  FINAL CLEANUP
   *
   ***************************************************************************/
  for ( i=0; i<nconfig; i++ )
    dstCloseUnit(unit[i]);

  closeAtmosDB();
  closeScatterDB();

  if (par[0]->flag_eventlist)
    fclose(eventlistfile);

  detachTACalibration();  // closes all of the calibration file channels

  freeFADCBuffer();   // free up the memory mapping for the FADC buffer
#ifdef HISTOGRAM_MODE
  endHistograms();
#endif

#ifdef QUIET_PRINT
  fclose(rtsfp);
#endif

  t1 = times(&tbuff);
  total_time = (double)(t1 - t0) / (double)(CLOCKS_PER_SEC/10000);
  total_hour = (int)total_time / 3600;
  total_min  = ((int)total_time % 3600) / 60;
  total_sec  = ((int)total_time % 3600) % 60;
  total_dec  = (int)(100.0*(total_time - floor(total_time)));

  printf("\n***** END CRITERIA MET : Trump will now exit "
       "[ %02d:%02d:%02d.%02d ]\n\n",
       total_hour, total_min, total_sec, total_dec);

  printf("  RUN SUMMARY:\n");
  printf("  -number EAS thrown: %d\n", n[0]);
  printf("  -number EAS passing E/Rp cut: ");
  for ( i=0; i<nconfig; i++ )
    printf(" %d:%d", i, verp[i]);
  printf("\n");

  printf("  -number EAS triggering FD site: ");
  for ( i=0; i<nconfig; i++ )
    printf(" %d:%d", i, v[i]);
  printf("\n\n");

  //  iomon_dumpStreamStatus(stdout);

  return 0;
}

int initTRUMP(int nc, RuntimeParameters *par[]) {
  int rc = 0;
  
  /* 
   *  This should be handled in the configuration-loading routine.  I have 
   *    already started a function designed to reconcile inconsistencies in 
   *    the configuration.
   */
  if ( par[0]->flag_eventlist && par[0]->flag_ontime ) {
    perr("Cannot specify BOTH an on-time file and an eventlist!\n");
    return -1;
  }

  /* Open ontimes files, if they exist */
  if ( par[0]->flag_ontime )
    rc = loadOnTimes(nc, par);

  /**  uncomment this if you want to revert to traditional 'dtime' mode.*
  if ( par[0]->flag_ontime ) {
    int i;
    for ( i=0; i<nc; i++ ) {
      par[i]->nevt = 1;
      par[i]->ntry = 1;
    }
  }
  **/

//  setSeeds(0,par->seed);
#ifdef _OPENMP
  int rs;
  for (rs=0; rs<omp_get_max_threads(); rs++)
    ranq2setup(rs,par[0]->seed+rs);
#else
  ranq2setup(par[0]->seed);
#endif

  rc += initAirShower(par[0]); /* initialize consts for Air Shower Sim */

  /*
   *  load the shower libraries
   */
  rc += loadShowerLibraries(par[0]);

  /*
   *  initialize atmospheric models
   */
  rc += loadAtmosDB(par[0]->atmosdbFile);
  rc += loadScatterDB(par[0]->scatterFile);

  return rc;
}

fdraw_dst_common *initFDSite(RuntimeParameters *par, UCalibration *calib, 
                       FDSiteGeometry *fdsg, TGeom *tgeom,
                       int banklist, fdraw_dst_common *fdraw) {
  int rc = 0;

  int outunit, outmode=MODE_WRITE_DST;
  char outfile[MAX_STRLEN];

  if ( strcmp(par->geometryFile,"(none)") == 0) {
    switch (par->siteid) {
      case BLACK_ROCK_SITEID:
        sprintf(par->geometryFile,"%s/fdgeom/geobr.dst.gz",RTDATA);
        break;
      case LONG_RIDGE_SITEID:
        sprintf(par->geometryFile,"%s/fdgeom/geolr.dst.gz",RTDATA);
        break;
      case MIDDLE_DRUM_SITEID:
        sprintf(par->geometryFile,"%s/fdgeom/geomd.dst.gz",RTDATA);
        break;
      case TALE_SITEID:
        sprintf(par->geometryFile,"%s/fdgeom/geotl.dst.gz",RTDATA);
        break;
      default:
        perr("Warning: unknown site; using BR geometry\n");
        sprintf(par->geometryFile,"%s/fdgeom/geobr.dst.gz",RTDATA);
    }
    perr("warning: using default geometry for site %d: %s\n",
         par->siteid,par->geometryFile);
  }
  
  
  loadFDSiteGeometry(par->geometryFile, fdsg);

  calib->nmir = fdsg->nmir;
  
  switch (par->siteid) {
    case BLACK_ROCK_SITEID:
    case LONG_RIDGE_SITEID:
      
      (calib->tc)->nmir = calib->nmir;
      
      /* Load defaults for all calibration DST banks */
      getTADefaultCalibration(par, calib->tc);

      /* Load time-independent DST banks */
      getTATimeIndependentCalibration(par, calib->tc);
      
      break;
      
    case MIDDLE_DRUM_SITEID:
      break;
    case TALE_SITEID:
      getTLDefaultCalibration(par, calib->tl);
      break;
  }

  /* Initialize RayTrace with calibration values */
  initRayTrace (fdsg, tgeom, calib);

  
  switch (par->siteid) {
    case BLACK_ROCK_SITEID:
    case LONG_RIDGE_SITEID:      
  
      /*
      *  trigger initialization
      */
      initTAElectronics();
      break;
    case MIDDLE_DRUM_SITEID:
      break;
    case TALE_SITEID:
      initTLElectronics(fdsg/*,calib->tl*/);
      break;
  }

  /* 
   *  open the output DST files and initialize them for writing
   */
  
  switch (par->siteid) {
    case BLACK_ROCK_SITEID:
      outunit = FDDSTOUT0_FILEDES;
      addBankList(banklist, BRRAW_BANKID);
      fdraw = &brraw_;
      break;
      
    case LONG_RIDGE_SITEID:
      outunit = FDDSTOUT1_FILEDES;
      addBankList(banklist, LRRAW_BANKID);
      fdraw = &lrraw_;
      break;
      
    case MIDDLE_DRUM_SITEID:
      outunit = FDDSTOUT2_FILEDES;
      break;
      
    case TALE_SITEID:
      outunit = FDDSTOUT3_FILEDES;
      addBankList(banklist, FRAW1_BANKID);
      addBankList(banklist, FTRG1_BANKID);
      fdraw = &fdraw_;
      break;
      
    default:
      outunit = FDDSTOUT4_FILEDES;
      addBankList(banklist, FDRAW_BANKID);
      fdraw = &fdraw_;
  }


  sprintf(outfile, "%s/%s", par->pdirn, par->outfn);
  rc = dstOpenUnit(outunit, outfile, outmode);
  if ( rc != 0 ) {
    fprintf(stderr, "\nunable to open output file for writing... abort!\n");
    fprintf(stderr, "  parent directory: %s\n", par->pdirn);
    fprintf(stderr, "  output file name: %s\n\n", par->outfn);
    return NULL;
  }

  return fdraw;
}
