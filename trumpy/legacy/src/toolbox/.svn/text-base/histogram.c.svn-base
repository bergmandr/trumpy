#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "constants.h"
#include "control.h"
#include "event.h"
#include "fdconstants.h"
#include "histogram.h"

static int hidx = 0;
static char histfn[MAX_STRLEN];

static boolean isOpen=FALSE;

void openNewHistogramFile(void) {
  char newfn[MAX_STRLEN];

  endHistograms();

  /* 
   *  have to handle renaming the hist file after first call, 
   *     calls to endHistograms will handle the rest.
   */
  if ( hidx == 0 ) {
    strcpy(newfn, histfn);
    char *p = strstr(newfn, ".hst");
    sprintf(p, ".%d.hst", hidx);
  }

  rename(histfn, newfn);

  perr("renamed histogram file to: %s\n", newfn);

  initHistograms(histfn);

  hidx++;
}

void initHistograms(char *path) {

  if ( isOpen )
    return;

  pout("Opening histogram file '%s'\n", path);
  openHistogramFile(path);

  strcpy(histfn, path);

  /*
   *  Define the row tags for the ntuples, since they're all the same, but 
   *    don't book them yet.
   */
  rowtag(1, "loge");
  rowtag(2, "pspec");
  rowtag(3, "xmax");
  rowtag(4, "nmax");
  rowtag(5, "zen");
  rowtag(6, "ximp");
  rowtag(7, "yimp");
  rowtag(8, "rp");
  rowtag(9, "psi");
  rowtag(10, "istrig");

/*   hbprof(NFL_ATT_PROF, "Fluorescence Attennuation", 80, 0.0, 40.0, 0.0, 1.0); // track.c */
/*   hbprof(NCVDIR_ATT_PROF, "Direct Cherenkov Attennuation", 80, 0.0, 40.0, 0.0, 1.0); // track.c */
/*   hbprof(NCVMIE_ATT_PROF, "Mie-Scattered Cherenkov Attennuation", 80, 0.0, 40.0, 0.0, 1.0); // track.c */
/*   hbprof(NCVRAY_ATT_PROF, "Rayleigh-Scattered Cherenkov Attennuation", 80, 0.0, 40.0, 0.0, 1.0); // track.c */

/*   hbook1(202, "[Dm]", 100, -50.0, 50.0); // trigger-simulator.c */
/*   hbook1(203, "[D]V", 100, -50.0, 50.0); // trigger-simulator.c */

/*   hbook1(211, "CL (all ro)", 100, -6.0, 6.0); // bank-manager.c */
/*   hbook1(212, "[Dm] (all ro)", 100, -50.0, 50.0); // bank-manager.c */
/*   hbook1(213, "[D]V (all ro)", 100, -50.0, 50.0); // bank-manager.c */
/*   hbook1(221, "CL (ro)", 100, -6.0, 6.0); // bank-manager.c */
/*   hbook1(222, "[Dm] (ro)", 100, -50.0, 50.0); // bank-manager.c */
/*   hbook1(223, "[D]V (ro)", 100, -50.0, 50.0); // bank-manager.c */
/*   hbook1(231, "CL (failed)", 100, -6.0, 6.0); // bank-manager.c */
/*   hbook1(232, "[Dm] (failed)", 100, -50.0, 50.0); // bank-manager.c */
/*   hbook1(233, "[D]V (failed)", 100, -50.0, 50.0); // bank-manager.c */

  /* in electronics.c */
  hbook1(NSBG_H1D, "[m]?NSBG!", 200, 0.0, 20.0);
  hbook1(NNOISE_TUBES_H1D, "Num. Noise Tubes per Camera", 257, -0.5, 256.5);
  hbook1(CL_H1D, "Tube Significance", 120, 0.0, 6.0); 
  hbook1(CL2_H1D, "Tube Significance", 120, 0.0, 6.0);
  hbook1(PEDRMS_H1D, "log( RMS )", 100, 0.0, 2.0);

  hbook1(INPUT_MEAN_H1D, "Input Mean FADC", 100, 200.0, 400.0);
  hbook1(INPUT_SDEV_H1D, "Input St. Dev. of FADC (16-bin sample)", 100, 5.0, 10.0);
  hbook1(OUTPUT_MEAN_H1D, "Output Mean FADC", 100, 200.0, 400.0);
  hbook1(OUTPUT_SDEV_H1D, "Output St. Dev. of FADC (16-bin sample)", 100, 5.0, 10.0);
  hbook1(DELTA_MEAN_H1D, "[D]Mean FADC (out-in)", 100, -5.0, 5.0);
  hbook1(DELTA_SDEV_H1D, "[D]St. Dev. of FADC (out-in, 16-bin sample)", 100, -2.5, 2.5);

  hbook1(RA_H1D, "V/16g^2!", 200, 0.0, 20.0);
  hbook1(TIME_H1D, "PE Time", NFADC_BINS, -0.5, (float)(NFADC_BINS-0.5));

  hbook1(NPE_H1D, "log(num)", 100, 0.0, 5.0);

  hbprof(402, "dEdX/N?CH! vs Shower Age", 100, 0.0, 2.0, 0.0, 10000.0);
  hbprof(403, "dEdX/N?CH! vs Height (km)", 100, 0.0, 50.0, 0.0, 10000.0);
  hbprof(404, "dEdX/N?CH! vs Pressure (torr)", 100, 0.0, 800.0, 0.0, 10000.0);
  hbprof(405, "y?fl!/N?CH! vs Shower Age", 100, 0.0, 2.0, 0.0, 10000.0);
  hbprof(406, "y?fl!/N?CH! vs Height (km)", 100, 0.0, 50.0, 0.0, 10000.0);
  hbprof(407, "y?fl!/N?CH! vs Pressure (torr)", 100, 0.0, 800.0, 0.0, 10000.0);

  isOpen = TRUE;
  return;
}

void addConfigHists(int iconfig) {
  char chtitl[80];

  sprintf(chtitl, "Number of Trials Since Last Success (config %d)", iconfig);
  hbook1(TRIAL_DIST_BASE+iconfig, chtitl, 50, 0.5, 50.5);  

  sprintf(chtitl, "configuration %d", iconfig);
  hbkrow(10*iconfig, chtitl, NTUPLE_NENT);

  return;
}

void endHistograms(void) {
  char newfn[MAX_STRLEN];

  if ( isOpen )
    closeHistogramFile();

  isOpen = FALSE;

  if ( hidx > 0 ) {
    strcpy(newfn, histfn);
    char *p = strstr(newfn, ".hst");
    sprintf(p, ".%d.hst", hidx);

    rename(histfn, newfn);

    perr("renamed histogram file to: %s\n", newfn);
  }

  return;
}
