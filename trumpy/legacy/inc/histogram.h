#ifndef _HISTOGRAM_H_
#define _HISTOGRAM_H_

#include "sshbook.h"

#define NTUPLE_NENT        10

/*
 *  Define some labels to stand for histogram ID numbers
 */
#define NFL_ATT_PROF      101
#define NCVDIR_ATT_PROF   102
#define NCVMIE_ATT_PROF   103
#define NCVRAY_ATT_PROF   104

#define NSBG_H1D          301
#define NNOISE_TUBES_H1D  302
#define CL_H1D            303
#define PEDRMS_H1D        304
#define INPUT_MEAN_H1D    305
#define INPUT_SDEV_H1D    306
#define OUTPUT_MEAN_H1D   307
#define OUTPUT_SDEV_H1D   308
#define DELTA_MEAN_H1D    309
#define DELTA_SDEV_H1D    310
#define RA_H1D            311
#define TIME_H1D          312
#define CL2_H1D           313

#define NPE_H1D           401

#define TRIAL_DIST_BASE   900

void initHistograms(char *path);

void addConfigHists(int iconfig);

void openNewHistogramFile(void);

void endHistograms(void);

#endif
