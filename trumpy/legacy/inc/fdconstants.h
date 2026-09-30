#ifndef _FDCONSTANTS_H_
#define _FDCONSTANTS_H_

#define NCAMERAS_SITE  (GEOFD_MAXMIR)
#define NTUBES_CAMERA  (GEOFD_MIRTUBE)
#define MAX_NADJCAM     5

/*
 *  Constants specific to FD simulation
 */
#define PMT_NUMDIV	80
#define PMT_BIN		 1

#define NFADC_BINS fdraw_nt_chan_max
#define FRAME_SIZE 256
#define HALF_FRAME_SIZE 128

#define CALIB_LAMBDA   337.1    /* nm */
#define CLF_WAVELENGTH 355.0e-9 /* meters  */
// #define CLF_WAVELENGTH 337.1e-9 /* meters */

#ifdef JPCOMPARE
#  warning "compiling with JPCOMPARE!"
#  define RHO_ALPHA 1.18
#  define RXF_FACTOR 1.00
#  define VOLTAGE_OFFSET 0.282
#  define NOISE_AMPLITUDE 0.9
#  define NSBG_FUDGE 0.982  // night-sky background fudge factor
#else
/* the nominal values... */
//#  define RXF_FACTOR (1./1.15) // CRAYS correction for calibration versions before 1.4
#  define RXF_FACTOR 1.0  // for calibration 1.4
#  define RHO_ALPHA 1.057
//#  define VOLTAGE_OFFSET 0.0
#  define VOLTAGE_OFFSET -0.175  // tuned for good D/MC comparison
#  define NOISE_AMPLITUDE 1.0
#  define NSBG_FUDGE 1.0

/* values that make for good D/MC comparison */
#endif

#define NIGHT_SKY_BG    9.0  /* #pe/100ns */

#ifdef NO_NOISE
#undef NIGHT_SKY_BG
#undef NOISE_AMPLITUDE
#define NIGHT_SKY_BG 0.0
#define NOISE_AMPLITUDE 0.0
#endif


#define BLACK_ROCK_SITEID   0
#define LONG_RIDGE_SITEID   1
#define MIDDLE_DRUM_SITEID  2
#define TALE_SITEID         3
/* 
 *  latitudes & longitudes are in radians, altitudes are in meters 
 */
/***
#define CLF_LATITUDE      0.685860814
#define CLF_LONGITUDE    -1.970629442
#define CLF_ALTITUDE      1382.0
***/

// #define BR_LATITUDE       0.683964863
// #define BR_LONGITUDE     -1.967190271
// #define BR_ALTITUDE       1404.0
// 
// #define LR_LATITUDE       0.684307296
// #define LR_LONGITUDE     -1.974342106
// #define LR_ALTITUDE       1554.0
// 
// #define MD_LATITUDE       0.688930674
// #define MD_LONGITUDE     -1.972111401
// #define MD_ALTITUDE       1600.0

// #include "parameters.h"
#define nhit_max 3584
#define nhit_mir_max 320
// #define max_tim 1024
#define max_tim 2048
#define npe_max 500
#define ntube_eye 3584
#define nchan_mir 320
// #define nt_chan_max 1024
#define nt_chan_max 2048
#define nhit_dsp_max 100


// formerly in structure glb_prm_
#define DT_ADC 100.00
#define THRESH_1 32.0
#define THRESH_2 22.0

#define DB_TRIG_VERSION 3
#define DB_ADJ_OPT 2

#define TAU_FADC 35.
#define TAU_16 500.
#define TAU_NEG 500000.

#endif
