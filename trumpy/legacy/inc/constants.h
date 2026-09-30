#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

/* 
 *  physical constants 
 */
#define SPEED_OF_LIGHT       2.99792458e-01  /* m/ns */
#define FINE_STRUCTURE_CONST 7.297352570e-03
#define ELECTRON_MASS        5.10991918e+05  /* eV/c^2 */
#define PLANCK_CONSTANT      4.13566733e-15  /* eV s */
#define TWO_PI_ALPHA         4.5850618e-2
#define CELSIUS_TO_KELVIN    273.15
#define TOP_OF_ATMOSPHERE    6.0e+04         /* meters */
#define JOULE_PER_EV         1.60217646e-19
#define HC                   1239.84187e-09  /* Planck's const * speed of light (eV*m) */

/* earth's mean radius in m at 39.29693N latitude, 1382m altitude. */
// #define EARTH_RADIUS 6370984.3

/* Earth's mean radius in m at 39.2969179 N latitude, 1370.046 m altitude */
#define EARTH_RADIUS 6370972.3
/* 
 *  operation constants 
 */
enum {
  NULL_FILEDES, 
  FDDSTOUT0_FILEDES,
  FDDSTOUT1_FILEDES, 
  FDDSTOUT2_FILEDES, 
  FDDSTOUT3_FILEDES,
  FDDSTOUT4_FILEDES,
  FDGEOM_FILEDES, 
  SHOWERLIB_FILEDES, 
  PARAGLAS_FILEDES, 
  BG3_FILEDES, 
  PMTQE_FILEDES, 
  PMTUNIF_FILEDES, 
  MIRROR_FILEDES, 
  PMTGAIN_FILEDES, 
  PMTCAL0_FILEDES,
  PMTCAL1_FILEDES,
  PMTCAL2_FILEDES,
  ONTIME_FILEDES, 
  ATMOSDB_FILEDES, 
  FDSCAT_FILEDES
};

#define PSPEC_PROTON     14
#define PSPEC_IRON     5626

#define DEDX_MODEL  1   // this is default
// #define DEDX_MODEL  2     // this is temporary, for studying Java shower lib
// #define DEDX_MODEL 14   // TAS QGSJet01c
// #define DEDX_MODEL 15   // TAS QGSJetII-04


#define NWAVELEN_BANDS 50
#define LAMBDA0 250.0e-09
#define DLAMBDA   5.0e-09

/* computation constants */
typedef int boolean;
#define TRUE   1
#define FALSE  0

#define DOUBLE_FP_PRECISION 0.1110223024625157e-15 /* 1/2^53 */
#define MAX_DOUBLE 1.7976931348623157e+308   /* largest finite number from 0 */

#define HALF_DEGREE  8.7266462599716479e-03
#define LN02         6.9314718055994531e-01
#define LN10         M_LN10                  // 2.3025850929940457e+00
#define LOGEXP       M_LOG10E                // 4.3429448190325183e-01
#define COS15        9.659258263e-01         /* cos(15) in degrees */

#endif
