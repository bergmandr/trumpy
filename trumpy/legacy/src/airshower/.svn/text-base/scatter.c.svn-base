#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "constants.h"
#include "control.h"
#include "event.h"
#include "toolbox.h"

#include "airshower.h"
#include "atmosphere-withdb.h"
#include "scatter.h"

//#define H_MIE 1.2

//#define MIE_ATTEN_LENGTH 25000.0        // meters
//#define SCALE_HEIGHT      1000.0        // meters
/* a little clunky, but should work for now... */
#ifdef ALTVAOD
#  define MIE_ATTEN_LENGTH 29400.0        // meters
// #  define MIE_ATTEN_LENGTH 14700.0        // meters
#  define SCALE_HEIGHT      1000.0        // meters
#else
#  define MIE_ATTEN_LENGTH (fdscat_.hzalen)        // meters
#  define SCALE_HEIGHT     (fdscat_.schght)        // meters
#endif
#define MIXING_HEIGHT        0.0        // meters

static boolean initialized = FALSE;
//static const double bb=0.05, xz=860.0, atscl=7.5, hMie=1.2;
static double ozab[NWAVELEN_BANDS];
static double btkp[NWAVELEN_BANDS];
//static double lMie[NWAVELEN_BANDS];
//static double fnu[NWAVELEN_BANDS];
static double mie_atten_length[NWAVELEN_BANDS];
static const double ozone[40] = {
  0.00, 3.20, 6.10, 8.60, 10.7, 12.7, 14.7, 16.6, 18.8, 21.2, 
  24.2, 28.3, 33.9, 41.4, 50.4, 60.0, 69.9, 80.4, 92.0, 105.,          
  121., 138., 157., 176., 196., 214., 231., 246., 260., 271., 
  281., 289., 296., 302., 307., 312., 315., 318., 321., 323.
};

static void init(void) {
  int i;
  double y1, y2;
  double wl;
  double s0;

  wl = LAMBDA0 + DLAMBDA / 2.0;
  wl *= 1.e9;
  s0 = 0.1;		// btkp at 334 nm (from mcru)
  for ( i=0; i<NWAVELEN_BANDS; i++ ) {
    ozab[i] = exp( 5.25 - 11.0*(wl-280.0)/80.0 );
    //    fnu[i] = exp( 2.4 - 0.56*(wl-280.0)/100.0 );
    y1 = 0.095 + 1.0e-6 * ( wl - 405.0 ) * ( wl - 405.0 );
    y2 = 0.080 + 9.0e-6 * ( wl - 417.5 ) * ( wl - 417.5 );
    btkp[i] = min(y1, y2);
    mie_atten_length[i] = MIE_ATTEN_LENGTH * (s0 / btkp[i]);
    //    lMie[i] = (0.04/(0.86*0.05)) * 1.2/btkp[i];

    wl += DLAMBDA*1.e9;
  }
  
  initialized = TRUE;
}

/*
 * Computes the attenuation factors (0<u<1) for each wavelength band between 
 *   altitudes 'ha' & 'hb', for a segment of length 'dl'.  This code is taken 
 *   directly from MCRU (attenn.f).  The constants are extrapolated from those 
 *   provided in MCRU.  This code is functional, but will need to be 
 *   re-written for clarity (SS).
 *
 * Inputs:
 *   double ha - altitude of one edge of shower segment in meters above sea 
 *                 level.
 *   double hb - atlitude of the opposite edge in meters above sea level.
 *   double dl - segment length in meters.
 *
 * Output:
 *   double *att - attenuation factors (0-1) for each wavelength band.
 */

/*
 * (DRB)
 * The btkp values computed below are based on a quadratic fit to the
 * BTKP values in ATTENN.F from MCRU.  The code specifying these
 * values was:
 *       DATA BTKP/.110,.108,.106,.105,.103,.101,.100,.099,.098,
 *      +      .097,.096,.093,.088,.084,.081,.080/
 * where one would access BTKP for a given wavelength, WL, in nm:
 *       IWL = 9*INT(WL-280.0)
 *       FOO = BTKP(IWL)
 *
 * The ozab values computed below are based on a exponential fit to
 * the BTKP values in ATTENN.F from MCRU.  The code specifying these
 * values was:
 *      DATA OZAB/119.,66.2,12.94,4.05,1.26,.47,.163,.025,
 *     +     .0065,.0013,6*0./
 * where one would access OZAB for a given wavelength, WL, in nm:
 *       IWL = 9*INT(WL-280.0)
 *       FOO = OZAB(IWL)
 */

//void getAttenuation(int season, double ra[3], double rb[3], double *att) {
void getAttenuation(fdatmos_param_dst_common *season, 
		    double ra[3], double rb[3], double *att) {
  int i;
  double dl, dr[3], dx;
  double lat, lon, ha, hb;
  double miefrac[NWAVELEN_BANDS];
  double ozonefrac[NWAVELEN_BANDS];
  double rayleighfrac[NWAVELEN_BANDS];

  if ( ! initialized ) init();

  subvec(ra, rb, dr);
  dl = sqrt( dotprod(dr, dr) );

  /* initialize attenuation factors */
  for ( i=0; i<NWAVELEN_BANDS; i++ )
    att[i] = 1.0;  // no need to get fancy (SS 10.4.11)
  
  if ( dl < 1.0 )
    return;

  r2lla(ra, &lat, &lon, &ha);
  r2lla(rb, &lat, &lon, &hb);

  getMieScatterFraction(ha, hb, dl, miefrac);

  /* find amount of atmosphere in path length */
  dx = getSlantDepth(season, ra, rb);
  getRayleighScatterFraction(dx, rayleighfrac);

  getOzoneAbsorbtionFraction(ha, hb, dl, ozonefrac);

  for ( i=0; i<NWAVELEN_BANDS; i++ )
    att[i] = ( 1.0 - miefrac[i] ) * ( 1.0 - rayleighfrac[i] ) 
#ifndef JPCOMPARE
      * ( 1.0 - ozonefrac[i] )
#endif
      ;
}

/* Taken from mcru's aero_trans.f */
void getMieScatterFraction(double ha, double hb, double dl, double *frac) {
  int i;
  double h1, h2, e1, e2;
  double lambda, arg;

  if ( !initialized) init();

  h1 = min(ha, hb) - CLF_ALTITUDE;
  h2 = max(ha, hb) - CLF_ALTITUDE;

  if (h1 < 0.)
    h1 = 0.;
  if (h2 < 0.)
    h2 = 0.;
  
  e1 = exp( -(h1 - MIXING_HEIGHT) / SCALE_HEIGHT );
  e2 = exp( -(h2 - MIXING_HEIGHT) / SCALE_HEIGHT );

  for (i=0; i<NWAVELEN_BANDS; i++) {
    lambda = mie_atten_length[i];

    /* horizontal case */
    if ( (h2 - h1) < 1.) {
      if (h1 > MIXING_HEIGHT)
        lambda *= exp ((h1-MIXING_HEIGHT) / SCALE_HEIGHT);
      arg = -dl / lambda;
      frac[i] = 1. - exp(arg);
    }
    /* within mixing layer */
    else if (h2 <= MIXING_HEIGHT) {
      arg = -dl / lambda;
      frac[i] = 1. - exp(arg);
    }
    /* above mixing layer (this is the most common) */
    else if (h1 >= MIXING_HEIGHT) {
      arg = (e2 - e1) * dl * SCALE_HEIGHT / ((h2 - h1) * lambda);
      frac[i] = 1. - exp(arg);
    }
    /* crossing mixing layer */
    else {
      arg = (e2 - 1.) * dl * SCALE_HEIGHT / ((h2 - h1) * lambda);
      frac[i] = exp( -(MIXING_HEIGHT - h1) * dl / ((h2 - h1) * lambda) );
      frac[i] *= exp(arg);
      frac[i] = 1. - frac[i];
    }
  }
}


/*
 * 2010/01/19 ELB
 * The mie scattering transmission function from the Java code.
 *
 * This is essentially copied from the Java function 
 * USAtmosphere.getTransmissionOfMie( ), however I made changes
 * so that it takes the same parameters and returns the fraction 
 * in the same way as the above function getMieScatterFraction
 * 
 * Inputs:
 * ha -- height of starting point (for example shower segment)
 * hb -- height of ending point (for example FD station)
 * dl -- distance between segments
 * Returns:
 * frac -- Mie scattering fraction in wavelength bins
 */
void getMieScatterFractionJava(double ha, double hb, double dl,
			       double *frac){
  
  /*
   * Will leave the Java variable names even though the style
   * is not relevant, so this method can be followed in the Java 
   * easily.  
   * We can always change this easily...
   */
  double fMieHeight = 1.0e3; // meters
  double fMieLength = 2.94e4; // meters 
  /*
   * Of course, if this is used in the long run we should
   * stay consistance with the rest of this file
  double fMieHeight = MIE_ATTEN_LENGTH;
  double fMieHeight = SCALE_HEIGHT;
   */
  
  double fMieBaseWaveLength = 360.0; // nm
  double fGroundHeight = 1430.0; // meters

  double h = ha - fGroundHeight;
  double hd = hb - fGroundHeight;
  /*
   * Calculate theta this way.  This add a slight difference
   * between the Java code and TRUMP.  The Java code uses a 
   * function in the J3Vector class
   */
  double costheta = (ha - hb)/dl;
  double wl = LAMBDA0;
  int i;
  for (i=0; i<NWAVELEN_BANDS; i++) {

    double lam = fMieBaseWaveLength / wl / 1.0e9;
    double logT = (exp(-h / fMieHeight) - exp(-hd / fMieHeight))
      * fMieHeight / fMieLength / costheta * lam;

    if (logT > 0.0) {
      logT = 0.0;
    }
    frac[i] = 1.0 - exp(logT);
    wl += DLAMBDA;
  }
}

/*
 *   Taken from mcru's aerp.f using the Longtin desert aerosol phase 
 *   function (mcru really used the "Jerry Elbert" model).
 */
#define BINSIZ (M_PI/90.0)
static const double phs[91] = {
  8.00e+02, 2.85e+00, 9.69e-01, 7.51e-01, 6.65e-01, 6.06e-01, 5.56e-01, 
  5.11e-01, 4.69e-01, 4.30e-01, 3.94e-01, 3.60e-01, 3.29e-01, 3.00e-01, 
  2.73e-01, 2.49e-01, 2.27e-01, 2.06e-01, 1.88e-01, 1.70e-01, 1.55e-01, 
  1.43e-01, 1.32e-01, 1.20e-01, 1.08e-01, 9.65e-02, 8.95e-02, 8.24e-02, 
  7.54e-02, 6.83e-02, 6.13e-02, 5.70e-02, 5.27e-02, 4.84e-02, 4.41e-02, 
  3.98e-02, 3.72e-02, 3.46e-02, 3.21e-02, 2.95e-02, 2.69e-02, 2.53e-02, 
  2.37e-02, 2.20e-02, 2.04e-02, 1.88e-02, 1.79e-02, 1.71e-02, 1.62e-02, 
  1.54e-02, 1.45e-02, 1.40e-02, 1.36e-02, 1.31e-02, 1.27e-02, 1.22e-02, 
  1.20e-02, 1.18e-02, 1.16e-02, 1.14e-02, 1.12e-02, 1.12e-02, 1.12e-02, 
  1.12e-02, 1.12e-02, 1.12e-02, 1.13e-02, 1.15e-02, 1.16e-02, 1.18e-02, 
  1.19e-02, 1.22e-02, 1.25e-02, 1.27e-02, 1.30e-02, 1.33e-02, 1.42e-02, 
  1.52e-02, 1.61e-02, 1.71e-02, 1.80e-02, 2.01e-02, 2.21e-02, 2.42e-02, 
  2.62e-02, 2.83e-02, 3.07e-02, 3.30e-02, 3.54e-02, 3.79e-02, 4.01e-02
};

double mieScatterFunction(double q) {
  int k;
  double q0, lphs;

  k = (int) floor( q / BINSIZ );
  q0 = (double)k * BINSIZ;

  lphs = log( phs[k] );
  lphs += ( q - q0 )*log( phs[k+1] / phs[k] )/BINSIZ;

  return exp( lphs );
}
#undef  BINSIZ

/* Taken from mcru's raylp.f */
void getRayleighScatterFraction(double dx, double *frac) {
  int i;
  double wl, xr;

  wl = LAMBDA0 * 1.0e+10;	// get wavelength in Angstroms
  for ( i=0; i<NWAVELEN_BANDS; i++ ) {
    /* 1.16e-11 = 2969.6 / 4000^4
     * 2969.6 g/cm2 = rayleigh attenuation length
     * 4000^4 is the wavelength scale factor (angstroms)
     * see page 101 in Andreas Zech's thesis
     */
    xr = 1.16e-11 * pow4(wl);
    frac[i] = 1.0 - exp( -dx / xr );

    wl += DLAMBDA * 1.0e+10;
  }
}

/* Taken from mcru's raylp.f */
#define NORM 5.96831037e-02  /* = 1 / (16 * pi / 3) */
double rayleighScatterFunction(double q) {
  double cq = cos(q);

  return ( 1.0 + cq*cq ) * NORM;
}
#undef  NORM

/* Taken from mcru's atten_rayl_ozon.f */
void getOzoneAbsorbtionFraction(double ha, double hb, double dl, double *fr) {
  int i, k1, k2;
  double h1, h2, dh, coslp;
  double del, oz1, oz2;
  double dodh, dodh1, dodh2;
  double arg, ozfac;

  if ( ! initialized ) init();

  dl *= 1.0e-3;

  h1 = max(ha, hb) * 1.0e-3;
  h2 = min(ha, hb) * 1.0e-3;
  dh = h1 - h2;

  coslp = dh / dl;

  /* 
   * Find ozone penetration depth in units of centimeters at atmospheric 
   *   pressure. */
  k1 = (int)( floor(h1) );
  if ( k1 < 39 ) {
    del = h1 - (double)k1;
    oz1 = ( 1.0 - del )*ozone[k1] + del*ozone[k1+1];
    dodh1 = ( ozone[k1+1] - ozone[k1] ) * 1.0e-03;
  }
  else {
    oz1 = 325.0;
    dodh1 = 0.002;
  }

  k2 = (int)( floor(h2) );
  if ( k2 < 39 ) {
    del = h2 - (double)k2;
    oz2 = ( 1.0 - del )*ozone[k2] + del*ozone[k2+1];
    dodh2 = ( ozone[k2+1] - ozone[k2] ) * 1.0e-03;
  }
  else {
    oz2 = 325.0;
    dodh2 = 0.002;
  }

  ozfac = ( oz1 - oz2 ) * 1.0e-03;

  for ( i=0; i<NWAVELEN_BANDS; i++ ) {
    dodh = ( dodh1 + dodh2 ) / 2.0;
    if ( dh < 0.1 )
      arg = dodh * ozab[i] * dl;
    else
      arg = ozfac * ozab[i] / coslp;
    
    fr[i] = 1.0 - exp( -arg );
  }
}

#undef MIE_ATTEN_LENGTH
#undef SCALE_HEIGHT
#undef MIXING_HEIGHT
