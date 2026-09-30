#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "constants.h"
#include "control.h"
#include "airshower.h"
#include "event.h"
#include "atmosphere-withdb.h"
#include "nerling.h"
#include "flash.h"

#define KPA_TO_TORR 7.5006

/*
 * FLASH spectrum at 155 Torr, recombined into 5 nm bins from 300-420 nm
 * This should be normalized to unity, but with the recombination (and the
 * parts on the ends) this may no longer be strictly true.
 * Taken from Abbasi et al., Astropart. Phys. 29 (2008) 77, Fig 7.
 */
/*
static const double flashSpectrumBins[NWAVELEN_BANDS] = { 
  0.0012, 0.0033, 0.0569, 0.0702, 0.0017, 0.0078, 
  0.0139, 0.2494, 0.0017, 0.0038, 0.0597, 0.1759, 
  0.0019, 0.0081, 0.0251, 0.0933, 0.0461, 0.0361, 
  0.0846, 0.0330, 0.0089, 0.0157, 0.0000,-0.0010 };
*/

/*
 * FLASH spectrum at 155 Torr, by lines at 155 Torr
 * Taken from Abbasi et al., Astropart. Phys. 29 (2008) 77, Fig 9.
 */
/*
static const double flashSpectrumLines[][2] = {
  {311.7,0.0051},
  {313.6,0.0430},
  {315.9,0.0890},
  {328.5,0.0111},
  {330.9,0.0014},
  {333.9,0.0023},
  {337.1,0.2617},
  {346.9,0.0020},
  {350.0,0.0030},
  {353.7,0.0579},
  {357.7,0.1641},
  {367.2,0.0081},
  {371.1,0.0173},
  {375.6,0.0544},
  {380.5,0.0977},
  {391.4,0.1167},
  {394.3,0.0053},
  {399.8,0.0400},
  {405.9,0.0207},
  {414.1,0.0004},
  {420.1,0.0000}
};
*/

/* 
 *  Sean's extraction of Figure 9 from Abbasi et al. (2008) 
 */
static const double flashSpectrumLines[21][2] = {
  {311.7, 0.005133},
  {313.6, 0.042940}, 
  {315.9, 0.088875}, 
  {328.5, 0.011127}, 
  {330.9, 0.001427}, 
  {333.9, 0.002283}, 
  {337.1, 0.261339}, 
  {346.9, 0.001998}, 
  {350.0, 0.002996}, 
  {353.7, 0.057774}, 
  {357.7, 0.163908}, 
  {367.2, 0.008133}, 
  {371.1, 0.017262}, 
  {375.6, 0.054351}, 
  {380.5, 0.097574}, 
  {391.4, 0.116546}, 
  {394.3, 0.005279}, 
  {399.8, 0.039942}, 
  {405.9, 0.020684}, 
  {414.1, 0.000429}, 
  {420.1, 0.000000}
};

/* Binned version of the above line spectrum
static const double flashSpectrum[NWAVELEN_BANDS] = { 
  0.0000, 0.0000, 0.0481, 0.0890, 0.0000, 0.0111, 
  0.0037, 0.2617, 0.0000, 0.0020, 0.0611, 0.1641, 
  0.0000, 0.0081, 0.0173, 0.0544, 0.0977, 0.0000, 
  0.1220, 0.0400, 0.0000, 0.0207, 0.0004, 0.0000
};
*/

/* Line spectrum, set up for auto-binning */
static double flashNormalization = 0.0;
static double flashSpectrum[NWAVELEN_BANDS] = { 0.0 };

/* Sean's calculation of the binned line spectrum
static const double flashSpectrum[NWAVELEN_BANDS] = {
  0.000000, 0.000000, 0.048073, 0.088875, 0.000000, 0.011127, 
  0.003710, 0.261339, 0.000000, 0.001998, 0.060770, 0.163908, 
  0.000000, 0.008133, 0.017262, 0.054351, 0.097574, 0.000000, 
  0.121825, 0.039942, 0.000000, 0.020684, 0.000429, 0.000000
};
*/

/*
 * (DRB) I fit the flash FY/m to the form A.P/(1+B.P).  With the HR
 * filter (one extra point), I found B = 0.0844(60).  I then fit the
 * FY/MeV (no filter) to A'.P/(1+B.P)/P = A'/(1+B.P) with the B from
 * above.  I found A' = 1386(10).  If I fit the FY/m (no filter) with
 * B fixed, I find A = 0.4352(31).
 */
static const double flashAp = 1386.;
/* static const double flashA  = 0.4352; */
static const double flashB  = 0.0844;
static const double flashT = 304.;

//static double flashFYByMeter(double P, double T) {return (flashA*P)/(1.+flashB*P*sqrt(T/flashT));}
inline double flashFYByMeV(double P, double T) {
  return flashAp/(1.+flashB*P*sqrt(T/flashT));
}

/* From Kakimoto et al., NIMA 372 (1996) 527. */
static const double kA1 = 89.0;
static const double kB1 = 1.85;
static const double kA2 = 55.0;
static const double kB2 = 6.50;
static const double kdEdX = 1.668; // From MCRU::newlight
double kakimotoFYByMeter(double rho, double T) {
  double sqrtT = sqrt(T);
  return rho*(kA1/(1.+rho*kB1*sqrtT)+kA2/(1.+rho*kB2*sqrtT));
}

/*
 *  This function really only dependent on altitude z, since rho and T depend 
 *    only on z.  For z->0, rho ~1.25 kg/m^3, T ~300K, fy/MeV = 15.3 photons 
 *    per MeV.  For z->infinity, rho ~0, T won't matter, fy/MeV = 863 photons 
 *    per MeV.
 */
double kakimotoFYByMeV(double rho, double T) {
  double sqrtT = sqrt(T);
  double a = kA1 / ( 1.0 + rho*kB1*sqrtT );
  double b = kA2 / ( 1.0 + rho*kB2*sqrtT );
  return 10.0 * ( a + b ) / kdEdX;
}

/*
 * The FYbyMeter (fl) and FYbyMeV (fe) are related like this:
 *   fl = fe * dE/dX * rho
 * Thus to find the dE/dX that FLASH must have used we solve for it
 *   dE/dX = fl/fe * 1/rho
 * Using 760 Torr, 304 K, and 1.205 mg/cc we have
 *    dE/dX = 5.077/21.28 / 12.05 = 1.980
 */
/* static const double flashdEdX = 1.980; */

/* Returns the flourescence yield in photons/MeV */
//void getFluorescenceYield(int season, double h, double *yfl) {
void getFluorescenceYield(fdatmos_param_dst_common *season,
				 double h, double *yfl) {
  int i;
  double P = KPA_TO_TORR*getPressureByAltitude(season,h);
  double T = getTemperatureByAltitude(season,h);
  double fy = flashFYByMeV(P,T)/1.e6;

  for ( i=0; i<NWAVELEN_BANDS; i++ )
    yfl[i] = fy * flashSpectrum[i];
}

/* Returns the flourescence yield in photons/m */
//void getFluorescenceYieldByMeter(int season, double age, double h, double *yfl) {
inline void getFluorescenceYieldByMeter(fdatmos_param_dst_common *season, 
					double age, double h, double *yfl) {
  int i;

  // double P = KPA_TO_TORR*getPressureByAltitude(season,h);
  double rho = getDensityByAltitude(season,h);
  double T = getTemperatureByAltitude(season,h);
  // double fy = flashFYByMeter(P,T);
  double fy = kakimotoFYByMeter(1000.*rho,T);
  // double dEdXRatio = getAlphaEff(age)/flashdEdX/1.e6;
  double dEdXRatio = getAlphaEff(age)/kdEdX/1.e6;
  fy *= dEdXRatio;

  for ( i=0; i<NWAVELEN_BANDS; i++ )
    yfl[i] = fy * flashSpectrum[i];
}
/* Returns the flourescence yield in photons/m */
//void getFluorescenceYieldByMeter(int season, double age, double h, double *yfl) {

void getFluorescenceYieldByMeV(fdatmos_param_dst_common *season, double h, double *yfl) {
  int i, k;
  double lambda;
  if ( fabs(flashNormalization) < DOUBLE_FP_PRECISION ) {
    for ( i=0; i<21; i++ ) {
      lambda = flashSpectrumLines[i][0] * 1.0e-9;
      k = (int)((lambda-LAMBDA0)/DLAMBDA);

      if ( k >= 0 && k < NWAVELEN_BANDS )
	flashSpectrum[k] += flashSpectrumLines[i][1];

      if ( lambda >= 300.0e-9 && lambda < 400.0e-9 )  // the values '300' and '400' are unrelated to TRUMP, they are specific to the Kakimoto FY normalization
	flashNormalization += flashSpectrumLines[i][1];
    }

    /**
    pout("flashNormalization=%lf\n", flashNormalization);
    for( i=0; i<NWAVELEN_BANDS; i++ ) {
      double l = LAMBDA0 + (double)i*DLAMBDA + DLAMBDA/2.0;
      pout("%.1lf  %e\n", 1.0e9*l, flashSpectrum[i]/flashNormalization);
    }
    **/
  }

  double p, t, rho;
  getPTDByAltitude(season, h, &p, &t, &rho);

  double fy = kakimotoFYByMeV(1000.0*rho, t) / flashNormalization;
  for ( i=0; i<NWAVELEN_BANDS; i++ )
    yfl[i] = fy * flashSpectrum[i];
}

#undef KPA_TO_TORR 
