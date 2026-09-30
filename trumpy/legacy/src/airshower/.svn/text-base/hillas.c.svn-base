#include <math.h>
#include "constants.h"
#include "hillas.h"

/*
 * Taken from Hillas J. Phys. G: Nucl. Phys. 8 (1982), p. 1466
 * and mcru's hillas_T.f
 * T refers to the fraction of electrons above threshold
 * at a given shower age.
 */

double getHillasE0 (double s) {
  /*
   * In the paper, the if statement is s >= 0.4
   * mcru has done it this way, so we are following suit.
   */
  double E0;
  if (s >= 0.0)
    E0 = 44.0 - 17.0 * (s-1.46)*(s-1.46);
  else
    E0 = 26.;

  return 1.0e6 * E0;	// E0 in eV
}

double getHillasT (double energy, double s) {
  double E  = energy * 1.0e-6;		// E in MeV
  double E0 = getHillasE0(s) * 1.0e-6;	// E0 in MeV

  double a = (0.89*E0 - 1.2)/(E0 + E);
  double b = 1. + 1.0e-4*s*E;

  return pow(a, s)/(b*b);
}

double getHillasdTdE (double energy, double s) {
  double E  = energy * 1.0e-6;		// E in MeV
  double E0 = getHillasE0(s) * 1.0e-6;	// E0 in MeV

  double a = 1.0 / (E0 + E);
  double b = 2.0e-4 / (1.0 + 1.0e-4*s*E);
  double dTdE = getHillasT (energy, s) * s * (a + b);

  return dTdE;
}

/* Taken from mcru's simshwf1.f */
double getHillasIntegrand (double energy, double s, double delta) {
  double ceff, dTdE;

  if (energy > 1.0e6) {
    ceff = 2.0 * delta - ELECTRON_MASS * ELECTRON_MASS / 
           (energy * energy);
    dTdE = getHillasdTdE (energy, s);
    return ceff * dTdE;
  }
  else
    return 0.;
}
