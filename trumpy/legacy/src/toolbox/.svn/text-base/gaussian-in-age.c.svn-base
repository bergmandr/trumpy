#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "constants.h"
#include "control.h"
#include "event.h"

#include "airshower.h"
#include "showerlib.h"

/*
 *  double gaussianInAgeFunction(GaisserHillasParameters, double)
 *
 *  Author:  Tom Stroman, University of Utah, 
 *      adapted from gaisserHillasFunction by Sean R. Stratton
 *           Rutgers University, Dept. of Physics & Astronomy
 *
 *  Compute the Gaussian-in-age function using the given parameter set at 
 *    atmospheric depth 'x'.  The Gaussian-in-age function (via Bill Hanlon) is:
 *
 *  GIA(x) = N_{max}\exp\left( -\frac{2}{\sigma_s^2} \left({\frac{x-x_m}{x+2x_m}\right)^2\right)
 *
 *  Input:
 *    GaisserHillasParameters *gh -- 3 parameters for GIA, packaged
 *    in a set of Gaisser-Hillas fit parameters, 
 *      including 'X0' (g/cm^2) as a stand-in for 'sigma_s', 
 *      'Xmax' (g/cm^2), 
 *      'Nmax', and ignores 'lambda' (g/cm^2)
 *
 *    double x -- Depth of atmosphere to evaluate GIA function in g/cm^2
 *
 *  Returns:
 *    The value of the Gaussian-in-age function with the given parameters at 
 *      atmospheric depth 'x'.
 */
double gaussianInAgeFunction(GaisserHillasParameters *gh, double x) {
  double a, b;
  a = (x - gh->xmax)/(x + 2*gh->xmax);
  a *= a;
  
  b = -2*a/(gh->x0*gh->x0);
  return gh->nmax * exp( b );
}
