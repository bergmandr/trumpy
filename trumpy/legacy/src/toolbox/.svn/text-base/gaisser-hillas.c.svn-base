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
 *  double gaisserHillasFunction(GaisserHillasParameters, double)
 *
 *  Author:  Sean R. Stratton
 *           Rutgers University, Dept. of Physics & Astronomy
 *
 *  Compute the Gaisser-Hillas function using the given parameter set at 
 *    atmospheric depth 'x'.  The Gaisser Hillas Function is:
 *
 *  GH(x) = N_{max} \left( \frac{x-x_0}{x_m-x_0} \right)^{\frac{x_m-x_0}{\lambda}} e^{\frac{x_m-x}{\lambda}}
 *
 *  Input:
 *    GaisserHillasParameters *gh -- Set of Gaisser-Hillas fit parameters, 
 *      including 'X0' (g/cm^2), 'Xmax' (g/cm^2), 'Nmax', and 'lambda' 
 *      (g/cm^2).
 *
 *    double x -- Depth of atmosphere to evaluate GH function in g/cm^2
 *
 *  Returns:
 *    The value of the Gaisser-Hillas function with the given parameters at 
 *      atmospheric depth 'x'.
 */
double gaisserHillasFunction(GaisserHillasParameters *gh, double x) {
  double a, b;

  a = ( x - gh->x0 ) / ( gh->xmax - gh->x0 );
  if ( a < 0.00 )
    return 0.00;

  b = gh->xmax - x + ( gh->xmax - gh->x0 ) * log( a );
  b = b / gh->lambda;

  return gh->nmax * exp( b );
}
