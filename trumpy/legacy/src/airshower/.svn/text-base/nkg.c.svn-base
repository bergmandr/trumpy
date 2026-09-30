#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "constants.h"
#include "control.h"
#include "event.h"
#include "atmosphere-withdb.h"
#include "random.h"
#include "nkg.h"

/*
 * Gets moliere radius in air (meters).
 * Taken from mcru's kmolu.f
 * Mysterious constant c = (20/76) * 36.4
 */
double getMoliereRadius (ATMP_t *season, double altitude) {

  double c = 9.578947368;
  double density = getDensityByAltitude (season, altitude);

  return c / (density * 100.);
}

/*
 * Given the age, rp unit vector and SDP normal
 * of the shower, returns the vector v[] that displaces
 * the light from the shower axis.
 */
void getNKGdr (double age, double molrad, 
	       double rpuv[], double npln[], double v[]) {
  int i;

  double r1, r2;
  double rm, phi, x, y;

  if ( fabs(age) < DOUBLE_FP_PRECISION ) {
    for (i=0; i<3; i++)
      v[i] = 0.;
    return;
  }

  r2 = 99.;

  do {

    r1 = RANDOM_NUMBER;
    r2 = RANDOM_NUMBER;

    rm = 3. * molrad * pow(r1, 1./age);
    phi = 2. * M_PI * RANDOM_NUMBER;

    x =	rm * sin(phi);
    y =	rm * cos(phi);

    for (i=0; i<3; i++)
      v[i] = x * npln[i] + y * rpuv[i];

  }  while ( r2 > pow(1+rm/molrad, age-4.5) );

}
