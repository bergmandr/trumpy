#include <math.h>

#include "random.h"

/*
 * RNG control functions.
 */
static int hiseed=0x02468ACE;
static int loseed=0x13579BDF;

void getSeeds(int *hi, int *lo) {
  *hi = hiseed;
  *lo = loseed;
}

void setSeeds(int his, int los) {
  hiseed = his;
  loseed = los;
}

double getRandomNumber() {
  return ranss(&hiseed, &loseed);
}

/*
 *  int prand(double nu) -
 *  Created 7.08.08 by Sean R. Stratton
 *
 *  Returns a random number from 0 to inf with Poisson distribution of mean 
 *  'nu'.  May put in a safety later in case this function is called with an 
 *  unreasonably large value for 'nu'.
 */
int prand(double nu) {
  int k=0;
  double u;

  u = RANDOM_NUMBER;
  while ( u >= exp( -nu ) ) {
    k++;
    u *= RANDOM_NUMBER;
  }

  return k;
}

/*
 *  double grand(void) - 
 *  Created on 11.17.04 by Sean R Stratton - Rutgers Dept. of Physics & Astron.
 *
 *  Returns a random # on -inf to inf with Gauss Dist. mean=0 and sigma=1.
 *  Method comes from Physics Review 'D' p240 under "Monte Carlo Techniques."
 */
double grand(void) {

  double rr, v1, v2, gx;

  do {
    v1 = 2.00*RANDOM_NUMBER - 1.00;
    v2 = 2.00*RANDOM_NUMBER - 1.00;

    rr = v1*v1 + v2*v2;
  }  while ( rr > 1.00 );

  gx = v1 * sqrt( -2.00*log(rr) / rr );

  return gx;
}

/*
 *  double trand(void) -
 *  Created 7.08.08 by Sean R. Stratton.
 *
 *  Returns a random number from -1 to 1 using a triangular distribution.  
 *  This is useful for scattering PE as seen by a PMT.
 */
double trand(void) {
  double u;

  /*
   * note:  this is NOT the same as: u = RANDOM_NUMBER + RANDOM_NUMBER
   */
  u = 2.0 * RANDOM_NUMBER;

  if ( u < 1.00 )
    return sqrt( u ) - 1.0;
  else
    return 1.00 - sqrt( u - 1.0 );

}
