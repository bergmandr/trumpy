#include "random.h"

#define DENOM 2147483648.0
#define SPLIT      65536
#define A          31513
#define C       73773773

/*
 * double ranss(int *, int *)
 *
 * Author: Sean R. Stratton
 *         Rutgers University, Dept. of Physics & Astronomy
 *
 * Generates a pseudo-random number in the range [ 0.0, 1.0 ).  The 
 *   progression will repeat with a period of T = 2^62 = 4611686018427387904 
 *   iterations.  The precision will be no worse than 2^(-53) (the limit of 
 *   double-precision floating point numbers relative to 1.00) and no better 
 *   than 2^(-62).
 *
 * Inputs:
 *   int *hiseed - The most significant 31 bits of the pseudo-random number.
 *   int *loseed - The least significant 31 bits of the pseudo-random number.
 *
 * Outputs:
 *   Both seed values are updated according to the progression:
 *
 *                 X_{i+1} = ( A X_i + C ) mod( 2^62 )
 *
 *     with 'hiseed' set to the high 31 bits and 'loseed' to the low 31 bits.
 *
 * Returns:
 *   A pseudo-random number evenly distributed in the range [ 0.0, 1.0 ).  
 *     That is, return value can be equal to zero but will never equal one.
 *
 */
double ranss(int *hiseed, int *loseed) {
  int hi, lo, over, temp;
  double u;

  if ( *loseed < 0 )
    *loseed = -(*loseed);

  if ( *hiseed < 0 )
    *hiseed = -(*hiseed);

  lo = *loseed % SPLIT;
  hi = *loseed / SPLIT;

  temp = A*lo + C;
  lo   = temp % SPLIT;
  over = temp / SPLIT;

  temp = A*hi + over;
  hi   = temp % (SPLIT/2);
  over = temp / (SPLIT/2);

  *loseed = lo + SPLIT*hi;

  u = (double)(*loseed)/DENOM;

  lo = *hiseed % SPLIT;
  hi = *hiseed / SPLIT;

  temp = A*lo + over;
  lo   = temp % SPLIT;
  over = temp / SPLIT;

  temp = A*hi + over;
  hi   = temp % (SPLIT/2);

  *hiseed = lo + SPLIT*hi;

  return ( (double)(*hiseed) + u )/DENOM;
}

#undef DENOM
#undef SPLIT
#undef A
#undef C
