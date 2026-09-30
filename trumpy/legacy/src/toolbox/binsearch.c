#include "toolbox.h"

/*
 *  Use binary search method to find the greatest value in the (sorted) array 
 *    'xi' less than 'x'.  Returns that value's index number in the array.
 *
 *  Input:
 *    double x   --- target.
 *    int 'j'    --- the number of elements in the array 'xi'
 *    double *xi --- general array of doubles.
 *
 *  Returns:
 *    Index number of the greatest value in the array less than 'x'.
 */
int binsearch(double x, int j, double *xi) {
  int i, k;

  i=0;
  while ( j > (i+1) ) {
    k = ( i + j ) / 2;
    if ( x < xi[k] )
      j = k;
    else
      i = k;
  }

  return i;
}
