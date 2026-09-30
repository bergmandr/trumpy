#include <stdio.h>

/**
 *  heapsortp - fill the pointer array 'p' with the addresses of the elements
 *    of array 'x' in ascending order using the heap sort algorithm.  The order
 *    of the entries in 'x' is left unchanged.
 *
 *  author: Sean R. Stratton, Rutgers Univ. dept of Physics & Astronomy
 *
 *  arguments:
 *    n    - number of elements in arrays 'x' and 'p'
 *    x[]  - array of numbers to be sorted
 *    *p[] - array of pointers to the entries in x in ascending order
 *
 *  returns
 *    void
 */
void heapsortp(int n, double x[], double *p[]) {
  int i, j, k;
  double *pa;

  // organize the heap (hiring phase)
  for ( i=0; i<n; i++ ) {
    pa = x + i;
    j = i;

    while ( j > 0 ) {
      k = (j-1) >> 1;
      if ( *pa < *p[k] ) break;
      p[j] = p[k];
      j = k;
    }
    p[j] = pa;
  }

  // sort the heap (promotion phase)
  for ( i=i-1; i>=0; i-- ) {
    pa = p[i];    // save the top element
    p[i] = p[0];  // retire the CEO into its place

    j = 0;
    while ( j < i ) {
      k = (j<<1) + 1;
      if ( k >= i ) break;  // there are no candidates
      if ( (k+1) < i )
	k = ( *p[k+1] > *p[k] ) ? k+1 : k;
      p[j] = p[k];
      j = k;
    }

    while ( j > 0 ) {
      k = (j-1) >> 1;
      if ( *pa < *p[k] ) break;
      p[j] = p[k];
      j = k;
    }
    p[j] = pa;
  }
}

/**
 *  heapsort - reorder the entries of array 'x' into ascending order using the 
 *    heap sort algorithm.
 * 
 *  author: Sean R. Stratton, Rutgers Univ. dept of Physics & Astronomy
 *
 *  arguments:
 *    n    - number of elements in 'x[]'
 *    x[]  - array of numbers to be sorted
 *
 *  returns
 *    void
 */
void heapsort(int n, double x[]) {
  int i, j, k;
  double xa;

  // organize the heap (hiring phase)
  for ( i=0; i<n; i++ ) {
    xa = x[i];
    j = i;

    while ( j > 0 ) {
      k = (j-1) >> 1;
      if ( xa < x[k] ) break;
      x[j] = x[k];
      j = k;
    }
    x[j] = xa;
  }

  // sort the heap (promotion phase)
  for ( i=i-1; i>=0; i-- ) {
    xa = x[i];    // save the top element
    x[i] = x[0];  // retire the CEO into its place

    j = 0;
    while ( j < i ) {
      k = (j<<1) + 1;
      if ( k >= i ) break;  // there are no candidates
      if ( (k+1) < i )
	k = ( x[k+1] > x[k] ) ? k+1 : k;
      x[j] = x[k];
      j = k;
    }

    while ( j > 0 ) {
      k = (j-1) >> 1;
      if ( xa < x[k] ) break;
      x[j] = x[k];
      j = k;
    }
    x[j] = xa;
  }
}
