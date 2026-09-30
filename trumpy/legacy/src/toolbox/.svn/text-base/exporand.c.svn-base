#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void) {
  int i;
  double rplo=10.00, rphi=40.0, rp;
  double u;

  rplo = rplo*rplo*rplo;
  rphi = rphi*rphi*rphi;

  for ( i=0; i<100000; i++ ) {
    u = drand48();

    rp = cbrt( rplo + u*(rphi-rplo) );

    printf("%e\n", rp);
  }

  return 0;
}
