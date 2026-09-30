#include <math.h>
#include "toolbox.h"

void linefit (int num, double x[], double y[], double ey[], 
              double *m, double *em, double *b, double *eb, double *chi2) {

  int i;

  double s   = 0.;
  double sx  = 0.;
  double sy  = 0.;
  double sxx = 0.;
  double sxy = 0.;

  double delta;

  double r;

  for (i=0; i<num; i++) {
    s   += 1.        / (ey[i]*ey[i]);
    sx  += x[i]      / (ey[i]*ey[i]);
    sy  += y[i]      / (ey[i]*ey[i]);
    sxx += x[i]*x[i] / (ey[i]*ey[i]);
    sxy += x[i]*y[i] / (ey[i]*ey[i]);
  }

  delta = s * sxx - (sx * sx);

  *b = (sxx * sy - sx * sxy) / delta;
  *m = (s * sxy  - sx *  sy) / delta;

  *eb = sqrt(sxx / delta);
  *em = sqrt(s / delta);

  *chi2 = 0.;
  for (i=0; i<num; i++) {
    r = y[i] - (*m * x[i] + *b);
    *chi2 += (r / ey[i]) * (r / ey[i]);
  }
}
