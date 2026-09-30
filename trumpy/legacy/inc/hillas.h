#ifndef _HILLAS_H_
#define _HILLAS_H_

/* in "hillas.c" */
double getHillasE0(double s);
double getHillasT(double energy, double s);
double getHillasdTdE(double energy, double s);
double getHillasIntegrand (double energy, double s, double delta);

#endif
