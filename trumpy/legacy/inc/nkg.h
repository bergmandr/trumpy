#ifndef _NKG_H_
#define _NKG_H_

/* in "nkg.c" */
double getMoliereRadius (ATMP_t *season, double altitude);
void getNKGdr (double age, double molrad,
               double rpuv[], double npln[], double v[]);

#endif
