#ifndef _SCATTER_H_
#define _SCATTER_H_

/* in "scatter.c" */
//void getAttenuation(double ha, double hb, double dl, double *att);
//void getAttenuation(int season, double ra[3], double rb[3], double *att);
void getAttenuation(fdatmos_param_dst_common *season,
		    double ra[3], double rb[3], double *att);
void getMieScatterFraction(double ha, double hb, double dl, double *frac);
void getMieScatterFractionJava(double ha, double hb, double dl, double *frac);
void getRayleighScatterFraction(double dx, double *frac);
void getOzoneAbsorbtionFraction(double ha, double hb, double dl, double *fr);

double mieScatterFunction(double q);
double rayleighScatterFunction(double q);

#endif
