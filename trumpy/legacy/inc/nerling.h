#ifndef _NERLING_H_
#define _NERLING_H_

/* in "nerling.c" */
#ifdef __cplusplus
extern "C"{
#endif
double getAlphaEff(double s);
double getDifferentialElectronSpectrum(double s, double e);
//void getCvPhotonIncrRate(int season, double s, double h, double dl, double *dncv);
//double cvPhaseFunction(int season, double s, double h, double q);
void getCvPhotonIncrRate(fdatmos_param_dst_common *season,
			 double s, double h, double *dncv);
double cvPhaseFunction(fdatmos_param_dst_common *season,
		       double s, double h, double q);
#ifdef __cplusplus
}
#endif

#endif
