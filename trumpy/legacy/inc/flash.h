#ifndef _FLASH_H_
#define _FLASH_H_

/* in "flash.c" */
//void getFluorescenceYield(int seanson, double h, double *yfl);
//void getFluorescenceYieldByMeter(int seanson, double age, double h, double *yfl);
#ifdef __cplusplus
extern "C"{
#endif
double kakimotoFYByMeV(double rho, double T);
double kakimotoFYByMeter(double rho, double T);
void getFluorescenceYield(fdatmos_param_dst_common *season, 
			  double h, double *yfl);
void getFluorescenceYieldByMeter(fdatmos_param_dst_common *season, 
				 double age, double h, double *yfl);

void getFluorescenceYieldByMeV(fdatmos_param_dst_common *season, 
			       double h, double *yfl);
#ifdef __cplusplus
}
#endif
#endif
