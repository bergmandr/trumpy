#ifndef _WAVEFORM_GENERATOR_H_
#define _WAVEFORM_GENERATOR_H_

int initWFG(void);

int generateWaveform(double *npe, double *tpe, double tref, int len, int ped, 
		     double nsbg, double gain, int *fadc);

#endif
