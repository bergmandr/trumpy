#ifndef _TAELECTRONICS_H_
#define _TAELECTRONICS_H_

#define MIN_READOUT_CL  3.00

#define TAU   56.00   // gives a power spectrum contour closer to what we 
#define ALPHA  5.00   //   see in the data

typedef struct {
  int it;
  double va[5];
  double vt[2];
} WaveformState;

/* in "electronics.c" */
#ifdef __cplusplus
extern "C"{
#endif
void initTAElectronics(void);
int initWaveformGenerator(double tau, double alpha);
int generateWaveform(double *npe, double *tpe, double tref, int len,
                     int ped,
                     double nsbg,
                     double gain,
                     int nfadc, short *fadc, WaveformState *st);
int simTATriggerResponse(TACalibration *cb, PETimes *pt, fdraw_dst_common *fdraw);

int initSimpleWaveformGenerator(double tau);
int generateSimpleWaveform(double *npe, double *tpe, double tref, int len,
                     int ped,
                     double nsbg,
                     double gain,
                     int nfadc, short *fadc, WaveformState *st);

void freeFADCBuffer(void);                     
                     
/* in "fadcSig.c" */
double fadcSig(int mean, int var, int nwf, int *wf);
#ifdef __cplusplus
}
#endif
#endif
