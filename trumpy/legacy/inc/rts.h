#ifndef _RTS_H_
#define _RTS_H_

typedef struct {
  unsigned int itry      : 24;  // trial number
  unsigned int triglev   :  2;  // 2=successful trigger, 1=passed E/Rp cut, 0=failed
  unsigned int configID  :  2;  // configuration index (0=master, 1=slave1, 2=slave2, etc)
  unsigned int pspec     :  4;  // particle species, 1=proton, 2=iron, etc.

  unsigned int rtday     : 14;  // The date the simulation was actually run is indicated
  unsigned int rtmonth   :  6;  //   by the prefix 'rt' (for run time).  'sec' is in 
  unsigned int rtyear    : 12;  //   multiples of 100, corresponding to the CPU clock

  unsigned int rtsec100  : 14;
  unsigned int rtminute  :  6;
  unsigned int rthour    : 12;

/*   unsigned int ttsec100  : 14;  // Accumulated simulation run time indicated by the  */
/*   unsigned int ttminute  :  6;  //   prefix 'tt'. */
/*   unsigned int tthour    : 12; */

  unsigned int stday     : 14;  // Simulation time used for deriving calibration data
  unsigned int stmonth   :  6;  //   indicated by the prefix 'st'.  'sec' is still in 
  unsigned int styear    : 12;  //   multiples of 100 for name-space consistency, but
  unsigned int stsec100  : 14;  //   stored values are all multiples of 100.  Fraction
  unsigned int stminute  :  6;  //   of a second is expressed entirely in 'stnano'
  unsigned int sthour    : 12;

  unsigned int stnano    : 32;

  int          loge1E4   : 32;  // 10^4*Log of thrown energy in log_10(E/eV)
  int          xcore     : 32;  // x-coordinate of core (cm, CLF, projected to z=0)
  int          ycore     : 32;  // y-coordinate of core (cm, CLF, projected to z=0)

  int          logn1E4   : 32;  // 10^4*Log of MC Nmax
  int          xmax1E4   : 32;  // 10^4*Xmax in g/cm^2
  int          lambda10  : 16;  // 10*profile width parameter in g/cm^2
  int          x010      : 16;  // 10*effective depth of 1st interaction in g/cm^2

  int          rpcm      : 32;  // MC Rp (cm)
  int          psirad1E4 : 32;  // 10^4*MC Psi (radians)
  int          zenrad1E4 : 32;  // 10^4*MC zenith angle of shower axis (radians)
  int          azmrad1E4 : 32;  // 10^4*MC azimuth angle of shower axis (radians)
} RuntimeSpecsBitField;

#endif
