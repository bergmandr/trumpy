#ifndef _CALIBRATION_H_
#define _CALIBRATION_H_


#include "tacalibration.h"
#include "tlcalibration.h"

typedef struct {
  int nmir;
  TACalibration *tc;   // traditional calibration structure
//   MDCalibraton *mc;    // new calibration structure for use with MD
//   TLCalibration *tlc;  // new calibration structure for use with TL
  TLCalibration *tl;
  
} UCalibration;


#endif
