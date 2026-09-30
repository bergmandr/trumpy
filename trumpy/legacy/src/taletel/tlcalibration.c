#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "control.h"
#include "event.h"
#include "fdconstants.h"
#include "tlcalibration.h"

void getTLDefaultCalibration(RuntimeParameters *par, TLCalibration *tes) {
  int mir, j;
  double gain_adc,ped_adc,ped_16,gain_16;
  
  gain_adc = 1.0;
  ped_adc = 10.50;
  ped_16 = 20.50;
  gain_16 = 1.89;
  
  for (mir=0; mir<GEOTL_MAXMIR; mir++) {
//     gain_16 = (geo->ring[mir]==2)?1.50:1.89;
    
    for (j=0; j<GEOFD_MIRTUBE; j++) {
      tes->ped_fadc[mir][j] = ped_adc;
      tes->gain_fadc[mir][j] = gain_adc;
    }
    for (j=0; j<16; j++) {
      tes->ped_fadc[mir][256+j] = ped_16;
      tes->ped_fadc[mir][288+j] = ped_16;
      tes->gain_fadc[mir][256+j] = gain_16;
      tes->gain_fadc[mir][288+j] = gain_16;
      tes->ped_fadc[mir][272+j] = ped_adc;
      tes->ped_fadc[mir][304+j] = ped_adc;
      tes->gain_fadc[mir][272+j] = gain_adc/9.8;
      tes->gain_fadc[mir][304+j] = gain_adc/9.8;
    }
  }
     
  tes->calib_corr = 1.0;
  return;
}