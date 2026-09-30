#ifndef _TL_CALIBRATION_H_
#define _TL_CALIBRATION_H_



typedef struct {
  double xmu_noise;
//   double tau_
  double ped_fadc[GEOTL_MAXMIR][nchan_mir];
  double gain_fadc[GEOTL_MAXMIR][nchan_mir];

//   int db_prescale;
  int db_deadmir;
  int db_dead[GEOTL_MAXMIR];
  
  double calib_corr;
} TLCalibration;

#endif

void getTLDefaultCalibration(RuntimeParameters *par, TLCalibration *tes);
