#ifndef _BANK_MANAGER_H_
#define _BANK_MANAGER_H_

void compactFDRawBank(fdraw_dst_common* fdraw);
void clearTrumpMCBank(trumpmc_dst_common *tmc);
void buildTrumpMCBank(int siteid, const AirShower* as,
		      const geofd_dst_common* fdsg, 
		      const GaisserHillasParameters* gh, const Track* t, 
		      const ObservedTrack* ot, const PETimes* pet,  
		      const TACalibration* calib, const RuntimeParameters* par, trumpmc_dst_common* tmc);


#endif
