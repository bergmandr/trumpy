#ifndef _EVENTLIST_H_
#define _EVENTLIST_H_

#define MAXEVENTS       1000
#define NUMPARAMS       18
#define MAXLENGTH       1024


#ifdef __cplusplus
extern "C"{
#endif

int getNextEvent (int nc, FILE *infile, RuntimeParameters *par, AirShower *as, 
		  GaisserHillasParameters *gh, fdraw_dst_common *f[]);
int getEventListNum (char *infilename);
int isOnTimeEvent(RuntimeParameters *par);

/* from "nextEventTime.c" */
int loadOnTimes(int nc, RuntimeParameters *par[]);
int nextTandemEventTime(int nc, RuntimeParameters *p[], fdraw_dst_common *f[]);
int nextEventTime(int nc, RuntimeParameters *p[], fdraw_dst_common *f[]);
void fdraw2date(fdraw_dst_common *fdr, RuntimeParameters *par);
void fraw12date(fraw1_dst_common *fr1, RuntimeParameters *par);
void cal2jul(RuntimeParameters *par);
void jul2cal(RuntimeParameters *par);
void eposec2jul(int eposec, RuntimeParameters *par);
#ifdef __cplusplus
}
#endif

#endif
