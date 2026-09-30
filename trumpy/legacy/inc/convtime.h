/*
 *  convtime.h
 *
 *  Author:  Sean R. Stratton
 *           Rutgers University, dept. of Physics & Astronomy
 *
 *  This file contains prototypes for the set of functions that convert dates 
 *    and times between Gregorian, Julian, Modified Julian, and "Computer" 
 *    dates.  At least in the scope of this file, the following commonly used 
 *    variables carry the meanings described in this table:
 *
 *                |          |
 *      variable  |  type    |  meaning
 *                |          |
 *    ------------+----------+------------------------------------------------
 *                |          |
 *       julday   |  int     |  Number of days passed since noon, 
 *                |          |    1 Jan 4713 BCE
 *                |          |
 *    ------------+----------+------------------------------------------------
 *                |          |
 *       mjlday   |  double  |  Time elapsed during current Julian 
 *                |          |    period in days.
 *                |          |
 *    ------------+----------+------------------------------------------------
 *                |          |
 *                |  int     |  (Pair of arguments), the Gregorian calendar 
 *       ymdsec   |          |    year as integer in YYYYMMDD format.
 *                |  double  |  Time elapsed into day in seconds.
 *                |          |
 *    ------------+----------+------------------------------------------------
 *                |          |
 *                |          |  Time elapsed since the Epoch (midnight, 
 *       eposec   |  double  |    1 Jan. 1970) in seconds, commonly used 
 *                |          |    by computers.
 *                |          |
 *
 *  The functions are explained in greater detail in 'convtime.c', where they 
 *    are defined.
 *
 *
 *  Revisions:
 *
 *     5.05.2009  :  Initial revision.
 *
 */
#ifndef _CONVTIME_H_
#define _CONVTIME_H_

#define SECPDAY  8.6400000e+04  /* Number of seconds in one day */
#define NSECPDAY 8.6400000e+13  /* Number of nanoseconds in one day */
#define MJLDOFF  2.4000005e+06  /* Offset between Modified and Real Julian dates */
#define JDAYEPO  2.4405875e+06  /* Julian day of the Epoch (midnight, 1 Jan 1970) */

/*
 *  NOTE: All routines use UT standard
 */
#ifdef __cplusplus
extern "C"{
#endif
void   julday2caldat(int julday, int *year, int *month, int *day);
int    caldat2julday(int yr, int mo, int da);

void   mjlday2ymdsec(double mjlday, int *ymd, double *sec);
double ymdsec2mjlday(int ymd, double sec);

double mjlday2eposec(double mjlday);
double eposec2mjlday(double eposec);

void   eposec2ymdsec(double eposec, int *ymd, double *sec);
double ymdsec2eposec(int ymd, double sec);
#ifdef __cplusplus
}
#endif
/*
 *  Macros for compatability with older time-conversion functions.
 */
#define sec_to_jday(T)        (eposec2mjlday((double)(T))+MJLDOFF)
#define ymd_to_jday(Y,M,D,S)  (ymdsec2mjlday((10000*(Y)+100*(M)+(D)),(S))+MJLDOFF)

#endif
