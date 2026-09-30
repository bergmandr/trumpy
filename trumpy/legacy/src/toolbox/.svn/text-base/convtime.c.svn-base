/*
 *  convtime.c
 *
 *  Author:  Sean R. Stratton
 *           Rutgers University, dept. of Physics & Astronomy
 *
 *  This file contains definitions for the set of functions that convert dates 
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
 *
 *  Revisions:
 *
 *     5.05.2009  :  Initial revision.
 */
#include <stdio.h>
#include <math.h>
#include "convtime.h"

/*
 *  Converts a _real_ Julian day (number of days since noon, 1 Jan 4713 BCE) 
 *    into the Gregorian calendar year, month, and day.  Note, the year month 
 *    and day represent the date of the beginning of the new Julian day.  That 
 *    is, a Julian day of '1' corresponds to noon, 2 Jan 4713 BCE.  The 
 *    comments are instructions given for computing the Gregorian calendar 
 *    date from a Julian day found on Wikipedia.
 *
 *  Input:
 *    int julday --- integer number of days since noon, 1 Jan 4713 BCE.
 *
 *  Output:
 *    int *year --- Gregorian year where the given Julian day begins.
 *    int *month --- Gregorian month where the given Julian day begins.
 *    int *day --- Gregorian day where the given Julian day begins.
 *
 *  Returns:
 *    void
 */
void julday2caldat(int julday, int *year, int *month, int *day) {
  int j, g, dg, c, dc, b, db, a, da, y, m, d;

  /*
   *  From J, compute a relative Julian day number j from a Gregorian epoch 
   *    starting on March 1, -4800 (i.e. March 1, 4801 BCE in the proleptic 
   *    Gregorian Calendar), the beginning of the Gregorian quadricentennial 
   *    32,044 days before the epoch of the Julian Period.
   */
  j = julday + 32044;

  /*
   *  From j, compute the number g of Gregorian quadircentennial cycles 
   *    elapsed (there are exactly 146,097 days per cycle) since the epoch; 
   *    subtract the days for this number of cycles, it leaved dg days since 
   *    the beginning of the current cycle.
   */
  g = j / 146097;
  dg = j % 146097;

  /*
   *  From dg, compute the number c (from 0 to 4) of Gregorian centennial 
   *    cycles (there are exactly 36,524 days per Gregorian centennial cycle) 
   *    elapsed since the beginning of the current Gregorian quadricentennial 
   *    cycle, number reduced to a maximum of 3 (this reduction occurs for the 
   *    last day of a leap centennial year where c would be 4 if it were not 
   *    reduced); subtract the number of days for this number of Gregorian 
   *    centennial cycles, it leaves dc days since the beginning of a 
   *    Gregorian century.
   */
  c = ( dg/36524 + 1 ) * 3 / 4;
  dc = dg - c*36524;

  /*
   *  From dc, compute the number b (from 0 to 24) of Julian quadrennial 
   *    cycles (there are exactly 1,461 days in 4 years, except for the last 
   *    cycle which may be incomplete by 1 day) since the beginning of the 
   *    Gregorian century; subtract the number of days for this number of 
   *    Julian cycles, it leaves db days in the Gregorian century.
   */
  b = dc / 1461;
  db = dc % 1461;

  /*
   *  From db, compute the number a (from 0 to 4) of Roman annual cycles 
   *    (there are exactly 365 days per Roman annual cycle) since the 
   *    beginning of the Julian quadrennial cycle, number reduced to a maximum 
   *    of 3 (this reduction occurs for the leap day, if any, where a would be 
   *    4 if it were not reduced); subtract the number of days for this number 
   *    of annual cycles, it leaves da days in the Julian year (that begins on 
   *    March 1).
   */
  a = ( db/365 + 1 ) * 3 / 4;
  da = db - a*365;

  /*
   *  Convert the four components g, c, b, a into the number y of years since 
   *    the epoch, by summing their values weighted by the number of years 
   *    that each component represents (respectively 400 years, 100 years, 4 
   *    years, and 1 year).
   */
  y = g*400 + c*100 + b*4 + a;

  /*
   *  With da, compute the number m (from 0 to 11) of months since March 
   *    (there are exactly 153 days per 5-month cycle; however, these 5-month 
   *    cycles are offset by 2 months within the year, i.e. the cycles start 
   *    in May, and so the year starts with an initial fixed number of days on 
   *    March 1, the month can be computed from this cycle by a Euclidian 
   *    division by 5); subtract the number of days for this number of months 
   *    (using the formula above), it leaves d days past since the beginning 
   *    of the month.
   */
  m = ( da*5 + 308 ) / 153 - 2;
  d = da - (m+4)*153/5 + 122;

  /*
   *  The Gregorian date (Y, M, D) can then be deduced by simple shifts from 
   *    (y, m, d).
   */
  *year  = y - 4800 + (m+2)/12;
  *month = (m+2)%12 + 1;
  *day   = d + 1;
  return;
}

/*
 *  This is _real_ Julian day (from noon, 1 Jan 4713 BCE)
 *
 *  Input:
 *    int yr --- Gregorian calendar year.
 *    int mo --- Gregorian calendar month.
 *    int da --- Gregorian calendar day.
 *
 *  Returns:
 *    int --- The number of days that have passes since 1 Jan 4713 BCE at noon.
 */
int caldat2julday(int yr, int mo, int da) {
  int a, y, m, julday;
 
  a = ( 14 - mo ) / 12;
  y = yr + 4800 - a;
  m = mo + 12*a - 3;
  julday = da + (153*m+2)/5 + 365*y + y/4 - y/100 + y/400 - 32045;
  return julday;
}

/*
 *  Converts a Modified Julian day to a YMD and seconds into the day.
 *
 *  Input:
 *    double mjlday --- The exact number of days (as a decimal fraction) 
 *          passed during the current Julian period, with 1/2-day offset.
 *
 *  Output:
 *    int *ymd --- The Gregorian calendar year, month, and day expressed in 
 *          integer form as YYYYMMDD.
 *    double *sec --- The number of seconds into the Gregorian calendar day.
 *
 *  Returns:
 *    void
 */
void mjlday2ymdsec(double mjlday, int *ymd, double *sec) {
  int yr, mo, da, julday;
  double imjlday, fmjlday;

  imjlday = floor(mjlday);
  fmjlday = mjlday - imjlday;
  *sec = SECPDAY * fmjlday;

  /* convert to real julian date */
  julday = (int)floor( imjlday + MJLDOFF + 0.50 );
  julday2caldat(julday, &yr, &mo, &da);
  *ymd = 10000*yr + 100*mo + da;

  return;
}

/*
 *  Converts a YMD and seconds into the day into a Modified Julian day.
 *
 *  Input:
 *    int ymd --- The Gregorian calendar date expressed in integer form 
 *          as YYYYMMDD.
 *    double sec --- The decimal fraction of seconds passed since midnight.
 *
 *  Returns:
 *    double --- The exact number of days (as a decimal fraction) that have 
 *          passed during the current Julian period, with 1/2-day offset.
 */
double ymdsec2mjlday(int ymd, double sec) {
  int yr, mo, da;
  double mjlday;

  sec = sec/SECPDAY - 0.500;

  yr = ymd / 10000;
  mo = ( ymd / 100 ) % 100;
  da = ymd % 100;

  mjlday = (double)caldat2julday(yr, mo, da) - MJLDOFF;

  return mjlday + sec;
}

/*
 *  Converts a Modified Julian day into the number of seconds passed since 
 *    the Epoch (midnight, 1 Jan. 1970), commonly used by computers.
 *
 *  Input:
 *    double mjlday --- Modified Julian date.
 *
 *  Returns:
 *    double --- The time elapsed since the Epoch (midnight, 1 Jan. 1970) in 
 *          seconds.
 */
double mjlday2eposec(double mjlday) {
  return SECPDAY * ( mjlday + MJLDOFF - JDAYEPO );
}

/*
 *  Converts the number of seconds since the Epoch into a Modified Julian day.
 *
 *  Input:
 *    double eposec --- Time elapsed since midnight, 1 Jan. 1970 in seconds.
 *
 *  Returns:
 *    double --- Modified Julian date corresponding to the given time.
 */
double eposec2mjlday(double eposec) {
  return eposec/SECPDAY + JDAYEPO - MJLDOFF;
}

/*
 *  Converts the number of seconds since the Epoch into a YMD and number of 
 *    seconds from midnight.
 *
 *  Input:
 *    double eposec --- Time elapsed since midnight, 1 Jan. 1970 in seconds.
 *
 *  Output:
 *    int *ymd --- The Gregorian calendar year, month, and day, expressed as 
 *          an integer in YYYYMMDD format.
 *    double *sec --- The time elapsed since midnight in seconds.
 *
 *  Returns:
 *    void
 */
void eposec2ymdsec(double eposec, int *ymd, double *sec) {
  int yr, mo, da, julday;
  double epoday, iepoday, fepoday;

  epoday = eposec / SECPDAY;
  iepoday = floor(epoday);
  fepoday = epoday - iepoday;
  *sec = SECPDAY * fepoday;

  julday = (int)floor( iepoday + JDAYEPO + 0.5 );
  julday2caldat(julday, &yr, &mo, &da);
  *ymd = 10000*yr + 100*mo + da;

  return;
}

/*
 *  Converts a YMD and seconds from midnight into the number of seconds from 
 *    the Epoch.
 *
 *  Input:
 *    int ymd --- The Gregorian calendar date expressed as an integer in 
 *          YYYYMMDD format.
 *
 *  Returns:
 *    double --- The time elapsed since midnight, 1 Jan. 1970, in seconds.
 */
double ymdsec2eposec(int ymd, double sec) {
  int yr, mo, da;
  double epoday;

  yr = ymd / 10000;
  mo = ( ymd / 100 ) % 100;
  da = ymd % 100;

  epoday = (double)caldat2julday(yr, mo, da) - JDAYEPO - 0.50;
  return SECPDAY*epoday + sec;
}


/*
 *  To compile in test mode, run the command:
 *
 *     gcc -o convtime.run convtime.c -lm -I../../inc -DCONVTIME_TEST_MODE
 */
#ifdef CONVTIME_TEST_MODE
#include <stdlib.h>

/*
 *  Tests each of the function defined above for consistency.  When run w/o 
 *    arguments (or the wrong number of arguments), the program will use Sun, 
 *    14 Jan 2007 @ 13:18:59.9 as reference for testing (JDay=2454115.05486).  
 *    Otherwise, it accepts the arguments: ${ymdday} ${hr} ${mn} ${sec} 
 *    (in that order), where ${ymdday} is the 8-digit calendar date in 
 *    YYYYMMDD format, ${hr} is the UT hour, ${mn} is the UT minute, and 
 *    ${sec} is time elapsed into the minute in seconds (is a decimal 
 *    fraction).
 */
int main(int argc, char *argv[]) {
  int ymd = 20070114;
  int yr, mo, da;
  double hr = 13;
  double mn = 18;
  double sec = 59.9;
  double mjlday = 2454115.05486 - 2400000.5;
  double eposec;

  if ( argc == 5 ) {
    ymd = atoi(argv[1]);
    hr = atof(argv[2]);
    mn = atof(argv[3]);
    sec = atof(argv[4]);
  }

  yr = ymd / 10000;
  mo = ( ymd / 100 ) % 100;
  da = ymd % 100;

  printf("input: %2d/%2d/%4d  %2.0lf:%02.0lf:%012.9lf\n\n", 
	 mo, da, yr, hr, mn, sec);

  sec += 3600.0*hr + 60.0*mn;

  int julday = caldat2julday(yr, mo, da);
  printf("julday = %d\n", julday);
  julday2caldat(julday, &yr, &mo, &da);
  printf("    to caldat: %2d/%2d/%4d\n\n", mo, da, yr);

  mjlday = ymdsec2mjlday(ymd, sec);
  printf("mjlday = %lf\n", mjlday);
  mjlday2ymdsec(mjlday, &ymd, &sec);
  yr = ymd / 10000;
  mo = ( ymd / 100 ) % 100;
  da = ymd %100;
  hr = floor(sec/3600.0);
  mn = floor( ( sec - 3600.0*hr ) / 60.0 );
  sec = sec - 3600.0*hr - 60.0*mn;
  printf("    to date, time: %2d/%2d/%4d %2.0lf:%02.0lf:%012.9lf\n\n", 
	 mo, da, yr, hr, mn, sec);

  printf("eposec = %lf\n", eposec=mjlday2eposec(mjlday));
  printf("    to mjlday: %lf\n\n", mjlday=eposec2mjlday(eposec));

  sec += 3600.0*hr + 60.0*mn;
  eposec = ymdsec2eposec(ymd, sec);
  printf("eposec = %lf\n", eposec);
  eposec2ymdsec(eposec, &ymd, &sec);
  yr = ymd / 10000;
  mo = ( ymd / 100 ) % 100;
  da = ymd %100;
  hr = floor(sec/3600.0);
  mn = floor( ( sec - 3600.0*hr ) / 60.0 );
  sec = sec - 3600.0*hr - 60.0*mn;
  printf("    to date, time: %2d/%2d/%4d %2.0lf:%02.0lf:%012.9lf\n\n", 
	 mo, da, yr, hr, mn, sec);  

  return 0;
}

#endif
