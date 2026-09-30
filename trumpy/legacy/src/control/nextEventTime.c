#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "constants.h"
#include "control.h"
#include "event.h"
#include "random.h"
#include "toolbox.h"
#include "fdconstants.h"

#include "airshower.h"
#include "showerlib.h"
#include "eventlist.h"

int loadOnTimes(int nc, RuntimeParameters *par[]) {
  int i, n, rc, s=100, e;
  int jd, js, jn;
  double dt, diff;
  fdped_dst_common *ped;

  if ( par[0]->dt > 0.00 ) {
    dt = par[0]->dt;
  }
  else {
    perr("no dt specified.  Using default (%.1f sec)", DEFAULT_DT);
    dt = DEFAULT_DT;
  }

  pout("loading on-times from fdped files:\n");
  
  int wb = newBankList(s);
  int gb = newBankList(s);

  jd = 9999999;
  js = (int)SECPDAY;
  jn = 999999999;

  for ( i=0; i<nc; i++ ) {
    pout("  %s\n", par[i]->ontime);

//     if (      par[i]->siteid == BLACK_ROCK_SITEID)
//       ped = &brped_;
//     else if ( par[i]->siteid == LONG_RIDGE_SITEID)
//       ped = &lrped_;
//     else {
//       perr("invalid site ID (%d)\n", par[i]->siteid);
//       return -1;
//     }
    switch (par[i]->siteid) {
      case BLACK_ROCK_SITEID:
        ped = &brped_;
        break;
      case LONG_RIDGE_SITEID:
        ped = &lrped_;
        break;
      default:
        perr("site ID (%d) not compatible with fdped\n", par[i]->siteid);
        return -1;
    }
    par[i]->dt = dt;

    rc = dstOpenUnit(ONTIME_FILEDES, par[i]->ontime, MODE_READ_DST);
    if ( rc != 0 ) {
      perr("unable to open '%s' for reading (%d).", par[i]->ontime, rc);
      continue;
    }

    n = 0;
    while ( eventRead(ONTIME_FILEDES, wb, gb, &e) >= 0 ) {
      par[i]->part[n]   = ped->part;

      par[i]->t0day[n]  = ped->julian_start;
      par[i]->t0sec[n]  = ped->jsecond_start;
      par[i]->t0nsec[n] = ped->jsecfrac_start;
      
      par[i]->t1day[n]  = ped->julian_end;
      par[i]->t1sec[n]  = ped->jsecond_end;
      par[i]->t1nsec[n] = ped->jsecfrac_end;

      diff  = (double)(jd - par[i]->t0day[n]) * SECPDAY;
      diff += (double)(js - par[i]->t0sec[n]);
      diff += (double)(jn - par[i]->t0nsec[n]) / 1.0e9;

      if ( diff > 0.) {
      jd = par[i]->t0day[n];
      js = par[i]->t0sec[n];
      jn = par[i]->t0nsec[n];
      }

      n++;
    }

    par[i]->numparts = n;

    dstCloseUnit(ONTIME_FILEDES);
  }

  delBankList(wb);
  delBankList(gb);

  /*
   *  Reconcile runtime parameters with 'jd'
   */
  par[0]->jday  = jd;
  par[0]->jsec  = js;
  par[0]->nsec  = jn;

  /* Get calendar date from julian */
  jul2cal(par[0]);

  for ( i=0; i<nc; i++ ) { 
    par[i]->jday  = par[0]->jday;
    par[i]->jsec  = par[0]->jsec;
    par[i]->year  = par[0]->year;
    par[i]->month = par[0]->month;
    par[i]->day   = par[0]->day;
    par[i]->hour  = par[0]->hour;
    par[i]->min   = par[0]->min;
    par[i]->sec   = par[0]->sec;
    par[i]->nsec  = par[0]->nsec;
  }

  return 0;
}

/* this function could and should be updated, but I'm leaving it alone for 
 *   the sake of ensuring the program will behave as expected when reverting 
 *   to single-detector mode. */
int nextEventTime(int nc, RuntimeParameters *p[], fdraw_dst_common *f[]) {
  static boolean ifirst = TRUE;
  static double jd, js, jn;
  static int part;

  boolean scan;
  int i;

  double diff0, diff1;

  if ( ! p[0]->flag_ontime ) {  /* use date from right now */
    time_t now = time(&now);

    eposec2jul((int)now, p[0]);
    jul2cal(p[0]);

    jd = (double)p[0]->jday;
    js = (double)p[0]->jsec;
    jn = (double)p[0]->nsec;

    p[0]->flag_mirref  = FALSE;
    p[0]->flag_pmtGain = FALSE;
    p[0]->flag_pmtCal  = FALSE;
  }
  else {
    if ( ifirst ) {
      jd = 1.0e+137;
      js = SECPDAY;
      jn = NSECPDAY;
      for ( i=0; i<nc; i++ ) {
      p[i]->ipart = 0;
      p[i]->jday = p[i]->t0day[0];
      p[i]->jsec = p[i]->t0sec[0];
      p[i]->nsec = p[i]->t0nsec[0];
      f[i]->part = (integer2)p[i]->part[0];
      f[i]->event_num = -1;

      diff0  = ( jd - (double)p[i]->jday );            // [sic]
      diff0 += ( js - (double)p[i]->jsec ) * SECPDAY;  // [sic]
      diff0 += ( jn - (double)p[i]->nsec ) / 1.0e9;

      if ( diff0 > 0. ) {
        jd = (double)p[i]->jday;
        js = (double)p[i]->jsec;
        jn = (double)p[i]->nsec;
      }
      }
    }

    do {  /* loop until I'm satisfied */
      jn += floor( -p[0]->dt * log(RANDOM_NUMBER) * 1.0e9 );
      while ( jn >= 1.0e+9 ) {
        jn -= 1.0e+9;
        js += 1.0;
      }
      while ( js >= SECPDAY ) {
        js -= SECPDAY;
        jd += 1.0;
      }

      scan = FALSE;
      for ( i=0; i<nc; i++ ) {
        diff0  = ( jd - (double)p[i]->t0day[p[i]->ipart] ) * SECPDAY;
        diff0 += ( js - (double)p[i]->t0sec[p[i]->ipart] );
        diff0 += ( jn - (double)p[i]->t0nsec[p[i]->ipart] ) / 1.0e9;

        diff1  = ( (double)p[i]->t1day[p[i]->ipart] - jd ) * SECPDAY;
        diff1 += ( (double)p[i]->t1sec[p[i]->ipart] - js );
        diff1 += ( (double)p[i]->t1nsec[p[i]->ipart] - jn ) / 1.0e9;

        if ( diff1 < 0. ) {
          /* part needs to advance */
          if ( ++p[i]->ipart >= p[i]->numparts )
            return -1;

          scan = TRUE;
        }
        else if ( diff0 < 0. ) {
          /* jn needs to advance */
          scan = TRUE;
        }
        else
          part = p[i]->part[p[i]->ipart];
      }
    } while ( scan );
  }

  /*
   *  Reconcile runtime parameters and fdraw banks with 'jd'
   */
  p[0]->jday  = (int)jd;
  p[0]->jsec  = (int)js;
  p[0]->nsec  = (int)jn;

  /* Get calendar date from julian */
  jul2cal(p[0]);

  for ( i=0; i<nc; i++ ) { 
    p[i]->jday  = p[0]->jday;
    p[i]->jsec  = p[0]->jsec;
    p[i]->year  = p[0]->year;
    p[i]->month = p[0]->month;
    p[i]->day   = p[0]->day;
    p[i]->hour  = p[0]->hour;
    p[i]->min   = p[0]->min;
    p[i]->sec   = p[0]->sec;
    p[i]->nsec  = p[0]->nsec;

    f[i]->part = p[i]->part[p[i]->ipart];
    f[i]->julian = p[i]->jday;
    f[i]->jsecond = p[i]->jsec;
    f[i]->ctdclock = p[i]->nsec / 25;
    f[i]->gps1pps_tick = 0;
  }

  ifirst = FALSE;
  return 0;
}

int nextTandemEventTime(int nc, RuntimeParameters *p[], fdraw_dst_common *f[]) {
  static boolean ifirst = TRUE;
  static double jd, js, jn;
  int i;
  double diff0, diff1;

  int flag_ontime = 0x1;
  for ( i=0; i<nc; i++ )
    flag_ontime &= p[i]->flag_ontime;

  if ( ! flag_ontime ) {  /* use date from right now */
    time_t now = time(&now);

    eposec2jul((int)now, p[0]);
    jul2cal(p[0]);

    jd = (double)p[0]->jday;
    js = (double)p[0]->jsec;
    jn = (double)p[0]->nsec;

    // there will certainly be no calibration data from this time
    p[0]->flag_mirref  = FALSE;
    p[0]->flag_pmtGain = FALSE;
    p[0]->flag_pmtCal  = FALSE;

    return 0xFF;
  }

  if ( ifirst ) {
    jd = 1.0e+137;
    js = SECPDAY;
    jn = NSECPDAY;
    for ( i=0; i<nc; i++ ) {
      p[i]->ipart = 0;
      p[i]->jday = p[i]->t0day[0];
      p[i]->jsec = p[i]->t0sec[0];
      p[i]->nsec = p[i]->t0nsec[0];
      f[i]->part = (integer2)p[i]->part[0];
      f[i]->event_num = -1;

      diff0  = ( jd - (double)p[i]->jday ) * SECPDAY;
      diff0 += ( js - (double)p[i]->jsec );
      diff0 += ( jn - (double)p[i]->nsec ) * 1.0e-9;

      if ( diff0 > 0. ) {
      jd = (double)p[i]->jday;
      js = (double)p[i]->jsec;
      jn = (double)p[i]->nsec;
      }
    }

    ifirst = FALSE;
  }

  // this version lets TRUMP main do the loop, because in tandem mode, part 
  //   times will be staggered in time.  Total on-time for a night will be the 
  //   difference between the end of the later part and the beginning of the 
  //   earliest part.  Dead time may need to be simulated as well.
  jn += floor( -p[0]->dt * log(RANDOM_NUMBER) * 1.0e9 );
  while ( jn >= 1.0e+9 ) {
    jn -= 1.0e+9;
    js += 1.0;
  }
  while ( js >= SECPDAY ) {
    js -= SECPDAY;
    jd += 1.0;
  }

  int site_mask = 0x0;
  boolean fullStop = TRUE;
  for ( i=0; i<nc; i++ ) {
    // make sure the run IDs are up to date
    diff1  = ( (double)p[i]->t1day[p[i]->ipart] - jd ) * SECPDAY;
    diff1 += ( (double)p[i]->t1sec[p[i]->ipart] - js );
    diff1 += ( (double)p[i]->t1nsec[p[i]->ipart] - jn ) / 1.0e9;
    while ( diff1 < 0. ) {
      if ( ++p[i]->ipart >= p[i]->numparts )
      break;

      diff1  = ( (double)p[i]->t1day[p[i]->ipart] - jd ) * SECPDAY;
      diff1 += ( (double)p[i]->t1sec[p[i]->ipart] - js );
      diff1 += ( (double)p[i]->t1nsec[p[i]->ipart] - jn ) / 1.0e9;
    }

    if ( p[i]->ipart >= p[i]->numparts )
      continue;
    else
      fullStop = FALSE;

    diff0  = ( jd - (double)p[i]->t0day[p[i]->ipart] ) * SECPDAY;
    diff0 += ( js - (double)p[i]->t0sec[p[i]->ipart] );
    diff0 += ( jn - (double)p[i]->t0nsec[p[i]->ipart] ) / 1.0e9;

    // diff1>0 same as saying part end time after sim time, 
    // diff2>0 same as saying sim time after part start time
    if ( (diff1>0.) && (diff0>=0.) )
      site_mask |= 0x1 << i;

  }

  if ( fullStop )
    return -1;

  /*
   *  Reconcile runtime parameters and fdraw banks with 'jd'
   */
  p[0]->jday  = (int)jd;
  p[0]->jsec  = (int)js;
  p[0]->nsec  = (int)jn;

  /* Get calendar date from julian */
  jul2cal(p[0]);

  for ( i=0; i<nc; i++ ) { 
    p[i]->jday  = p[0]->jday;
    p[i]->jsec  = p[0]->jsec;
    p[i]->year  = p[0]->year;
    p[i]->month = p[0]->month;
    p[i]->day   = p[0]->day;
    p[i]->hour  = p[0]->hour;
    p[i]->min   = p[0]->min;
    p[i]->sec   = p[0]->sec;
    p[i]->nsec  = p[0]->nsec;

    f[i]->part = p[i]->part[p[i]->ipart];
    f[i]->julian = p[i]->jday;
    f[i]->jsecond = p[i]->jsec;
    f[i]->ctdclock = p[i]->nsec / 25;
    f[i]->gps1pps_tick = 0;
  }

  return site_mask;
}

void fdraw2date(fdraw_dst_common *fdr, RuntimeParameters *par) {
  par->jday = fdr->julian;
  par->jsec = fdr->jsecond;
  par->nsec = (fdr->ctdclock - fdr->gps1pps_tick) * 25;

  jul2cal(par);
}

void fraw12date(fraw1_dst_common *fr1, RuntimeParameters *par) {
  par->jday = fr1->julian;
  par->jsec = fr1->jsecond;
  par->nsec = fr1->jclkcnt;
  
  jul2cal(par);
}

/* Get jday and jsec from year, month, day, hour, min, sec */
void cal2jul(RuntimeParameters *par) {
  /* ymd_to_jday returns real julian day (from noon) */
  par->jday  = (int)ymd_to_jday(par->year, par->month, par->day, 0.);
  par->jsec  = 43200 + 3600*par->hour + 60*par->min + par->sec;

  while (par->jsec >= 86400) {
    par->jsec -= 86400;
    par->jday++;
  }
}

/* Get year, month, day, hour, min, sec, from jday and jsec */
void jul2cal(RuntimeParameters *par) {
  par->hour = par->jsec / 3600 + 12;
  if (par->hour >= 24) {
    julday2caldat(par->jday+1, &par->year, &par->month, &par->day);   
    par->hour -= 24;
  }
  else
    julday2caldat(par->jday,   &par->year, &par->month, &par->day); 

  par->min = (par->jsec / 60) % 60;
  par->sec = par->jsec % 60;
}

void eposec2jul(int eposec, RuntimeParameters *par) {
  par->jday = eposec / 86400 + (int)JDAYEPO;
  par->jsec = 43200 + eposec % 86400;
  par->nsec = 0;

  while (par->jsec >= 86400) {
    par->jsec -= 86400;
    par->jday++;
  }
}
