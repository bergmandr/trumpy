/*
 *  atmosphere-withdb.c
 *
 *  Author:  Sean R. Stratton
 *           Rutgers University, dept. of Physics & Astronomy
 *
 *  This file defines the set of functions for accessing and interpolating 
 *    data from the atmospheric data base (DST).
 *
 *  Revisions:
 *     4.24.09  Base Revision
 *     4.28.09  Added Atmospheric transmission functions
 *     5.06.09  Did some code optimization...
 *                * removed extra call to _getLayerIndex() from _evaluatePrTp()
 *                * changed ATMOS_EPS from 10^-6 to 10^-4
 *                * replaced calls to R_to_latlonalt with r2lla
 *                      -- SS
 *
 *     6.29.10  Made a change in the getAtmosphere function below.  This caused the
 *                 code to seg fault if dates were not time ordered.  ELB
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "constants.h"
#include "control.h"
#include "airshower.h"
#include "event.h"
#include "atmosphere-withdb.h"
#include "toolbox.h"

/*
 *  Constants used by the hydrostatic formulae
 */
#define REARTH 6.369e+03  /* km */
#define G  9.80665  /* m/s^2 */
#define M  28.9644  /* g/mol */
#define R  8.31432  /* J/mol*K */
//#define GMR  (G*M/R)  /* for computing geopotential altitude */
#define GMR  .341631947363103657e+02


/*
 *  Constants used for defining the US standard atmosphere.
 */
#define NUM_US1976_LAYERS 5
#define P0 1.01325e+03  /* Standard sea level air pressure, (hPa) */
const double US1976_GEOPOTENTIAL[NUM_US1976_LAYERS] = { 
  0.00, 11.0, 20.0, 32.0, 47.0 
};
const double US1976_TEMPERATURE[NUM_US1976_LAYERS] = { 
  288.15, 216.65, 216.65, 228.65, 270.65 
};


/*
 *  variables for handling the atmos DST files
 */
static int lunatmp = -1;
static int pwb, phb, pev;
static char *pdbname;

static int lunatmt = -1;
static int twb, thb, tev;
static char *tdbname;

/*
 *  "fall-back" banks so functions can evaluate thermodynamic quantities for 
 *    altitudes beyond the limits of the data points in the DST bank.
 */
static ATMP_t _ufparam_;
static ATMP_t *ufparam = &_ufparam_;

static ATMP_t _ofparam_;
static ATMP_t *ofparam = &_ofparam_;


/**
 *
 *  Finds the corresponding layer index for the given height. (private)
 *
 *  Input:
 *    double h --- geopotential
 *    int j --- the number of entries in the array 'hb'
 *    double *hb --- array of base-level geopotential altitudes.
 *
 *  Returns:
 *    int --- the index number of the array 'hb' with the base geopotential 
 *          altitude corresponding to 'h'.
 *
 */
static int _getLayerIndex(double h, int j, float *hb) {
  int i, k;

  /*
   *  Use a binary search to find the index 'i' of array 'hb' corresponding to 
   *    the layer containing height 'h'.
   */
  i = 0;
  while ( j > (i+1) ) {
    k = ( i + j ) / 2;
    if ( h < hb[k] )
      j = k;
    else
      i = k;
  }

  return i;
}

/**
 *
 *  Computes the values of 'P' & 'T' by interpolating the data in the 'aparam' 
 *    bank with the hydrostatic formulae.  (private)
 *
 *  Input:
 *    fdatmos_param_t *aparam --- A DST Bank, likely a default-type version.
 *    double h --- The geopotential for which the computation should be made.
 *    int mode --- Tells the evaluator to use either linear extrapolation or 
 *          interpolate using the equations of hystrostatic equilibrium.
 *
 *  Output:
 *    double *p --- The address of the variable representing pressure (units 
 *          are automatically consistent with bank data).
 *    double *t --- The address of the variable representing temperature 
 *          (units are consistent with bank data).
 *
 *  Returns:
 *    void.
 */
static void _evaluatePrTp(ATMP_t *aparam, int idx, double h, double *p, double *t) {
  int idxp1;       /* =idx+1, done for efficiency */
  double hb, dh;   /* base height, altitude above base height */
  double tb, gt;   /* base temp., temp. gradient */
  double pb, lnp;  /* base pressure, ln(pressure/base pressure) */

  pb = aparam->pressure[idx];
  tb = aparam->temperature[idx];

  hb = aparam->height[idx];
  dh = h - hb;

  idxp1 = idx + 1;
  if ( idxp1 < aparam->nItem )
    gt = ( aparam->temperature[idxp1] - tb ) / ( aparam->height[idxp1] - hb );
  else
    gt = 0.00;

  if ( fabs(gt) < DOUBLE_FP_PRECISION )
    lnp = -GMR * dh / tb;
  else
    lnp = -GMR * log( (tb+gt*dh) / tb ) / gt;

  *p = pb * exp( lnp );
  *t = tb + gt*dh;
}

/**
 *
 *  Checks that bank 'j' agrees with bank 'i' in their evaluation of T & P at 
 *    the data point 'k' of bank 'i'.  (private)
 *
 *  Input:
 *    fdatmos_param_t *iparam --- A DST bank, likely filled with real data.
 *    int k --- The data point from bank 'i' to be used in the comparison.
 *    fdatmos_param_t *jparam --- A DST bank, likely an over- or under-flow 
 *          state.
 *
 *  Returns:
 *    int --- '0' if neither T nor P agree, '1' if only P agrees, '2' if only 
 *          T agrees, '3' if both T and P agree.
 *
 */
/* static int _compareBoundaryValues(ATMP_t *iparam, int k, ATMP_t *jparam) { */
/*   int istat; */
/*   double hanc, tanc, panc, tex, pex; */

/*   /\* */
/*    *  alias the data targeted for comparison */
/*    *\/ */
/*   panc = iparam->pressure[k]; */
/*   tanc = iparam->temperature[k]; */
/*   hanc = iparam->height[k]; */

/*   int idx = _getLayerIndex(hanc, jparam->nItem, jparam->height); */
/*   _evaluatePrTp(jparam, idx, hanc, &pex, &tex); */

/*   istat = 0; */
/*   if ( fabs(tanc-tex) < 1.0e-03*tex ) */
/*     istat += 1; */

/*   if ( fabs(panc-pex) < 1.0e-03*pex ) */
/*     istat += 2; */

/*   return istat; */
/* } */

/**
 *
 *  Sets the given bank to default PT values, using the US1976 standard 
 *    atmosphere model.  (public)
 *
 *  Input:
 *    none
 *
 *  Output:
 *    fdatmos_param_t *aparam --- A DST bank initialized to US standard 
 *          atmosphere.
 *
 *  Returns:
 *    void
 *
 */
void getDefaultParam(ATMP_t *aparam) {
  int i;
  double pb, hb, tb, dh, dt, gt, lnp;

  aparam->dateFrom = 0;
  aparam->dateTo = 0x7FFFFFFF;
  aparam->nItem = NUM_US1976_LAYERS;

  aparam->height[0] = (real4)US1976_GEOPOTENTIAL[0];
  aparam->pressure[0] = P0;  /* bank pressures are in units of hPa */
  aparam->pressureError[0] = 0.00;
  aparam->temperature[0] = (real4)US1976_TEMPERATURE[0];
  aparam->temperatureError[0] = 0.00;
  aparam->dewPoint[0] = 0.00;
  aparam->dewPointError[0] = 0.00;

  hb = aparam->height[0];
  pb = aparam->pressure[0];
  tb = aparam->temperature[0];
  for ( i=1; i<aparam->nItem; i++ ) {
    aparam->height[i] = (real4)US1976_GEOPOTENTIAL[i];
    aparam->temperature[i] = (real4)US1976_TEMPERATURE[i];

    dt = US1976_TEMPERATURE[i] - tb;
    dh = US1976_GEOPOTENTIAL[i] - hb;
    gt = dt / dh;

    if ( fabs(gt) < DOUBLE_FP_PRECISION )
      lnp = -GMR * dh / tb;
    else
      lnp = -GMR * log( US1976_TEMPERATURE[i] / tb ) / gt;

    aparam->pressure[i] = (real4)( pb * exp( lnp ) );

    aparam->pressureError[i] = 0.00;
    aparam->temperatureError[i] = 0.00;
    aparam->dewPoint[i] = 0.00;
    aparam->dewPointError[i] = 0.00;

    hb = aparam->height[i];
    pb = aparam->pressure[i];
    tb = aparam->temperature[i];
  }

  for ( i=aparam->nItem; i<FDATMOS_PARAM_MAXITEM; i++ ) {
    aparam->height[i] = 0.00;
    aparam->pressure[i] = 0.00;
    aparam->pressureError[i] = 0.00;
    aparam->temperature[i] = 0.00;
    aparam->temperatureError[i] = 0.00;
    aparam->dewPoint[i] = 0.00;
    aparam->dewPointError[i] = 0.00;
  }
}

/**
 *
 *  Sets the default values for the atrans bank using mean HiRes parameters.
 *
 *  Input:
 *    none.
 *
 *  Output:
 *    fdscat_dst_common *atrans --- DST bank containing the light transmission
 *          properties for the given night.
 *
 *  Returns:
 *    void.
 *
 */
void getDefaultTrans(ATMT_t *atrans) {

  /*
   *  set default Rayleigh & Mie scattering params
   */
  atrans->startTime = (double)0x0;
  atrans->endTime = (double)0x7FFFFFFF;

  atrans->hzalen = 25000.0;
  atrans->vaodep = 0.04;
  atrans->schght = 1000.0;

}

/**
 *
 *  Generates an 'out-of-bounds' structure to handle P & T evaluation at 
 *    altitudes that do not fall within the boundaries of the data points.
 *    (private)
 *
 *  Input:
 *    fdatmos_param_t *aparam --- A DST bank, likely real data.
 *    int idx --- The data point where the fall-back bank will be anchored.
 *
 *  Output:
 *    fdatmos_param_t *xparam --- The fall-back DST bank.
 *
 *  Returns:
 *    void
 *
 */
static void _getOutOfBoundState(ATMP_t *aparam, int ianchor, ATMP_t *xparam) {
  int i, idx;
  double hanc, tanc, panc;
  double dh, hb, gt, tb, pb, lnp;
  double texpect, pexpect;
  double dt, prat;

  /*
   *  start by setting 'xparam' to standard atmosphere
   */
  getDefaultParam(xparam);

  /*
   *  save the anchor points
   */
  hanc = aparam->height[ianchor];
  tanc = aparam->temperature[ianchor];
  panc = aparam->pressure[ianchor];

  /* 
   *  find the corresponding layer for the anchor point 
   */
  idx = _getLayerIndex(hanc, xparam->nItem, xparam->height);

  /*
   *  Alter the default atmosphere by shifting the boundary layer temperatures 
   *    and pressures to align with the anchor point T & P.
   */

  /*
   *  compute expected T & adjust the breakpoints accordingly.
   */
  _evaluatePrTp(xparam, idx, hanc, &pexpect, &texpect);

  dt = tanc - texpect;
  for ( i=0; i<xparam->nItem; i++ )
    xparam->temperature[i] += (real4)dt;  /* T values are shifted vertically */

  /*
   *  compute expected P & adjust the breakpoints accordingly.
   *  need to re-evaluate, since base level pressures depend on T values.
   */
  _evaluatePrTp(xparam, idx, hanc, &pexpect, &texpect);

  prat = panc / pexpect;

  /* update the pressure values for higher levels, starting with this one */
  xparam->pressure[idx] = xparam->pressure[idx] * (real4)prat;
  pb = xparam->pressure[idx];
  hb = xparam->height[idx];
  tb = xparam->temperature[idx];
  for ( i=(idx+1); i<xparam->nItem; i++ ) {
    dh = xparam->height[i] - hb;
    dt = xparam->temperature[i] - tb;
    gt = dt / dh;

    if ( fabs(gt) < DOUBLE_FP_PRECISION )
      lnp = -GMR * dh / tb;
    else
      lnp = -GMR * log( (tb+dt) / tb ) / gt;

    xparam->pressure[i] = (real4)( pb * exp( lnp ) );

    hb = xparam->height[i];
    tb = xparam->temperature[i];
    pb = xparam->pressure[i];
  }

  /* now work backwards through the lower levels */
  pb = xparam->pressure[idx] * prat;
  hb = xparam->height[idx];
  tb = xparam->temperature[idx];
  for ( i=(idx-1); i>=0; i-- ) {
    dh = xparam->height[i] - hb;
    dt = xparam->temperature[i] - tb;
    gt = dt / dh;

    if ( fabs(gt) < DOUBLE_FP_PRECISION )
      lnp = -GMR * dh / tb;
    else
      lnp = -GMR * log( (tb+dt) / tb ) / gt;

    xparam->pressure[i] = (real4)( pb * exp( lnp ) );

    hb = xparam->height[i];
    tb = xparam->temperature[i];
    pb = xparam->pressure[i];
  }
}

/**
 *
 *  Reads the next DB Record from the attached DST file. (private)
 *
 *  Input:
 *    none.
 *
 *  Returns:
 *    int --- the return status from the call to 'eventRead()' (either number 
 *          of events read or '-1' for EOF.
 *
 */
inline int _nextAtmosDBRecord(void) {
  if ( lunatmp < 0 ) {
    perr("no DB file loaded!\n");
    return 0;
  }

//   return eventRead(lunatmp, pwb, phb, &pev);

// modify to read GDAS banks and populate fdatmos_param_ accordingly. Keep GDAS a secret.
  int i,rc;

  rc = eventRead(lunatmp, pwb, phb, &pev);
  if (tstBankList(phb,GDAS_BANKID)) {
    fdatmos_param_dst_common *f = &fdatmos_param_;
    memcpy(&fdatmos_param_,&gdasbank_,sizeof(fdatmos_param_));    
    for (i=0; i<f->nItem; i++) {
      f->temperature[i] += CELSIUS_TO_KELVIN;
      f->dewPoint[i] += CELSIUS_TO_KELVIN;
    }
  }
  return rc;
}

/**
 *
 *  Reads the next DB Record from the attached DST file. (private)
 *
 *  Input:
 *    none.
 *
 *  Returns:
 *    int --- the return status from the call to 'eventRead()' (either number 
 *          of events read or '-1' for EOF.
 *
 */
inline int _nextScatterDBRecord(void) {
  if ( lunatmt < 0 ) {
    perr("no DB file loaded!\n");
    return 0;
  }
  
  return eventRead(lunatmt, twb, thb, &tev);
}

/**
 *
 *  Opens an atmosphere DST file.  (public)
 *
 *  Input:
 *    char *path --- The full path of a DST file.
 *
 *  Returns:
 *    '0' if file was opened successfully, '-1' if an error occurred.
 *
 */
int loadAtmosDB(char *path) {
#ifndef JPCOMPARE
  int rc;

  rc = dstOpenUnit(ATMOSDB_FILEDES, path, MODE_READ_DST);

  if ( rc == 0 ) {
    pwb = newBankList(200);
    phb = newBankList(200);
    pdbname = path;
    lunatmp = ATMOSDB_FILEDES;

    _nextAtmosDBRecord();
    perr("loading atmosphere data : %s\n", path);
    return 0;
  }

  perr("error opening '%s' for reading!\n", path);
  perr("reverting to mean US1976 atmosphere.\n");

  lunatmp = -1;
  return -1;
#else
  perr("entered JPCOMPARE mode.  Defaulting to US1976 atmosphere.\n");
  lunatmp = -1;
  return 0;
#endif
}

/**
 *
 *  Opens an atmospheric transmission DST file.  (public)
 *
 *  Input:
 *    char *path --- The full path of a DST file.
 *
 *  Returns:
 *    '0' if file was opened successfully, '-1' if an error occurred.
 *
 */
int loadScatterDB(char *path) {
  int rc;

  rc = dstOpenUnit(FDSCAT_FILEDES, path, MODE_READ_DST);

  if ( rc == 0 ) {
    twb = newBankList(100);
    thb = newBankList(100);
    tdbname = path;
    lunatmt = FDSCAT_FILEDES;

    _nextScatterDBRecord();
    perr("loading fdscat data : %s\n", path);
    return 0;
  }

  perr("error opening '%s' for reading!\n", path);
  perr("reverting to mean values (25,0.04,1).\n");

  lunatmt = -1;
  return -1;
}

/**
 *
 *  Closes the atmosphere DST file.  (public)
 *
 *  Input:
 *    none.
 *
 *  Returns:
 *    The return status from the call to 'dstCloseUnit()', or '-1' if no DST 
 *          file was attached.
 *
 */
int closeAtmosDB(void) {
  int istat;

  if ( lunatmp < 0 )
    return -1;

  delBankList(pwb);
  delBankList(phb);

  istat = dstCloseUnit(lunatmp);

  return istat;
}

/**
 *
 *  Closes the atmospheric transmission DST file.  (public)
 *
 *  Input:
 *    none.
 *
 *  Returns:
 *    The return status from the call to 'dstCloseUnit()', or '-1' if no DST 
 *          file was attached.
 *
 */
int closeScatterDB(void) {
  int istat;

  if ( lunatmt < 0 )
    return -1;

  delBankList(twb);
  delBankList(thb);

  istat = dstCloseUnit(lunatmt);

  return istat;
}

/**
 *
 *  Gets the physical parameters that define the atmosphere at the time given 
 *    by 'jsec'.   Parameters for the US1976 standard atmosphere model are 
 *    used if a valid record is not found.  (public)
 *
 *  Input:
 *    double jsec --- A time in number of seconds after epoch (midn. 1/1/1970)
 *
 *  Output:
 *    fdatmos_param_t *aparam --- A DST bank initialized with data from the 
 *          attached DST file.
 *    fdatmos_trans_t *atrans --- A DST bank initialized with data from the 
 *          attached DST file.
 *
 *  Returns:
 *    '0' if a record was found that includes the given jsec, '-1' otherwise.
 *
 * *** ELB 2010/06/29 Made one change which will allow this to work with
 *     dates that are not time ordered.  Otherwise this shouldn't have any effect
 */

int getAtmosphere(const RuntimeParameters *par, ATMP_t *aparam, ATMT_t *atrans) {
  int istat, loop, i;
  double jsec = mjlday2eposec(par->jday-MJLDOFF);

  /*
   *  Seek the next atmosphere (thermo) for the given jsec
   */
  loop = 0;

/*
 * ***ELB changed the below to while and added a condition on
 * on lunatm so that this doesn't explode it the USStdAtmos is 
 * used. 
 */
  istat = lunatmp; 
    // do {
  while (istat > 0){

    if ( jsec >= (double)fdatmos_param_.dateTo ) {
      /* 
       *  need to seek forward 
       */

      istat = _nextAtmosDBRecord();

    }
    else if ( jsec >= (double)fdatmos_param_.dateFrom ) {
      /* 
       *  we're in the right time frame.  note that we don't want to simply 
       *    repoint structures here, but copy their contents.  this is why we 
       *    chose the lines "*a = f" instead of "a = &f".
       */

      *aparam = fdatmos_param_;

      /* these lines are necessary to correct for DST files produced by non-TRUMP
       * sources (i.e., not those produced by "raob" software). Sometimes the height
       * is expressed in km above the CLF plane, instead of the geopotential height.
       * -TAS 20150228
       */
      if (aparam->height[0] < 0.6 // this is underground everywhere in Utah
          && aparam->pressure[0] < 900) // this seems adequate to catch a mismatch
      {
//         fprintf(stderr,"WARNING: Adding %.3f km to atmosphere reference heights. Adding %f to temperatures and dewpoints.\n",0.001*CLF_ALTITUDE,CELSIUS_TO_KELVIN);
        for (i=0; i<aparam->nItem; i++) {
          aparam->height[i] += 0.001*CLF_ALTITUDE;
          aparam->temperature[i] += CELSIUS_TO_KELVIN;
          aparam->dewPoint[i] += CELSIUS_TO_KELVIN;
        }
      }
      
      
      
      istat = 0;
    }
    else {
      /* re-load the file and try once more */
      //***ELB Think there's a bug here istat will be zero and will fall through
      if(closeAtmosDB() <= 0){
        fprintf(stderr, "Warning: failed to close previous atmospheric DB\n");
      }
      if(loadAtmosDB(pdbname) <=0){
        fprintf(stderr, "Error: failed to reload atmospheric DB\n");
        istat = -1;
      }
      else{
        istat = 1; // ***ELB will go back around
      }
      // ***ELB if it's already gone around once and not found the
      // atmosphere it never will.  Stop and set istat = -1
      if ( (++loop) > 1 )
        istat = -1;
    }
  }
  /*
   * ***ELB Also made this change to be consistant with the change from 
   * do -> while
   */
//  } while ( istat > 0 );

  /*
   *  Revert to US Standard atmosphere if no appropriate record is found.
   */
  if ( istat < 0 ) {
    perr("unable to find consistent data record.\n");
    perr("reverting to US1976 model.\n");
    perr("  jsecond = %lf\n", jsec);
    perr("  DB file = '%s'\n", pdbname);

    getDefaultParam(aparam);
    fdatmos_param_ = *aparam;
  }

  /* get the corresponding over- & under- flow states */
  _getOutOfBoundState(aparam, 0, ufparam);
  _getOutOfBoundState(aparam, aparam->nItem-1, ofparam);


  /*
   *  Seek the next atmosphere (scatter) for the given jsec
   */
  loop = 0;
  do {
    if ( jsec >= fdscat_.endTime ) {
      /* 
       *  need to seek forward 
       */

      istat = _nextScatterDBRecord();

    }
    else if ( jsec >= fdscat_.startTime ) {
      /* 
       *  we're in the right time frame.  note that we don't want to simply 
       *    repoint structures here, but copy their contents.  this is why we 
       *    chose the lines "*a = f" instead of "a = &f".
       */

      *atrans = fdscat_;
      istat = 0;
    }
    else {
      /* re-load the file and try once more */
      istat = closeScatterDB();
      istat += loadScatterDB(tdbname);

      if ( (++loop) > 1 )
        istat = -1;
    }
  } while ( istat > 0 );

  /*
   *  Revert to mean HiRes parameters if no appropriate record is found.
   */
  if ( istat < 0 ) {
    perr("unable to find consistent data record.\n");
    perr("reverting to mean scatter params.\n");
    perr("  jsecond = %lf\n", jsec);
    perr("  DB file = '%s'\n", tdbname);

    getDefaultTrans(atrans);
  }

  return istat;
}
/**
 *
 *  Evaluates the air pressure, temperature, and density at altitude z above 
 *    sea level. (public)
 *
 *  Input:
 *    fdatmos_param_t *aparam --- A DST bank that defines the thermodynamic 
 *          quantities of the atmosphere.
 *    double z --- The real sea-level altitude (meters)
 *
 *  Output:
 *    double *p --- Air pressure at altitude 'z', in the context of the given 
 *          atmophere. (Pa)
 *    double *t --- Temperature at altitude 'z', in the context of the given 
 *          atmophere. (Kelvin)
 *    double *d --- Density of air at altitude 'z', in the context of the 
 *          given atmophere. (m^-3)
 *
 */
void getPTDByAltitude(ATMP_t *aparam, double z, double *p, double *t, double *d) {
  int idx;
  double h, dh, x, lnp;

  if ( z < 0.00 ) {
    fprintf(stderr, "getPDTByAltitude(): Not equipped to evaluate PDT below sea level!\n");
    *p = 0.00;
    *d = 0.00;
    *t = 0.00;
    return;
  }

  /* convert 'z' to geopotential altitude in km */
  z *= 1.0e-03;
  h = z * REARTH / ( z + REARTH );

  if ( h < aparam->height[0] ) {
    /* 
     *  underflow state.  revert to US standard atmosphere, anchoring to T[0] 
     *    and P[0] 
     */
    idx = _getLayerIndex(h, ufparam->nItem, ufparam->height);
    _evaluatePrTp(ufparam, idx, h, p, t);
  }
  else if ( h < aparam->height[(aparam->nItem-1)] ) {
    /*
     *  normal range.  interpolate between points if distance between layer 
     *    heights is sufficiently small, like 1km.  Otherwise, use the 
     *    hydrostatic formulae.
     */

    idx = _getLayerIndex(h, aparam->nItem, aparam->height);

    /* expect that ( idx < nItem-1 ), or this code may seg-fault! */

    dh = aparam->height[idx+1] - aparam->height[idx];

    if ( dh <= 1.0 ) {
      x = ( h - aparam->height[idx] ) / dh;

      /* temperature is linearly extrapolated between points */
      *t = (1.0-x)*aparam->temperature[idx] + x*aparam->temperature[idx+1];

      /* pressure is log-linearly extrapolated between points */
      lnp = (1.0-x)*log(aparam->pressure[idx]);
      lnp += x*log(aparam->pressure[idx+1]);

      *p = exp( lnp );
    }
    else {
      /* try using hydrostatic formulae instead */
      _evaluatePrTp(aparam, idx, h, p, t);
    }
  }
  else {
    /* 
     *  overflow state.  revert to US standard atmosphere, anchoring to T[n-1] 
     *    and P[n-1] 
     */

    idx = _getLayerIndex(h, ofparam->nItem, ofparam->height);
    _evaluatePrTp(ofparam, idx, h, p, t);
  }

  *p *= 1.0e+02;  /* convert from hPa to Pa */
  *d = M/R * (*p)/(*t) * 1.0e-6;  /* ideal gas law */
  *p *= 1.0e-03;  /* convert from Pa to kPa */
}

/**
 *
 *  Evaluates the air pressure at altitude z above sea level. (public)
 *
 *  Input:
 *    fdatmos_param_t *aparam --- A DST bank that defines the thermodynamic 
 *          quantities of the atmosphere.
 *    double z --- The real sea-level altitude (meters)
 *
 *  Returns:
 *    double --- Air pressure at altitude 'z', in the context of the given 
 *          atmophere. (kPa)
 *
 */
double getPressureByAltitude(ATMP_t *aparam, double z) {
  double p, d, t;

  getPTDByAltitude(aparam, z, &p, &t, &d);

  return p;
}

/*
 *  Evaluates the temperature at altitude z above sea level. (public)
 *
 *  Input:
 *    fdatmos_param_t *aparam --- A DST bank that defines the thermodynamic 
 *          quantities of the atmosphere.
 *    double z --- The real sea-level altitude (meters)
 *
 *  Returns:
 *    double --- Temperature at altitude 'z', in the context of the given 
 *          atmophere. (Kelvin)
 */
double getTemperatureByAltitude(ATMP_t *aparam, double z) {
  double p, d, t;

  getPTDByAltitude(aparam, z, &p, &t, &d);

  return t;
}

/*
 *  Evaluates the density of the air at altitude z above sea level. (public)
 *
 *  Input:
 *    fdatmos_param_t *aparam --- A DST bank that defines the thermodynamic 
 *          quantities of the atmosphere.
 *    double z --- The real sea-level altitude (meters)
 *
 *  Returns:
 *    double --- Density of air at altitude 'z', in the context of the given 
 *          atmophere. (g/cm^3)  (comment was wrong before. SS-10.18.11)
 */
double getDensityByAltitude(ATMP_t *aparam, double z) {
  double p, d, t;

  getPTDByAltitude(aparam, z, &p, &t, &d);

  return d;
}

#undef REARTH
#undef G
#undef M
#undef R
#undef GMR

#undef NUM_US1976_LAYERS 
#undef P0

/*
 *  The following functions are meant to reconcile this code with 
 *    'atmosphere.c', which will eventually be retired.
 */
#define ATMOS_EPS 1.0e-5
#define UNIV_GAS_CONST 8.31432   // N.m/(mol.K)
#define AIR_MOLAR_MASS 0.0289644 // kg/mol
#define KPA_TO_TORR 7.5006

/*
 * This taken from http://emtoolbox.nist.gov, in the documentation
 * section for a "Shop-floor Formula for the Index of Refraction of
 * Air," and using 0% relative humidity
 */
#define DELTA_COEFFICIENT 7.86e-4
/*
 * The grammage rate of change is just the density 
 * but in units of g/cm2/m (or 100*g/cm3).
 */
double getGrammageRateOfChange(ATMP_t *aparam, double altitude) {
  return -100.*getDensityByAltitude(aparam, altitude);
}

/* The grammage is just the pressure in hPa */
double getGrammageByAltitude(ATMP_t *aparam, double altitude) {
  return 10.0*getPressureByAltitude(aparam, altitude);
}

/*
 * This static function is meant to be passed to an integrator
 * The static variables are set in getSlantDepth
 */
static double pointA[3];
static double fromAtoBnorm[3];
static ATMP_t *seasonAtoB;
double setupAtoB(double ptA[3], double ptB[3], ATMP_t *season) {
  double fromAtoB[3];
  cpyvec(ptA, pointA);
  subvec(ptB, ptA, fromAtoB);
  unitVector(fromAtoB, fromAtoBnorm);
  double distAtoB = dotprod(fromAtoB, fromAtoBnorm);
  seasonAtoB = season;
  return distAtoB;
}

double getDensityAlongVector(double dl) {
  double point[3];
  double lat,lon,alt;

  addvec(pointA, dl*fromAtoBnorm, point);

  r2lla(point, &lat, &lon, &alt);
  if (alt > TOP_OF_ATMOSPHERE) 
    return 0.;

  else if (alt < 0.)
    return getDensityByAltitude(seasonAtoB,0.);
  else 
    return getDensityByAltitude(seasonAtoB,alt);
}

/*
double getCurrentDensityByAltitude(double alt) {
  if ( alt >= 0.00 ) 
    return getDensityByAltitude(seasonAtoB, alt);
  else
    return getDensityByAltitude(seasonAtoB, 0.00);
  }
*/

double getLRhoAlongVector(double dl) {
  return dl*getDensityAlongVector(dl);
}

double getLRhoIntegralBetweenPoints(ATMP_t *season, 
				    double ptA[3], double ptB[3]) {
  double lat,lon,altA,altB;

  r2lla(ptA, &lat, &lon, &altA);
  r2lla(ptB, &lat, &lon, &altB);

  if (altA > TOP_OF_ATMOSPHERE && altB > TOP_OF_ATMOSPHERE)
    return 0.;

  double distAtoB = setupAtoB(ptA,ptB,season);
  double dl = 0.;

  if (altA > TOP_OF_ATMOSPHERE)
    dl = getDistanceToAtmosphere(ptA, fromAtoBnorm, 0.1);
  else if (altB > TOP_OF_ATMOSPHERE)
    distAtoB += getDistanceToAtmosphere(ptB, fromAtoBnorm, 0.1); // should be negative

  double lRhoAtoB = 1.e2*qsimp(getLRhoAlongVector,dl,distAtoB,ATMOS_EPS);

  return lRhoAtoB;
}

/*
 * (DRB) getGrammageBetweenPoints is almost a clone of getSlantDepth.
 * The difference being that the points in getGrammageBetweenPoints
 * are in Earth center coordinates.
 */
double getGrammageBetweenPoints(ATMP_t *season, double ptA[3], double ptB[3]) {
  /* Do sanity check */
  double lat, lon, altA, altB;

  r2lla(ptA, &lat, &lon, &altA);
  r2lla(ptB, &lat, &lon, &altB);

  if (altA > TOP_OF_ATMOSPHERE && altB > TOP_OF_ATMOSPHERE)
    return 0.;
  else if (altA < 0. || altB < 0.)
    perr("attempt to evaluate grammage below sea level.\n");

  double distAtoB = setupAtoB(ptA,ptB,season);
  double dl = 0.;

  if (altA > TOP_OF_ATMOSPHERE)
    dl = getDistanceToAtmosphere(ptA, fromAtoBnorm, 0.1);
  else if (altB > TOP_OF_ATMOSPHERE)
    distAtoB += getDistanceToAtmosphere(ptB, fromAtoBnorm, 0.1); // should be negative
  double gramAtoB = 1.e2*qsimp(getDensityAlongVector, dl, distAtoB, ATMOS_EPS);

  return gramAtoB;
}
  
/*
 * It seems that this should return the actual amount of material
 *   between points rt and rb
 * getSlantDepth should return the same number as
 *   getGrammageBetweenPoints
 *
 * rt & rb should be in CLF coordinates. 
 */
double getSlantDepth(ATMP_t *season, double rt[3], double rb[3]) {
  double dr[3], dtop, len, lat, lon, alta, altb;

  r2lla(rt, &lat, &lon, &alta);
  r2lla(rb, &lat, &lon, &altb);

  if (alta > TOP_OF_ATMOSPHERE && altb > TOP_OF_ATMOSPHERE) 
    return 0.00;

  seasonAtoB = season;

  cpyvec(rt, pointA);

  subvec(rb, rt, dr);
  unitVector(dr, fromAtoBnorm);

  if (alta > TOP_OF_ATMOSPHERE)
    dtop = getDistanceToAtmosphere(pointA, fromAtoBnorm, 0.1);
  else
    dtop = 0.00;

  len = dotprod(dr, fromAtoBnorm);
  if (altb > TOP_OF_ATMOSPHERE)
    len += getDistanceToAtmosphere(rb, fromAtoBnorm, 0.1); // should be negative
  
  /* We should do this by real numerical integration ... */
  return 100.*qsimp(getDensityAlongVector,dtop,len,ATMOS_EPS);
}

/*
 * getDistanceByGrammage returns the distance one must move along a
 * given vector from a given point to have a given grammage traversed.
 *   direction must be a UNIT vector.
 */
double getDistanceByGrammage(ATMP_t *season, 
			     double ptA[3], 
			     double direction[3], 
			     double grammage,
			     double eps) {
  double r[3];
  double s[3];
  double dr;
  double dl = 0., dls = 0.;
  double dx = 0.;
  double dxLeft = grammage;
  double dxdr;
  double lat, lon, alt, alts;
  cpyvec(ptA, r);
  r2lla(r, &lat, &lon, &alt);		// new

  /* Get r to top of atmosphere, no change to x */
  if (alt > TOP_OF_ATMOSPHERE) {
    dl = getDistanceToAtmosphere(r, direction, 0.1); // Get 10 cm into atm
    addvec(r, dl*direction, s);
    cpyvec(s, r);
    r2lla(r, &lat, &lon, &alt);
  }

  while (fabs(dxLeft) > eps) {
    assert ( alt < TOP_OF_ATMOSPHERE );
    dxdr = 100.*getDensityByAltitude(season, alt);
    dr = dxLeft / dxdr;
    dl += dr;
    addvec(r, dr*direction, s);

    /* recalculate s if it is above the atmosphere */
    r2lla (s, &lat, &lon, &alts);
    if (alts > TOP_OF_ATMOSPHERE) {
      dls = getDistanceToAtmosphere(s, direction, 0.1);
      addvec (r, dls*direction, s);
    }

    double gbp = getGrammageBetweenPoints(season, r, s);
    dx += (dr/fabs(dr)) * gbp;
    dxLeft = grammage - dx;
    cpyvec(s, r);
    r2lla(r, &lat, &lon, &alt);
  }

  return dl;
}

double getLargeDistanceByGrammage(ATMP_t *season, 
				  double ptA[3], 
				  double direction[3], 
				  double grammage,
				  double eps) {
  double r[3];
  double dl;
  double distance = 0.;
  double dxLeft = grammage;
  double dX = 1.0;
  cpyvec(ptA, r);

  double intgram = (double)((int)grammage);
  double extragram = grammage - intgram;

  do {
    dl = getDistanceByGrammage(season, r, direction, dX, eps);
    addvec(r, dl*direction, r);
    dxLeft -= dX;
    distance += dl;
  } while (dxLeft > dX);

  distance += getDistanceByGrammage(season, r, direction, extragram, eps);

  return distance;
}

double getDistanceToAtmosphere(const double ptA[3], 
			       const double direction[3], 
			       double eps) {
  double dir[3];
  double r[3];
  double s[3];
  double dr;
  double dl = 0.;
  double lat,lon,alt;
  double dh, dhdr;
  cpyvec(direction, dir);
  cpyvec(ptA, r);
  r2lla(r, &lat,&lon,&alt);
  if (alt <= TOP_OF_ATMOSPHERE) return 0.;

  dh = alt - TOP_OF_ATMOSPHERE;
  while ( fabs(dh) > eps || dh > DOUBLE_FP_PRECISION ) {
    dhdr = -dotprod(r, dir)/magvec(r);
    dr = dh / dhdr;
    addvec(r, dr*direction, s);
    dl += dr;
    cpyvec(s, r);
    r2lla(r, &lat,&lon,&alt);
    dh = alt - TOP_OF_ATMOSPHERE;
  }
  return dl;
}

/*
 * This taken from http://emtoolbox.nist.gov, in the documentation
 * section for a "Shop-floor Folmula for the Index of Refraction of
 * Air, and using 0% relative humidity
 */
double getDeltaByPressTemp(double P, double T) {
  return DELTA_COEFFICIENT * P/T;
}

#undef ATMOS_EPS
#undef UNIV_GAS_CONST
#undef AIR_MOLAR_MASS
#undef KPA_TO_TORR
#undef DELTA_COEFFICIENT

