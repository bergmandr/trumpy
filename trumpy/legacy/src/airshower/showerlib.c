#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "constants.h"
#include "control.h"
#include "event.h"
#include "random.h"
#include "toolbox.h"

#include "airshower.h"
#include "showerlib.h"

static int numLibraries;
static ShowerLibrary *library[MAX_LIBRARIES];

int loadShowerLibraries(RuntimeParameters *par) {
  int i;
  int rc = 0;

  numLibraries = par->nshowlib;
  for ( i=0; (i<par->nshowlib)&&(i<MAX_LIBRARIES); i++ ) {
    library[i] = newInstanceOf(ShowerLibrary);
    rc += loadShowerLibrary(par->showlibFile[i], library[i]);
    library[i]->pspec = par->species[i];
  }

  return rc;
}

int loadShowerLibrary(char *path, ShowerLibrary *lib) {
  int nent;
  int rc, u=SHOWERLIB_FILEDES, m=MODE_READ_DST;
  int w, g, s=200, e;
  char str[80];

  pout("loading shower library '%s'.\n", path);

  strcpy(str, path);

  w = newBankList(s);
  g = newBankList(s);

  rc = dstOpenUnit(u, str, m);
  if ( rc != 0 ) {
    perr("unable to open '%s' for reading (%d)", path, rc);
    return rc;
  }

  /* get shower scale object */
  eventRead(u, w, g, &e);
  lib->scaledef = showscale_;

  /* When using Nerling's alpha_eff, don't scale shower */
  // TAS 20150127: disable this line! I need <0 for some applications.
//   lib->scaledef.nMaxScaleFactor = 1.0;

  nent = 0;
  while ( eventRead(u, w, g, &e) >= 0 ) {
    lib->entry[nent++] = showlib_;
  }

  rc = dstCloseUnit(u);
  delBankList(w);
  delBankList(g);

  return 0;
}

/* dummy function defs for shower library accessor routines */
//int getGaisserHillasParameters
//  (ShowerLibrary *lib, double loge, double secth, GaisserHillasParameters *gh) {
int getGaisserHillasParameters(AirShower *as, GaisserHillasParameters *gh) {
  int ie, iq, idr, ik, idx;
  int i, j;
  double width, ratio;
  double loge = as->loge;
  double secth = 1.0 / cos(as->zenith);
  ShowerLibrary *lib = NULL;

  for ( i=0; i<numLibraries; i++ ) {
    if ( library[i]->pspec == as->species ) {
      lib = library[i];
      break;
    }
  }

  if ( lib == NULL ) {
    perr("could not find matching library for part. spec. %d\n", as->species);
    return -1;
  }

  ie = binsearch(loge, lib->scaledef.nlE+1, lib->scaledef.lEEdge);
  iq = binsearch(secth, lib->scaledef.nth+1, lib->scaledef.thEdge);

/*   idr = ie*lib->scaledef.nth + iq; */
/*   ik = (int)( (double)NENTRIES_PER_DRAWER * RANDOM_NUMBER ); */
/*   idx = NENTRIES_PER_DRAWER*idr + ik; */
  ik = (int)( (double)lib->scaledef.numEntries[ie][iq] * RANDOM_NUMBER );
//   printf("selected index %d at random from %d entries (loge %f secth %f)\n",
//          ik,lib->scaledef.numEntries[ie][iq],loge,secth);
  idx = 0;
  for ( i=0; i<lib->scaledef.nlE; i++ ) {
    for ( j=0; j<lib->scaledef.nth; j++ ) {
      if ( i == ie && j == iq ) {
	i = lib->scaledef.nlE;
	break;
      }
      idx += lib->scaledef.numEntries[i][j];
    }
  }
  idx += ik;

  gh->nmax = showscale_getScaledNMax(&lib->scaledef, &lib->entry[idx], loge);
  gh->xmax = showscale_getScaledXMax(&lib->scaledef, &lib->entry[idx], loge);
  
  if (lib->scaledef.nMaxScaleFactor < 0) { // working with WIDTH as parameter to be scaled
//     width = getWidthFromLambdaX0(&lib->entry[idx]);
    ratio = getRatioFromLambdaX0(&lib->entry[idx]);
    width = showscale_getScaledWidth(&lib->scaledef, &lib->entry[idx], loge);
    gh->lambda = getLambdaFromScaled(width, ratio);
    gh->x0 = getX0FromScaled(width, ratio, gh->xmax) ;
  }
  else { // working with X0 as fit parameter to be scaled
    gh->x0 = showscale_getScaledX0(&lib->scaledef, &lib->entry[idx], loge);
    gh->lambda = lib->entry[idx].lambda;
  }

  



  

  if ( isnan(gh->nmax) ) {
    vperr("  ERROR!!!\n");
    vperr("  Nmax=%f  Xmax=%f  X0=%f  lambda=%f  (idx=%d)\n", 
	  gh->nmax, gh->xmax, gh->x0, gh->lambda, idx);
    vperr(" Entry: %e %e %e %e\n",lib->entry[idx].nmax,lib->entry[idx].xmax,
          lib->entry[idx].x0,lib->entry[idx].lambda);
    exit(66);
  }

  return idx;
}

