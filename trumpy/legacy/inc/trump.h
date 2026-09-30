#ifndef _TRUMP_H_
#define _TRUMP_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "event.h"
#include "random.h"
#include "constants.h"
#include "control.h"
#include "nerling.h"
#include "toolbox.h"

#include "atmosphere-withdb.h"
#include "airshower.h"
#include "showerlib.h"
#include "eventlist.h"
#include "track.h"
#include "nkg.h"

#include "fdconstants.h"
#include "fdsite.h"
#include "calibration.h"
#include "raytrace.h"
#include "acpttrack.h"
#include "taelectronics.h"
#include "tlelectronics.h"
//#include "pelist.h"

#include "bank-manager.h"

int initTRUMP(int nc, RuntimeParameters *par[]);
fdraw_dst_common *initFDSite(RuntimeParameters *par, UCalibration *calib,
			     FDSiteGeometry *fdsg, TGeom *tgeom,
			     int banklist, fdraw_dst_common *fdraw);

#endif
