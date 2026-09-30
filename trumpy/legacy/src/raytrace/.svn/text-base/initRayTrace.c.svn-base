#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "constants.h"
#include "control.h"
#include "event.h"
#include "fdconstants.h"
#include "toolbox.h"

#include "calibration.h"
#include "tacalibration.h"
#include "raytrace.h"
// #include "spotsize.h"

void initRayTrace (geofd_dst_common *geo, TGeom *tgm, UCalibration *calib) {
  int i, j, k, l;
  double segangle[GEOFD_MAXMIR];

  for (i=0; i<geo->nmir; i++) {
    segangle[i] = asin( (geo->seg_point2point / 2.) / geo->rcurve3[i] );
    for (j=0; j<geo->nseg[i]; j++)
      tgm->spot[i][j] = geo->seg_spot[i][j] * (0.01/0.195);
//     switch (geo->uniqID) {
//       case GEOFD_UNIQBRLRSCOTT:        
//         tgm->spot[i]  =  0.010;
//         break;
//       case GEOFD_UNIQBRLRTOKUNO:
//       case GEOFD_UNIQBRSTANSSP:
//       case GEOFD_UNIQLRSTANSSP:
//       case GEOFD_UNIQBRLRSSPTHOMAS:
//         tgm->spot[i]  =  ssp[geo->siteid][i]*(0.01/0.195);
//         break;
//       case GEOFD_UNIQBRSTAN:
//       case GEOFD_UNIQLRSTAN:
//       case GEOFD_UNIQBRLRTHOMAS:
//         tgm->spot[i]  = 0.0285*(0.01/0.195);
//         break;
//       default:
//         tgm->spot[i]  =  0.010;
//     }
  }
  tgm->mir_reflect    =  0.92;      // average "clean mirror" reflectivity
  tgm->camcover_trans =  0.910642;  // transmission at 337 nm
  tgm->filter_trans   =  0.891514;  // transmission at 337 nm
  tgm->n_air          =  1.000293;
  tgm->n_glass        =  1.50;
  tgm->n_filter       =  (geo->siteid<=1)?1.53:1.0;
  tgm->pmtQE          =  (geo->siteid<=1)?(calib->tc)->pmtmaxqe:0.25;
  tgm->thick_glass    =  0.005;
  tgm->thick_incam    =  0.020;
  tgm->thick_filter   =  0.004;

  for (i=0; i<geo->nmir; i++) {
    tgm->sep[i] = geo->sep3[i] - tgm->thick_glass 
                                - tgm->thick_incam 
                                - tgm->thick_filter;
    tgm->h[i] = geo->rcurve3[i] - sqrt (geo->rcurve3[i]*geo->rcurve3[i] - 
                                        geo->diameters[i]*geo->diameters[i]/4.);
    tgm->maxseg[i] = 2. * geo->rcurve3[i] * sin (segangle[i] / 2.);

    tgm->meantime[i] =  geo->sep3[i] / SPEED_OF_LIGHT +
      sqrt(geo->sep3[i]*geo->sep3[i] +
           geo->diameters[i] * geo->diameters[i] / 4.) / SPEED_OF_LIGHT;
    tgm->meantime[i] /= 2.;
  }

  for (i=0; i<6; i++) {
    tgm->hexpt[i][0] = -cos ( (30+i*60)*D2R );
    tgm->hexpt[i][1] = -sin ( (30+i*60)*D2R );
    tgm->hexpt[i][2] = 0.;
  }

  /* Get rotation angle for each segment and get vectors
   * from points to centers of mirror segments
   */
  for (l=0; l<geo->nmir; l++) {
    if (geo->camtype[l] == GEOFD_TA) {
      for (i=0; i<geo->nseg[l]; i++) {
        tgm->seg_the[l][i] = M_PI - acos (geo->vseg3[l][i][2]);
        tgm->seg_phi[l][i] = atan2 ( geo->vseg3[l][i][1], geo->vseg3[l][i][0] );
        for (j=0; j<6; j++) {
          for (k=0; k<3; k++)
            tgm->hexmirpt[l][i][j][k] = tgm->hexpt[j][k];

          rotz (tgm->hexmirpt[l][i][j],  tgm->seg_phi[l][i], tgm->hexmirpt[l][i][j]);
          roty (tgm->hexmirpt[l][i][j],  tgm->seg_the[l][i], tgm->hexmirpt[l][i][j]);
          rotz (tgm->hexmirpt[l][i][j], -tgm->seg_phi[l][i], tgm->hexmirpt[l][i][j]);
        }
      }
    }
    else if (geo->camtype[l] == GEOFD_MD || geo->camtype[l] == GEOFD_TALE) {
      tgm->cosang[l] = sqrt(geo->rcurve3[l]*geo->rcurve3[l] -
                          geo->seg_flat2flat *
                          geo->seg_flat2flat / 4.) /
                          geo->rcurve3[l];
    }
  }

  switch (geo->siteid) {
    case BLACK_ROCK_SITEID:
    case LONG_RIDGE_SITEID:
      tgm->max_sensitivity = (calib->tc)->pmtmaxunif;
      for (i=0; i<PMT_NUMDIV; i++)
        for (j=0; j<PMT_NUMDIV; j++)
          tgm->pmt_map[i][j] = (calib->tc)->pmtrelunif[i][j];

      buildTACameraPoles(geo, tgm);
      break;
      
    case MIDDLE_DRUM_SITEID:
      
//       break;
    case TALE_SITEID:
      tgm->max_sensitivity = 1.0;
      for (i=0; i<PMT_NUMDIV; i++)
        for (j=0; j<PMT_NUMDIV; j++)
          tgm->pmt_map[i][j] = 1.0;
        
      buildMDCameraPoles(geo, tgm);
      break;
  }
}

void clearRayTrace (RayTrace *ray) {
  free(ray->pmt);                 
  free(ray->tpe);
}
