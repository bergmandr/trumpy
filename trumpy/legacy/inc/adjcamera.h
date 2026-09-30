#ifndef _ADJCAMERA_H_
#define _ADJCAMERA_H_

#include "geofd_dst.h"

static const int nAdjCameras[GEOFD_MAXMIR] = { 
  3, 3, 5, 5, 5, 5, 5, 5, 5, 5, 3, 3
};

static const int adjCameraList[GEOFD_MAXMIR][MAX_NADJCAM] = {
  {1, 2, 3, -1, -1}, {0, 2, 3, -1, -1}, {0, 1, 3, 4, 5}, {0, 1, 2, 4, 5},
  {2, 3, 5, 6, 7}, {2, 3, 4, 6, 7}, {4, 5, 7, 8, 9}, {4, 5, 6, 8, 9},
  {6, 7, 9, 10, 11}, {6, 7, 8, 10, 11}, {8, 9, 11, -1, -1},
  {8, 9, 10, -1, -1}
};


#endif
