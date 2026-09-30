#ifndef  _TRACK_FINDER_H_
#define  _TRACK_FINDER_H_

#include "constants.h"

#define  L1_CLUSTER_SIZE  5
#define  L2_CLUSTER_SIZE  3

int initTFSim(void);

int isLevel1Trigger(short *hit_pt);
int isLevel2Trigger(short *hit_pt);

#endif
