#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "trump.h"  // TODO: narrow this down a bit

void loadAcptMap(int icam, float acpt[256][200*200]);

int acptmap(TGeom *tgm, RayTrace *ray) {
  static int lcam = -1;
  static float acpt[256][200*200];

  int i;
  double pcum[256] = { 0.0 };
  
  if ( ray->cam != lcam ) {
    loadAcptMap(ray->cam, acpt);
    lcam = ray->cam;
  }
  
  double r = magvec(ray->vsite);
  double x = ray->vsite[0] / r;
  double y = ray->vsite[1] / r;

  double alt = asin(-y);
  double azm = asin(x/sqrt(1-y*y));

  int altbin = (int)( ( alt - (-10.*M_PI/180.) ) / (0.1*M_PI/180.) );
  if ( (altbin<0) || (altbin>=200) )
    return -1;

  int azmbin = (int)( ( azm - (-10.*M_PI/180.) ) / (0.1*M_PI/180.) );
  if ( (azmbin<0) || (azmbin>=200) )
    return -1;

  int k = 200*altbin + azmbin;

  // find cumulative probability for any tube to accept PE from this direction
  pcum[0] = acpt[0][k];
  for ( i=1; i<256; i++ )
    pcum[i] = pcum[i-1] + acpt[i][k];
  
  double u = RANDOM_NUMBER;  // in range [0,1)

  if ( u >= pcum[255] )
    return -1;  // not accepted by any tube

  for ( i=0; i<256; i++ ) {
    if ( u < pcum[i] ) {
      ray->time = r/SPEED_OF_LIGHT + tgm->meantime[ray->cam];
      return i;
    }
  }

  return -66;  // an error must have occurred!
}

void loadAcptMap(int icam, float acpt[256][200*200]) {
  int i;
  char path[256];

  sprintf(path, "%s/dat/acptmaps/acptmap%02d.dat", TRUMP_HOME, icam%4);
  
  FILE *fp = fopen(path, "r");
  if ( fp == NULL ) {
    perr("unable to open %s for reading.\n", path);
    return;
  }
  
  pout("loaded %s.\n", path);

  for ( i=0; i<256; i++ )
    fread(acpt[i], sizeof(float), 200*200, fp);
  fclose(fp);

  return;
}
