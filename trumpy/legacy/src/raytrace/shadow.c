#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "control.h"
#include "event.h"
#include "fdconstants.h"
#include "toolbox.h"

#include "calibration.h"
#include "tacalibration.h"
#include "raytrace.h"

#define SQR(X) ((X)*(X))
#define unitVector(A,B) ({\
      double den=sqrt(SQR(A[0])+SQR(A[1])+SQR(A[2]));  \
      B[0]=A[0]/den; B[1]=A[1]/den; B[2]=A[2]/den;     \
    })


#define DUMPX3D     // Create X3D model files for the mirrors (poles only) - SBT 

void buildTACameraPoles (geofd_dst_common *geo, TGeom *tgm) {
  int cam, i, j;
  double vertbar[3];               // direction of vert. bars on cluster stand
  double horzbar[3];               // direction of horz. bars on cluster stand
  double deepbar[3];               // direction of deep  bars on cluster stand
  double zclust[3];                // location of front of cluster box
  double zcam[3];                  // location of camera (mirror?) center;    
  double ybar[3];                  // Direction to front of bay - SBT
  double xbar[3];                  // Direction to right side of bay - SBT

  const double prad = 0.02;        // pole radius (meters)
  const double pmid = 0.01;        // small pole radius (meters)
  const double pwid = 1.32;        // separation between legs of
                                   // main support poles across face of
                                   // camera (meters)
  const double pbig = 0.076;       // big support pole radius (meters)
  const double pdep = 0.40;        // separation between legs of
                                   // main support poles along depth of
                                   // camera (meters)
  const double phgt = 6.745;       // height of main support poles (meters)
// double hring[2] = {6.11, 2.94}; // height of rings 1 and 2 camera center to ground (meters)
  double hring[2];                 // Must be computed - SBT 20150430

  const double wrad1 = 0.0175;     // Pole radii for walkway bannister
  const double wrad2 = 0.0048;

  const double dR1 = 1.5267;       // Height above floor of center of ring 2 mirrors (from drawings)
  const double dR2 = 3.9727;       // Height of center of ring 1 mirrors above center of ring 2 mirrors (from drawings)

  double sinang, cosang;

  double theta;                    // Standard rotation angle from mirror to bay coordinates (should be good enough) SBT 20150430

  for (cam=0; cam<geo->nmir; cam++) {
    sinang = sin(geo->mir_the[cam]);
    cosang = cos(geo->mir_the[cam]);

    vertbar[0] = 0.;               
    vertbar[1] =  sinang;
    vertbar[2] =  cosang;
    
    horzbar[0] = 1.;
    horzbar[1] = 0.;
    horzbar[2] = 0.;

    deepbar[0] = 0.;
    deepbar[1] = -cosang;
    deepbar[2] = sinang;

    zclust[0]  = 0.;
    zclust[1]  = 0.;
    zclust[2]  = -geo->rcurve3[cam] + geo->sep3[cam];

    zcam[0] = 0.;
    zcam[1] = 0.;
    zcam[2] = -geo->rcurve3[cam];
    
    // Directions relative to the telescope bay - SBT 20150430

    if (cam == 0 || cam == 1 ||
	cam == 4 || cam == 5 ||
	cam == 8 || cam == 9)
      theta = +9.0*M_PI/180.0;
    else
      theta = -9.0*M_PI/180.0;
    
    for(i=0;i<3;i++)
      {
	xbar[i] =  cos(theta)*horzbar[i] + sin(theta)*deepbar[i];
	ybar[i] = -sin(theta)*horzbar[i] + cos(theta)*deepbar[i];
      }
    
    // Added calculation of camera center height above the ground (floor of building) - SBT 20150430
    // This are unique for each camera because the camera separations vary.

    if( (geo->ring[cam]-1) == 0)
	hring[geo->ring[cam]-1] = dR1 + dR2 + geo->sep3[cam]*cos(geo->mir_the[cam]);
    else
	hring[geo->ring[cam]-1] = dR1 + geo->sep3[cam]*cos(geo->mir_the[cam]);

    // main support poles, qty 6 (three per side), approximating actual side panels
    for (i=0; i<6; i++) {
      tgm->length[i] = phgt;
      tgm->radius[i] = prad;
      for (j=0; j<3; j++)
	tgm->pole[cam][i][j] = vertbar[j];
    }

    // horizontal poles, qty 14, in 7 front-and-back pairs
    for (i=6; i<20; i++) {
      tgm->length[i] = pwid;
      tgm->radius[i] = prad;
      for (j=0; j<3; j++)
	tgm->pole[cam][i][j] = horzbar[j];
    }
    
    // diagonal bottom: 20, 23
    for (i=20; i<26; i+=3) {
      /* 1.497414 = sqrt(0.707 * 0.707 + 1.320 * 1.320) */
      tgm->length[i] = 1.497414;
      tgm->radius[i] = prad;
      for (j=0; j<3; j++)
	tgm->pole[cam][i][j] = 0.707 * vertbar[j] + 1.320 * horzbar[j];
      unitVector(tgm->pole[cam][i], tgm->pole[cam][i]);
    }
    
    // diagonal middle: 21, 24
    for (i=21; i<26; i+=3) {                  
      /* 1.783928 = sqrt(1.200 * 1.200 + 1.320 * 1.320) */
      tgm->length[i] = 1.783928;
      tgm->radius[i] = pmid;
      for (j=0; j<3; j++)
      tgm->pole[cam][i][j] = 1.200 * vertbar[j] + 1.320 * horzbar[j];
      unitVector(tgm->pole[cam][i], tgm->pole[cam][i]);
    }

    // diagonal top: 22, 25
    for (i=22; i<26; i+=3) {
      /* 1.626315 = sqrt(0.950 * 0.950 + 1.320 * 1.320) */
      tgm->length[i] = 1.626315;
      tgm->radius[i] = pmid;
      for (j=0; j<3; j++)
	tgm->pole[cam][i][j] = 0.950 * vertbar[j] + 1.320 * horzbar[j];
      unitVector(tgm->pole[cam][i], tgm->pole[cam][i]);
    }

    // reverse diagonal bottom: 26, 29
    for (i=26; i<32; i+=3) {
      /* 1.497414 = sqrt(0.707 * 0.707 + 1.320 * 1.320) */
      tgm->length[i] = 1.497414;
      tgm->radius[i] = prad;
      for (j=0; j<3; j++)
      tgm->pole[cam][i][j] = -0.707 * vertbar[j] + 1.320 * horzbar[j];
      unitVector(tgm->pole[cam][i], tgm->pole[cam][i]);
    }

    // reverse diagonal middle: 27, 30
    for (i=27; i<32; i+=3) {                  
      /* 1.783928 = sqrt(1.200 * 1.200 + 1.320 * 1.320) */
      tgm->length[i] = 1.783928;
      tgm->radius[i] = pmid;
      for (j=0; j<3; j++)
      tgm->pole[cam][i][j] = -1.200 * vertbar[j] + 1.320 * horzbar[j];
      unitVector(tgm->pole[cam][i], tgm->pole[cam][i]);
    }

    // reverse diagonal top: 28, 31
    for (i=28; i<32; i+=3) {
      /* 1.626315 = sqrt(0.950 * 0.950 + 1.320 * 1.320) */
      tgm->length[i] = 1.626315;
      tgm->radius[i] = pmid;
      for (j=0; j<3; j++)
      tgm->pole[cam][i][j] = -0.950 * vertbar[j] + 1.320 * horzbar[j];
      unitVector(tgm->pole[cam][i], tgm->pole[cam][i]);
    }

    // flanking poles
    for (i=32; i<34; i++) { 
      tgm->length[i] = 4.552;
      tgm->radius[i] = pbig;    
      for (j=0; j<3; j++)
      tgm->pole[cam][i][j] = vertbar[j];
    }

    // horizontal bars to flanking pole
    for (i=34; i<38; i+=2) {
      /* 2.366479 = sqrt(1.265 * 1.265 + 2.000 * 2.000) */
      tgm->length[i] = 2.366479; 
      tgm->radius[i] = prad;
      for (j=0; j<3; j++)  
      tgm->pole[cam][i][j] = 1.265 * horzbar[j] + 2.000 * deepbar[j];
      unitVector(tgm->pole[cam][i], tgm->pole[cam][i]);
    }

    // horizontal bars to other flanking pole
    for (i=35; i<38; i+=2) {
      /* 2.366479 = sqrt(1.265 * 1.265 + 2.000 * 2.000) */
      tgm->length[i] = 2.366479; 
      tgm->radius[i] = prad;
      for (j=0; j<3; j++)  
      tgm->pole[cam][i][j] = -1.265 * horzbar[j] + 2.000 * deepbar[j];
      unitVector(tgm->pole[cam][i], tgm->pole[cam][i]);
    }

    // top flank connector
    i = 38;
    /* 3.195955 = sqrt(1.265 * 1.265 + 2.000 * 2.000 + 2.148 * 2.148) */
    tgm->length[i] = 3.195955;   
    tgm->radius[i] = prad;
    for (j=0; j<3; j++)    
      tgm->pole[cam][i][j] = 
      1.265 * horzbar[j] + 
      2.000 * deepbar[j] + 
      2.148 * vertbar[j];
    unitVector(tgm->pole[cam][i], tgm->pole[cam][i]);

    // other top flank connector
    i = 39;
    tgm->length[i] = tgm->length[i-1];
    tgm->radius[i] = prad;
    for (j=0; j<3; j++)
      tgm->pole[cam][i][j] = 
      -1.265 * horzbar[j] + 
       2.000 * deepbar[j] + 
       2.148 * vertbar[j];
    unitVector(tgm->pole[cam][i], tgm->pole[cam][i]);

    // opposite flank connector
    i = 40;              
    /* 3.332375 = sqrt(2.585 * 2.585 + 2.000 * 2.000 + 0.650 * 0.650) */
    tgm->length[i] = 3.332375;          
    tgm->radius[i] = prad;
    for (j=0; j<3; j++)
      tgm->pole[cam][i][j] = 
      2.585 * horzbar[j] + 
      2.000 * deepbar[j] + 
      0.650 * vertbar[j];
    unitVector(tgm->pole[cam][i], tgm->pole[cam][i]);

    // other opposite flank connector
    i = 41;              
    tgm->length[i] = tgm->length[i-1];    
    tgm->radius[i] = prad;
    for (j=0; j<3; j++)
      tgm->pole[cam][i][j] = 
      -2.585 * horzbar[j] + 
       2.000 * deepbar[j] + 
       0.650 * vertbar[j];
    unitVector(tgm->pole[cam][i], tgm->pole[cam][i]);

    /* New bar for wind baffle sheet */
    i = 42;
    tgm->length[i] = 4.714;
    tgm->radius[i] = 0.053;
    for (j=0; j<3; j++)tgm->pole[cam][i][j] = vertbar[j];

    /* Second new bar for wind baffle sheet - SBT 20140430 */ 

    i = 54;
    tgm->length[i] = 4.714;
    tgm->radius[i] = 0.053;
    for (j=0; j<3; j++)tgm->pole[cam][i][j] = vertbar[j];

    /* New outside bar for wind baffle sheet  - SBT 20140430*/

    i = 55;
    tgm->length[i] = 8.000;
    tgm->radius[i] = 0.051;
    for (j=0; j<3; j++)tgm->pole[cam][i][j] = vertbar[j];

    /* Ring 1 walkway */
    if (geo->ring[cam] == 1) {
      for (i=43; i<45; i++) {
      tgm->length[i] = 1.5785;
      tgm->radius[i] = wrad1;
      for (j=0; j<3; j++)
        tgm->pole[cam][i][j] = horzbar[j];
      }

      for (i=45; i<47; i++) {
      tgm->length[i] = 1.5785;
      tgm->radius[i] = wrad2;
      for (j=0; j<3; j++)
        tgm->pole[cam][i][j] = horzbar[j];
      }

      for (i=47; i<53; i++) {
      tgm->length[i] = 1.1049;
      tgm->radius[i] = wrad1;
      for (j=0; j<3; j++)
        tgm->pole[cam][i][j] = vertbar[j];
      }
    }
    else {
      for (i=43; i<53; i++) {
      tgm->length[i] = 0.;
      tgm->radius[i] = 0.;
      for (j=0; j<3; j++)
        tgm->pole[cam][i][j] = 0.;
      }
    }

    /* Reinforcement column stay (page MBF-40-BD-1106, p. 44 on wiki)
     * http://www-ta.icrr.u-tokyo.ac.jp/private/Personal/mtakeda/Docs/20060301MES_BRM.pdf
     */

    i=53;
    //    tgm->length[i] = 3.314;
    tgm->length[i] = 3.850;  // SBT 20150430 (was too short)
    tgm->radius[i] = 0.05;
    for (j=0; j<3; j++) 
      tgm->pole[cam][i][j] = horzbar[j];      
      
    for (i=0; i<3; i++) {

    // The reference point is the base of the right main support pole (#00) - SBT 20150430

    if(geo->ring[cam] == 1)
      tgm->base[cam][0][i]  = zcam[i] + 3.0838*deepbar[i] - (pwid/2.)*horzbar[i] - (dR1 + dR2)*vertbar[i];
    else
      tgm->base[cam][0][i]  = zcam[i] + 2.7153*deepbar[i] - (pwid/2.)*horzbar[i] - dR1*vertbar[i];

      tgm->base[cam][1][i]  = tgm->base[cam][0][i]  + (pdep/2.) * deepbar[i];
      tgm->base[cam][2][i]  = tgm->base[cam][1][i]  + (pdep/2.) * deepbar[i];

      tgm->base[cam][3][i]  = tgm->base[cam][0][i]  + pwid * horzbar[i];  // Left main support pole (#03)
      tgm->base[cam][4][i]  = tgm->base[cam][1][i]  + pwid * horzbar[i];
      tgm->base[cam][5][i]  = tgm->base[cam][2][i]  + pwid * horzbar[i];

      // see p. 38
      tgm->base[cam][6][i]  = tgm->base[cam][0][i]  + 0.707 * vertbar[i];
      tgm->base[cam][7][i]  = tgm->base[cam][6][i]  + 1.200 * vertbar[i];
      tgm->base[cam][8][i]  = tgm->base[cam][7][i]  + 1.570 * vertbar[i];
      tgm->base[cam][9][i]  = tgm->base[cam][8][i]  + 0.200 * vertbar[i];
      tgm->base[cam][10][i] = tgm->base[cam][9][i]  + 0.705 * vertbar[i];
      tgm->base[cam][11][i] = tgm->base[cam][10][i] + 0.950 * vertbar[i];
      tgm->base[cam][12][i] = tgm->base[cam][11][i] + 1.413 * vertbar[i];

      tgm->base[cam][13][i] = tgm->base[cam][6][i]  + pdep * deepbar[i];
      tgm->base[cam][14][i] = tgm->base[cam][7][i]  + pdep * deepbar[i];
      tgm->base[cam][15][i] = tgm->base[cam][8][i]  + pdep * deepbar[i];
      tgm->base[cam][16][i] = tgm->base[cam][9][i]  + pdep * deepbar[i];
      tgm->base[cam][17][i] = tgm->base[cam][10][i] + pdep * deepbar[i];
      tgm->base[cam][18][i] = tgm->base[cam][11][i] + pdep * deepbar[i];
      tgm->base[cam][19][i] = tgm->base[cam][12][i] + pdep * deepbar[i];

      tgm->base[cam][20][i] = tgm->base[cam][0][i];
      tgm->base[cam][21][i] = tgm->base[cam][6][i];
      tgm->base[cam][22][i] = tgm->base[cam][10][i];

      tgm->base[cam][23][i] = tgm->base[cam][20][i] + pdep * deepbar[i];
      tgm->base[cam][24][i] = tgm->base[cam][21][i] + pdep * deepbar[i];
      tgm->base[cam][25][i] = tgm->base[cam][22][i] + pdep * deepbar[i];

      tgm->base[cam][26][i] = tgm->base[cam][6][i];
      tgm->base[cam][27][i] = tgm->base[cam][7][i];
      tgm->base[cam][28][i] = tgm->base[cam][11][i]; 
      
      tgm->base[cam][29][i] = tgm->base[cam][26][i] + pdep * deepbar[i];
      tgm->base[cam][30][i] = tgm->base[cam][27][i] + pdep * deepbar[i];
      tgm->base[cam][31][i] = tgm->base[cam][28][i] + pdep * deepbar[i];

      tgm->base[cam][32][i] = tgm->base[cam][0][i]  - 1.265 * horzbar[i] 
                                                    - 2.000 * deepbar[i];
      tgm->base[cam][33][i] = tgm->base[cam][3][i]  + 1.265 * horzbar[i] 
                                                    - 2.000 * deepbar[i];

      tgm->base[cam][34][i] = tgm->base[cam][32][i] + 2.600 * vertbar[i]; 
      tgm->base[cam][35][i] = tgm->base[cam][33][i] + 2.600 * vertbar[i];

      tgm->base[cam][36][i] = tgm->base[cam][34][i] + 1.727 * vertbar[i];    
      tgm->base[cam][37][i] = tgm->base[cam][35][i] + 1.727 * vertbar[i];

      tgm->base[cam][38][i] = tgm->base[cam][36][i];
      tgm->base[cam][39][i] = tgm->base[cam][37][i];

      tgm->base[cam][40][i] = tgm->base[cam][38][i] - 0.180 * vertbar[i];
      tgm->base[cam][41][i] = tgm->base[cam][39][i] - 0.180 * vertbar[i];

      /* New bars for wind baffle sheet - updated SBT 20150430 */ 

      if (   cam == 0 || cam == 1 
          || cam == 4 || cam == 5
          || cam == 8 || cam == 9)
	{
	  tgm->base[cam][42][i] = zcam[i] + 4.7228*ybar[i] + 2.2043*xbar[i];
	  tgm->base[cam][54][i] = zcam[i] + 4.8307*ybar[i] + 2.2043*xbar[i];
	  tgm->base[cam][55][i] = zcam[i] + 5.2857*ybar[i] + 2.2043*xbar[i];
	}
      else
	{
	  tgm->base[cam][42][i] = zcam[i] + 4.7228*ybar[i] - 2.2043*xbar[i];
	  tgm->base[cam][54][i] = zcam[i] + 4.8307*ybar[i] - 2.2043*xbar[i];
	  tgm->base[cam][55][i] = zcam[i] + 5.2857*ybar[i] - 2.2043*xbar[i];
	}

      if(geo->ring[cam] == 1)
	{
	  tgm->base[cam][42][i] -= (dR1 + dR2)*vertbar[i];
	  tgm->base[cam][54][i] -= (dR1 + dR2)*vertbar[i];
	  tgm->base[cam][55][i] -= (dR1 + dR2)*vertbar[i];
	}
      else
	{
	  tgm->base[cam][42][i] -= dR1*vertbar[i];
	  tgm->base[cam][54][i] -= dR1*vertbar[i];
	  tgm->base[cam][55][i] -= dR1*vertbar[i];
	}
                                
      tgm->base[cam][43][i] = zcam[i] + 0.9271 * deepbar[i] - 
      0.3683 * vertbar[i] - 1.6574 * horzbar[i];
      tgm->base[cam][44][i] = zcam[i] + 0.9271 * deepbar[i] - 
        0.3683 * vertbar[i] + 0.0699 * horzbar[i];

      tgm->base[cam][45][i] = tgm->base[cam][43][i] - 0.5525 * vertbar[i];
      tgm->base[cam][46][i] = tgm->base[cam][44][i] - 0.5525 * vertbar[i];

      tgm->base[cam][47][i] = tgm->base[cam][43][i] - 1.1049 * vertbar[i];
      tgm->base[cam][48][i] = tgm->base[cam][47][i] + 0.7938 * horzbar[i];
      tgm->base[cam][49][i] = tgm->base[cam][48][i] + 0.7938 * horzbar[i];

      tgm->base[cam][50][i] = tgm->base[cam][44][i] - 1.1049 * vertbar[i];
      tgm->base[cam][51][i] = tgm->base[cam][50][i] + 0.7938 * horzbar[i];
      tgm->base[cam][52][i] = tgm->base[cam][51][i] + 0.7938 * horzbar[i];
      
      tgm->base[cam][53][i] = tgm->base[cam][32][i] + 3.8900 * vertbar[i];  
    }
  }
}
  

void buildMDCameraPoles (geofd_dst_common *geo, TGeom *tgm) {
  int cam, ring, i, j;

  double vertbar[3];                  // direction of vert. bars on clust stand
  double horzbar[3];                  // direction of horz. bars on clust stand
  double deepbar[3];                  // direction of deep  bars on clust stand
  double zclust[3];                   // location of front of cluster box

  double prad = 0.0222;               // support pole radius (meters)
  double pdep = 0.1524;               // support pole depth (meters)
  double crad = 0.0127;               // cross bar radius (meters)
  double clen = 1.4478;               // cross bar length (meters)
  double pwid[2] = {0.7906, 0.7938};  // separation between legs of
                                      // main support poles across face of
                                      // camera (meters)

  double phgt[2] = {1.7844, 2.2225};  // height of main support poles (meters)
  double poff[2] = {0.1264, 0.0833};  // offset of support poles from camera
                                      // center

  double c1off[2] = {0.4261, 0.9073}; // offset of back cross beams from camera
                                      // center
  double c2off[2] = {0.2991, 0.7803}; // offset of front cross beams from 
                                      // camera center

  double wrad = 0.0508;               // cable radius
  double wlen = 2.0000;               // cable length

  double sinang, cosang;

  for (cam=0; cam<geo->nmir; cam++) {
    ring = geo->ring[cam] - 1;

    sinang = sin(geo->mir_the[cam]);
    cosang = cos(geo->mir_the[cam]);

    vertbar[0] = 0.;
//     vertbar[1] = -sinang;  //LMS
    vertbar[1] =  sinang;  //TAS
    vertbar[2] =  cosang;

    horzbar[0] = 1.;
    horzbar[1] = 0.;
    horzbar[2] = 0.;

    deepbar[0] = 0.;
//     deepbar[1] =  cosang;  //LMS
    deepbar[1] = -cosang;  //TAS
    deepbar[2] =  sinang;

    zclust[0]  = 0.;
    zclust[1]  = 0.;
    zclust[2]  = -geo->rcurve[cam] + geo->sep[cam] + geo->cam_depth/2.;

    /* main support poles */
    for (i=0; i<6; i++) {
      tgm->length[i] = phgt[ring];
      tgm->radius[i] = prad;
      for (j=0; j<3; j++)
      tgm->pole[cam][i][j] = -vertbar[j];
    }

    /* back cross beams */
    i = 6;
    tgm->length[i] = clen;
    tgm->radius[i] = crad;
    for (j=0; j<3; j++)
      tgm->pole[cam][i][j] = -pwid[ring] * horzbar[j] - 1.2319 * vertbar[j];
    unitVector(tgm->pole[cam][i], tgm->pole[cam][i]);

    i = 7;
    tgm->length[i] = clen;
    tgm->radius[i] = crad;
    for (j=0; j<3; j++)
      tgm->pole[cam][i][j] =  pwid[ring] * horzbar[j] - 1.2319 * vertbar[j];
    unitVector(tgm->pole[cam][i], tgm->pole[cam][i]);

    /* front cross beams */
    i = 8;
    tgm->length[i] = clen;
    tgm->radius[i] = crad;
    for (j=0; j<3; j++)
      tgm->pole[cam][i][j] = -0.3366 * horzbar[j] - 
      1.3589 * vertbar[j] + 0.3556 * deepbar[j];
    unitVector(tgm->pole[cam][i], tgm->pole[cam][i]);

    i = 9;
    tgm->length[i] = clen;
    tgm->radius[i] = crad;
    for (j=0; j<3; j++)
      tgm->pole[cam][i][j] =  0.3366 * horzbar[j] - 
      1.3589 * vertbar[j] + 0.3556 * deepbar[j];
    unitVector(tgm->pole[cam][i], tgm->pole[cam][i]);

    i = 10;
    tgm->length[i] = wlen;
    tgm->radius[i] = wrad;
    for (j=0; j<3; j++)
      tgm->pole[cam][i][j] = vertbar[j];

    for (i=0; i<3; i++) {
      tgm->base[cam][0][i] = zclust[i] - (pwid[ring]/2.) * horzbar[i]
      + poff[ring] * vertbar[i];
      tgm->base[cam][1][i] = tgm->base[cam][0][i] +
      (0.0095 + 2.0*prad) * deepbar[i];
      tgm->base[cam][2][i] = tgm->base[cam][0][i] - 
      (0.0095 + 2.0*prad) * deepbar[i];
      tgm->base[cam][3][i] = tgm->base[cam][0][i] + pwid[ring] * horzbar[i];
      tgm->base[cam][4][i] = tgm->base[cam][3][i] + 
      (0.0095 + 2.0*prad) * deepbar[i];
      tgm->base[cam][5][i] = tgm->base[cam][3][i] - 
      (0.0095 + 2.0*prad) * deepbar[i];

      tgm->base[cam][6][i] = zclust[i] + (pwid[ring]/2.) * horzbar[i]
      - (pdep/2.) * deepbar[i]
      - c1off[ring] * vertbar[i];
      tgm->base[cam][7][i] = tgm->base[cam][6][i] - pwid[ring] * horzbar[i];

      tgm->base[cam][8][i] = zclust[i] + (pwid[ring]/2.) * horzbar[i]
      + (pdep/2.) * deepbar[i]
      - c2off[ring] * vertbar[i];
      tgm->base[cam][9][i] = tgm->base[cam][8][i] - pwid[ring] * horzbar[i];
    }
    tgm->base[cam][10][0] = 0.;
    tgm->base[cam][10][1] = 0.;
    tgm->base[cam][10][2] = -geo->rcurve[cam] + geo->sep[cam] + 
      geo->cam_depth;
  }
}

int hitCameraBox (geofd_dst_common *geo, TGeom *tgm, RayTrace *ray,
                  double v[], double w[]) {
  int i;
  double back[3];
  double front[3];

  /* Checking back and front planes of camera box */
  for (i=0; i<3; i++) {
    back[i]  = w[i] + ((geo->sep[ray->cam] - tgm->h[ray->cam]) / v[2]) * v[i];
    front[i] = w[i] + ((geo->sep[ray->cam] + 
                  geo->cam_depth - tgm->h[ray->cam]) / v[2]) * v[i];
  }

  if ( ( fabs(back[0])  < (geo->cam_width/2. ) &&
         fabs(back[1])  < (geo->cam_height/2.) ) ||
       ( fabs(front[0]) < (geo->cam_width/2. ) &&
         fabs(front[1]) < (geo->cam_height/2.) )
     )
    return TRUE;
  else
    return FALSE;
}


int hitFlasherBox (geofd_dst_common *geo, TGeom *tgm, RayTrace *ray,
                   double v[], double w[]) {
  int i;
  
  double uvled_sep = 0.056;
  double uvled_width = 0.117;
  double uvled_height = 0.155;
  double uvled_depth = 0.060;
  
  double back[3];
  double front[3];
  
  for (i=0; i<3; i++) {
    back[i]  = w[i] + ( (uvled_sep - tgm->h[ray->cam]) / v[2]) * v[i];
    front[i] = w[i] + ( (uvled_sep + uvled_depth - tgm->h[ray->cam]) / v[2]) * v[i];
  }

  if ( ( fabs(back[0])  < (uvled_width/2. ) &&
         fabs(back[1])  < (uvled_height/2.) ) ||
       ( fabs(front[0]) < (uvled_width/2.) &&
         fabs(front[1]) < (uvled_height/2.) )
     )
    return TRUE;
  else
    return FALSE;
}


int hitCameraPoles (int num, TGeom *tgm, RayTrace *ray, double u[], double v[]) {
  int i, j;
  double pu[3];                    // ray axis vector
  double pr[3];                    // vector from point on pole to
                                   // intersection with mirror     
  double pd[3];                    // unit vector pointing in direction
                                   // of closest approach of ray and pole

  double pudotpole, prdotpole, prdotpd, prdotpu, len;

  FILE *fp;
  static int PASS0[200];
  char fname[128];
  double val;

  for (i=0; i<3; i++)
    pu[i] = -v[i];
  unitVector (pu, pu);

#ifdef DUMPX3D
  if(PASS0[ray->cam] == 0)
    {
      sprintf(fname,"tacam%.2d.x3d", ray->cam);
      fp = fopen(fname,"w");

      fprintf(fp,"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
      fprintf(fp,"<!DOCTYPE X3D PUBLIC \"ISO//Web3D//DTD X3D 3.0//EN\"\n");
      fprintf(fp,"\"http://www.web3d.org/specifications/x3d-3.0.dtd\">\n");

      fprintf(fp,"<X3D version=\'3.0\' profile=\'Immersive\'>\n");
      fprintf(fp,"<head>\n");
      fprintf(fp,"<meta name=\'filename\' content=\'tacam%.2d.x3d\'/>\n", ray->cam);
      fprintf(fp,"</head>\n");
      fprintf(fp,"<Scene>\n");
      fprintf(fp,"<ProtoDeclare name=\'TwoColorTable\'>\n");
      fprintf(fp,"<ProtoInterface>\n");
      fprintf(fp,"<field name=\'textColor\' type=\'SFColor\' value=\'.1 .8 .3\' accessType=\'initializeOnly\'/>\n");
      fprintf(fp,"<field name=\'legColor\' type=\'SFColor\' value=\'.8 .4 .7\' accessType=\'initializeOnly\'/>\n");
      fprintf(fp,"</ProtoInterface>\n");
      fprintf(fp,"<ProtoBody>\n");
      fprintf(fp,"<Transform>\n");
      for (i=0; i<num; i++) 
	{
	  // Text
	  fprintf(fp,"<Transform ");
	  fprintf(fp,"rotation=\'");
	  fprintf(fp,"%lf 0.0 %lf ", tgm->pole[ray->cam][i][2], -tgm->pole[ray->cam][i][0]); // Cross product of y-axis and pole
	  fprintf(fp,"%lf\' ", acos(tgm->pole[ray->cam][i][1]));  // Rotation angle
	  fprintf(fp,"translation=\'");
	  for (j=0; j<3; j++)
	    {
	      val = tgm->base[ray->cam][i][j] + tgm->pole[ray->cam][i][j]*tgm->length[i]/2.0e0;
	      if(j==2) val += 2.0*tgm->radius[i];
	      fprintf(fp,"%lf ", val);
	    }
	  fprintf(fp,"\'>\n");
	  fprintf(fp,"<Shape DEF=\'Element%d_0\'>\n", i);
	  fprintf(fp,"<Appearance>\n");
	  fprintf(fp,"<Material DEF=\'LegMaterial%d_0\' diffuseColor=\'0.0 0.0 1.0\'>\n", i);
	  fprintf(fp,"<IS>\n");
	  fprintf(fp,"<connect nodeField=\'diffuseColor\' protoField=\'textColor\'/>\n");
	  fprintf(fp,"</IS>\n");
	  fprintf(fp,"</Material>\n");
	  fprintf(fp,"</Appearance>\n");
	  fprintf(fp,"<Text string=\'\"%.2d\"\'>\n", i);
	  fprintf(fp,"<FontStyle DEF=\'testFontStyle%d\' justify=\'\"MIDDLE\" \"MIDDLE\"\' size=\'%lf\'/>\n", i, 5.0*tgm->radius[i]);
	  fprintf(fp,"</Text>\n"); 
	  fprintf(fp,"</Shape>\n");
	  fprintf(fp,"</Transform>\n");

	  // Pole
	  fprintf(fp,"<Transform ");
	  fprintf(fp,"rotation=\'");
	  fprintf(fp,"%lf 0.0 %lf ", tgm->pole[ray->cam][i][2], -tgm->pole[ray->cam][i][0]); // Cross product of y-axis and pole
	  fprintf(fp,"%lf\' ", acos(tgm->pole[ray->cam][i][1]));  // Rotation angle
	  fprintf(fp,"translation=\'");
	  for (j=0; j<3; j++)fprintf(fp,"%lf ", tgm->base[ray->cam][i][j] + tgm->pole[ray->cam][i][j]*tgm->length[i]/2.0e0);
	  fprintf(fp,"\'>\n");
	  fprintf(fp,"<Shape DEF=\'Element%d_1\'>\n", i);
	  fprintf(fp,"<Appearance>\n");
	  fprintf(fp,"<Material DEF=\'LegMaterial%d_1\' diffuseColor=\'1.0 0.0 0.0\'>\n", i);
	  fprintf(fp,"<IS>\n");
	  fprintf(fp,"<connect nodeField=\'diffuseColor\' protoField=\'legColor\'/>\n");
	  fprintf(fp,"</IS>\n");
	  fprintf(fp,"</Material>\n");
	  fprintf(fp,"</Appearance>\n");
	  fprintf(fp,"<Cylinder height=\'%lf\' radius=\'%lf\' solid=\'true\'/>\n", tgm->length[i], tgm->radius[i]); 
	  fprintf(fp,"</Shape>\n");
	  fprintf(fp,"</Transform>\n");
	  
	}
      
      fprintf(fp,"</Transform>\n");
      fprintf(fp,"</ProtoBody>\n");
      fprintf(fp,"</ProtoDeclare>\n");
      fprintf(fp,"<ProtoInstance name=\'TwoColorTable\'>\n");
      fprintf(fp,"<fieldValue name=\'legColor\' value=\'1 0 0\'/>\n");
      fprintf(fp,"</ProtoInstance>\n");
      fprintf(fp,"<NavigationInfo type=\'\"EXAMINE\"\'/>\n");
      fprintf(fp,"</Scene>\n");
      fprintf(fp,"</X3D>\n");

      fclose(fp);
      PASS0[ray->cam] = 1;
    }
#endif

  for (i=0; i<num; i++) {
    for (j=0; j<3; j++)
      pr[j] = u[j] - tgm->base[ray->cam][i][j];

//     crossProduct (pu, tgm->pole[ray->cam][i], pd);
    crsprod (pu, tgm->pole[ray->cam][i], pd);
    unitVector (pd, pd);

//     pudotpole = dotProduct (pu, tgm->pole[ray->cam][i]);
//     prdotpole = dotProduct (pr, tgm->pole[ray->cam][i]);
//     prdotpd   = dotProduct (pr, pd);
//     prdotpu   = dotProduct (pr, pu);
    pudotpole = dotprod (pu, tgm->pole[ray->cam][i]);
    prdotpole = dotprod (pr, tgm->pole[ray->cam][i]);
    prdotpd   = dotprod (pr, pd);
    prdotpu   = dotprod (pr, pu);
    
    len = prdotpole/pudotpole - prdotpu;
    len /= (1./pudotpole - pudotpole);

    if ( fabs(prdotpd) < tgm->radius[i] && 
       len > 0. && len < tgm->length[i] )
      return TRUE;
  }

  return FALSE;
}
