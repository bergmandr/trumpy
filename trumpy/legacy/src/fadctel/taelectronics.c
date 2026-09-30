#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/mman.h>
#include "constants.h"
#include "control.h"
#include "event.h"
#include "fdconstants.h"
#include "random.h"
#include "toolbox.h"

#include "atmosphere-withdb.h"
#include "airshower.h"
#include "showerlib.h"
#include "eventlist.h"
#include "track.h"

#include "fdconstants.h"
#include "adjcamera.h"
#include "fdsite.h"
#include "tacalibration.h"
#include "calibration.h"
#include "raytrace.h"
#include "acpttrack.h"
#include "taelectronics.h"

#include "histogram.h"

#define SQR(X) ((X)*(X))

/*
 *  TF Trigger Type:
 *
 *  * Type '1' tells the simulation to scan time slices 0-255 for fired tubes
 *      to check for the triggering requirement.
 *  * Type '2' tells the simulation to scan time slices 128-383 for fired tubes
 *      to check for the triggering requirement.  In addition, at least one 
 *      tube must have fired before time slice 256.
 *
 */
#ifndef TFTRIG_TYPE
#define TFTRIG_TYPE 1
#endif

/*
 *  TF Trigger Version:
 *
 *  * Version '1' uses TAMA algorithm for including adjacent camera data w/o 
 *      trigger verification.  i.e. Only cameras with non-zero TF trigger code 
 *      and those adjacent to them are read out.
 *  * Version '2' uses TAMA algorithm for including adjacent camera data with 
 *      trigger verification.  i.e. Cameras with non-zero TF trigger code are 
 *      re-evaluated for trigger condition using data from subsequent frames 
 *      across the waveform.  The trigger code read out is 10*code + additional
 *      frame steps needed to get trigger.
 *
 */
#ifndef TFTRIG_VER
#define TFTRIG_VER  1
#endif

/*
 *  The SDF firmware has the capability of including 4- & 8- bin sums when 
 *    determining when a tube fires.  These are both typically disabled.  The 
 *    code for the 8-bin sums is included here, but is commented out with the 
 *    following compiler switch.  This may be overridden by including the 
 *    argument '-DENABLE_8BIN_SUM' when compiling.
 *
 */
//#define ENABLE_8BIN_SUM

#if ( TFTRIG_VER == 1 )
#define TFCODE1     1
#define TFCODE2     2
#elif ( TFTRIG_VER == 2 )
#define TFCODE1    10
#define TFCODE2    20
#endif

#if ( TFTRIG_TYPE == 1 )
#define TFTRIG_MASK 0xFFFFFFFF
#elif ( TFTRIG_TYPE == 2 )
#define TFTRIG_MASK 0xFFFF0000
#endif

#define PE_MIN 100
#define TRIGGER_CL 6.0

#define BIN_DT       25.0         /* nanoseconds (see constants.h)           */
#define RES          25           /* number of elements of time const. array */
#define DRES         25.0         /* double-precision version of resolution  */

static unsigned int trigpat5[1364];
static unsigned int trigpat3[82];
enum {
  LEV2_TOP, LEV2_BOT, LEV2_LT, LEV2_RT
};

#ifndef VOLTAGE_OFFSET
#  define VOLTAGE_OFFSET 0.282  /* in units of FADC */
#endif

#ifndef NOISE_AMPLITUDE
#  define NOISE_AMPLITUDE 0.9
#endif

/*
 * Best-fit values of the two parameters of the approximate one PE transfer 
 *   function.  The 'true' transfer function comes from the convolution of the 
 *   transfer functions from the shaper circuit (two-pole, 50ns RC const.), 
 *   the preamp (single pole, 1ns), and the output of the PMT pulse (5ns wide 
 *   Gaussian).  The function used here ignores the preamp TF and approximates 
 *   the Gaussian with a 5-pole filter (time constant 'alpha').  To 
 *   reconstruct this formula, convolude the functions,
 *
 *     f(t) = t*exp(-t/tau)/(tau)^2
 *
 *     g(t) = t^4*exp(-t/alpha)/(24*(alpha)^5)
 *
 * The result is a hot mess of an equation with the form,
 *
 *     F(t) = exp(-t/alpha)*\sum_{i} A_i * t^i + exp(-t/tau)*\sum_{i} B_i*t^i
 *
 *   Since this is a linear combination of exponentials with simple arguments, 
 *   a recursion relation can be constructed to propagate the TF values over 
 *   time (even though they are multiplied by order 'n' polynomials).
 *
 */
/*
#define TAU   50.92   // best fit results to Taketa's data file
#define ALPHA  5.39
*/
/*
#define TAU   56.00   // gives a power spectrum contour closer to what we 
#define ALPHA  5.00   //   see in the data
*/

/* the values of the one PE TF in BIN_DT/RES ns resolution */
static double tfa1[RES], tfa2[RES], tfa3[RES], tfa4[RES], tfa5[RES];
static double tft1[RES], tft2[RES];

/* constants used to propagate TF's using recursion relation */
static double tca1, tca2, tca3, tca4, tca5, tct1, tct2;

#ifdef HISTOGRAM_MODE
static double _getPedRMS(short *wf) {
  int i;
  double x = 0.00;
  double xx = 0.00;
  for ( i=0; i<NFADC_BINS; i++ ) {
    x += (double)wf[i];
    xx += (double)wf[i] * (double)wf[i];
  }

  x /= (double)NFADC_BINS;
  xx = sqrt(xx/(double)NFADC_BINS - x*x);

  return xx;
}

static int _getSampleMeanSdev(int ssize, short *wf, double *mean, double *sdev) {
  int i, j, x;
  double sumx, sumxx, nx;

  /* sample size must be a power of 2 */
  int shift = 0;
  while ( ssize > 1 ) {
    shift++;
    ssize /= 2;
  }
  ssize = 0x1 << shift;

  sumx = 0.00;
  sumxx = 0.00;
  for ( i=0; i<NFADC_BINS; i+=ssize ) {
    x = 0;
    for ( j=0; j<ssize; j++ )
      x += wf[i+j];

    sumx += (double)x;
    sumxx += (double)x * (double)x;
  }

  nx = (double)( NFADC_BINS / ssize );

  *mean = sumx / nx;
  *sdev = sqrt( (nx*sumxx - SQR(sumx))/(nx*(nx-1)) );

  return ssize;
}
#endif  /* #ifdef HISTOGRAM_MODE */

#define NOISE_ARY_LEN 512
static const double sampled_noise[NOISE_ARY_LEN] = {
  -0.32992554,  0.40359497, -0.53295898,  0.20819092, -0.17852783,  0.06814575, -0.51992798, -0.12612915,  
  -0.84960938,  0.18539429, -0.12292480,  0.48901367,  0.07690430,  0.64721680,  0.51708984,  0.66149902, 
  -0.13198853,  0.88992310, -0.05374146,  0.47933960, -0.01422119,  0.29129028, -0.06585693,  0.18765259, 
  -0.55706787,  0.13488770, -0.12832642,  0.12341309,  0.04449463, -0.20074463,  0.29931641,  0.49993896, 
  -0.03781128,  0.46688843, -0.15713501,  0.42956543, -0.05242920,  0.41293335, -0.04858398,  0.29449463,  
  -0.63302612,  0.15231323, -0.29190063,  0.12936401, -0.03640747,  0.09667969,  0.30114746,  0.41021729, 
  -0.03826904,  0.50457764, -0.08670044,  0.47583008, -0.15814209,  0.26156616, -0.28771973,  0.09243774, 
  -0.70980835, -0.06707764, -0.26351929,  0.00274658, -0.26828003, -0.31420898,  0.17324829,  0.42263794, 
  -0.22229004,  0.40539551, -0.21194458,  0.02835083, -0.20513916,  0.28707886, -0.17984009,  0.25558472, 
  -0.52694702,  0.06817627, -0.28952026,  0.08825684, -0.07141113, -0.11590576,  0.26477051,  0.43377686, 
  -0.26336670,  0.40695190, -0.10159302,  0.25158691, -0.15872192,  0.21035767, -0.31597900,  0.04702759, 
  -0.69015503, -0.06622314, -0.30859375,  0.08834839, -0.08065796, -0.16763306,  0.16705322,  0.32958984, 
  -0.36791992,  0.33267212, -0.22146606,  0.16677856, -0.18322754,  0.32217407, -0.24539185,  0.07617188, 
  -0.65124512, -0.01202393, -0.28457642,  0.08157349, -0.13540649, -0.23565674,  0.14462280,  0.35293579, 
  -0.29995728,  0.36297607, -0.17749023,  0.13690186, -0.23123169,  0.16445923, -0.30471802,  0.05029297, 
  -0.70816040, -0.11050415, -0.43716431, -0.06015015, -0.21807861, -0.22991943,  0.19042969,  0.35415649, 
  -0.34149170,  0.33993530, -0.59893799,  0.23220825, -0.09603882,  0.13769531, -0.52529907, -0.20285034, 
  -0.93505859,  0.13723755, -0.08334351,  0.52148438,  0.06912231,  0.59167480,  0.48736572,  0.67544556, 
  -0.06225586,  0.95294189, -0.03555298,  0.45117188, -0.05334473,  0.27691650, -0.04837036,  0.21493530, 
  -0.57150269,  0.06893921, -0.18136597,  0.10443115,  0.09979248, -0.10745239,  0.36038208,  0.46661377,  
  -0.09304810,  0.46075439, -0.12973022,  0.48028564, -0.02557373,  0.39129639, -0.09646606,  0.26843262, 
  -0.61917114,  0.20257568, -0.20501709,  0.21047974, -0.00421143,  0.09158325,  0.30007935,  0.38717651, 
  -0.03234863,  0.50021362, -0.11621094,  0.44674683, -0.17449951,  0.26400757, -0.28924561,  0.09594727, 
  -0.69021606, -0.04736328, -0.26977539,  0.00845337, -0.24374390, -0.31158447,  0.18283081,  0.36727905, 
  -0.31243896,  0.37869263, -0.14639282,  0.09912109, -0.21057129,  0.23245239, -0.24276733,  0.26162720,  
  -0.45578003,  0.11831665, -0.24884033,  0.07916260, -0.18133545, -0.21652222,  0.25366211,  0.48336792, 
  -0.20080566,  0.39816284, -0.18295288,  0.20608521, -0.12283325,  0.28115845, -0.26962280,  0.01654053, 
  -0.75289917, -0.10839844, -0.33624268,  0.08807373, -0.07827759, -0.18951416,  0.13021851,  0.30993652, 
  -0.33892822,  0.40090942, -0.15005493,  0.17626953, -0.23489380,  0.25875854, -0.27850342,  0.09939575, 
  -0.62341309, -0.03048706, -0.35415649,  0.00463867, -0.15280151, -0.18933105,  0.23074341,  0.42779541,  
  -0.32104492,  0.26455688, -0.26821899,  0.13507080, -0.17880249,  0.21627808, -0.31091309, -0.01187134, 
  -0.79043579, -0.14703369, -0.40460205,  0.00741577, -0.13092041, -0.18118286,  0.15258789,  0.27999878, 
  -0.38848877,  0.34796143, -0.56518555,  0.22323608, -0.14514160,  0.06619263, -0.56185913, -0.18106079, 
  -0.88342285,  0.20111084, -0.06970215,  0.48703003,  0.01965332,  0.60452271,  0.55566406,  0.70462036, 
  -0.12869263,  0.84902954, -0.09082031,  0.49224854,  0.01721191,  0.30181885, -0.07403564,  0.15573120,  
  -0.59844971,  0.10543823, -0.12103271,  0.15869141,  0.07736206, -0.18066406,  0.27908325,  0.43911743, 
  -0.05349731,  0.51254272, -0.11407471,  0.43994141, -0.06002808,  0.39773560, -0.04696655,  0.30902100, 
  -0.61593628,  0.15463257, -0.26116943,  0.17291260, -0.02276611,  0.10543823,  0.29174805,  0.37799072, 
  -0.05938721,  0.52090454, -0.06140137,  0.50845337, -0.14141846,  0.24819946, -0.30987549,  0.10803223, 
  -0.67657471, -0.02929688, -0.28524780, -0.04296875, -0.29071045, -0.30758667,  0.22427368,  0.44558716, 
  -0.25659180,  0.36755371, -0.23852539,  0.01950073, -0.18856812,  0.28640747, -0.19326782,  0.25881958, 
  -0.50787354,  0.08108521, -0.24533081,  0.12557983, -0.30422974,  0.00177002,  0.37750244,  0.50216675, 
  -0.22830200,  0.42385864, -0.07543945,  0.29998779, -0.04940796,  0.32241821, -0.31759644,  0.00161743, 
  -0.68737793, -0.02224731, -0.28820801,  0.12017822, -0.09490967, -0.19296265,  0.17675781,  0.38534546, 
  -0.27954102,  0.36895752, -0.20715332,  0.15283203, -0.15997314,  0.36102295, -0.20840454,  0.08676147, 
  -0.65051270, -0.01809692, -0.29360962,  0.07864380, -0.11337280, -0.16970825,  0.19876099,  0.38073730, 
  -0.31427002,  0.31796265, -0.20532227,  0.16778564, -0.18246460,  0.22488403, -0.28555298,  0.03747559, 
  -0.72891235, -0.12719727, -0.40081787,  0.01504517, -0.12734985, -0.17059326,  0.16763306,  0.29647827, 
  -0.37329102,  0.36648560, -0.51760864,  0.27194214, -0.05499268,  0.12948608, -0.54501343, -0.18695068, 
  -0.87322998,  0.19076538, -0.08056641,  0.48770142,  0.05368042,  0.61349487,  0.59210205,  0.76339722,  
  -0.03485107,  0.91589355, -0.06265259,  0.45300293, -0.01605225,  0.32745361, -0.02325439,  0.20971680, 
  -0.55761719,  0.10281372, -0.12390137,  0.17523193,  0.14291382, -0.04769897,  0.44934082,  0.55603027, 
   0.01062012,  0.56066895, -0.06127930,  0.57443237,  0.11337280,  0.55462646,  0.04113770,  0.34594727, 
  -0.53018188,  0.33270264, -0.05487061,  0.35018921,  0.07940674,  0.17724609,  0.38742065,  0.53393555, 
   0.11331177,  0.64901733,  0.00848389,  0.53027344, -0.09857178,  0.36413574, -0.14413452,  0.27239990,  
  -0.56698608,  0.02490234, -0.20895386,  0.09130859, -0.13500977, -0.18026733,  0.31533813,  0.51019287, 
  -0.19033813,  0.46615601, -0.12277222,  0.15097046, -0.07250977,  0.40524292, -0.11535645,  0.33782959, 
  -0.43612671,  0.17315674, -0.14889526,  0.19937134,  0.21582031, -0.05426025,  0.33303833,  0.51324463, 
  -0.19757080,  0.43542480, -0.12438965,  0.28500366, -0.08166504,  0.27487183, -0.31723022, -0.11419678, 
  -0.86175537, -0.20080566, -0.51019287, -0.11468506, -0.29913330, -0.42660522,  0.07046509,  0.28182983, 
  -0.42050171,  0.29476929, -0.24371338,  0.08990479, -0.23928833,  0.27444458, -0.26263428,  0.03933716, 
  -0.72299194, -0.10058594, -0.36721802,  0.01641846, -0.15319824, -0.20455933,  0.19952393,  0.37567139, 
  -0.35717773,  0.21630859, -0.36093140,  0.04632568, -0.24548340,  0.17843628, -0.36199951, -0.09338379, 
  -0.85858154, -0.20367432, -0.42864990, -0.00274658, -0.20028687, -0.30871582,  0.02239990,  0.02239990
};

static int _cluster(int nt, int *tl0, int neighbor[NTUBES_CAMERA][NTUBES_CAMERA]) {
  int i, j, k, it, jt;

  int nl[NTUBES_CAMERA], tl[NTUBES_CAMERA];
  int nmax = 0;

  for ( i=0; i<nt; i++ )
    tl[i] = tl0[i];

  int nn = 0;
  while ( nt > 0 ) {
    nl[nn++] = tl[0];
    nt--;
    for ( i=0; i<nt; i++ )
      tl[i] = tl[i+1];

    for ( i=0; i<nn; i++ ) {
      it = nl[i];
      for ( j=0; j<nt; j++ ) {
	jt = tl[j];
	if ( neighbor[it][jt] ) {
	  nl[nn++] = tl[j];
	  nt--;
	  for ( k=j--; k<nt; k++ )
	    tl[k] = tl[k+1];
	}
      }
    }

    if ( nn > nmax )
      nmax = nn;

    nn = 0;
  }

  return nmax;
}

static void _initTriggerPatterns(void) {
  int i, j, r, c;
  int neighbor[NTUBES_CAMERA][NTUBES_CAMERA] = {{ 0 }};

  for (c=0; c<16; c++) {
    for (r=0; r<16; r++) {
      i = 16*c + r;

      /* Next in column */
      if (r != 15) {
        j = i+1;
        neighbor[i][j] = 1;
        neighbor[j][i] = 1;
      }

      /* Next in row */
      if (i < 240) {
        j = i+16;
        neighbor[i][j] = 1;
        neighbor[j][i] = 1;
      }

      /* Add neighbors in next column for only odd rows */
      if (r%2 == 1) {
        if (r != 0 && i<240) { // Not necessary, make easy to switch odd/even
          j = i+15;
          neighbor[i][j] = 1;
          neighbor[j][i] = 1;
        }
        if (r != 15 && i<239) {
          j = i+17;
          neighbor[i][j] = 1;
          neighbor[j][i] = 1;
        }
      }

    }
  }

  int p1, p2, p3, p4, p5, n, icnt;
  int tlist[NTUBES_CAMERA];

  /*
   *  do 5x5 combos first
   */
  icnt = 0;
  for ( p1=0; p1<21; p1++ ) {
    tlist[0] = 16*(p1%5) + p1/5;
    for ( p2=p1+1; p2<22; p2++ ) {
      tlist[1] = 16*(p2%5) + p2/5;
      for ( p3=p2+1; p3<23; p3++ ) {
	tlist[2] = 16*(p3%5) + p3/5;
	for ( p4=p3+1; p4<24; p4++ ) {
	  tlist[3] = 16*(p4%5) + p4/5;
	  for ( p5=p4+1; p5<25; p5++ ) {
	    tlist[4] = 16*(p5%5) + p5/5;
	    n = _cluster(5, tlist, neighbor);
	    if ( n == 5 )
	      trigpat5[icnt++] = (0x1<<p1) + (0x1<<p2) + (0x1<<p3) + 
		(0x1<<p4) + (0x1<<p5);
	  }
	}
      }
    }
  }
  pout("initialized %d 5-tube trigger patterns.\n", icnt);

  /*
   *  now do 4x4 combos
   */
  icnt = 0;
  for ( p1=0; p1<14; p1++ ) {
    tlist[0] = 16*(p1%4) + p1/4;
    for ( p2=p1+1; p2<15; p2++ ) {
      tlist[1] = 16*(p2%4) + p2/4;
      for ( p3=p2+1; p3<16; p3++ ) {
	tlist[2] = 16*(p3%4) + p3/4;
	n = _cluster(3, tlist, neighbor);
	if ( n == 3 )
	  trigpat3[icnt++] = (0x1<<p1) + (0x1<<p2) + (0x1<<p3);
      }
    }
  }
  pout("initialized %d 3-tube trigger patterns.\n", icnt);

}

static unsigned int _getTriggerMask(int mean, int vari, short *fadc) {
  int i, j, k, hsum, sum, sig;
  unsigned int mask;
  int thld;
  int fend = 3 * HALF_FRAME_SIZE;

#ifdef ENABLE_8BIN_SUM
  /* this part specifically for checking 8-bin set */
  thld = 36 * vari / 2;
  mask = 0xFFFFFFFF;
  hsum = 0;
  for ( j=HALF_FRAME_SIZE-4; j<fstart; j++ )
    hsum += fadc[j];

  sum = hsum;
  for ( j=HALF_FRAME_SIZE; j<fend; j+=4 ) {
    for ( k=0; k<4; k++ )
      sum += fadc[j+k];

    sig = ( sum - mean/2 ) * ( sum - mean/2 );
    if ( (sig>=thld) && (sum>mean/2) )
      return mask;

    mask = 0xFFFFFFFF >> (j/8);
    sum -= hsum;
    hsum = sum;
  }
#endif

  /* this part scans the 16, 32, 64, & 128-bin sums */
  for ( i=8; i<FRAME_SIZE; i*=2 ) {
    thld = 36*vari;
    mask = 0xFFFFFFFF;
    hsum = 0;
    for ( j=HALF_FRAME_SIZE-i; j<HALF_FRAME_SIZE; j++ )
      hsum += fadc[j];

    sum = hsum;
    for ( j=HALF_FRAME_SIZE; j<fend; j+=i ) {
      for ( k=0; k<i; k++ )
        sum += fadc[j+k];

      sig = ( sum - mean ) * ( sum - mean );
      if ( (sig>=thld) && (sum>mean) )
        return mask;

      mask >>= i/8;
      sum -= hsum;
      hsum = sum;
    }

    mean *= 2;
    vari *= 2;
  }

  return 0x0;
}

static int _isLevel1Trigger(unsigned int *pm) {
  int i, j, it, tc, ir, ic, id;
  unsigned int wm, tm=0x80000000;

  /*
   *  Scan 5x5 windows for clusters of 5 adjacent, coincident tubes
   */
  for ( i=0; i<32; i++ ) {
    for ( tc=0; tc<188; tc++ ) {
      if ( tc%16 == 12 )
	tc += 4;

      wm = 0x00000000;
      for ( it=0; it<25; it++ ) {
	ir = tc%16 + it/5;
	ic = tc/16 + it%5;
	id = 16*ic + ir;
	if ( (pm[id]&tm) != 0 ) {
	  if ( tc%2 == 1 && ir%2 == 0 )
	    wm |= 0x1<<(it-1);
	  else
	    wm |= 0x1<<it;
	}
      }

      for ( j=0; j<1364; j++ )
	if ( (wm&trigpat5[j]) == trigpat5[j] )
	  return 1;
    }

    tm >>= 1;
  }

  return 0;
}

static int _isLevel2Trigger(unsigned int *pm, int side) {
  int i, j, it, tc, ir, ic, id;
  int ts, te, dt;
  unsigned int wm, tm=0x80000000;

  switch ( side ) {
  case LEV2_TOP:
    ts = 16;
    te = 192;
    dt = 16;
    break;

  case LEV2_BOT:
    ts = 28;
    te = 204;
    dt = 16;
    break;

  case LEV2_LT:
    ts = 192;
    te = 205;
    dt = 1;
    break;

  case LEV2_RT:
    ts = 0;
    te = 13;
    dt = 1;
    break;

  default:
    perr("invalid side (%d), use one of: LEV2_TOP, LEV2_BOT, LEV2_RT, or "
	 "LEV2_LT\n", side);
    return -1;
    ;;
  }

  /*
   *  Scan 4x4 windows for clusters of 3 adjacent, coincident tubes
   */
  for ( i=0; i<32; i++ ) {
    for ( tc=ts; tc<te; tc+=dt ) {
      wm = 0x00000000;
      for ( it=0; it<16; it++ ) {
	ir = tc%16 + it/4;
	ic = tc/16 + it%4;
	id = 16*ic + ir;
	if ( (pm[id]&tm) != 0 ) {
	  if ( tc%2 == 1 && ir%2 == 0 )
	    wm |= 0x1<<(it-1);
	  else
	    wm |= 0x1<<it;
	}
      }

      for ( j=0; j<82; j++ )
	if ( (wm&trigpat3[j]) == trigpat3[j] )
	  return 1;
    }

    tm >>= 1;
  }

  return 0;
}

#ifdef HISTOGRAM_MODE
static double _fadcSig(int mean, int var, int nwf, short *wf) {
  int i, sum16, sum32, smax;

  /* initialize the 16-bin sum */
  sum16 = 0;
  for ( i=0; i<16; i++ )
    sum16 += wf[i];

  /* initialize the 32-bin sum, while scanning the first of the 16-bin sums */
  smax = sum16;
  sum32 = sum16;
  for ( i=16; i<32; i++ ) {
    sum32 += wf[i];

    sum16 += wf[i] - wf[i-16];
    if ( sum16 > smax )
      smax = sum16;
  }

  /* scan the rest of the fadc trace for the largest sum */
  for ( i=32; i<nwf; i++ ) {
    sum16 += wf[i] - wf[i-16];
    sum32 += wf[i] - wf[i-32];

    if ( sum16 > smax )
      smax = sum16;

    if ( sum32 > 2*smax )
      smax = sum32 / 2;
  }

  return (double)( smax - mean ) / sqrt( (double)var );
}
#endif

/*
 * Initializes the waveform generator.
 *
 * Author: Sean R. Stratton
 *         Rutgers University, Dept. of Physics & Astronomy
 *
 * no arguments.
 *
 * returns:
 *   zero upon successful completion.
 *
 */
int initWaveformGenerator(double tau, double alpha) {
  int i;
  double dt, t, t2, t3, t4, kappa, gamma;
  double C1, C2, C3, C4, C5, C6, C7;

  /* construct the constants of the one PE transfer function */
  kappa = 24.0 * alpha*alpha*alpha * ( tau - alpha ) * ( tau - alpha );
  gamma = alpha * tau / ( tau - alpha );
  C5 = 1.0 / kappa;
  C4 = C5 * 4.0 * gamma * 2.00 / 1.00;
  C3 = C4 * 3.0 * gamma * 3.00 / 2.00;
  C2 = C3 * 2.0 * gamma * 4.00 / 3.00;
  C1 = C2 * 1.0 * gamma * 5.00 / 4.00;

  C7 = 24.0 * gamma*gamma*gamma / kappa;
  C6 = -C7 * 5.0 * gamma; 

  dt = BIN_DT/DRES;
  /* t is time from pe to measurement with BIN_DT/RES ns resolution */
  t = dt/2.00;
  for ( i=0; i<RES; i++ ) {
    t2 = t * t;
    t3 = t * t2;
    t4 = t * t3;

    tfa1[i] =   C5                                             *exp(-t/alpha);
    tfa2[i] = ( C4 + 4.0*C5*t )                                *exp(-t/alpha);
    tfa3[i] = ( C3 + 3.0*C4*t + 6.0*C5*t2 )                    *exp(-t/alpha);
    tfa4[i] = ( C2 + 2.0*C3*t + 3.0*C4*t2 + 4.0*C5*t3 )        *exp(-t/alpha);
    tfa5[i] = ( C1 +     C2*t +     C3*t2 +     C4*t3 + C5*t4 )*exp(-t/alpha);

    tft1[i] =   C7          * exp(-t/tau);
    tft2[i] = ( C6 + C7*t ) * exp(-t/tau);

    t += dt;
  }

  /* tc's are used to propogate the transfer function over time */
  tca1 = exp(-BIN_DT/alpha);
  tca2 = BIN_DT * tca1;
  tca3 = BIN_DT * tca2;
  tca4 = BIN_DT * tca3;
  tca5 = BIN_DT * tca4;

  tct1 = exp(-BIN_DT/tau);
  tct2 = BIN_DT * tct1;

  return 0;
}

/*
 * Simulates the FADC response to a given signal.  This incorporates a 
 *   recursion relation method used by 'new_pe' from HiRes2 MC.  It can be 
 *   shown that a transfer function comprised of a linear combination of 
 *   exponentials with simple arguments (i.e. exp(-t/RC)) multiplied by any 
 *   polynomial can be propagated in time by multiplying itself by a constant 
 *   and adding its lower order terms multiplied by constants.  This is 
 *   explained in much detail on the TA Wiki.
 *
 * Author: Sean R. Stratton
 *         Rutgers University, Dept. of Physics & Astronomy
 *
 * arguments:
 *   double *npe -- array of length 'len' containing numbers of npe weights
 *   double *tpe -- array of length 'len' containing npe arrival times (ns)
 *   double tref -- reference time to start of FADC bin '0'.
 *   int len     -- number of elements of the given arrays
 *   int ped     -- tube FADC pedestal
 *   double nsbg -- night sky background level in npe/25ns
 *   double gain -- tube gain in units of FADC counts
 *   int *fadc   -- the simulated FADC response to the given signal
 *
 * returns:
 *   zero upon sucessful completion.  Fills the array 'fadc' with the 
 *     generated waveform.
 *
 * CAVEAT:
 *   In reality, the FADC digitizes the voltage every 25ns, and the total 
 *     number of FADC is the sum over four digitizations.  To account for 
 *     intrinsic noise, the digitization here is performed after summing the 
 *     voltages every 25ns.  That is, if $\hat{D}$ is the "digitization 
 *     operator", then this function does,
 *
 *       FADC = \hat{D}\left[ \sum_{i=1}^{4} V_i \right]
 *
 *     instead of,
 *
 *       FADC = \sum_{i=1}^{4} \hat{D}\left[ V_i \right]
 *
 *     This is an unresolved issue as of 7.17.08 which should be addressed.
 */
int generateWaveform(double *npe, double *tpe, double tref, int len, 
		     int ped, 
		     double nsbg, 
		     double gain, 
		     int nfadc, short *fadc, WaveformState *st) {
  //  const int ixx=532521516;  /* VERY important constant */

  int i, j, it, idx, ifadc;
  double u, expnsbg;
  double t, tnext;
  double va1, va2, va3, va4, va5, vt1, vt2, vout;

  /*
   * This prevents the while loops from having to constantly execute the 
   *   'exp' math function.
   */
#ifdef NO_NOISE  
  nsbg = 0.; // remove this line
#endif
  
  expnsbg = exp(-nsbg/4.0);

  va1 = st->va[0];
  va2 = st->va[1];
  va3 = st->va[2];
  va4 = st->va[3];
  va5 = st->va[4];

  vt1 = st->vt[0];
  vt2 = st->vt[1];

  /* simulate the waveform */
  idx = 0;
  tnext = BIN_DT;
  for ( i=0; i<nfadc; i++ ) {
    vout = 0.00;

    for ( j=0; j<4; j++ ) {
      /* update readout voltage terms */
      va5 = tca1*va5 +     tca2*va4 +     tca3*va3 +     tca4*va2 + tca5*va1;
      va4 = tca1*va4 + 2.0*tca2*va3 + 3.0*tca3*va2 + 4.0*tca4*va1;
      va3 = tca1*va3 + 3.0*tca2*va2 + 6.0*tca3*va1;
      va2 = tca1*va2 + 4.0*tca2*va1;
      va1 = tca1*va1;

      vt2 = tct1*vt2 + tct2*vt1;
      vt1 = tct1*vt1;

      /* add noise pe */
      /* this loop simulates Poisson statistics */
      u = RANDOM_NUMBER;
      while ( u >= expnsbg ) {
        it = (int)( RANDOM_NUMBER * DRES );
        va1 += tfa1[it];
        va2 += tfa2[it];
        va3 += tfa3[it];
        va4 += tfa4[it];
        va5 += tfa5[it];

        vt1 += tft1[it];
        vt2 += tft2[it];

        u *= RANDOM_NUMBER;
      }

      /* add signal pe */
      for ( ; idx<len; idx++ ) {  // idx initialized to 0 at start of function
        t = tpe[idx] - tref;
        if ( t > tnext )
          break;

        it = (int)( ( tnext - t ) * DRES/BIN_DT );
        if ( it < RES ) {
          va1 += npe[idx] * tfa1[it];
          va2 += npe[idx] * tfa2[it];
          va3 += npe[idx] * tfa3[it];
          va4 += npe[idx] * tfa4[it];
          va5 += npe[idx] * tfa5[it];

          vt1 += npe[idx] * tft1[it];
          vt2 += npe[idx] * tft2[it];
        }
      } 

      vout += va5 + vt2;
      tnext += BIN_DT;
    }

    //    ifadc = (int)( gain*BIN_DT*vout + 0.9*sampled_noise[st->it] + VOLTAGE_OFFSET ) + ped;
    ifadc = (int)( gain*BIN_DT*vout + NOISE_AMPLITUDE*sampled_noise[st->it] + VOLTAGE_OFFSET ) + ped;
    fadc[i] = (short)( (ifadc<0x3FFF) ? ifadc : 0x3FFF );

    st->it = ( st->it + 1 ) % NOISE_ARY_LEN;
  }

  st->va[0] = va1;
  st->va[1] = va2;
  st->va[2] = va3;
  st->va[3] = va4;
  st->va[4] = va5;

  st->vt[0] = vt1;
  st->vt[1] = vt2;

  return 0;
}

static short *mbuf = NULL;
static short *fadcbuf[fdraw_nmir_max][NTUBES_CAMERA];

void freeFADCBuffer() {
  size_t size = fdraw_nmir_max * NTUBES_CAMERA *
    (NFADC_BINS+HALF_FRAME_SIZE) * sizeof(short);
  munmap(mbuf, size);
}


void initTAElectronics(void) {
  pout("initializing electronics simulation package.\n");
  _initTriggerPatterns();
  initWaveformGenerator(TAU, ALPHA);
  size_t size = fdraw_nmir_max * NTUBES_CAMERA *
    (NFADC_BINS+HALF_FRAME_SIZE) * sizeof(short);
  mbuf = mmap(NULL, size, (PROT_READ|PROT_WRITE),
             (MAP_PRIVATE|MAP_ANONYMOUS), -1, 0);
  int i, j;
  for ( i=0; i<fdraw_nmir_max; i++ )
    for ( j=0; j<NTUBES_CAMERA; j++ )
      fadcbuf[i][j] = mbuf + (NFADC_BINS+HALF_FRAME_SIZE)*(NTUBES_CAMERA*i+j);
}

int simTATriggerResponse(TACalibration *cb, PETimes *pt, fdraw_dst_common *fdraw) {
  int i, j, k, f, nf, ic, nc;
  int petot, pemax;
  double t0, tref, dt;

  WaveformState st[fdraw_nmir_max][NTUBES_CAMERA];
//   short fadcbuf[fdraw_nmir_max][NTUBES_CAMERA][NFADC_BINS+HALF_FRAME_SIZE] = {
//     {{0}}
//   }; /* 4MB */
  size_t size = fdraw_nmir_max * NTUBES_CAMERA *
    (NFADC_BINS+HALF_FRAME_SIZE) * sizeof(short);
  memset(mbuf, 0x00, size);
  unsigned int sdfmask[fdraw_nmir_max][NTUBES_CAMERA], tfmask;
  double nsbg[fdraw_nmir_max][NTUBES_CAMERA];

  int trigmir[NCAMERAS_SITE] = { 0 };
  int tftime = 0;

#ifdef HISTOGRAM_MODE
  short *wf;
  float nnoise_tubes = 0.0f;
  double omean, osdev, cl, lrms;
#endif  /* #ifdef HISTOGRAM_MODE */

  /*
   *  Select reference time at random between 3/2 & 2 frame lengths before 
   *    arrival of 1st PE.
   */
#ifndef NO_NOISE
  t0 = pt->t1 - ceil((double)HALF_FRAME_SIZE*100.0*(3.0+RANDOM_NUMBER));
#else
  t0 = pt->t1 - HALF_FRAME_SIZE*300.0;
#endif
  nf = (int)( ( pt->t2 - t0 ) / ((double)HALF_FRAME_SIZE*100.0) ) + 1;

  /*
   *  Count PE in each frame to see if we could possibly have a trigger
   */
  tref = t0;
  pemax = 0;
  for ( f=0; f<nf; f++ ) {
    petot = 0;

    for ( i=0; i<pt->nmir; i++ ) {
      ic = pt->mir[i];
      for ( j=0; j<NTUBES_CAMERA; j++ ) {
        for ( k=0; k<pt->len[i][j]; k++ ) {
          dt = pt->t[i][j][k] - tref;
          if ( dt >= 0.00 && dt < ((double)FRAME_SIZE*100.0) )
            petot += (int)pt->n[i][j][k];
        }
      }
    }

    if ( petot > pemax )
      pemax = petot;

    tref += (double)HALF_FRAME_SIZE*100.0;
  }

  if ( pemax < PE_MIN ) {
    pout("trigger FAILED.  Too few PE : %d\n", pemax);
    return 0;
  }

  /*
   *  For mirrors with any PE, initialize the sky noise and 1st waveform
   */
  tref = t0 - (double)HALF_FRAME_SIZE*100.0;
  for ( i=0; i<pt->nmir; i++ ) {
    ic = pt->mir[i];
    for ( j=0; j<NTUBES_CAMERA; j++ ) {
      if ( cb->liveflag[ic][j] ) {
      /*	nsbg[i][j] = NSBG_ALPHA*SQR(cb->mean[ic][j]-16*cb->pedestal[ic][j]) / */
      /*	  (16.*cb->vari[ic][j]); */
        nsbg[i][j] = (double)cb->vari[ic][j] * NSBG_FUDGE * RHO_ALPHA /
          ( 16.0 * cb->pmtgain[ic][j] * cb->pmtgain[ic][j] );

        /* initialize the waveform state */
        st[i][j].it = 3 * HALF_FRAME_SIZE;
        for ( k=0; k<5; k++ )
          st[i][j].va[k] = 0.00;
        for ( k=0; k<2; k++ )
          st[i][j].vt[k] = 0.00;

        /* initialize the sky noise and first waveform */
        if ( cb->pedestal[ic][j] < 0 ) {
          perr("pedestal for %d %d < 0!\n", ic, j);
          cb->liveflag[ic][j] = FALSE;
          continue;
        }
        if ( cb->pmtgain[ic][j] < 0.25 ) {
          perr("excessvely low gain in %d %d (%lf) -- will produce too much NSBG.\n", ic, j, cb->pmtgain[ic][j]);
          cb->liveflag[ic][j] = FALSE;
          continue;
        }
        generateWaveform(pt->n[i][j], pt->t[i][j], tref, pt->len[i][j], 
                  cb->pedestal[ic][j], nsbg[i][j], cb->pmtgain[ic][j], 
            NFADC_BINS, fadcbuf[i][j]+HALF_FRAME_SIZE, &st[i][j]);
      }
    }
  }

  /*
   *  Loop over frames and scan for triggers
   */
  tref = t0 + (double)(NFADC_BINS-HALF_FRAME_SIZE)*100.0;  
  for ( f=1; f<nf; f++ ) {
    for ( i=0; i<pt->nmir; i++ ) {
      ic = pt->mir[i];
      trigmir[ic] = 0;
      tfmask = 0x0;
      for ( j=0; j<NTUBES_CAMERA; j++ ) {
        if ( cb->liveflag[ic][j] ) {
          /* shift current fadc trace by one half frame */
          for ( k=0; k<NFADC_BINS; k++ )
            fadcbuf[i][j][k] = fadcbuf[i][j][k+HALF_FRAME_SIZE];

          /* simulate the next half-frame fadc */
          generateWaveform(pt->n[i][j], pt->t[i][j], tref, pt->len[i][j], 
                    cb->pedestal[ic][j], nsbg[i][j], cb->pmtgain[ic][j], 
              HALF_FRAME_SIZE, fadcbuf[i][j]+NFADC_BINS, &st[i][j]);

          /* compute the trigger mask */
#if ( TFTRIG_TYPE == 1 )
          sdfmask[i][j] = _getTriggerMask(cb->mean[ic][j], cb->vari[ic][j], 
              fadcbuf[i][j]);
#elif ( TFTRIG_TYPE == 2 )
          sdfmask[i][j] = _getTriggerMask(cb->mean[ic][j], cb->vari[ic][j], 
              fadcbuf[i][j]+HALF_FRAME_SIZE);
#endif  /* #if ( TFTRIG_TYPE == ? ) */
          
          tfmask |= sdfmask[i][j];
        }
        else {
          sdfmask[i][j] = 0x0;
        }
      }

      if ( (tfmask&TFTRIG_MASK) != 0x0 ) {
        /* Check for level 1 trigger. */
        int tlev = _isLevel1Trigger(sdfmask[i]);
        if ( tlev == 1 ) {
          trigmir[ic] = TFCODE1;
        }
        else if ( tlev == 0 ) {
          /* Level 1 trigger failed... check edges for level 2 trigger. */
          int tflag = _isLevel2Trigger(sdfmask[i], LEV2_TOP);
          tflag += _isLevel2Trigger(sdfmask[i], LEV2_BOT);
          tflag += _isLevel2Trigger(sdfmask[i], LEV2_LT);
          tflag += _isLevel2Trigger(sdfmask[i], LEV2_RT);

          if ( tflag > 0 ) {
            trigmir[ic] = TFCODE2;
          }
        }
        else {
          vperr("error occurred returning from _isLevel1Trigger().\n");
        }
      }
    }

    /* 
     *  Check if station triggered 
     */
    boolean triggered = FALSE;
    for ( i=0; (i<pt->nmir)&&(triggered==0); i++ ) {
      ic = pt->mir[i];
      if ( trigmir[ic] == TFCODE1 ) {
        triggered = TRUE;
      }
      else if ( trigmir[ic] == TFCODE2 ) {
        /* check adjacent mirrors for level 2 trigger */
        for ( j=0; j<nAdjCameras[ic]; j++ ) {
          if ( trigmir[adjCameraList[ic][j]] == TFCODE2 )
            triggered = TRUE;
        }
      }
    }

    if ( triggered ) {
      pout("triggered frame %d/%d\n", f, nf);

      tftime = (int)((tref-HALF_FRAME_SIZE*100.)/25. + 1.);

#if ( TFTRIG_VER == 2 )
      /*
       * check subsequent frames for potential triggers.
       */
#if ( TFTRIG_TYPE == 1 )
      int exframe = HALF_FRAME_SIZE;
#elif ( TFTRIG_TYPE == 2 )
      int exframe = 2*HALF_FRAME_SIZE;
#endif  /* ( TFTRIG_TYPE == ? ) */
      while ( exframe < (NFADC_BINS-HALF_FRAME_SIZE) ) {
        for ( i=0; i<pt->nmir; i++ ) {
          ic = pt->mir[i];
          for ( j=0; j<NTUBES_CAMERA; j++ ) {
            if ( cb->liveflag[ic][j] )
              sdfmask[i][j] = _getTriggerMask(cb->mean[ic][j], cb->vari[ic][j], 
            wf[i][j]+exframe-HALF_FRAME_SIZE);
            else
              sdfmask[i][j] = 0x0;
          }

          /* Check for level 1 trigger. */
          int tlev = _isLevel1Trigger(sdfmask[i]);
          if ( tlev == 1 ) {
            trigmir[ic] = TFCODE1 + exframe/HALF_FRAME_SIZE;
          }
          else if ( tlev == 0 ) {
            /* Level 1 trigger failed... check edges for level 2 trigger. */
            int tflag = _isLevel2Trigger(sdfmask[i], LEV2_TOP);
            tflag += _isLevel2Trigger(sdfmask[i], LEV2_BOT);
            tflag += _isLevel2Trigger(sdfmask[i], LEV2_LT);
            tflag += _isLevel2Trigger(sdfmask[i], LEV2_RT);

            if ( tflag > 0 ) {
              trigmir[ic] = TFCODE2 + exframe/HALF_FRAME_SIZE;
            }
          }
          else {
            vperr("error occurred returning from _isLevel1Trigger().\n");
          }
        }
        exframe += HALF_FRAME_SIZE;
      }
#endif  /* ( TFTRIG_VER == 2 ) */

#ifdef HISTOGRAM_MODE
      for ( i=0; i<pt->nmir; i++ ) {
        ic = pt->mir[i];
        for ( j=0; j<NTUBES_CAMERA; j++ ) {
          if ( cb->liveflag[ic][j] ) {
            hf1(NSBG_H1D, (float)nsbg[i][j], 1.0);

            /*
            *  let's define a "noise tube" as one that received no signal 
            *    PE in the current time window.
            */
            wf = &fadcbuf[i][j][HALF_FRAME_SIZE];
            cl = fadccl(cb->mean[ic][j], cb->vari[ic][j], NFADC_BINS, wf);
            hf1(CL2_H1D, (float)cl, 1.0);
            cl = _fadcSig(cb->mean[ic][j], cb->vari[ic][j], NFADC_BINS, wf);
            hf1(CL_H1D, (float)cl, 1.0);
            lrms = log10( _getPedRMS(wf) );
            hf1(PEDRMS_H1D, (float)lrms, 1.0);
            petot = 0;
            for ( k=0; k<pt->len[i][j]; k++ ) {
              dt = pt->t[i][j][k] - tref + 
                (double)(NFADC_BINS-HALF_FRAME_SIZE)*100.0;
              if ( dt >= 0.00 && dt < ((double)NFADC_BINS*100.0) )
                petot += (int)pt->n[i][j][k];
            }
            if ( petot == 0 ) {
              hf1(INPUT_MEAN_H1D, (float)(cb->mean[ic][j]/16.0), 1.0);
              hf1(INPUT_SDEV_H1D, (float)(sqrt(cb->vari[ic][j])/4.0), 1.0);

              _getSampleMeanSdev(16, wf, &omean, &osdev);
              hf1(OUTPUT_MEAN_H1D, (float)(omean/16.0), 1.0);
              hf1(OUTPUT_SDEV_H1D, (float)osdev, 1.0);

              hf1(DELTA_MEAN_H1D, (float)(omean-cb->mean[ic][j]), 1.0);
              hf1(DELTA_SDEV_H1D, (float)((osdev-sqrt(cb->vari[ic][j]))/4.0), 1.0);

              hf1(RA_H1D, (float)(cb->vari[ic][j]/
            (16.0*cb->pmtgain[ic][j]*cb->pmtgain[ic][j])), 1.0);

              if ( cl >= 3.00 )
                nnoise_tubes += 1.0f;
            }
            else {
              for ( k=0; k<NFADC_BINS; k++ ) {
                float pswf = (float)( wf[k] - cb->mean[ic][j]/16 );
                hf1(TIME_H1D, (float)k, pswf);
              }
            }
          }
        }
      }
#endif /* #ifdef HISTOGRAM_MODE */

      /* 
       *  Fill adjacent mirrors with sky noise.
       */
      int keepmir[NCAMERAS_SITE] = { 0 };
      for ( i=0; i<pt->nmir; i++ ) {
        ic = pt->mir[i];
        if ( trigmir[ic] > 0 ) {
          keepmir[ic] = 1;
          for ( j=0; j<nAdjCameras[ic]; j++ )
            keepmir[adjCameraList[ic][j]] = 1;
        }
      }

      int ncameras = 0, icamera[fdraw_nmir_max];
      for ( i=0; i<NCAMERAS_SITE; i++ ) {
        int newcam = TRUE;
        for ( j=0; j<pt->nmir; j++ ) {
          if ( pt->mir[j] == i )
            newcam = FALSE;
        }
        if ( newcam && keepmir[i] )
          icamera[ncameras++] = i;

        if ( (pt->nmir+ncameras) > fdraw_nmir_max ) {
          vperr("Maximum number of mirrors reached!\n");
          abort();
        }
      }

      for ( i=0; i<ncameras; i++ ) {
        ic = icamera[i];
        nc = pt->nmir + i;
        for ( j=0; j<NTUBES_CAMERA; j++ ) {
          if ( cb->liveflag[ic][j] ) {
      /*	    nsbg[nc][j] = NSBG_ALPHA*SQR(cb->mean[ic][j]-16*cb->pedestal[ic][j])/(16.*cb->vari[ic][j]); */
            nsbg[nc][j] = (double)cb->vari[ic][j] * NSBG_FUDGE * RHO_ALPHA /
                      ( 16.0 * cb->pmtgain[ic][j] * cb->pmtgain[ic][j] );

#ifdef HISTOGRAM_MODE
            hf1(NSBG_H1D, (float)nsbg[nc][j], 1.0);
#endif  /* #ifdef HISTOGRAM_MODE */

            st[nc][j].it = 3 * HALF_FRAME_SIZE;
            for ( k=0; k<5; k++ )
              st[nc][j].va[k] = 0.00;
            for ( k=0; k<2; k++ )
              st[nc][j].vt[k] = 0.00;

            /* initialize the sky noise and simulate a full waveform */
            if ( cb->pedestal[ic][j] < 0 ) {
              perr("pedestal for %d %d < 0!\n", ic, j);
              cb->liveflag[ic][j] = FALSE;
              continue;
            }
            if ( cb->pmtgain[ic][j] < 0.25 ) {
              perr("excessvely low gain in %d %d (%lf) -- will produce too much NSBG.\n", ic, j, cb->pmtgain[ic][j]);
              cb->liveflag[ic][j] = FALSE;
              continue;
            }

            generateWaveform(NULL, NULL, 0.0, 0, cb->pedestal[ic][j], 
                      nsbg[nc][j], cb->pmtgain[ic][j], HALF_FRAME_SIZE+NFADC_BINS, 
                      fadcbuf[nc][j], &st[nc][j]);

            sdfmask[nc][j] = _getTriggerMask(cb->mean[ic][j], cb->vari[ic][j], 
              fadcbuf[nc][j]);

#ifdef HISTOGRAM_MODE
            wf = &fadcbuf[nc][j][HALF_FRAME_SIZE];
            cl = fadccl(cb->mean[ic][j], cb->vari[ic][j], NFADC_BINS, wf);
            hf1(CL2_H1D, (float)cl, 1.0);

            cl = _fadcSig(cb->mean[ic][j], cb->vari[ic][j], NFADC_BINS, wf);
            hf1(CL_H1D, (float)cl, 1.0);

            lrms = log10( _getPedRMS(wf) );
            hf1(PEDRMS_H1D, (float)lrms, 1.0);
            hf1(INPUT_MEAN_H1D, (float)cb->mean[ic][j], 1.0);
            hf1(INPUT_SDEV_H1D, (float)sqrt(cb->vari[ic][j])/4.0f, 1.0);

            _getSampleMeanSdev(16, wf, &omean, &osdev);
            hf1(OUTPUT_MEAN_H1D, (float)(omean/16.0), 1.0);
            hf1(OUTPUT_SDEV_H1D, (float)(osdev/4.0), 1.0);

            hf1(DELTA_MEAN_H1D, (float)(omean-cb->mean[ic][j])/16.0f, 1.0);
            hf1(DELTA_SDEV_H1D, (float)(osdev-sqrt(cb->vari[ic][j]))/4.f, 1.0);

            hf1(RA_H1D, (float)(cb->vari[ic][j]/
              (16.0*cb->pmtgain[ic][j]*cb->pmtgain[ic][j])), 1.0);

            if ( cl >= 3.00 )
              nnoise_tubes += 1.0f;
#endif  /* #ifdef HISTOGRAM_MODE */
          }
          else {
            sdfmask[nc][j] = 0x0;
          }
        }
      }

      /* 
       *  Fill the fdraw data 
       */
      fdraw->event_code = 0;
      fdraw->ctd_version = 0;
      fdraw->tf_version = 0;
      fdraw->sdf_version = 0;

      /* Adjust ctdclock for trigger time */
      fdraw->ctdclock += tftime;
      while ( (fdraw->ctdclock - fdraw->gps1pps_tick) >= 40000000 ) {
        fdraw->ctdclock -= 40000000;
        fdraw->jsecond++;
      }
      while ( (fdraw->ctdclock - fdraw->gps1pps_tick) < 0 ) {
        fdraw->ctdclock += 40000000;
        fdraw->jsecond--;
      }

      while ( fdraw->jsecond >= 86400 ) {
        fdraw->jsecond -= 86400;
        fdraw->julian++;
      }
      while ( fdraw->jsecond < 0 ) {
        fdraw->jsecond += 86400;
        fdraw->julian--;
      }

      for ( i=0; i<fdraw_nmir_max; i++ ) {
        fdraw->trig_code[i] = trigmir[i];
        fdraw->second[i] = 0;
        fdraw->microsec[i] = 0;
        fdraw->clkcnt[i] = 0;
        fdraw->mir_num[i] = -1;
        fdraw->num_chan[i] = 0;
        fdraw->tf_mode[i] = 0;
        fdraw->tf_mode2[i] = 0;
      }

      int nmir = 0;
      for ( i=0; i<(pt->nmir+ncameras); i++ ) {
        if ( i < pt->nmir )
          ic = pt->mir[i];
        else
          ic = icamera[i-pt->nmir];

        if ( ! keepmir[ic] )
          continue;

        fdraw->trig_code[nmir] = trigmir[ic];
        fdraw->second[nmir] = 0;
        fdraw->microsec[nmir] = 0;
        fdraw->clkcnt[nmir] = fdraw->ctdclock;
        fdraw->mir_num[nmir] = (short)ic;
        fdraw->num_chan[nmir] = NTUBES_CAMERA;  /* to be compacted later */
        fdraw->tf_mode[nmir] = 0;
        fdraw->tf_mode2[nmir] = 0;

        for ( j=0; j<NTUBES_CAMERA; j++ ) {
          if ( cb->liveflag[ic][j] ) {
            fdraw->hit_pt[nmir][j] = (short)( ( sdfmask[i][j] != 0 ) ? 1 : 0 );
            fdraw->channel[nmir][j] = (short)j;
            fdraw->sdf_peak[nmir][j] = 0;
            fdraw->sdf_tmphit[nmir][j] = 0;
            fdraw->sdf_mode[nmir][j] = 0;
            fdraw->sdf_ctrl[nmir][j] = 0;
            fdraw->sdf_thre[nmir][j] = (short)6;  /* this should be a macro */
            for ( k=0; k<4; k++ )
              fdraw->mean[nmir][j][k] = (unsigned short)cb->mean[ic][j];
            for ( k=0; k<4; k++ )
              fdraw->disp[nmir][j][k] = (unsigned short)cb->vari[ic][j];
            for ( k=0; k<NFADC_BINS; k++ )
              fdraw->m_fadc[nmir][j][k] = 
                (short)fadcbuf[i][j][k+HALF_FRAME_SIZE];
          }
          else {
            fdraw->hit_pt[nmir][j] = 0;
            fdraw->channel[nmir][j] = -1;
            fdraw->sdf_peak[nmir][j] = 0;
            fdraw->sdf_tmphit[nmir][j] = 0;
            fdraw->sdf_mode[nmir][j] = 0;
            fdraw->sdf_ctrl[nmir][j] = 0;
            fdraw->sdf_thre[nmir][j] = 0x7FFF;
            for ( k=0; k<4; k++ )
              fdraw->mean[nmir][j][k] = 0;
            for ( k=0; k<4; k++ )
              fdraw->disp[nmir][j][k] = 0x7FFF;
            for ( k=0; k<NFADC_BINS; k++ )
              fdraw->m_fadc[nmir][j][k] = 0;
          }
        }
        fdraw->hit_pt[nmir][j] = 0;
        nmir++;
      }
      fdraw->num_mir = nmir;

#ifdef HISTOGRAM_MODE
      nnoise_tubes /= (float)nmir;
      hf1(NNOISE_TUBES_H1D, nnoise_tubes, 1.0);
#endif  /* #ifdef HISTOGRAM_MODE */

      pout("trigger SUCCEEDED. (t1=%.1lf t2=%.1lf)\n", 
	   pt->t1-tref+(double)(NFADC_BINS-HALF_FRAME_SIZE)*100.0, 
	   pt->t2-tref+(double)(NFADC_BINS-HALF_FRAME_SIZE)*100.0);

      return 1;
    }

    tref += (double)HALF_FRAME_SIZE*100.0;
  }

  pout("trigger FAILED.\n");
  return 0;
}
