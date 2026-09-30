#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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
#include "calibration.h"
#include "raytrace.h"
#include "acpttrack.h"
#include "taelectronics.h"

#include "histogram.h"

#define SQR(X) ((X)*(X))

#define BIN_DT       25.0         /* nanoseconds (see constants.h)           */
#define RES          25           /* number of elements of time const. array */
#define DRES         25.0         /* double-precision version of resolution  */

#ifndef VOLTAGE_OFFSET
#  define VOLTAGE_OFFSET 0.282  /* in units of FADC */
#endif

#ifndef NOISE_AMPLITUDE
#  define NOISE_AMPLITUDE 0.9
#endif

/* the values of the one PE TF in BIN_DT/RES ns resolution */
static double tft1[RES], tft2[RES];

/* constants used to propagate TF's using recursion relation */
static double tct1, tct2;

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
int initSimpleWaveformGenerator(double tau) {
  int i;
  double dt, t;

  dt = BIN_DT/DRES;
  /* t is time from pe to measurement with BIN_DT/RES ns resolution */
  t = dt/2.00;
  for ( i=0; i<RES; i++ ) {
    tft1[i] = exp(-t/tau) / tau;
    tft2[i] = t * tft1[i] / tau;

    t += dt;
  }

  /* tc's are used to propogate the transfer function over time */
  tct1 = exp(-BIN_DT/tau);
  tct2 = BIN_DT * tct1 / tau;

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
int generateSimpleWaveform(double *npe, double *tpe, double tref, int len, 
		     int ped, 
		     double nsbg, 
		     double gain, 
		     int nfadc, short *fadc, WaveformState *st) {
  //  const int ixx=532521516;  /* VERY important constant */

  int i, j, it, idx, ifadc;
  double u, expnsbg;
  double t, tnext;
  double vt1, vt2, vout;

  /*
   * This prevents the while loops from having to constantly execute the 
   *   'exp' math function.
   */
  expnsbg = exp(-nsbg/4.0);

  vt1 = st->vt[0];
  vt2 = st->vt[1];

  /* simulate the waveform */
  idx = 0;
  tnext = BIN_DT;
  for ( i=0; i<nfadc; i++ ) {
    vout = 0.00;

    for ( j=0; j<4; j++ ) {
      /* update readout voltage terms */
      vt2 = tct1*vt2 + tct2*vt1;  // add some noise to shaper
      vt1 = tct1*vt1;

      /* add noise pe */
      /* this loop simulates Poisson statistics */
      u = RANDOM_NUMBER;
      while ( u >= expnsbg ) {
	it = (int)( RANDOM_NUMBER * DRES );

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
	  vt1 += npe[idx] * tft1[it];
	  vt2 += npe[idx] * tft2[it];
	}
      } 

      vout += vt2;
      tnext += BIN_DT;
    }

    //    ifadc = (int)( gain*BIN_DT*vout + 0.9*sampled_noise[st->it] + VOLTAGE_OFFSET ) + ped;
    ifadc = (int)( gain*BIN_DT*vout + NOISE_AMPLITUDE*sampled_noise[st->it] + VOLTAGE_OFFSET ) + ped;
    fadc[i] = (short)( (ifadc<0x3FFF) ? ifadc : 0x3FFF );

    st->it = ( st->it + 1 ) % NOISE_ARY_LEN;
  }

  st->vt[0] = vt1;
  st->vt[1] = vt2;

  return 0;
}

