#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "constants.h"
#include "control.h"
#include "airshower.h"
#include "random.h"
#include "toolbox.h"

#define IMPACT_POINT_GENERATION

static int nb;
static double bp[MAX_NBREAK+1];
static double ea[MAX_NBREAK+1], eb[MAX_NBREAK+1], ec[MAX_NBREAK+1];

int initAirShower(RuntimeParameters *par) {
  int i, j;
  double el[MAX_NBREAK+1], eg[MAX_NBREAK];
  double a, n, g, dE;

  if ( fabs(par->logehi-par->logelo) < DOUBLE_FP_PRECISION ) { 
    /* 
     *  monoenergetic mode 
     */
    pout("monoenergetic mode (logE=%lf)\n", par->logelo);

    ea[0] = par->logelo;
    eb[0] = 0.0;
    ec[0] = 0.0;
    bp[0] = 1.00;

    return 0;
  }

  /*
   *  find where the breakpoints are wrt the energy bounds
   */
  el[0] = par->logelo;
  eg[0] = 0.00; /* not used */
  nb = 1;
  if ( par->nbreak > 0 ) {
    for ( i=0; i<par->nbreak; i++ ) {
      if ( par->logelo < par->ebreak[i] ) {
	for ( j=i; j<par->nbreak; j++ ) {
	  if ( par->logehi < par->ebreak[j] )
	    break;

	  eg[nb] = par->eslope[j];
	  el[nb] = par->ebreak[j];
	  nb++;
	}

	eg[nb] = par->eslope[j];
	el[nb] = par->logehi;
	nb++;
	break;
      }
    }
  }
  else {
    eg[1] = par->eslope[0];
    el[1] = par->logehi;
    nb = 2;
  }

  /*
   *  initialize energy distribution constants
   */
  n = 0.00;
  bp[0] = n;
  a = 1.00;
  for ( i=1; i<nb; i++ ) {
    g = 1.0 - eg[i];
    if ( fabs(g) < DOUBLE_FP_PRECISION ) {
      pout("distribution is flat in the range %lf<logE<%lf\n", el[i-1], el[i]);
      ea[i] = el[i-1];
      eb[i] = 1.0 / a;
      ec[i] = 0.00;

      n += a * ( el[i] - el[i-1] );
    }
    else {
      dE = pow(10.0, g*el[i]);
      dE -= pow(10.0, g*el[i-1]);

      ea[i] = pow(10.0, g*el[i-1]) - g*n/a;
      eb[i] = g / a;
      ec[i] = 1.0 / g;

      n += a * dE / g;
    }

    bp[i] = n;

    if ( i == nb-1 )
      break;

    a *= pow(10.0, el[i]*(eg[i+1]-eg[i]));
  }

  for ( i=0; i<nb; i++ ) {
    bp[i] /= n;
    eb[i] *= n;
  }

  /*
  printf("nb = %d\n", nb);
  printf("eg = %e  %e  %e  %e\n", eg[0], eg[1], eg[2], eg[3]);
  printf("el = %e  %e  %e  %e\n", el[0], el[1], el[2], el[3]);
  printf("bp = %e  %e  %e  %e\n", bp[0], bp[1], bp[2], bp[3]);
  printf("ea = %e  %e  %e  %e\n", ea[0], ea[1], ea[2], ea[3]);
  printf("eb = %e  %e  %e  %e\n", eb[0], eb[1], eb[2], eb[3]);
  printf("ec = %e  %e  %e  %e\n", ec[0], ec[1], ec[2], ec[3]);
  */

  return 0;
}

int createAirShower(const RuntimeParameters *par, AirShower *as) {
  int i;
  double u;
  double zup, qup, pup, upv[3];
  double zdn, qdn, pdn, dnv[3];
  double rm[3][3];

  // pout("createAirShower()", "creating extensive air shower.");

#ifdef IMPACT_POINT_GENERATION
  /* find the impact point on the Earth's surface in terms of CLF coords. */
  pup = 2.0*M_PI * RANDOM_NUMBER;
  //  qup = 1.0 * M_PI / 180.0;
  //  qup = 0.5 * M_PI / 180.0;
  qup = par->lat;
  u = RANDOM_NUMBER;
  zup = u * ( 1.0 - cos(qup) ) + cos(qup);
#  ifdef ORIGINAL_IPG
  upv[0] = cos(pup) * sqrt(1.0 - zup*zup);
  upv[1] = sin(pup) * sqrt(1.0 - zup*zup);
  upv[2] = zup;
#  else
  upv[0] = cos(pup) * tan(qup) * sqrt(u);
  upv[1] = sin(pup) * tan(qup) * sqrt(u);
  upv[2] = 1.0;
#  endif

  /*
   * find the direction vector of the shower track, then rotate to it's 
   *   impact point.
   */
  pdn = 2.0*M_PI * RANDOM_NUMBER;
  qdn = par->thetahi;
  /*
   * This assumes that area is independant of angle,
   *   zdn = RANDOM_NUMBER * ( 1.0 - cos(qdn) ) - 1.00;
   * where 
   *   p(q)dq = sin(q)dq = -d(cos(q))
   * However, the area is effectively smaller by cos(q)
   *   p(q)dq = sin(q)cos(q)dq = -cos(q)d(cos(q))
   */
  double cosqdn = cos(qdn);
  double cos2qdn = cosqdn*cosqdn;
#  ifndef ZENITH_KLUDGE
  double cos2q = cos2qdn + RANDOM_NUMBER * (1.-cos2qdn);
  zdn = -sqrt(cos2q);
  dnv[0] = cos(pdn) * sqrt(1.0 - zdn*zdn);
  dnv[1] = sin(pdn) * sqrt(1.0 - zdn*zdn);
  dnv[2] = zdn;
#  else
  double select;
  do {
    double cos2q = cos2qdn + RANDOM_NUMBER * (1.-cos2qdn);
    zdn = -sqrt(cos2q);
    dnv[0] = cos(pdn) * sqrt(1.0 - zdn*zdn);
    dnv[1] = sin(pdn) * sqrt(1.0 - zdn*zdn);
    dnv[2] = zdn;
    double zenith = 180. * acos(-zdn) / M_PI;
    // select = 1. - 0.01*zenith;
    select = 0.5*cos2q + 0.5;
  } while ( RANDOM_NUMBER > select );
#  endif
#else
  for (;;) {// Repeat until successful
    /* Find direction of shower; flat on 2PI, but out to thetahi */
    pdn = 2.0*M_PI * RANDOM_NUMBER;
    qdn = par->thetahi;
    zdn = RANDOM_NUMBER * ( 1.0 - cos(qdn) ) - 1.00;
    dnv[0] = cos(pdn) * sqrt(1.0 - zdn*zdn);
    dnv[1] = sin(pdn) * sqrt(1.0 - zdn*zdn);
    dnv[2] = zdn;
    
    /* Find position on disk perpendicular to dnv */
    prp = 2.0*M_PI * RANDOM_NUMBER;
    rpm = RANDOM_NUMBER * (par->rphi*par->rphi);
    rpm = sqrt(rpm);
    
    /* Find Rp point wrt CLF */
    rpv[0] = rpm * cos(prp);
    rpv[1] = rpm * sin(prp);
    rpv[2] = 0;
    zmatrix(dnv, rm);
    applyRotation(rm, rpv, rpv);
    
    /* Find convert Rp point to unit earth sphere with CLF at (0,0,1) */
    rpv[0] /= EARTH_RADIUS;
    rpv[1] /= EARTH_RADIUS;
    rpv[2] /= EARTH_RADIUS;
    rpv[2] += 1.;

    /*
     * Find intersection point between track through Rp point
     * and unit sphere
     */
    double dnvdotrpv = dotprod(dnv, rpv);
    double rpv2 = dotprod(rpv, rpv);
    double discriminant = dnvdotrpv*dnvdotrpv - (rpv2-1);
    if (discriminant <= 0.)
      continue; // Try again, no solution or tangent
    else {
      int i;
      double s1[3], s2[3], d1, d2;
      double sqrtDisc = sqrt(discriminant);
      d1 = -dnvdotrpv-sqrtDisc;
      d2 = -dnvdotrpv+sqrtDisc;
      for (i=0; i<3; i++) {
	s1[i] = rpv[i] + d1*dnv[i];
	s2[i] = rpv[i] + d2*dnv[i];
      }
      assert( dotprod(s1,dnv)*dotprod(s2,dnv) <= 0. );
      if (dotprod(s1,dnv) < 0.)
	for (i=0; i<3; i++)
	  upv[i] = s1[i];
      else 
	for (i=0; i<3; i++)
	  upv[i] = s2[i];
      break; // We found the solution, go on to next step
    }
  }
#endif

  zmatrix(upv, rm);
  /*
   * zmatrix gets rotation from (0, 0, 1) to (z0, z1, z2)
   * Need to invert this matrix to rotate TO the CLF coordinate system 
   */
  matrixInverse (rm, rm);
  applyRotation(rm, dnv, as->trackuv);

  as->zenith = acos(-zdn);

  /* give impact coordinate units of km */
  as->impactv[0] = upv[0] * EARTH_RADIUS;
  as->impactv[1] = upv[1] * EARTH_RADIUS;
  as->impactv[2] = ( upv[2] - 1.00 ) * EARTH_RADIUS;

  /* get shower energy */
  u = RANDOM_NUMBER;

  i = 0;
  while ( (u>bp[i]) && (i<nb) )
    i++;

  if ( fabs(ec[i]) < DOUBLE_FP_PRECISION )
    as->loge = ea[i] + eb[i]*u;
  else
    as->loge = ec[i] * log10( ea[i] + eb[i]*u );

  if ( par->nshowlib == 1 )
    as->species = par->species[0];
  else
    as->species = randomParticle(as->loge);

  return 0;
}

#define PROTC1    246.6893  /* The proton and Iron distribution functions were found by fitting */
#define PROTC2     52.4137  /*   the showers of the Nerling shower library.  The showers come   */
#define IRONC1     73.8595  /*   in 5 bins of sec(q), so a weighted average was used, assuming  */
#define IRONC2     61.3975  /*   zenith angles will follow a sin(2q) distribution.              */
#define DATA1C1  -127.97
#define DATA1C2    93.33
#define DATA2C1   209.22
#define DATA2C2    55.23
#define PSPECBP    8.85

int randomParticle(double loge) {
  double fprot, firon, fdist, x, piron;
  double legev = loge - 9.0;  // convert to GeV

  fprot = PROTC2 * legev + PROTC1;
  firon = IRONC2 * legev + IRONC1;
  if ( legev < PSPECBP )
    fdist = DATA1C2 * legev + DATA1C1;
  else
    fdist = DATA2C2 * legev + DATA2C1;

  x = ( fprot - fdist ) / ( fdist - firon );
  piron = x / ( 1.0 + x );

  if ( RANDOM_NUMBER <= piron )
    return PSPEC_IRON;
  else
    return PSPEC_PROTON;

}

#undef PROTC1
#undef PROTC2 
#undef IRONC1
#undef IRONC2
#undef DATA1C1
#undef DATA1C2
#undef DATA2C1
#undef DATA2C2
#undef PSPECBP

