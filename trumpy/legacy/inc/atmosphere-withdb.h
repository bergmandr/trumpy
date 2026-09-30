/*
 *  atmosphere-withdb.h
 *
 *  Author:  Sean R. Stratton
 *           Rutgers University, dept. of Physics & Astronomy
 *
 *  Set of functions for accessing the atmospheric database (DST) and 
 *    computing thermodynamic and scattering parameters by run night.
 *
 *  Revisions:
 *     4.24.09  Base revision
 *     4.28.09  Added atmospheric transmission functions
 *
 */
#ifndef _ATMOSPHERE_WITHDB_H_
#define _ATMOSPHERE_WITHDB_H_

/*
 *  because 'fdatmos_xxxxx_dst_common' is a pain in the ass to keep writing.
 */
#define ATMP_t fdatmos_param_dst_common
#define ATMT_t fdscat_dst_common

/*
 *  Public functions
 */
#ifdef __cplusplus
extern "C"{
#endif
void getDefaultParam(ATMP_t *);
void getDefaultTrans(ATMT_t *);
int loadAtmosDB(char *);
int loadScatterDB(char *);
int closeAtmosDB(void);
int closeScatterDB(void);
int getAtmosphere(const RuntimeParameters *, ATMP_t *, ATMT_t *);
void getPTDByAltitude(ATMP_t *, double, double *, double *, double *);
double getPressureByAltitude(ATMP_t *, double);
double getTemperatureByAltitude(ATMP_t *, double);
double getDensityByAltitude(ATMP_t *, double);
#ifdef __cplusplus
}
#endif

//static double getCurrentDensityByAltitude(double alt);

/*
 *  Private functions
 */
#ifdef MAKE_PRIVATE_PUBLIC
int _getLayerIndex(double, int, float *);
void _evaluatePrTp(ATMP_t *, int, double, double *, double *);
int _compareBoundaryValues(ATMP_t *, int, ATMP_t *);
void _getOutOfBoundState(ATMP_t *, int, ATMP_t *);
int _nextAtmosDBRecord(void);
#endif

#if 0
/*
 *  The following arrays come from the 1966 supplement to the 1976 US standard
 *    atmosphere.  They are currently unused, but should be kept for reference.
 */
static double GWINTER[] = { 0.00, 3.00, 10.0, 19.0, 27.0, 32.0 };
static double TWINTER[] = { 272.15, 261.65, 219.65, 215.15, 215.15, 219.15 };

static double GSUMMER[] = { 0.00, 2.00, 6.00, 13.0, 17.0, 27.0, 32.0 };
static double TSUMMER[] = { 294.15, 285.15, 261.15, 215.65, 215.65, 227.65, 238.15 };

static double GUSMEAN[] = { 0.00, 11.0, 20.0, 32.0 };
static double TUSMEAN[] = { 288.15, 216.65, 216.65, 228.65 };
#endif

/*
 *  Functions for reconciling with soon-to-be-retired 'atmosphere.h'
 */
#ifdef __cplusplus
extern "C"{
#endif
double getGrammageRateOfChange(ATMP_t *, double);
double getGrammageByAltitude(ATMP_t *, double);
double setupAtoB(double ptA[3], double ptB[3], ATMP_t *season);
double getDensityAlongVector(double dl);
double getLRhoAlongVector(double dl);
double getLRhoIntegralBetweenPoints(ATMP_t *season, 
				    double ptA[3], double ptB[3]);
double getGrammageBetweenPoints(ATMP_t *season, double ptA[3], double ptB[3]);
double getSlantDepth(ATMP_t *season, double rt[3], double rb[3]);
double getDistanceByGrammage(ATMP_t *season, 
			     double ptA[3], 
			     double direction[3], 
			     double grammage,
			     double eps);
double getLargeDistanceByGrammage(ATMP_t *season, 
				  double ptA[3], 
				  double direction[3], 
				  double grammage,
				  double eps);
double getDistanceToAtmosphere(const double ptA[3], 
			       const double direction[3], 
			       double eps);
double getDeltaByPressTemp(double P, double T);
#ifdef __cplusplus
}
#endif

#endif
