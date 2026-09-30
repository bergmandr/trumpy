#ifndef _ATMOSPHERE_H_
#define _ATMOSPHERE_H_

#include "constants.h"

#define ATMOS_EPS 1.0e-6
#define UNIV_GAS_CONST 8.31432 // N.m/(mol.K)
#define AIR_MOLAR_MASS 0.0289644 // kg/mol
#define KPA_TO_TORR 7.5006

/* This taken from http://emtoolbox.nist.gov, in the documentation
   section for a "Shop-floor Folmula for the Index of Refraction of
   Air, and using 0% relative humidity */
#define DELTA_COEFFICIENT 7.86e-4

enum {
  SPRING, 
  SUMMER, 
  FALL, 
  WINTER
};

static const double SUMMER_ATMOSPHERIC_DATA[31*4] = {
     0.0,  1.19E-03, 1013.5,  296.22,
     1.0,  1.08E-03,  902.2,  291.14,
     2.0,  9.76E-04,  801.5,  286.19,
     3.0,  8.85E-04,  710.3,  279.78,
     4.0,  8.00E-04,  627.8,  273.55,
     5.0,  7.21E-04,  553.3,  267.43,
     6.0,  6.48E-04,  486.3,  261.30,
     7.0,  5.83E-04,  425.9,  254.75,
     8.0,  5.22E-04,  371.8,  248.21,
     9.0,  4.66E-04,  323.4,  241.69,
    10.0,  4.15E-04,  280.3,  235.17,
    11.0,  3.69E-04,  241.9,  228.65,
    12.0,  3.26E-04,  207.8,  222.15,
    13.0,  2.87E-04,  177.8,  215.65,
    14.0,  2.45E-04,  151.8,  215.65,
    15.0,  2.09E-04,  129.5,  215.65,
    16.0,  1.79E-04,  110.5,  215.65,
    17.0,  1.52E-04,  94.35,  215.65,
    18.0,  1.29E-04,  80.56,  216.85,
    19.0,  1.10E-04,  68.85,  218.05,
    20.0,  9.36E-05,  58.89,  219.25,
    21.0,  7.97E-05,  50.42,  220.45,
    22.0,  6.79E-05,  43.20,  221.65,
    23.0,  5.79E-05,  37.04,  222.85,
    24.0,  4.94E-05,  31.79,  224.05,
    25.0,  4.22E-05,  27.31,  225.25,
    26.0,  3.61E-05,  23.47,  226.45,
    27.0,  3.09E-05,  20.19,  227.65,
    28.0,  2.64E-05,  17.39,  229.75,
    29.0,  2.25E-05,  15.00,  231.85,
    30.0,  1.93E-05,  12.95,  233.95
};

static const double SPRFAL_ATMOSPHERIC_DATA[31*4] = {
     0.0,  1.23E-03, 1013.3,  288.15,
     1.0,  1.11E-03,  898.7,  281.65,
     2.0,  1.01E-03,  795.0,  275.15,
     3.0,  9.09E-04,  701.1,  268.65,
     4.0,  8.19E-04,  616.4,  262.15,
     5.0,  7.36E-04,  540.2,  255.65,
     6.0,  6.60E-04,  471.8,  249.15,
     7.0,  5.90E-04,  410.6,  242.65,
     8.0,  5.25E-04,  356.0,  236.15,
     9.0,  4.66E-04,  307.4,  229.65,
    10.0,  4.13E-04,  264.4,  223.15,
    11.0,  3.64E-04,  226.3,  216.65,
    12.0,  3.11E-04,  193.3,  216.65,
    13.0,  2.66E-04,  165.1,  216.65,
    14.0,  2.27E-04,  141.0,  216.65,
    15.0,  1.94E-04,  120.4,  216.65,
    16.0,  1.65E-04,  102.9,  216.65,
    17.0,  1.41E-04,  87.87,  216.65,
    18.0,  1.21E-04,  75.05,  216.65,
    19.0,  1.03E-04,  64.10,  216.65,
    20.0,  8.80E-05,  54.75,  216.65,
    21.0,  7.49E-05,  46.78,  217.65,
    22.0,  6.37E-05,  40.00,  218.65,
    23.0,  5.43E-05,  34.22,  219.65,
    24.0,  4.63E-05,  29.30,  220.65,
    25.0,  3.95E-05,  25.11,  221.65,
    26.0,  3.37E-05,  21.53,  222.65,
    27.0,  2.88E-05,  18.47,  223.65,
    28.0,  2.46E-05,  15.86,  224.65,
    29.0,  2.10E-05,  13.63,  225.65,
    30.0,  1.80E-05,  11.72,  226.65
};

static const double WINTER_ATMOSPHERIC_DATA[31*4] = {
     0.0,  1.30E-03, 1018.0,  272.59,
     1.0,  1.16E-03,  897.3,  269.00,
     2.0,  1.04E-03,  789.6,  265.43,
     3.0,  9.23E-04,  693.7,  261.85,
     4.0,  8.28E-04,  607.9,  255.77,
     5.0,  7.41E-04,  531.0,  249.74,
     6.0,  6.61E-04,  462.4,  243.70,
     7.0,  5.88E-04,  401.2,  237.68,
     8.0,  5.22E-04,  346.8,  231.66,
     9.0,  4.61E-04,  298.7,  225.66,
    10.0,  4.06E-04,  256.2,  219.65,
    11.0,  3.49E-04,  219.3,  219.15,
    12.0,  2.99E-04,  187.6,  218.65,
    13.0,  2.56E-04,  160.4,  218.15,
    14.0,  2.20E-04,  137.1,  217.65,
    15.0,  1.88E-04,  117.2,  217.15,
    16.0,  1.61E-04,  100.1,  216.15,
    17.0,  1.38E-04,  85.49,  216.15,
    18.0,  1.18E-04,  72.98,  215.65,
    19.0,  1.01E-04,  62.28,  215.15,
    20.0,  8.60E-05,  53.13,  215.15,
    21.0,  7.34E-05,  45.33,  215.15,
    22.0,  6.26E-05,  38.68,  215.15,
    23.0,  5.34E-05,  33.00,  215.15,
    24.0,  4.56E-05,  28.15,  215.15,
    25.0,  3.89E-05,  24.02,  215.15,
    26.0,  3.32E-05,  20.49,  215.15,
    27.0,  2.83E-05,  17.48,  215.15,
    28.0,  2.41E-05,  14.92,  215.95,
    29.0,  2.05E-05,  12.74,  216.75,
    30.0,  1.74E-05,  10.89,  217.55
};

/* in "atmosphere.c" */
double getDensityByAltitude(int season, double altitude);
double getGrammageByAltitude(int season, double altitude);
double getTemperatureByAltitude(int season, double altitude);
double getPressureByAltitude(int season, double altitude);
double getAltitudeByGrammage(int season, double x);
double getGrammageRateOfChange(int season, double altitude);
double getLRhoIntegralBetweenPoints(int season, double ptA[3], double ptB[3]);
double getGrammageBetweenPoints(int season, double ptA[3], double ptB[3]);
double getSlantDepth(int season, double rtop[3], double rbot[3]);
double getDistanceByGrammage(int season, double ptA[3], double direction[3], double grammage, double eps);
double getDistanceToAtmosphere(const double ptA[3], 
			       const double direction[3], 
			       double eps);
double getDeltaByPressTemp(double P, double t);
#endif
