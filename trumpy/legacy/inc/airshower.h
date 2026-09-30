#ifndef _AIRSHOWER_H_
#define _AIRSHOWER_H_

typedef struct {
  int species;        // primary particle species
  double loge;        // shower energy in log(E/eV)
  double zenith;      // shower zenith angle at impact coordinate
  double impactv[3];  // impact point in CLF coordinates
  double trackuv[3];  // direction vector along shower's path
} AirShower;

/* in "airshower.c" */
#ifdef __cplusplus
extern "C"{
#endif
int initAirShower(RuntimeParameters *par);
int createAirShower(const RuntimeParameters *par, AirShower *as);
int randomParticle(double loge);
#ifdef __cplusplus
}
#endif

#endif
