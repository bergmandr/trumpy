#ifndef _SHOWERLIB_H_
#define _SHOWERLIB_H_

// #define NENTRIES_PER_DRAWER  500
// #define NSHOWLIB_ENTRIES    22500
#define MAX_SHOWLIB_ENTRIES 100000

typedef struct {
  int pspec;
  showscale_dst_common scaledef;
//   showlib_dst_common entry[NSHOWLIB_ENTRIES];
  showlib_dst_common entry[MAX_SHOWLIB_ENTRIES];
} ShowerLibrary;

/*
 * these functions should remain the same no matter what shower library is 
 *   being used!!
 */

typedef struct {
  double x0;
  double xmax;
  double nmax;
  double lambda;
} GaisserHillasParameters;

/* in "showerlib.c" */
#ifdef __cplusplus
extern "C"{
#endif
int loadShowerLibraries(RuntimeParameters *par);
int loadShowerLibrary(char *path, ShowerLibrary *lib);
//int getGaisserHillasParameters(ShowerLibrary *lib, double loge, double secth, GaisserHillasParameters *gh);
int getGaisserHillasParameters(AirShower *as, GaisserHillasParameters *gh);

/* in "gaisser-hillas.c" */
double gaisserHillasFunction(GaisserHillasParameters *gh, double x);

/* in "gaussian-in-age.c" */
double gaussianInAgeFunction(GaisserHillasParameters *gh, double x);
#ifdef __cplusplus
}
#endif
#endif
