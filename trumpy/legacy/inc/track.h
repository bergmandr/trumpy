#ifndef _TRACK_H_
#define _TRACK_H_

#define DEPTH_STEPSIZE   1.0
#define DEPTH_RESOLUTION 1.0e-3

#include "event.h"

typedef struct {
  /* Segment of a Shower Track */
  double position;   /* distance from top of first segment */
  double sdepth;     /* slant depth (at start of segment) */
  double age;        /* age at middle of segment */
  double dlseg;      /* distance (m) to start of next shower segment */
  double dlmid;      /* density weighted distance (m) to middle of
			 segment  */
  double dxseg;      /* total grammage in shower segment */
  double height;     /* altitude of top of shower segment (m) */
  double dedep;      /* energy deposition (summed over shower) rate
			(eV/(g/cm2)) at top of shower segment */
  double nch;        /* no. of charged part. in this segment from GH */
  double molrad;     /*  Moliere radius (meters) */
  double nfl[NWAVELEN_BANDS]; /* fluorescence photons generated in 5
				 nm wavelength bands within the
				 (following) segment */
  double pcv[NWAVELEN_BANDS]; /* Cherenkov photons added to beam in
			         the previous segment */
  double ncv[NWAVELEN_BANDS]; /* Cherenkov photons in beam at top of
				 segment */
} TrackSegment;

typedef struct {
  /* Shower Track */
  int nseg;               /* number of air shower segments */
  double xoffset;         /* offset in slant depth (g/cm^2) */
  double length;          /* total length of air shower, i.e., the
			     distance from the impact point to the top of
			     the first segment */
  double t0;              /* Time origin for event */
  TrackSegment *segment;   /* Array of individual segments */
  fdatmos_param_dst_common aparam;         /* current molecular atmosphere */
  fdscat_dst_common atrans;        /* current aerosol parameters */
} Track;

/* in "track.c" */
#ifdef __cplusplus
extern "C"{
#endif
int buildTrack(RuntimeParameters *par, AirShower *as, 
	       GaisserHillasParameters *gh, Track *track);
int buildTrackGeom(RuntimeParameters *par, AirShower *as, Track *track);
void fillTrackLight(RuntimeParameters *par, AirShower *as, 
		    GaisserHillasParameters *gh, Track *track);
void clearTrack(Track *track);
double findTrackLength(double zenith, double altitude);
int getSegmentsInView(AirShower *as, Track *track, double *cvec,
		      double *zvec, double tolerance, int *segnum);
int getTelescopesInView(double r[3], int nc, double cvec[][3],
			double zvec[][3], double tolerance, int *camera);

/* in "track_up.c" */
int buildLaserTrack(double lambda, RuntimeParameters *par, AirShower *as, 
		     GaisserHillasParameters *gh, Track *track);
int buildUpwardTrackGeom(RuntimeParameters *par, AirShower *as, Track *track);
void fillLaserTrackLight(double lambda, RuntimeParameters *par, AirShower *as, 
			 GaisserHillasParameters *gh, Track *track);
#ifdef __cplusplus
}
#endif
#endif
