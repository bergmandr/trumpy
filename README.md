# TRUMP in python

This is an implementation of the TRUMP air-fluorescence detector simulation MC program. The fundamental C `structs` have been implemented as python `dataclass`es. 

An example: `Track`
```python
@dataclass
class Track:
    """Intrinsic physical properties of the shower or laser."""
    # Global parameters
    species: int            # e.g., 1 for proton, 5626 for Fe
    log_e: float            # log10(E/eV)
    zenith: float           # zenith angle in radians
    impact_v: np.ndarray    # shape (3,) impact point in CLF coordinates
    track_uv: np.ndarray    # shape (3,) unit vector of track direction

    # Vectorized Segments (Structure of Arrays)
    n_segments: int
    positions: np.ndarray   # shape (n_segments, 3)
    time_gen: np.ndarray    # shape (n_segments,) emission times
    altitude: np.ndarray    # shape (n_segments,) altitude at segment mid-points
    dl_seg: np.ndarray      # shape (n_segments,) segment lengths
    de_dep: np.ndarray      # shape (n_segments,) energy deposit profile
```
replaces `TrackSegment` and `Track` from the C code
```C
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
  fdatmos_param_dst_common aparam; /* current molecular atmosphere */
  fdscat_dst_common atrans;        /* current aerosol parameters */
} Track;
```

Gemini's expectation for pluggable components:
```
Simulation
 ├── Source (Track Generator)
 │    ├── AirShowerGenerator (Gaisser-Hillas, shower libraries)
 │    ├── LaserGenerator (Vertical/steered laser tracks)
 │    ├── BeamPointGenerator (FLUKA/ELS voxel-based deposits)
 │    └── EventListReplay (Pre-generated CORSIKA / external lists)
 │
 ├── Atmosphere
 │    ├── Profile (Pressure, temperature, density via spline interpolation)
 │    ├── FluorescenceYield (FLASH, Kakimoto, AIRFLY models)
 │    └── ExtinctionEngine (Rayleigh scattering + Mie aerosol models)
 │
 ├── Experiment (Detector Hierarchy)
 │    └── Site (e.g., Black Rock, Long Ridge, Middle Drum, TALE)
 │         └── Telescope / Mirror (Center of curvature, focal length, segmented dish)
 │              ├── Optics (Paraglas filter, BG3, mirror reflectivity curves)
 │              ├── Camera (Hexagonal / rectangular PMT layout)
 │              │    └── PMTPointingTable (Vectorized directions for ~6000 tubes)
 │              └── Electronics & Trigger (FADC pipeline, cluster/track trigger)
 │
 └── RunSchedule
      ├── FixedTrialsSchedule (Trial-bounded or event-bounded)
      └── LiveTimeSchedule (Data-driven on-time, pedestal/noise tracking)
```