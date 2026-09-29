# TRUMP in python

This is an implementation of the TRUMP air-fluorescence detector simulation MC program. The fundamental C `structs` have been implemented as python `dataclass`es. 

## Simulation

The `Simulation` class is the top level organization in `trumpy` replacing much of what is in `trump.c`. It has a number of `Protocol`s for other class types needed to run the program
```python
class TrackSource(Protocol):
    def generate_event(self, rng: np.random.Generator, timestamp: float) -> Track:
class RunSchedule(Protocol):
    def get_trials(self) -> Iterator[Tuple[int, float]]:
class AtmosphereModel(Protocol):
    def get_fluorescence_yield(self, altitudes: np.ndarray, de_dep: np.ndarray, wavelengths: np.ndarray) -> np.ndarray:
    def get_transmission(self, emission_points: np.ndarray, mirror_centers: np.ndarray, wavelengths: np.ndarray) -> np.ndarray:
class Experiment(Protocol):
    wavelength_bands: np.ndarray
    mirror_centers: np.ndarray
    def passes_fast_cuts(self, track: Track) -> bool:
    def trace_photons(self, track: Track, photons_at_mirrors: np.ndarray, rng: np.random.Generator) -> ak.Array:
    def process_electronics(self, pe_times: ak.Array, rng: np.random.Generator) -> ak.Array:
```
basically saying that one needs a source of `Track`s, a run schuduler (to prelace the main loop in `trump.c`), an atmosphere, and the experimental setup in order to run the program. The `Simulation` class sets up everything and does the run.
```python
class Simulation:
    def __init__(self, experiment: Experiment, atmosphere: AtmosphereModel, source: TrackSource, schedule: RunSchedule):
    def run(self) -> ak.Array:
```
Note that the results are returned as an `awkward` array.

### Summary of Simulation components
Here's the original idea for pluggable components (now perhaps outdated):
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

## TrackSource

The `TrackSource` protocol requires a `generate_event` method which returns a `Track`. `Track` implements everything in the old C structs `Track` *and* `TrackSegment`.
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
replaces
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
The idea of removing `TrackSegment` as a class is that vectorization in numpy arrays greatly improves performance. (We may want to make the link between the C fields and the python fields more parallel to aid in assessing the veridity of trumpy.)

#### Implementation in `example_run.py`

For the test code `example_run.py` the `TrackSource` source protocol was implemented as such:
```python
class VerticalShowerGenerator:
    def __init__(self, log_e: float):
        self.log_e = log_e
        self.gh = GaisserHillasProfile(x0=-75.8, xmax=773.2, nmax=6.692e9, lambda_inv=59.9)
    def generate_event(self, rng: np.random.Generator, timestamp: float) -> Track:
        n_segments = 100
        slant_depths = np.linspace(0, 1200, n_segments)
        impact = np.array([0.0, 0.0, 1400.0])
        uv = np.array([0.0, 0.0, -1.0])
        distances = np.linspace(30000, 0, n_segments)
        segment_positions = impact + distances[:, np.newaxis] * (uv * -1.0)
        # Build the Track using our real dataclass
        return Track(species=1, log_e=self.log_e, zenith=0.0, impact_v=impact,
            track_uv=uv, positions=segment_positions, n_segments=n_segments,
            time_gen=np.linspace(0, 30000, n_segments),  # 30 us track time
            altitude=30000.0 - (slant_depths * 20.0),    # Dummy altitude conversion
            dl_seg=np.full(n_segments, 300.0),
            de_dep=self.gh.evaluate_dedep(slant_depths)  # Real energy deposit!
        )
```