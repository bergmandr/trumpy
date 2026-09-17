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