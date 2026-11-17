# `trumpy` Development Status

This document summarizes the current status of the `trumpy` Python port relative to the legacy C `trump` codebase, along with the roadmap for the next major subsystems to port.

---

## 1. Current Module Status

| Python Module | Legacy C Equivalent(s) | Current Status |
| :--- | :--- | :--- |
| [`trumpy/simulation.py`](trumpy/simulation.py) | [`trump.c`](trumpy/legacy/bin/trump.c), [`control.h`](trumpy/legacy/inc/control.h) | **Core loop implemented.** Defines protocols (`TrackSource`, `RunSchedule`, `AtmosphereModel`, `Experiment`) and orchestrates the 7-step loop returning an Awkward `ak.Array`. Currently assumes a single `Experiment` and skips Cherenkov production/scattering in the main loop. |
| [`trumpy/track.py`](trumpy/track.py) | [`track.h`](trumpy/legacy/inc/track.h), [`airshower.h`](trumpy/legacy/inc/airshower.h), [`fdsite.h`](trumpy/legacy/inc/fdsite.h) | **Data structures ported (SoA).** `Track` merges `AirShower`, `Track`, and `TrackSegment` into vectorized NumPy arrays (tested in [`tests/test_track.py`](tests/test_track.py)). `ObservedTrack` merges `ObservedTrack` and `ObservedTrackSegment`, though it is not yet wired into [`trumpy/simulation.py`](trumpy/simulation.py). |
| [`trumpy/shower.py`](trumpy/shower.py) | [`gaisser-hillas.c`](trumpy/legacy/src/toolbox/gaisser-hillas.c), [`nerling.c`](trumpy/legacy/src/airshower/nerling.c) | **Partial.** `GaisserHillasProfile` evaluates $N(X)$ and constant-rate $dE/dX$ (`avg_dedx_mev = 2.2`, tested in [`tests/test_track.py`](tests/test_track.py)). Missing Nerling age-dependent $dE/dX$, NKG lateral distribution ([`nkg.c`](trumpy/legacy/src/airshower/nkg.c)), and shower library sampling ([`showerlib.c`](trumpy/legacy/src/airshower/showerlib.c)). |
| [`trumpy/atmosphere.py`](trumpy/atmosphere.py) | [`atmosphere-withdb.c`](trumpy/legacy/src/airshower/atmosphere-withdb.c), [`scatter.c`](trumpy/legacy/src/airshower/scatter.c) | **Mostly implemented.** `AtmosphericProfile`, `GDASProfile`, `HydrostaticProfile`, `AerosolModel`, and `Atmosphere` are implemented with cubic splines and vectorized Mie/Rayleigh optical depth. *(Note: `Atmosphere.__init__` references `self.rayleigh_ref_grammage` and `self.rayleigh_ref_wl`, which still need to be initialized in `__init__`).* |
| [`trumpy/fluorescence.py`](trumpy/fluorescence.py) | [`flash.c`](trumpy/legacy/src/airshower/flash.c), [`flash.h`](trumpy/legacy/inc/flash.h) | **Stubbed.** Defines `FluorescenceModel` protocol and `KakimotoYield` skeleton (currently returns a dummy `5.0` photons/MeV/band). |
| [`trumpy/telescope.py`](trumpy/telescope.py) & [`trumpy/camera.py`](trumpy/camera.py) | [`fdsite.h`](trumpy/legacy/inc/fdsite.h), [`raytrace.h`](trumpy/legacy/inc/raytrace.h) (`TGeom`), [`tacalibration.h`](trumpy/legacy/inc/tacalibration.h) | **Data structures only.** Defines `TelescopeGeometry`, `MirrorTopology` (`IdealSphericalMirror`, `SegmentedMirror`), `OpticalLayer`, and `CameraAssembly`. Not yet connected to ray-tracing methods or site geometry loaders. |
| [`trumpy/optics.py`](trumpy/optics.py) | [`trace.c`](trumpy/legacy/src/raytrace/trace.c), [`shadow.c`](trumpy/legacy/src/raytrace/shadow.c), [`hitMirror.c`](trumpy/legacy/src/raytrace/hitMirror.c), [`camera.c`](trumpy/legacy/src/raytrace/camera.c) | **Standalone vectorized primitives.** Has `mirror_intersection`, `reflect_and_blur`, `map_to_pmts`, `is_within_aperture`, and shadow checks (`is_shadowed`, `is_shadowed_by_box`, `is_shadowed_by_cylinder_strut`). Missing Snell shifts (`snellShift`), hexagonal PMT boundary/crack checks (`pmtDetect`), and PMT QE/uniformity maps (`pmtQE`). |
| [`trumpy/electronics.py`](trumpy/electronics.py) | [`taelectronics.c`](trumpy/legacy/src/fadctel/taelectronics.c), [`tlelectronics.h`](trumpy/legacy/inc/tlelectronics.h) | **TA FADC digitization implemented; triggers & HiRes stubbed.** `TAFADCFrontEnd.digitize` and `generate_ta_impulse_response` work (used in [`example_run.py`](example_run.py)), but hardware trigger logic (`simTATriggerResponse`: tube threshold, 5-tube / 3-tube pattern recognition, adjacent camera readout) is not yet ported. |

---

## 2. Next Major Subsystems to Port

### Priority 1: Detector / `FDSite` (`Experiment`) & Ray-Tracing Orchestrator (`fdsite.c` + `acpttrack.c` + `doRayTrace.c`)
Currently, low-level optical primitives exist in [`trumpy/optics.py`](trumpy/optics.py) and detector dataclasses exist in [`trumpy/telescope.py`](trumpy/telescope.py), [`trumpy/camera.py`](trumpy/camera.py), and [`trumpy/electronics.py`](trumpy/electronics.py), but they are not yet connected into a concrete `Experiment` (`FDSite`) implementation:
1. **[`fdsite.c`](trumpy/legacy/src/raytrace/fdsite.c) (`getObservedTrack`, `applyEvsRpCut`):**
   - Transform `Track` from CLF coordinates into site/camera coordinates (`ObservedTrack`).
   - Compute $R_p$, $\psi$, viewing angles, mirror field-of-view culling ($\cos\theta > \cos 15^\circ$), solid-angle geometric factor $A \cos\theta / (4\pi r^2)$, and fast $E$ vs. $R_p$ cuts (`passes_fast_cuts`).
2. **[`acpttrack.c`](trumpy/legacy/src/raytrace/acpttrack.c) (`getAcceptedTrack`) & [`doRayTrace.c`](trumpy/legacy/src/raytrace/doRayTrace.c) (`doRayTrace`):**
   - Take `ObservedTrack`, apply wavelength reduction factors (`reductFactor`: mirror reflectivity $\times$ Paraglas $\times$ BG3 $\times$ PMT QE), Poisson-sample photon counts per segment/mirror, spread emission points along the segment and transversely via NKG ([`nkg.c`](trumpy/legacy/src/airshower/nkg.c)), invoke the vectorized ray tracer ([`trumpy/optics.py`](trumpy/optics.py)), compute arrival times (`tgen + ray.time`), and output the jagged `PETimes` (`ak.Array`) per camera and PMT.

### Priority 2: Concrete `TrackSource` / Air Shower Builder (`airshower.c` + `track.c` + `nerling.c` + `nkg.c`)
Replace the mock `VerticalShowerGenerator` in [`example_run.py`](example_run.py) with a physical track builder ([`track.c`](trumpy/legacy/src/airshower/track.c) `buildTrackGeom` and `fillTrackLight`, plus [`airshower.c`](trumpy/legacy/src/airshower/airshower.c) `createAirShower`):
- Randomly sample $E$, zenith, azimuth, and core impact position (`AirShowerGenerator`).
- Step along the track in grammage increments ($1\text{ g/cm}^2$), computing segment lengths `dlseg`, density-weighted midpoints `dlmid`, Molière radius `molrad`, Nerling energy deposit `dedep` ([`nerling.c`](trumpy/legacy/src/airshower/nerling.c)), and Cherenkov beam production/scattering (`pcv`, `ncv`).

### Priority 3: Real Fluorescence & Cherenkov Physics (`flash.c` + `scatter.c`)
- Port the FLASH / Kakimoto spectrum and collisional quenching from [`flash.c`](trumpy/legacy/src/airshower/flash.c) into [`trumpy/fluorescence.py`](trumpy/fluorescence.py) (replacing the dummy `5.0` yield in `KakimotoYield`), and initialize the Rayleigh reference attributes in [`trumpy/atmosphere.py`](trumpy/atmosphere.py).
- Add Direct, Rayleigh-scattered, and Mie-scattered Cherenkov light ([`track.c`](trumpy/legacy/src/airshower/track.c), [`fdsite.c`](trumpy/legacy/src/raytrace/fdsite.c), [`scatter.c`](trumpy/legacy/src/airshower/scatter.c)) to populate `ObservedTrack` (`n_cv_dir`, `n_cv_mie`, `n_cv_ray`).

### Priority 4: Hardware Trigger Simulator (`taelectronics.c` / `trigger-simulator.h`)
- Port the TA trigger recognition (`simTATriggerResponse` / `fadcSig` / 5-tube & 3-tube patterns in [`taelectronics.c`](trumpy/legacy/src/fadctel/taelectronics.c)) into [`trumpy/electronics.py`](trumpy/electronics.py) so `process_electronics` evaluates hardware trigger conditions rather than returning `is_triggered = True`.
