import numpy as np
import trumpy as tp

# 1. Atmosphere: Decoupled profile + aerosol extinction
atmos = tp.Atmosphere.from_radiosonde(
    sounding_file="raobslc_20141106.dst.gz",
    aerosol=tp.AerosolModel(hal=1.0, scale_height=1000.0),
    fluorescence_model="kakimoto",
)

# 2. Experiment: Geometry, optics, calibration, and electronics
detector = tp.Experiment.from_config(
    "detector_black_rock.yaml"
)  # Or programmatically assembled

# 3. Source & Track Generator: Polymorphic track producers
# Option A: Standard Cosmic Ray Air Showers
source = tp.sources.AirShowerGenerator(
    energy_range=(17.0, 20.0),
    spectral_index=3.0,
    primary="proton",
    zenith_range_deg=(0.0, 70.0),
    core_range_m=(( -20000, 20000 ), ( -20000, 20000 )),
)

# Option B: Upward Laser (CLF / FSL)
# source = tp.sources.UpwardLaser(
#     origin_clf=(0.0, 0.0, 0.0),
#     wavelength_nm=355.0,
#     energy_mj=2.0,
# )

# Option C: Direct Beam Deposit (ELS FLUKA input)
# source = tp.sources.BeamPointFile("els_deposit.txt")

# 4. Schedule & Timing: Controls trial count or on-time data matching
# Option A: Monte Carlo Trial Mode
schedule = tp.RunSchedule.fixed_trials(
    max_trials=10000, max_triggers=500, seed=42
)

# Option B: Real Run Replay (matching detector live-time and parts)
# schedule = tp.RunSchedule.from_ontime_db(
#     ontime_file="fdped_run.dst",
#     start_time="2020-05-15T04:00:00",
#     duration_sec=7200,
# )

# 5. Simulation Orchestrator
sim = tp.Simulation(
    experiment=detector,
    atmosphere=atmos,
    source=source,
    schedule=schedule,
    fast_cuts=[tp.cuts.EvsRpCut()],  # Skip raytracing for non-viable geometries
)

# Run the simulation -> returns an Awkward record array
results = sim.run(show_progress=True)

# 6. Output Handling
results.to_parquet("sim_output.parquet")
# or write to DST / awkward-compatible format