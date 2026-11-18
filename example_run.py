import numpy as np
import awkward as ak

from trumpy.track import Track
from trumpy.shower import GaisserHillasProfile, get_moliere_radius
from trumpy.camera import CameraAssembly
from trumpy.telescope import TelescopeGeometry, FDSite
from trumpy.electronics import (
    TAFADCFrontEnd,
    PMTState,
    generate_ta_impulse_response,
)
from trumpy.atmosphere import HydrostaticProfile, AerosolModel, Atmosphere
from trumpy.fluorescence import KakimotoYield
from trumpy.simulation import Simulation


# =====================================================================
# 1. Concrete Schedule (Run 2 trials)
# =====================================================================
class FixedTrialsSchedule:
    def __init__(self, n_trials: int):
        self.n_trials = n_trials

    def get_trials(self):
        for i in range(self.n_trials):
            yield i, float(i * 1000)


# =====================================================================
# 2. Concrete Track Source (Generates a vertical proton shower)
# =====================================================================
class VerticalShowerGenerator:
    def __init__(self, log_e: float, atmosphere: Atmosphere):
        self.log_e = log_e
        self.atmosphere = atmosphere
        self.gh = GaisserHillasProfile(
            x0=-75.8, xmax=773.2, nmax=6.692e9, lambda_inv=59.9
        )

    def generate_event(self, rng: np.random.Generator, timestamp: float) -> Track:
        nseg = 100
        impact = np.array([0.0, 0.0, 1400.0])
        uv = np.array([0.0, 0.0, -1.0])
        zenith = 0.0

        # Segment top positions tracing downward from 12 km altitude to ground (1.4 km)
        distances = np.linspace(10600.0, 0.0, nseg)
        segment_positions = impact + distances[:, np.newaxis] * (-uv)
        dlseg = np.full(nseg, 10600.0 / nseg)

        # Slant depth and density from the atmospheric model
        altitudes = segment_positions[:, 2]
        slant_depths = self.atmosphere.get_slant_depth(altitudes, zenith)
        safe_slant_depths = np.clip(slant_depths, 1e-3, None)
        dxseg = np.gradient(slant_depths)

        mid_altitudes = altitudes + 0.5 * dlseg * uv[2]
        densities = self.atmosphere.get_density(mid_altitudes)

        return Track(
            species=14,  # Proton (PSPEC_PROTON)
            log_e=self.log_e,
            zenith=zenith,
            impact_v=impact,
            track_uv=uv,
            positions=segment_positions,
            nseg=nseg,
            time_gen=np.linspace(0.0, 35358.0, nseg),
            dlseg=dlseg,
            dedep=self.gh.evaluate_dedep(slant_depths) * dxseg,
            age=3.0 / (1.0 + 2.0 * self.gh.xmax / safe_slant_depths),
            dlmid=0.5 * dlseg,
            dxseg=dxseg,
            nch=self.gh.evaluate_particles(slant_depths),
            molrad=get_moliere_radius(densities),
            nfl=np.zeros((nseg, 1)),
            pcv=np.zeros((nseg, 1)),
            ncv=np.zeros((nseg, 1)),
        )


# =====================================================================
# 3. Build a 256-PMT Telescope Array FD Site (FDSite)
# =====================================================================
def build_single_telescope_site() -> FDSite:
    # Standard 16x16 hexagonal TA camera with Paraglas window and BG3 filter
    camera = CameraAssembly.default_ta_camera(n_wavelengths=1)
    n_tubes = len(camera.pmt_vectors)

    # Place the FD station 10 km East of the shower core (x = 10,000 m, z = 1,400 m)
    # Point the mirror optical axis toward the shower (-x direction, elevated ~18 deg)
    axis_z = np.array([-1.0, 0.0, 0.32])
    axis_z /= np.linalg.norm(axis_z)
    axis_x = np.array([0.0, 1.0, 0.0])
    axis_y = np.cross(axis_z, axis_x)
    site2cam = np.vstack([axis_x, axis_y, axis_z])[np.newaxis, :, :]

    geometry = TelescopeGeometry(
        mirror_centers=np.array([[0.0, 0.0, 0.0]]),
        mirror_radii=np.array([6.067]),
        camera_rotations=site2cam,
        pmt_vectors=camera.pmt_vectors[np.newaxis, :, :],
        pmt_live_flags=np.ones((1, n_tubes), dtype=bool),
        mirror_pointing=axis_z[np.newaxis, :],
        mirror_diameters=np.array([3.3]),
        camera_separations=np.array([3.0]),
        spot_size_blur=np.array([0.005]),
        local_vsite=np.array([10000.0, 0.0, 1400.0]),
        site2clf=np.eye(3),
    )

    pmt_state = PMTState(
        gain=np.full(n_tubes, 1.0),
        pedestal=np.full(n_tubes, 100.0),
        pedestal_rms=np.full(n_tubes, 2.0),
        live_flag=np.ones(n_tubes, dtype=bool),
    )

    electronics = TAFADCFrontEnd(
        sample_rate_mhz=10.0,
        trace_length_bins=256,
        impulse_response=generate_ta_impulse_response(tau_ns=56.0, alpha=5.0),
        nsb_rate_mhz=9.0,
    )

    return FDSite(
        site_id=0,
        geometry=geometry,
        camera=camera,
        pmt_state=pmt_state,
        electronics=electronics,
        wavelength_bands=np.array([350.0]),
    )


# =====================================================================
# RUN THE SIMULATION
# =====================================================================
def main():
    print("Initializing Simulation with FDSite, CameraAssembly, and Atmosphere...")

    atmos = Atmosphere(
        profile=HydrostaticProfile(),
        aerosol=AerosolModel(
            hal=25000.0, scale_height=1000.0, ref_wavelength=355.0
        ),
        fy_model=KakimotoYield(),
    )

    site = build_single_telescope_site()

    sim = Simulation(
        experiment=site,
        atmosphere=atmos,
        source=VerticalShowerGenerator(log_e=19.0, atmosphere=atmos),
        schedule=FixedTrialsSchedule(n_trials=2),
        seed=1337,
    )

    print("Running...")
    results = sim.run()

    print("\n=== SIMULATION RESULTS ===")
    print(f"Awkward Array Schema:\n{results.type}")

    print("\n--- Event 0 Data ---")
    print(f"Log Energy: {results[0].log_e}")

    readout = results[0].readout
    waveforms = ak.to_numpy(readout.waveform)
    print(f"Total live PMTs read out: {len(readout.pmt_id)}")
    print(f"Waveform array shape: {waveforms.shape}")

    # Find the PMT with the highest peak FADC count above the 100-count pedestal
    peak_pmt_idx = int(np.argmax(np.max(waveforms, axis=1)))
    peak_bin_idx = int(np.argmax(waveforms[peak_pmt_idx]))
    peak_val = int(waveforms[peak_pmt_idx, peak_bin_idx])
    print(
        f"Brightest tube: PMT {peak_pmt_idx} (peak = {peak_val} FADC counts at time slice {peak_bin_idx})"
    )

    window_start = max(0, peak_bin_idx - 5)
    window_end = min(waveforms.shape[1], window_start + 20)
    print(
        f"FADC Trace for PMT {peak_pmt_idx} (bins {window_start}..{window_end - 1} around peak):"
    )
    print(waveforms[peak_pmt_idx, window_start:window_end].tolist())


if __name__ == "__main__":
    main()
