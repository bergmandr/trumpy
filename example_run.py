import numpy as np
import awkward as ak

# Import the modules we just built
from trumpy.track import Track
from trumpy.shower import GaisserHillasProfile
from trumpy.electronics import TAFADCFrontEnd, PMTState
from trumpy.simulation import Simulation

# =====================================================================
# 1. Concrete Schedule (Just run 2 trials)
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
    def __init__(self, log_e: float):
        self.log_e = log_e
        # Real GH Profile for a proton!
        self.gh = GaisserHillasProfile(x0=-75.8, xmax=773.2, nmax=6.692e9, lambda_inv=59.9)

    def generate_event(self, rng: np.random.Generator, timestamp: float) -> Track:
        n_segments = 100
        slant_depths = np.linspace(0, 1200, n_segments)
        impact = np.array([0.0, 0.0, 1400.0])
        uv = np.array([0.0, 0.0, -1.0])
        # Mock 3D segment positions (tracing backward from impact)
        distances = np.linspace(30000, 0, n_segments)
        segment_positions = impact + distances[:, np.newaxis] * (uv * -1.0)
        # Build the Track using our real dataclass
        return Track(
            species=1,
            log_e=self.log_e,
            zenith=0.0,
            impact_v=impact,
            track_uv=uv,
            positions=segment_positions,
            n_segments=n_segments,
            time_gen=np.linspace(0, 30000, n_segments),  # 30 us track time
            altitude=30000.0 - (slant_depths * 20.0),    # Dummy altitude conversion
            dl_seg=np.full(n_segments, 300.0),
            de_dep=self.gh.evaluate_dedep(slant_depths)  # Real energy deposit!
        )

# =====================================================================
# 3. Dummy Atmosphere (No extinction, flat yield)
# =====================================================================
class DummyAtmosphere:
    def get_fluorescence_yield(self, altitudes, de_dep, wavelengths):
        # Just assume 5 photons emitted per MeV deposited
        return de_dep[:, np.newaxis] * np.ones_like(wavelengths) * 5.0

    def get_transmission(self, emission_points, mirror_centers, wavelengths):
        N, M, W = len(emission_points), len(mirror_centers), len(wavelengths)
        transmission = np.ones((N, M, W))
        
        # Mirror area in square meters (approximate TA/HiRes mirror)
        mirror_area = 3.0 
        
        for m in range(M):
            # Calculate distance squared from each segment to this mirror
            diff = emission_points - mirror_centers[m]
            dist_sq = np.sum(diff**2, axis=1)
            
            # Solid angle fraction: A / (4 * pi * r^2)
            # Clip at 1.0 to prevent divide-by-zero if the track hits the mirror exactly
            dist_sq = np.clip(dist_sq, 1.0, None)
            solid_angle_fraction = mirror_area / (4.0 * np.pi * dist_sq)
            
            # Apply this scalar to all wavelength bands
            transmission[:, m, :] *= solid_angle_fraction[:, np.newaxis]
            
        return transmission

# =====================================================================
# 4. Dummy Experiment (1 PMT, real TA Electronics)
# =====================================================================
class SingleTubeExperiment:
    def __init__(self):
        self.wavelength_bands = np.array([350.0])  # 1 wavelength band
        self.mirror_centers = np.array([[10000.0, 0.0, 1400.0]]) # 1 mirror
        # Real PMT State (Gain = 1.0, Pedestal = 100 FADC counts)
        self.pmt_state = PMTState(
            gain=np.array([1.0]), 
            pedestal=np.array([100.0]), 
            pedestal_rms=np.array([2.0]), 
            live_flag=np.array([True])
        )
        # Real TA FADC Pipeline (10MHz sampling, 56ns shaping, 9MHz NSB)
        t = np.arange(0, 1000, 1.0)
        pulse = (t / 56.0)**5.0 * np.exp(-t / 56.0)
        pulse /= np.sum(pulse)
        self.electronics = TAFADCFrontEnd(
            sample_rate_mhz=10.0, 
            trace_length_bins=256, 
            impulse_response=pulse, 
            nsb_rate_mhz=9.0
        )

    def passes_fast_cuts(self, track: Track) -> bool:
        return True

    def trace_photons(self, track: Track, photons_at_mirrors: np.ndarray, rng: np.random.Generator):
        # photons_at_mirrors is now a realistic expected value (e.g., ~1500.4)
        expected_photons = np.sum(photons_at_mirrors)
        
        # Draw the actual integer number of photons using Poisson statistics
        total_photons = rng.poisson(expected_photons)
        
        print(f"Tracing {total_photons} photons hitting the mirror...")
        
        # Simulate photon arrival times spread over the 30,000 ns track duration
        arrival_times = rng.uniform(0, 30000, size=total_photons)
        
        # Return a jagged Awkward array: [ [t1, t2, t3...] ] for PMT 0
        return ak.Array([arrival_times])

    def process_electronics(self, pe_times: ak.Array, rng: np.random.Generator):
        # Run the real DSP pipeline
        readout = self.electronics.digitize(pe_times, self.pmt_state, rng)
        # Dummy trigger: always save the event
        return readout, True

# =====================================================================
# RUN THE SIMULATION
# =====================================================================
def main():
    print("Initializing Simulation...")
    sim = Simulation(
        experiment=SingleTubeExperiment(),
        atmosphere=DummyAtmosphere(),
        source=VerticalShowerGenerator(log_e=19.0),
        schedule=FixedTrialsSchedule(n_trials=2),
        seed=1337
    )
    
    print("Running...")
    results = sim.run()
    
    print("\n=== SIMULATION RESULTS ===")
    print(f"Awkward Array Schema:\n{results.type}")
    
    # Inspect the specific FADC readout for the first event
    print("\n--- Event 0 Data ---")
    print(f"Log Energy: {results[0].log_e}")
    
    readout = results[0].readout
    print(f"PMT IDs hit: {readout.pmt_id.tolist()}")
    
    # Print the first 20 bins of the FADC trace for PMT 0
    fadc_trace = readout.waveform[0].tolist()
    print(f"FADC Trace (First 20 bins, including 100 count pedestal and 9MHz NSB noise):")
    print(fadc_trace[:20])

if __name__ == "__main__":
    main()