import numpy as np
import awkward as ak
from typing import Protocol, Iterator, Tuple

# We import the Track dataclass we made earlier
from trumpy.track import Track

# --- Protocols defining the expected interfaces ---

class TrackSource(Protocol):
    def generate_event(self, rng: np.random.Generator, timestamp: float) -> Track:
        ...

class RunSchedule(Protocol):
    def get_trials(self) -> Iterator[Tuple[int, float]]:
        """Yields (trial_id, timestamp_seconds) until the run is complete."""
        ...

class AtmosphereModel(Protocol):
    def get_fluorescence_yield(self, altitudes: np.ndarray, de_dep: np.ndarray, wavelengths: np.ndarray) -> np.ndarray:
        ...
    def get_transmission(self, emission_points: np.ndarray, mirror_centers: np.ndarray, wavelengths: np.ndarray) -> np.ndarray:
        ...

class Experiment(Protocol):
    wavelength_bands: np.ndarray
    mirror_centers: np.ndarray
    
    def passes_fast_cuts(self, track: Track) -> bool:
        """Evaluates EvsRp and FOV cuts to quickly discard doomed tracks."""
        ...
    def trace_photons(self, track: Track, photons_at_mirrors: np.ndarray, rng: np.random.Generator) -> ak.Array:
        """Returns jagged PE arrival times grouped by PMT."""
        ...
    def process_electronics(self, pe_times: ak.Array, rng: np.random.Generator) -> ak.Array:
        """Digitizes waveforms and evaluates the hardware trigger."""
        ...


# --- The Main Orchestrator ---

class Simulation:
    """Orchestrates the Monte Carlo execution pipeline."""
    
    def __init__(
        self,
        experiment: Experiment,
        atmosphere: AtmosphereModel,
        source: TrackSource,
        schedule: RunSchedule,
        seed: int = 42
    ):
        self.experiment = experiment
        self.atmosphere = atmosphere
        self.source = source
        self.schedule = schedule
        
        # Initialize the global random number generator for this run
        self.rng = np.random.default_rng(seed)

    def run(self) -> ak.Array:
        """
        Executes the main event loop and returns an Awkward record array of triggered events.
        """
        triggered_events = []

        for trial_id, timestamp in self.schedule.get_trials():
            
            # 1. Generate the physical track
            track = self.source.generate_event(self.rng, timestamp)

            # 2. Fast geometric filters (e.g., E vs Rp cuts)
            if not self.experiment.passes_fast_cuts(track):
                continue

            # 3. Light Emission (Fluorescence)
            # Yield shape: (n_segments, n_wavelengths)
            photon_production = self.atmosphere.get_fluorescence_yield(
                track.altitude, track.de_dep, self.experiment.wavelength_bands
            )

            # 4. Atmospheric Transmission
            # Transmission shape: (n_segments, n_mirrors, n_wavelengths)
            transmission = self.atmosphere.get_transmission(
                track.positions, self.experiment.mirror_centers, self.experiment.wavelength_bands
            )
            
            # Broadcast to find actual photons reaching the apertures
            photons_at_mirrors = photon_production[:, np.newaxis, :] * transmission

            # 5. Ray Tracing (Optics)
            # pe_times is a jagged Awkward array grouped by PMT ID
            pe_times = self.experiment.trace_photons(track, photons_at_mirrors, self.rng)

            # 6. Electronics & Trigger
            # Applies analog shaping, PMT gains, noise, FADC sampling, and evaluates trigger
            readout, is_triggered = self.experiment.process_electronics(pe_times, self.rng)

            # 7. Data Collection
            if is_triggered:
                # We package the event into a dictionary that maps perfectly to Awkward records
                event_data = {
                    "trial_id": trial_id,
                    "timestamp": timestamp,
                    "primary_species": track.species,
                    "log_e": track.log_e,
                    "readout": readout  # Contains pmt_ids, FADC traces, etc.
                }
                triggered_events.append(event_data)

        # Convert the list of dicts into a unified Awkward array
        if not triggered_events:
            print("Run completed: 0 triggered events.")
            return ak.Array([])
            
        print(f"Run completed: {len(triggered_events)} triggered events.")
        return ak.from_iter(triggered_events)