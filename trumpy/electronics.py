from typing import Protocol
from dataclasses import dataclass
import numpy as np
import awkward as ak
import scipy.signal as signal

@dataclass
class PMTState:
    """Vectorized calibration data for the camera's PMTs."""
    gain: np.ndarray          # shape (n_tubes,)
    pedestal: np.ndarray      # shape (n_tubes,) mean pedestal 
    pedestal_rms: np.ndarray  # shape (n_tubes,) pedestal variance
    live_flag: np.ndarray     # shape (n_tubes,) boolean mask

class ElectronicsFrontEnd(Protocol):
    def digitize(self, pe_times: ak.Array, pmt_state: PMTState, rng: np.random.Generator):
        """
        Takes jagged PE arrival times and returns the digitized readout.
        """
        ...

class TAFADCFrontEnd:
    """
    Telescope Array FADC Electronics.
    Standard continuous digitization (e.g., 10 MHz sampling / 100 ns slices).
    """
    def __init__(
        self, 
        sample_rate_mhz: float, 
        trace_length_bins: int, 
        impulse_response: np.ndarray,
        nsb_rate_mhz: float  # Night Sky Background rate in photons per microsecond
    ):
        self.sample_rate_mhz = sample_rate_mhz
        self.dt_ns = 1000.0 / sample_rate_mhz
        self.trace_length = trace_length_bins
        self.impulse_response = impulse_response 
        
        # Convert NSB rate to expected photons per 1 ns bin
        self.nsb_rate_per_ns = nsb_rate_mhz / 1000.0

    def digitize(self, pe_times: ak.Array, pmt_state: PMTState, rng: np.random.Generator):
        n_tubes = len(pe_times)
        max_time_ns = int(self.trace_length * self.dt_ns)
        
        # 1. High-Res Analog Histogram (1 ns bins)
        flat_times = ak.to_numpy(ak.flatten(pe_times))
        
        pmt_indices = ak.Array(np.arange(n_tubes))
        broadcasted_indices = ak.broadcast_arrays(pmt_indices, pe_times)[0]
        tube_indices = ak.to_numpy(ak.flatten(broadcasted_indices))
        
        analog_signals, _, _ = np.histogram2d(
            tube_indices, flat_times, 
            bins=[n_tubes, max_time_ns], 
            range=[[0, n_tubes], [0, max_time_ns]]
        )
        
        # 2. Inject Night Sky Background (Photon Poisson Noise)
        if self.nsb_rate_per_ns > 0:
            analog_signals += rng.poisson(
                lam=self.nsb_rate_per_ns, 
                size=analog_signals.shape
            )
            
        # 3. Inject Electronic Noise (Gaussian White Noise)
        # We inject it here so the FFT correlates it bin-to-bin, matching the real data.
        # We calculate the power of the impulse response to ensure the final downsampled 
        # FADC bins have exactly the RMS requested by the calibration state.
        h_power = np.sum(self.impulse_response**2)
        noise_sigma_pe = (pmt_state.pedestal_rms / pmt_state.gain) / np.sqrt(h_power)
        
        analog_signals += rng.normal(
            scale=noise_sigma_pe[:, np.newaxis], 
            size=analog_signals.shape
        )
            
        # 4. Shape the signal (Impulse Response)
        shaped_signals = signal.fftconvolve(
            analog_signals, 
            self.impulse_response[np.newaxis, :], 
            mode='same', 
            axes=1
        )
        
        # 5. Downsample to the FADC clock speed (e.g., every 100 ns)
        step = int(self.dt_ns)
        fadc_traces = shaped_signals[:, ::step]
        
        # 6. Apply Calibration (Gains, Pedestals)
        # The noise was injected in PE units, so applying the gain here correctly 
        # scales the noise to FADC counts, alongside the signal.
        fadc_traces = (fadc_traces * pmt_state.gain[:, np.newaxis]) + pmt_state.pedestal[:, np.newaxis]
        
        # 7. ADC Clipping (14-bit) and integer cast
        fadc_traces = np.clip(np.round(fadc_traces), 0, 16383).astype(np.int32)
        
        # Return masked Awkward array with depth_limit=1 to prevent broadcasting errors
        live_mask = pmt_state.live_flag > 0
        return ak.zip({
            "pmt_id": np.arange(n_tubes)[live_mask],
            "waveform": fadc_traces[live_mask]
        }, depth_limit=1)

class HiRes2FADCFrontEnd:
    """
    HiRes-II FADC Electronics.
    Stubs for the 40 MHz clock domain downsampled to 100 ns slices.
    """
    def __init__(self, trace_length_bins: int, impulse_response: np.ndarray, nsb_rate_mhz: float):
        # We will need to design the impulse_response to explicitly mimic 
        # the 40 MHz cutoff and Fourier effects you mentioned.
        self.trace_length = trace_length_bins
        self.impulse_response = impulse_response
        self.nsb_rate_per_ns = nsb_rate_mhz / 1000.0

    def digitize(self, pe_times: ak.Array, pmt_state: PMTState, rng: np.random.Generator):
        # TODO: Implement 40MHz specific clock jitter / frequency cutoff behavior
        pass


class HiRes1QDCFrontEnd:
    """
    HiRes-I Sample-and-Hold Electronics.
    Integrates charge over a gate window using RC filtering rather than tracing FADC counts.
    """
    def __init__(self, rc_decay_time_ns: float, gate_width_ns: float, nsb_rate_mhz: float):
        self.tau = rc_decay_time_ns
        self.gate_width = gate_width_ns
        self.nsb_rate_per_ns = nsb_rate_mhz / 1000.0

    def digitize(self, pe_times: ak.Array, pmt_state: PMTState, rng: np.random.Generator):
        # TODO: Implement TDC triggers, gate windows, and RC-weighted charge integration
        pass

def generate_ta_impulse_response(
    tau_ns: float = 56.0, 
    alpha: float = 5.0, 
    length_ns: int = 1000
) -> np.ndarray:
    """
    Generates the single photoelectron (SPE) impulse response for the 
    Telescope Array FD shaping amplifier at 1 ns resolution.
    
    tau_ns: Time constant (56.0 ns matches TA data)
    alpha: Shaping parameter (5.0 matches TA data)
    length_ns: Total length of the response array to generate.
    """
    t = np.arange(0, length_ns, 1.0)
    
    # Calculate the shaping pulse
    # Add a tiny epsilon to t to avoid 0^0 if alpha = 0, though alpha is 5.0 here
    response = (t / tau_ns)**alpha * np.exp(-t / tau_ns)
    
    # Normalize the area to 1 so that convolving it with the analog signal
    # strictly conserves the total photoelectron count (charge).
    response /= np.sum(response)
    
    return response