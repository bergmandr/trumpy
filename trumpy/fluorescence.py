import numpy as np
from typing import Protocol

class FluorescenceModel(Protocol):
    def yield_per_mev(self, density: np.ndarray, temperature: np.ndarray, wavelengths: np.ndarray) -> np.ndarray:
        """
        Returns the number of photons produced per MeV of energy deposit.
        density: (N,) array of local atmospheric densities
        temperature: (N,) array of local atmospheric temperatures
        wavelengths: (W,) array of wavelength bands
        Returns: (N, W) array of photon yields
        """
        ...

class KakimotoYield:
    """Kakimoto fluorescence yield parameterization."""
    
    def __init__(self):
        # We would store the specific constants for the Kakimoto spectral lines,
        # collisional quenching cross-sections, and reference temperatures here.
        pass
        
    def yield_per_mev(self, density: np.ndarray, temperature: np.ndarray, wavelengths: np.ndarray) -> np.ndarray:
        # 1. Calculate pressure from density and temperature (ideal gas law or empirical)
        # 2. Calculate the quenching factors for each wavelength band
        # 3. Calculate the absolute yield
        
        # Placeholder for the actual vectorized math
        N = density.shape[0]
        W = wavelengths.shape[0]
        
        # returns an array of shape (N, W)
        return np.ones((N, W)) * 5.0 # Dummy value: 5 photons/MeV/band