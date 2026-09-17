import numpy as np
from dataclasses import dataclass

@dataclass
class GaisserHillasProfile:
    """Evaluates the longitudinal development of an air shower."""
    x0: float          # Depth of first interaction (g/cm^2)
    xmax: float        # Depth of maximum (g/cm^2)
    nmax: float        # Number of particles at maximum
    lambda_inv: float  # Attenuation length (g/cm^2)

    def evaluate_particles(self, x: np.ndarray) -> np.ndarray:
        """
        Vectorized evaluation of the Gaisser-Hillas profile.
        Returns 0 for depths prior to x0.
        """
        n_particles = np.zeros_like(x, dtype=np.float64)
        
        # Only evaluate where the shower has actually started
        valid = x > self.x0
        x_valid = x[valid]
        
        m = (self.xmax - self.x0) / self.lambda_inv
        
        term1 = (x_valid - self.x0) / (self.xmax - self.x0)
        term2 = np.exp((self.xmax - x_valid) / self.lambda_inv)
        
        n_particles[valid] = self.nmax * (term1 ** m) * term2
        
        return n_particles

    def evaluate_dedep(self, x: np.ndarray, avg_dedx_mev: float = 2.2) -> np.ndarray:
        """
        Approximates the energy deposit profile (MeV / g/cm^2).
        avg_dedx_mev scales the particle count by the average energy loss rate.
        """
        return self.evaluate_particles(x) * avg_dedx_mev