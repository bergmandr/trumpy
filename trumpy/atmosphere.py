import numpy as np
from scipy.interpolate import CubicSpline
from dataclasses import dataclass
from fluorescence import FluorescenceModel

class AtmosphericProfile:
    """Continuous atmospheric profile via cubic splines."""
    
    def __init__(self, altitudes, densities, pressures, temperatures):
        # Sort by altitude just in case
        idx = np.argsort(altitudes)
        alt, rho, P, T = altitudes[idx], densities[idx], pressures[idx], temperatures[idx]
        
        # Fit continuous splines
        self._rho_spline = CubicSpline(alt, rho, extrapolate=True)
        self._P_spline = CubicSpline(alt, P, extrapolate=True)
        self._T_spline = CubicSpline(alt, T, extrapolate=True)
        
        # Grammage is the integral of density from altitude h to infinity.
        # We can create a grammage spline by integrating the density spline!
        # SciPy splines have an .antiderivative() method.
        density_integral = self._rho_spline.antiderivative()
        
        # X(h) = integral from h to top_of_atmosphere
        top_of_atmos = 60000.0  # meters, from legacy constants.h
        total_depth = density_integral(top_of_atmos)
        
        # Evaluate cumulative grammage at our sample altitudes and fit a spline
        grammages = total_depth - density_integral(alt)
        self._grammage_spline = CubicSpline(alt, grammages, extrapolate=True)

    # Vectorized public accessors
    def density(self, altitude):
        return self._rho_spline(altitude)

    def pressure(self, altitude):
        return self._P_spline(altitude)

    def temperature(self, altitude):
        return self._T_spline(altitude)

    def grammage(self, altitude):
        """Returns vertical grammage (g/cm^2) at a given altitude."""
        return self._grammage_spline(altitude)

@dataclass
class AerosolModel:
    hal: float              # Horizontal Attenuation Length (meters) at reference wavelength
    scale_height: float     # Aerosol scale height (meters)
    ref_wavelength: float   # Reference wavelength (e.g., 355e-9 meters)
    gamma: float = 1.0      # Wavelength dependence exponent (typically ~1.0 for Mie)

    def mie_optical_depth(self, alt_a, alt_b, distances, wavelengths):
        """
        Calculates Mie optical depth for (N, M) paths and (W,) wavelengths.
        alt_a: shape (N, 1)
        alt_b: shape (1, M)
        distances: shape (N, M)
        wavelengths: shape (W,)
        """
        dh = alt_b - alt_a
        
        # Analytic integral of exp(-h / scale_height)
        # Handle perfectly horizontal rays (dh == 0) with np.where
        effective_path = np.where(
            np.abs(dh) > 1e-3,
            distances * (self.scale_height / dh) * (
                np.exp(-alt_a / self.scale_height) - np.exp(-alt_b / self.scale_height)
            ),
            distances * np.exp(-alt_a / self.scale_height) # Horizontal limit
        )
        
        # tau_base shape: (N, M)
        tau_base = effective_path / self.hal
        
        # Apply wavelength scaling (lambda_ref / lambda)^gamma
        # Broadcast (N, M, 1) * (W,) -> (N, M, W)
        wl_scaling = (self.ref_wavelength / wavelengths) ** self.gamma
        
        return tau_base[..., np.newaxis] * wl_scaling

class Atmosphere:
    """Manages the atmospheric state and computes optical transmission and fluorescence yield."""

    def __init__(self, profile, aerosol, fy_model: FluorescenceModel):
        self.profile = profile
        self.aerosol = aerosol
        self.fy_model = fy_model

    def get_transmission(self, emission_points, mirror_centers, wavelengths):
        """
        emission_points: (N, 3) 
        mirror_centers: (M, 3)
        wavelengths: (W,)
        
        Returns:
            transmission: (N, M, W)
        """
        # 1. Coordinate Extraction & Broadcasting
        # alt_a shape: (N, 1), alt_b shape: (1, M)
        alt_a = emission_points[:, 2][:, np.newaxis]
        alt_b = mirror_centers[:, 2][np.newaxis, :]
        
        # Distances shape: (N, M)
        diff_vectors = emission_points[:, np.newaxis, :] - mirror_centers[np.newaxis, :, :]
        distances = np.linalg.norm(diff_vectors, axis=2)
        
        # 2. Rayleigh Optical Depth (N, M)
        dh = alt_a - alt_b
        grammage_a = self.profile.grammage(alt_a)
        grammage_b = self.profile.grammage(alt_b)
        
        # Average density = dX / dh
        # Multiply by distance to get total slant grammage (g/cm^2)
        slant_grammage = np.where(
            np.abs(dh) > 1e-3,
            np.abs(grammage_a - grammage_b) / np.abs(dh) * distances,
            self.profile.density(alt_a) * distances * 0.1 # 0.1 converts kg/m^3 * m to g/cm^2
        )
        
        tau_ray_base = slant_grammage / self.rayleigh_ref_grammage
        wl_scaling_ray = (self.rayleigh_ref_wl / wavelengths) ** 4
        
        # Broadcast to (N, M, W)
        tau_ray = tau_ray_base[..., np.newaxis] * wl_scaling_ray
        
        # 3. Mie Optical Depth (N, M, W)
        tau_mie = self.aerosol.mie_optical_depth(alt_a, alt_b, distances, wavelengths)
        
        # 4. Total Transmission
        return np.exp(-(tau_ray + tau_mie))

    def get_fluorescence_yield(self, altitudes: np.ndarray, dedep: np.ndarray, wavelengths: np.ndarray) -> np.ndarray:
        """
        Calculates the total fluorescence photons produced across all segments.
        
        altitudes: (N,) array of segment altitudes
        dedep: (N,) array of energy deposited in MeV for each segment
        wavelengths: (W,) array of wavelength band centers
        
        Returns: (N, W) array of total photons emitted per segment per wavelength
        """
        # 1. Query the continuous splines for local conditions
        rho = self.profile.density(altitudes)
        temp = self.profile.temperature(altitudes)
        
        # 2. Get the specific yield (photons / MeV) from the injected model
        # Shape: (N, W)
        specific_yield = self.fy_model.yield_per_mev(rho, temp, wavelengths)
        
        # 3. Multiply by the actual energy deposit (N,) -> broadcast to (N, W)
        total_photons = specific_yield * dedep[:, np.newaxis]
        
        return total_photons
    