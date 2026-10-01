import numpy as np
import awkward as ak
from scipy.interpolate import CubicSpline
from dataclasses import dataclass
from trumpy.fluorescence import FluorescenceModel

# Physical constants
UNIV_GAS_CONST = 8.31432      # J/(mol K)
AIR_MOLAR_MASS = 0.0289644    # kg/mol
G0 = 9.80665                  # m/s^2 (standard gravity)

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
        density_integral = self._rho_spline.antiderivative()
        
        # X(h) = integral from h to top_of_atmosphere
        top_of_atmos = 60000.0  # meters, from legacy constants.h
        total_depth = density_integral(top_of_atmos)
        
        # Evaluate cumulative grammage at our sample altitudes and fit a spline
        # 0.1 factor converts kg/m^2 to g/cm^2
        grammages = (total_depth - density_integral(alt)) * 0.1
        self._grammage_spline = CubicSpline(alt, grammages, extrapolate=True)

    # Vectorized public accessors
    def get_density(self, altitude):
        return self._rho_spline(altitude)

    def get_pressure(self, altitude):
        return self._P_spline(altitude)

    def get_temperature(self, altitude):
        return self._T_spline(altitude)

    def get_grammage(self, altitude):
        """Returns vertical grammage (g/cm^2) at a given altitude."""
        return self._grammage_spline(altitude)
        
    def get_slant_depth(self, altitude, zenith_angle):
        """Flat earth approximation for tests."""
        return self.get_grammage(altitude) / np.cos(zenith_angle)

class GDASProfile(AtmosphericProfile):
    """Interpolates atmospheric profiles from awkward array GDAS records."""
    
    def __init__(self, gdas_records: ak.Array):
        # Record names align with legacy fdatmos_param_dst_common
        altitudes = ak.to_numpy(gdas_records["alt"])
        densities = ak.to_numpy(gdas_records["rho"])
        pressures = ak.to_numpy(gdas_records["pres"])
        temperatures = ak.to_numpy(gdas_records["temp"])
        
        # Inherit the spline generator for continuity
        super().__init__(altitudes, densities, pressures, temperatures)

class HydrostaticProfile(AtmosphericProfile):
    """Hydrostatic atmosphere model supporting layers with linear temperature gradients."""
    
    def __init__(self, layers=None):
        if layers is None:
            # Columns: [Base Altitude (m), Base Temp (K), Base Pressure (Pa), Lapse Rate (K/m)]
            self.layers = np.array([
                [0.0,     288.15, 101325.0,  -0.0065],
                [11000.0, 216.65, 22632.1,    0.0],
                [20000.0, 216.65, 5474.89,    0.001],
                [32000.0, 228.65, 868.019,    0.0028],
                [47000.0, 270.65, 110.906,    0.0],
                [51000.0, 270.65, 66.9389,   -0.0028],
                [71000.0, 214.65, 3.95642,   -0.002]
            ])
        else:
            self.layers = np.array(layers)
            
        # Sample over a fine grid to create the baseline arrays for the spline initialization
        altitudes = np.linspace(0, 100000.0, 1000)
        idx = np.searchsorted(self.layers[:, 0], altitudes, side='right') - 1
        idx = np.clip(idx, 0, len(self.layers) - 1)
        
        h_b = self.layers[idx, 0]
        T_b = self.layers[idx, 1]
        P_b = self.layers[idx, 2]
        L_b = self.layers[idx, 3]
        
        temperatures = T_b + L_b * (altitudes - h_b)
        
        pressures = np.empty_like(altitudes)
        isothermal = (L_b == 0)
        gradient = ~isothermal
        
        # Hydrostatic integration
        if np.any(gradient):
            exponent = (G0 * AIR_MOLAR_MASS) / (UNIV_GAS_CONST * L_b[gradient])
            pressures[gradient] = P_b[gradient] * (T_b[gradient] / temperatures[gradient]) ** exponent
            
        if np.any(isothermal):
            exponent = -(G0 * AIR_MOLAR_MASS * (altitudes[isothermal] - h_b[isothermal])) / (UNIV_GAS_CONST * T_b[isothermal])
            pressures[isothermal] = P_b[isothermal] * np.exp(exponent)
            
        densities = (pressures * AIR_MOLAR_MASS) / (UNIV_GAS_CONST * temperatures)
        
        super().__init__(altitudes, densities, pressures, temperatures)

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

    def __init__(self, profile: AtmosphericProfile, aerosol: AerosolModel, fy_model: FluorescenceModel):
        self.profile = profile
        self.aerosol = aerosol
        self.fy_model = fy_model

    # Expose profile getters to satisfy AtmosphereModel protocol
    def get_density(self, altitudes): return self.profile.get_density(altitudes)
    def get_pressure(self, altitudes): return self.profile.get_pressure(altitudes)
    def get_temperature(self, altitudes): return self.profile.get_temperature(altitudes)
    def get_grammage(self, altitudes): return self.profile.get_grammage(altitudes)
    def get_slant_depth(self, altitudes, zenith_angle): return self.profile.get_slant_depth(altitudes, zenith_angle)

    def get_transmission(self, emission_points, mirror_centers, wavelengths):
        """
        emission_points: (N, 3) 
        mirror_centers: (M, 3)
        wavelengths: (W,)
        
        Returns:
            transmission: (N, M, W)
        """
        alt_a = emission_points[:, 2][:, np.newaxis]
        alt_b = mirror_centers[:, 2][np.newaxis, :]
        
        diff_vectors = emission_points[:, np.newaxis, :] - mirror_centers[np.newaxis, :, :]
        distances = np.linalg.norm(diff_vectors, axis=2)
        
        dh = alt_a - alt_b
        grammage_a = self.profile.get_grammage(alt_a)
        grammage_b = self.profile.get_grammage(alt_b)
        
        # Average density = dX / dh
        # Multiply by distance to get total slant grammage (g/cm^2)
        slant_grammage = np.where(
            np.abs(dh) > 1e-3,
            np.abs(grammage_a - grammage_b) / np.abs(dh) * distances,
            self.profile.get_density(alt_a) * distances * 0.1 
        )
        
        tau_ray_base = slant_grammage / self.rayleigh_ref_grammage
        wl_scaling_ray = (self.rayleigh_ref_wl / wavelengths) ** 4
        
        tau_ray = tau_ray_base[..., np.newaxis] * wl_scaling_ray
        tau_mie = self.aerosol.mie_optical_depth(alt_a, alt_b, distances, wavelengths)
        
        return np.exp(-(tau_ray + tau_mie))

    def get_fluorescence_yield(self, altitudes: np.ndarray, dedep: np.ndarray, wavelengths: np.ndarray) -> np.ndarray:
        """
        Calculates the total fluorescence photons produced across all segments.
        
        altitudes: (N,) array of segment altitudes
        dedep: (N,) array of energy deposited in MeV for each segment
        wavelengths: (W,) array of wavelength band centers
        
        Returns: (N, W) array of total photons emitted per segment per wavelength
        """
        rho = self.profile.get_density(altitudes)
        temp = self.profile.get_temperature(altitudes)        
        specific_yield = self.fy_model.yield_per_mev(rho, temp, wavelengths)
        total_photons = specific_yield * dedep[:, np.newaxis]
        
        return total_photons
    