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


def get_moliere_radius(density_kg_m3: np.ndarray) -> np.ndarray:
    """
    Calculates Moliere radius in air (meters), ported from getMoliereRadius in nkg.c.
    Note: legacy getDensityByAltitude returns g/cm^3 and multiplies by 100;
    in SI units (kg/m^3), density_g_cm3 * 100 == density_kg_m3 * 0.1.
    """
    c = 9.578947368  # (20 / 76) * 36.4 g/cm^2
    return c / (np.asarray(density_kg_m3) * 0.1)


def get_nkg_dr(
    ages: np.ndarray,
    molrads: np.ndarray,
    rp_uv: np.ndarray,
    n_pln: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    """
    Vectorized NKG lateral displacement from the shower axis, ported from getNKGdr in nkg.c.

    Parameters
    ----------
    ages : (N,) array of shower ages for each photon
    molrads : (N,) array of Moliere radii (m) for each photon
    rp_uv : (3,) unit vector along Rp in site coordinates
    n_pln : (3,) unit vector normal to shower-detector plane in site coordinates
    rng : numpy.random.Generator

    Returns
    -------
    disp : (N, 3) array of lateral 3D displacement vectors
    """
    ages = np.asarray(ages, dtype=np.float64)
    molrads = np.asarray(molrads, dtype=np.float64)
    n_rays = len(ages)
    if n_rays == 0:
        return np.zeros((0, 3), dtype=np.float64)

    disp = np.zeros((n_rays, 3), dtype=np.float64)
    valid = np.abs(ages) > 1e-12
    if not np.any(valid):
        return disp

    # Vectorized rejection sampling matching legacy getNKGdr:
    # r1, r2 ~ U(0, 1), rm = 3 * molrad * r1^(1/age), accept while r2 <= (1 + rm/molrad)^(age - 4.5)
    rm = np.zeros(n_rays, dtype=np.float64)
    pending = valid.copy()
    while np.any(pending):
        idx = np.flatnonzero(pending)
        r1 = rng.random(len(idx))
        r2 = rng.random(len(idx))
        age_sub = ages[idx]
        mol_sub = molrads[idx]

        rm_cand = 3.0 * mol_sub * np.power(np.clip(r1, 1e-30, None), 1.0 / age_sub)
        accept = r2 <= np.power(1.0 + rm_cand / np.clip(mol_sub, 1e-12, None), age_sub - 4.5)

        accepted_idx = idx[accept]
        rm[accepted_idx] = rm_cand[accept]
        pending[accepted_idx] = False

    phi = 2.0 * np.pi * rng.random(n_rays)
    x = rm * np.sin(phi)
    y = rm * np.cos(phi)

    disp[valid] = x[valid, np.newaxis] * n_pln[np.newaxis, :] + y[valid, np.newaxis] * rp_uv[np.newaxis, :]
    return disp
