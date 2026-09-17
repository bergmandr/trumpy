import numpy as np
from dataclasses import dataclass

@dataclass
class Track:
    """Intrinsic physical properties of the shower or laser."""
    # Global parameters
    species: int            # e.g., 1 for proton, 5626 for Fe
    log_e: float            # log10(E/eV)
    zenith: float           # zenith angle in radians
    impact_v: np.ndarray    # shape (3,) impact point in CLF coordinates
    track_uv: np.ndarray    # shape (3,) unit vector of track direction

    # Vectorized Segments (Structure of Arrays)
    n_segments: int
    positions: np.ndarray   # shape (n_segments, 3)
    time_gen: np.ndarray    # shape (n_segments,) emission times
    altitude: np.ndarray    # shape (n_segments,) altitude at segment mid-points
    dl_seg: np.ndarray      # shape (n_segments,) segment lengths
    de_dep: np.ndarray      # shape (n_segments,) energy deposit profile

@dataclass
class ObservedTrack:
    """Site-specific projections, atmospheric transmission, and photon fluxes."""
    site_id: int
    rp: float               # Impact parameter distance
    psi: float              # Shower angle in the SD plane
    rp_uv: np.ndarray       # shape (3,) unit vector to point of closest approach
    n_pln: np.ndarray       # shape (3,) shower-detector plane normal

    # 1D Segment Arrays
    q_view: np.ndarray      # shape (n_segments,) viewing angles 

    # 2D Arrays: [n_mirrors, n_segments]
    distance: np.ndarray    # Distance from segment to each mirror

    # 3D Arrays: [n_mirrors, n_segments, n_wavelength_bands]
    # Wavelength-dependent light fluxes arriving at the mirrors
    n_fl: np.ndarray        # Fluorescence 
    n_cv_dir: np.ndarray    # Direct Cherenkov
    n_cv_mie: np.ndarray    # Mie-scattered Cherenkov
    n_cv_ray: np.ndarray    # Rayleigh-scattered Cherenkov    