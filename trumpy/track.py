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
    nseg: int               # Number of air shower segments **
    positions: np.ndarray   # shape (nseg, 3)
    time_gen: np.ndarray    # shape (nseg,) emission times
    altitude: np.ndarray    # shape (nseg,) altitude at segment mid-points
    dlseg: np.ndarray       # shape (nseg,) distance (m) to start of next shower segment **
    dedep: np.ndarray       # shape (nseg,) energy deposit profile **
                            # energy deposition (summed over shower) rate (eV/(g/cm2)) at top of shower segment

    # Copying TrackSegment **
    position: np.ndarray   # shape (nseg,) distance from top of first segment
    age: np.ndarray        # shape (nseg,) age at middle of each segment
    dlmid: np.ndarray      # shape (nseg,) density weighted distance (m) to middle of segment
    dxseg: np.ndarray      # shape (nseg,) total grammage in shower segment
    height: np.ndarray     # shape (nseg,) altitude of top of shower segment (m) 
    nch: np.ndarray        # shape (nseg,) number of charged particles in segment
    molrad: np.ndarray     # shape (nseg,) moliere radius at segment mid-point (m)

    nfl: np.ndarray        # shape (nseg,nwl) number of fluorescence photons produced in segment
    pcv: np.ndarray        # shape (nseg,nwl) Cherenkov photons added to beam in the previous segment
    ncv: np.ndarray        # shape (nseg,nwl) Cherenkov photons in beam at top of segment

@dataclass
class ObservedTrack:
    """Site-specific projections, atmospheric transmission, and photon fluxes."""
    site_id: int
    rp: float               # Impact parameter distance
    psi: float              # Shower angle in the SD plane
    rp_uv: np.ndarray       # shape (3,) unit vector to point of closest approach
    n_pln: np.ndarray       # shape (3,) shower-detector plane normal

    # 1D Segment Arrays
    q_view: np.ndarray      # shape (nseg,) viewing angles 

    # 2D Arrays: [n_mirrors, nseg]
    distance: np.ndarray    # Distance from segment to each mirror

    # 3D Arrays: [n_mirrors, nseg, n_wavelength_bands]
    # Wavelength-dependent light fluxes arriving at the mirrors
    n_fl: np.ndarray        # Fluorescence 
    n_cv_dir: np.ndarray    # Direct Cherenkov
    n_cv_mie: np.ndarray    # Mie-scattered Cherenkov
    n_cv_ray: np.ndarray    # Rayleigh-scattered Cherenkov    