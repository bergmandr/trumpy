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
    nseg: int               # Number of air shower segments
    positions: np.ndarray   # shape (nseg, 3) 3D coordinates at top of each segment
    time_gen: np.ndarray    # shape (nseg,) emission times
    dlseg: np.ndarray       # shape (nseg,) distance (m) to start of next shower segment
    dedep: np.ndarray       # shape (nseg,) energy deposit profile
    
    # Required physics variables
    age: np.ndarray         # shape (nseg,) age at middle of each segment
    dlmid: np.ndarray       # shape (nseg,) density weighted distance (m) to middle of segment
    dxseg: np.ndarray       # shape (nseg,) total grammage in shower segment
    nch: np.ndarray         # shape (nseg,) number of charged particles in segment
    molrad: np.ndarray      # shape (nseg,) moliere radius at segment mid-point (m)

    # Yield profiles
    nfl: np.ndarray         # shape (nseg, nwl) number of fluorescence photons produced in segment
    pcv: np.ndarray         # shape (nseg, nwl) Cherenkov photons added to beam in the previous segment
    ncv: np.ndarray         # shape (nseg, nwl) Cherenkov photons in beam at top of segment

    # Derived Geometrical Properties
    @property
    def height(self) -> np.ndarray:
        """shape (nseg,) altitude of top of shower segment (m)."""
        # Assuming Z-axis (index 2) represents elevation in CLF coordinates
        return self.positions[:, 2]

    @property
    def altitude(self) -> np.ndarray:
        """shape (nseg,) altitude at segment mid-points (m)."""
        # Step halfway down the segment along the track Z-axis.
        # track_uv[2] handles the zenith angle projection and sign automatically.
        return self.positions[:, 2] + 0.5 * self.dlseg * self.track_uv[2]

    @property
    def position(self) -> np.ndarray:
        """shape (nseg,) 1D slant distance from top of first segment (m)."""
        # Calculate the 3D displacement vector from the first segment top to all segments
        delta = self.positions - self.positions[0]
        # Project the displacement onto the track direction (dot product)
        # Using np.sum with axis=1 performs a broadcasted row-wise dot product
        return np.sum(delta * self.track_uv, axis=1)

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