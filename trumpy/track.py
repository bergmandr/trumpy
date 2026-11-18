import numpy as np
import awkward as ak
from dataclasses import dataclass, field
from typing import Optional

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

    @property
    def length(self) -> float:
        """Total length of air shower from top of first segment to impact point (m)."""
        return float(np.linalg.norm(self.impact_v - self.positions[0]))


@dataclass
class ObservedTrack:
    """Site-specific projections, atmospheric transmission, and photon fluxes (analogous to fdsite.h)."""
    site_id: int
    rp: float               # Impact parameter distance (m)
    psi: float              # Shower angle in the SD plane (radians)
    rp_uv: np.ndarray       # shape (3,) unit vector to point of closest approach in site coords
    n_pln: np.ndarray       # shape (3,) shower-detector plane normal in site coords

    # 1D Segment Arrays
    q_view: np.ndarray      # shape (nseg,) viewing angles (radians)

    # 2D Arrays: [n_mirrors, nseg]
    distance: np.ndarray    # Distance from segment midpoint to each mirror (m)

    # 3D Arrays: [n_mirrors, nseg, n_wavelength_bands]
    # Wavelength-dependent light fluxes arriving at the mirrors
    n_fl: np.ndarray        # Fluorescence
    n_cv_dir: np.ndarray    # Direct Cherenkov
    n_cv_mie: np.ndarray    # Mie-scattered Cherenkov
    n_cv_ray: np.ndarray    # Rayleigh-scattered Cherenkov

    # Optional site-coordinate geometry arrays mirroring fdsite.h ObservedTrack / ObservedTrackSegment
    uv: Optional[np.ndarray] = None           # shape (3,) shower direction unit vector in site coords
    rimp: Optional[np.ndarray] = None         # shape (3,) shower impact position in site coords
    local_vsite: Optional[np.ndarray] = None  # shape (3,) CLF to Rp origin in CLF coords
    local_vmir: Optional[np.ndarray] = None   # shape (n_mirrors, 3) mirror locations wrt Rp origin in site coords
    t_gen: Optional[np.ndarray] = None        # shape (nseg,) generation time at segment midpoint (ns)
    cos_the: Optional[np.ndarray] = None      # shape (n_mirrors, nseg) cosine of angle to mirror axis
    geo_factor: Optional[np.ndarray] = None   # shape (n_mirrors, nseg) solid-angle factor cq * A / (4*pi*r^2)
    dtheta: Optional[np.ndarray] = None       # shape (n_mirrors, nseg) angular extent of segment (radians)
    in_fov: Optional[np.ndarray] = None       # shape (n_mirrors, nseg) boolean mask (cos_the > cos(15 deg))


@dataclass
class PETimes:
    """
    Arrival times of photoelectrons indexed by camera and tube.
    Direct Python/Awkward analogue of PETimes in legacy/inc/acpttrack.h.
    """
    t1: float                       # Time of first PE across site (ns)
    t2: float                       # Time of last PE across site (ns)
    nmir: int                       # Number of participating cameras with >= 1 PE
    mir: np.ndarray                 # shape (nmir,) participating camera IDs
    t: ak.Array                     # Jagged array of PE arrival times: shape (n_cameras, n_tubes, var)
    n: Optional[ak.Array] = None    # Optional jagged array of PE weights/counts matching `t`

    def for_camera(self, cam_id: int) -> ak.Array:
        """Returns the (n_tubes, var) jagged array of PE arrival times for a single camera."""
        return self.t[cam_id]
