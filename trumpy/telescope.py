from dataclasses import dataclass
import numpy as np
import awkward as ak

@dataclass
class TelescopeGeometry:
    """Vectorized optics and camera layout for a given site."""
    mirror_centers: np.ndarray      # shape (n_mirrors, 3)
    mirror_radii: np.ndarray        # shape (n_mirrors,)
    camera_rotations: np.ndarray    # shape (n_mirrors, 3, 3)
    
    # PMT mapping
    pmt_vectors: np.ndarray         # shape (n_mirrors, 256, 3) pointing directions
    pmt_live_flags: np.ndarray      # shape (n_mirrors, 256) boolean mask

@dataclass
class MirrorTopology:
    """Base class for mirror geometries."""
    reflectance: np.ndarray         # Wavelength-dependent reflection curve

@dataclass
class IdealSphericalMirror(MirrorTopology):
    """A perfect continuous spherical cap with stochastic spot-size blurring."""
    radius_of_curvature: float
    center_of_curvature: np.ndarray # shape (3,)
    spot_size_blur: float           # Angular Gaussian blur applied to reflected rays

@dataclass
class SegmentedMirror(MirrorTopology):
    """A dish made of individual adjustable facets (e.g., for detailed CTA/TA alignment)."""
    # Vectorized arrays for N segments (e.g., 18 petals)
    segment_centers: np.ndarray     # shape (n_segments, 3)
    segment_normals: np.ndarray     # shape (n_segments, 3)
    segment_radii: np.ndarray       # shape (n_segments,)