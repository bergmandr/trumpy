from dataclasses import dataclass, field
import numpy as np

@dataclass
class OpticalLayer:
    """Represents a pane of glass, acrylic, or a UV filter in the optical path."""
    name: str
    thickness: float                # Thickness in meters
    refractive_index: float         # Can be a constant or a wavelength-dependent array
    transmittance: np.ndarray       # Wavelength-dependent transmission curve
    distance_from_pmt: float        # Position in the camera box

@dataclass
class CameraAssembly:
    """The focal plane and its protective/filtering layers."""
    pmt_vectors: np.ndarray         # shape (n_tubes, 3) pointing directions
    pmt_live_flags: np.ndarray      # shape (n_tubes,) boolean mask
    
    # List of optical layers (e.g., [ParaglasLayer, BG3FilterLayer])
    # HiRes might only have [UVFilterLayer]
    optical_layers: list[OpticalLayer] = field(default_factory=list)
