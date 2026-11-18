from dataclasses import dataclass, field
from typing import Optional
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

    # Optional physical camera dimensions and PMT plane coordinates (from geofd / TGeom)
    width: float = 1.152            # Camera box width (m)
    height: float = 1.000           # Camera box height (m)
    depth: float = 0.400            # Camera box depth (m)
    fov_limit_rad: float = np.deg2rad(0.55)  # Half-angle acceptance per tube
    pmt_positions: Optional[np.ndarray] = None  # shape (n_tubes, 2) (xtube, ytube) in meters
    pmt_flat2flat: float = 0.060    # Hexagonal PMT flat-to-flat size (m)
    pmt_point2point: float = 0.0693 # Hexagonal PMT point-to-point size (m)

    @classmethod
    def default_ta_camera(cls, n_wavelengths: int = 1) -> "CameraAssembly":
        """
        Builds a standard 256-PMT (16x16 hexagonal grid) Telescope Array camera assembly
        with Paraglas window and BG3 filter layers matching initRayTrace.c.
        """
        n_rows, n_cols = 16, 16
        flat2flat = 0.061
        point2point = flat2flat / np.cos(np.deg2rad(30.0))
        dx = flat2flat
        dy = flat2flat * np.cos(np.deg2rad(30.0))

        rows, cols = np.meshgrid(np.arange(n_rows), np.arange(n_cols), indexing="ij")
        x_tube = (cols - (n_cols - 1) / 2.0 + 0.25 * ((rows % 2) * 2 - 1)) * dx
        y_tube = (rows - (n_rows - 1) / 2.0) * dy
        pmt_positions = np.column_stack([x_tube.ravel(), y_tube.ravel()])

        # Approximate pointing vectors from focal length ~ 3.0 m
        focal_length = 3.0
        dirs = np.column_stack([
            -pmt_positions[:, 0],
            -pmt_positions[:, 1],
            np.full(len(pmt_positions), focal_length),
        ])
        dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)

        paraglas = OpticalLayer(
            name="Paraglas",
            thickness=0.005,
            refractive_index=1.50,
            transmittance=np.full(n_wavelengths, 0.910642),
            distance_from_pmt=0.024,
        )
        bg3 = OpticalLayer(
            name="BG3",
            thickness=0.004,
            refractive_index=1.53,
            transmittance=np.full(n_wavelengths, 0.891514),
            distance_from_pmt=0.0,
        )

        return cls(
            pmt_vectors=dirs,
            pmt_live_flags=np.ones(len(pmt_positions), dtype=bool),
            optical_layers=[paraglas, bg3],
            width=1.152,
            height=1.000,
            depth=0.400,
            fov_limit_rad=np.deg2rad(0.55),
            pmt_positions=pmt_positions,
            pmt_flat2flat=flat2flat,
            pmt_point2point=point2point,
        )

