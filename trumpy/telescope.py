from dataclasses import dataclass, field
from typing import Optional, Union
import numpy as np
import awkward as ak
from scipy.spatial.transform import Rotation

from trumpy.track import Track, ObservedTrack, PETimes
from trumpy.camera import CameraAssembly
from trumpy.electronics import PMTState, ElectronicsFrontEnd
from trumpy.shower import get_nkg_dr
from trumpy.optics import (
    mirror_intersection,
    reflect_and_blur,
    map_to_pmts,
    is_within_aperture,
    is_shadowed_by_box,
    snell_shift,
    pmt_detect_hex,
)

# Physical and FD constants from legacy/inc/constants.h and acpttrack.c
SPEED_OF_LIGHT = 0.299792458        # m / ns
COS15 = np.cos(np.deg2rad(15.0))    # 0.9659258263
COS13 = np.cos(np.deg2rad(13.0))    # 0.9743700648
SATURATION_RATE = 8192.0            # Accepted PE / 100 ns
MAX_SATURATED_TUBES = 16


@dataclass
class TelescopeGeometry:
    """Vectorized optics and camera layout for a given site (analogous to geofd_dst_common & TGeom)."""
    mirror_centers: np.ndarray      # shape (n_mirrors, 3) mirror positions in site coords (local_vmir)
    mirror_radii: np.ndarray        # shape (n_mirrors,) radii of curvature (rcurve3)
    camera_rotations: np.ndarray    # shape (n_mirrors, 3, 3) site-to-camera rotation matrices (site2cam)
    
    # PMT mapping
    pmt_vectors: np.ndarray         # shape (n_mirrors, n_tubes, 3) pointing directions
    pmt_live_flags: np.ndarray      # shape (n_mirrors, n_tubes) boolean mask

    # Optional site-level and mirror-level geometric parameters from geofd_dst_common / TGeom
    mirror_pointing: Optional[np.ndarray] = None    # shape (n_mirrors, 3) mirror optical axis unit vectors in site coords (vmir)
    mirror_diameters: Optional[np.ndarray] = None   # shape (n_mirrors,) physical aperture diameters (m)
    camera_separations: Optional[np.ndarray] = None # shape (n_mirrors,) mirror-to-camera distance (sep3, m)
    spot_size_blur: Optional[np.ndarray] = None     # shape (n_mirrors,) Gaussian angular blur (rad)
    local_vsite: Optional[np.ndarray] = None        # shape (3,) site origin in CLF coordinates (m)
    site2clf: Optional[np.ndarray] = None           # shape (3, 3) rotation matrix from site to CLF coordinates

    def __post_init__(self):
        n_mirrors = len(self.mirror_centers)
        if self.mirror_pointing is None:
            # Default mirror pointing vector is the +z axis of each camera coordinate frame
            # Since camera_rotations is site2cam, the 3rd row of site2cam is the mirror axis in site coords
            self.mirror_pointing = self.camera_rotations[:, 2, :].copy()
        if self.mirror_diameters is None:
            self.mirror_diameters = np.full(n_mirrors, 3.3, dtype=np.float64)
        if self.camera_separations is None:
            self.camera_separations = 0.5 * self.mirror_radii
        if self.spot_size_blur is None:
            self.spot_size_blur = np.full(n_mirrors, 0.005, dtype=np.float64)
        if self.local_vsite is None:
            self.local_vsite = np.zeros(3, dtype=np.float64)
        if self.site2clf is None:
            self.site2clf = np.eye(3, dtype=np.float64)


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
    segment_centers: np.ndarray     # shape (nseg, 3)
    segment_normals: np.ndarray     # shape (nseg, 3)
    segment_radii: np.ndarray       # shape (nseg,)


# =====================================================================
# Vectorized Site Geometry Utilities (ported from legacy/src/raytrace/fdsite.c)
# =====================================================================

def get_sdpn_vector(rp_uv: np.ndarray, track_uv: np.ndarray) -> np.ndarray:
    """
    Computes the Shower-Detector Plane (SDP) normal vector.
    Direct vectorized port of getSDPNVector() in fdsite.c.
    """
    n = np.cross(rp_uv, track_uv)
    norm = np.linalg.norm(n, axis=-1, keepdims=True)
    return np.where(norm > 0, n / np.clip(norm, 1e-15, None), 0.0)


def get_rp_origin(
    local_vsite: np.ndarray,
    local_vmir: np.ndarray,
    vmir: np.ndarray,
    site2clf: np.ndarray,
    n_pln: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, int]:
    """
    Finds the mirror whose pointing direction is closest to perpendicular to the SDP normal,
    and returns the Rp origin in CLF coordinates and mirror positions relative to that origin.
    Direct vectorized port of getRpOrigin() in fdsite.c.
    """
    if len(vmir) > 0:
        cos_angles = np.abs(vmir @ n_pln)
        best_mir = int(np.argmin(cos_angles))
        site_origin = local_vmir[best_mir]
        rot_site2clf = Rotation.from_matrix(site2clf)
        obt_local_vsite = rot_site2clf.apply(site_origin) + local_vsite
    else:
        best_mir = -1
        site_origin = np.zeros(3, dtype=np.float64)
        obt_local_vsite = local_vsite.copy()

    obt_local_vmir = local_vmir - site_origin[np.newaxis, :]
    return obt_local_vsite, obt_local_vmir, best_mir


def get_observed_track_origin(
    local_vsite: np.ndarray,
    local_vmir: np.ndarray,
    vmir: np.ndarray,
    site2clf: np.ndarray,
    track: Track,
) -> tuple[np.ndarray, np.ndarray, int]:
    """
    Computes the preliminary SDP normal wrt the site origin and finds the best mirror Rp origin.
    Direct vectorized port of getObservedTrackOrigin() in fdsite.c.
    """
    uv = track.track_uv
    rpv = track.impact_v - local_vsite
    rcosa = float(np.dot(rpv, uv))
    rpv = rpv - rcosa * uv

    rp_norm = np.linalg.norm(rpv)
    rp_uv = rpv / rp_norm if rp_norm > 1e-12 else np.array([1.0, 0.0, 0.0])
    n_pln_clf = get_sdpn_vector(rp_uv, uv)

    # Rotate n_pln from CLF into site coordinate system (site2clf^T)
    rot_clf2site = Rotation.from_matrix(site2clf).inv()
    n_pln_site = rot_clf2site.apply(n_pln_clf)

    return get_rp_origin(local_vsite, local_vmir, vmir, site2clf, n_pln_site)


def get_rp_vector(
    local_vsite: np.ndarray,
    local_vmir: np.ndarray,
    vmir: np.ndarray,
    site2clf: np.ndarray,
    track: Track,
) -> tuple[float, np.ndarray]:
    """
    Computes the Rp vector (from the detector Rp origin to the point of closest approach
    on the shower track) in site coordinates and its signed magnitude.
    Direct vectorized port of getRpVector() in fdsite.c.
    """
    obt_local_vsite, _, _ = get_observed_track_origin(
        local_vsite, local_vmir, vmir, site2clf, track
    )

    rpv_clf = track.impact_v - obt_local_vsite
    rcosa = float(np.dot(rpv_clf, track.track_uv))
    rpv_clf = rpv_clf - rcosa * track.track_uv

    rot_clf2site = Rotation.from_matrix(site2clf).inv()
    rpv_site = rot_clf2site.apply(rpv_clf)

    rp = float(np.linalg.norm(rpv_site))
    if track.track_uv[2] > 0.0:
        rp = -rp

    return rp, rpv_site


def get_psi_angle(
    local_vsite: np.ndarray,
    local_vmir: np.ndarray,
    vmir: np.ndarray,
    site2clf: np.ndarray,
    track: Track,
) -> float:
    """
    Computes the angle Psi (radians) between the shower track and the horizontal plane
    within the shower-detector plane at the detector site.
    Direct vectorized port of getPsiAngle() in fdsite.c.
    """
    _, rpv_site = get_rp_vector(local_vsite, local_vmir, vmir, site2clf, track)
    rot_clf2site = Rotation.from_matrix(site2clf).inv()
    uv_site = rot_clf2site.apply(track.track_uv)

    if np.abs(uv_site[2]) < 1e-12:
        return 0.0

    # Trace track from Rp to where it crosses the z = 0 plane in site coordinates
    rx = rpv_site[0] - rpv_site[2] * uv_site[0] / uv_site[2]
    ry = rpv_site[1] - rpv_site[2] * uv_site[1] / uv_site[2]
    rr = float(np.hypot(rx, ry))

    if rr < 1e-12:
        return float(np.pi / 2.0)

    cpsi = np.clip((rx * uv_site[0] + ry * uv_site[1]) / rr, -1.0, 1.0)
    psi = float(np.arccos(cpsi))
    if uv_site[2] > 0.0:
        psi = -psi

    return psi


def get_zenith_angle(site2clf: np.ndarray, track: Track) -> float:
    """
    Computes the shower zenith angle as viewed in the site coordinate frame.
    Direct port of getZenithAngle() in fdsite.c.
    """
    rot_clf2site = Rotation.from_matrix(site2clf).inv()
    uv_site = rot_clf2site.apply(track.track_uv)
    return float(np.arccos(np.clip(-uv_site[2], -1.0, 1.0)))


def get_viewing_angles(
    local_vsite: np.ndarray,
    site2clf: np.ndarray,
    track: Track,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Vectorized computation of segment viewing angles and chi angles.
    Direct vectorized port of getViewingAngles() in fdsite.c.
    """
    rot_clf2site = Rotation.from_matrix(site2clf).inv()
    R = rot_clf2site.as_matrix()

    x = track.impact_v - local_vsite
    dl = np.dot(x, R[2]) / np.dot(track.track_uv, R[2])

    ru = dl * track.track_uv - x
    ru /= np.clip(np.linalg.norm(ru), 1e-12, None)

    # Vector pointing from each segment top to site origin: s = -(positions - local_vsite)
    s = local_vsite[np.newaxis, :] - track.positions
    su = s / np.clip(np.linalg.norm(s, axis=1, keepdims=True), 1e-12, None)

    view_ang = np.arccos(np.clip(su @ track.track_uv, -1.0, 1.0))
    chi_ang = np.arccos(np.clip(su @ ru, -1.0, 1.0))
    return view_ang, chi_ang


def apply_e_vs_rp_cut(log_e: float, rp: float) -> bool:
    """
    Evaluates the fast E vs Rp geometric filter to discard non-viable tracks.
    Direct port of applyEvsRpCut() in fdsite.c.
    """
    rp_km = abs(rp) / 1000.0
    if rp_km <= 7.5:
        return True
    return bool(log_e > (16.0 + np.sqrt(rp_km - 7.5) / 2.5))


# =====================================================================
# Detector Site & Ray-Tracing Orchestrator (FDSite)
# =====================================================================

class FDSite:
    """
    Concrete implementation of the Experiment protocol representing a Fluorescence
    Detector site (e.g., Black Rock, Long Ridge, Middle Drum, TALE).
    Orchestrates site-frame projections (fdsite.c), ray-tracing (doRayTrace.c / trace.c),
    accepted photoelectron collection (acpttrack.c), and front-end electronics.
    """

    def __init__(
        self,
        site_id: int,
        geometry: TelescopeGeometry,
        camera: CameraAssembly,
        pmt_state: PMTState,
        electronics: ElectronicsFrontEnd,
        wavelength_bands: np.ndarray,
        reduct_factor: Optional[np.ndarray] = None,
        enable_evsrp_cut: bool = True,
        enable_fov_cut: bool = True,
        enable_nkg_spread: bool = True,
        enable_camera_shadow: bool = False,
    ):
        self.site_id = site_id
        self.geometry = geometry
        self.camera = camera
        self.pmt_state = pmt_state
        self.electronics = electronics
        self.wavelength_bands = np.asarray(wavelength_bands, dtype=np.float64)

        self.enable_evsrp_cut = enable_evsrp_cut
        self.enable_fov_cut = enable_fov_cut
        self.enable_nkg_spread = enable_nkg_spread
        self.enable_camera_shadow = enable_camera_shadow

        n_mirrors = self.n_mirrors
        n_wl = len(self.wavelength_bands)
        if reduct_factor is None:
            # Compute default optical efficiency (mirror reflectivity * optical layers * PMT QE)
            # Matches initRayTrace.c: mir_reflect (0.92) * camcover (0.9106) * filter (0.8915) * pmtQE (0.28)
            layer_trans = np.ones(n_wl, dtype=np.float64)
            for layer in self.camera.optical_layers:
                layer_trans = layer_trans * np.asarray(layer.transmittance)
            default_eff = 0.92 * 0.28 * layer_trans
            self.reduct_factor = np.broadcast_to(default_eff[np.newaxis, :], (n_mirrors, n_wl)).copy()
        else:
            self.reduct_factor = np.asarray(reduct_factor, dtype=np.float64)

        # Precompute mirror areas, chord target plane heights, and mean mirror-to-camera travel times (initRayTrace.c)
        radii = self.geometry.mirror_radii
        diams = self.geometry.mirror_diameters
        seps = self.geometry.camera_separations

        self.mirror_areas = np.pi * (diams ** 2) / 4.0
        self.target_h = radii - np.sqrt(np.clip(radii ** 2 - (diams ** 2) / 4.0, 0.0, None))
        self.mean_time_ns = 0.5 * (
            seps / SPEED_OF_LIGHT + np.sqrt(seps ** 2 + (diams ** 2) / 4.0) / SPEED_OF_LIGHT
        )

    @property
    def n_mirrors(self) -> int:
        return len(self.geometry.mirror_centers)

    @property
    def n_tubes(self) -> int:
        return self.geometry.pmt_vectors.shape[1]

    @property
    def mirror_centers(self) -> np.ndarray:
        """
        Returns mirror centers in global CLF coordinates: shape (n_mirrors, 3).
        Satisfies the Experiment protocol for Simulation.run().
        """
        rot_site2clf = Rotation.from_matrix(self.geometry.site2clf)
        return rot_site2clf.apply(self.geometry.mirror_centers) + self.geometry.local_vsite[np.newaxis, :]

    def passes_fast_cuts(self, track: Track) -> bool:
        """Evaluates E vs Rp geometric filter (applyEvsRpCut in fdsite.c)."""
        if not self.enable_evsrp_cut:
            return True
        rp, _ = get_rp_vector(
            self.geometry.local_vsite,
            self.geometry.mirror_centers,
            self.geometry.mirror_pointing,
            self.geometry.site2clf,
            track,
        )
        return apply_e_vs_rp_cut(track.log_e, rp)

    def get_observed_track(
        self,
        track: Track,
        attenuated_fluorescence: Optional[np.ndarray] = None,
    ) -> ObservedTrack:
        """
        Vectorized computation of site-frame geometry, solid angles, and arriving photon fluxes.
        Direct port of getObservedTrackGeom() and fillObservedTrackLight() in fdsite.c.

        Parameters
        ----------
        track : Track
        attenuated_fluorescence : optional (nseg, n_mirrors, n_wavelengths) array of
            fluorescence photons after atmospheric transmission (before solid angle geoFactor).

        Returns
        -------
        ObservedTrack
        """
        obt_local_vsite, obt_local_vmir, _ = get_observed_track_origin(
            self.geometry.local_vsite,
            self.geometry.mirror_centers,
            self.geometry.mirror_pointing,
            self.geometry.site2clf,
            track,
        )

        rot_clf2site = Rotation.from_matrix(self.geometry.site2clf).inv()
        uv_site = rot_clf2site.apply(track.track_uv)
        uv_norm = np.linalg.norm(uv_site)
        if uv_norm > 0:
            uv_site = uv_site / uv_norm

        # Impact point and segment positions in site coordinates wrt Rp origin
        rimp = rot_clf2site.apply(track.impact_v - obt_local_vsite)
        sv = rot_clf2site.apply(track.positions - obt_local_vsite)  # (nseg, 3) top of segments
        smv = sv + track.dlmid[:, np.newaxis] * uv_site[np.newaxis, :]  # (nseg, 3) middle of segments
        snv = sv + track.dlseg[:, np.newaxis] * uv_site[np.newaxis, :]  # (nseg, 3) bottom of segments

        # Rp, Psi, rp_uv, and SDP normal in site coordinates
        rp, rpv_site = get_rp_vector(
            self.geometry.local_vsite,
            self.geometry.mirror_centers,
            self.geometry.mirror_pointing,
            self.geometry.site2clf,
            track,
        )
        psi = get_psi_angle(
            self.geometry.local_vsite,
            self.geometry.mirror_centers,
            self.geometry.mirror_pointing,
            self.geometry.site2clf,
            track,
        )
        rp_norm = np.linalg.norm(rpv_site)
        rp_uv = rpv_site / rp_norm if rp_norm > 1e-12 else np.array([1.0, 0.0, 0.0])
        n_pln = get_sdpn_vector(rp_uv, uv_site)

        # Segment midpoint generation times (ns)
        t_gen = track.time_gen + track.dlmid / SPEED_OF_LIGHT

        # Viewing angle at middle of each segment wrt Rp origin
        ms_norm = np.clip(np.linalg.norm(smv, axis=1, keepdims=True), 1e-12, None)
        ms_uv = -smv / ms_norm
        q_view = np.arccos(np.clip(ms_uv @ uv_site, -1.0, 1.0))

        # Vectorized segment-to-mirror geometry: shape (n_mirrors, nseg, 3)
        rv = smv[np.newaxis, :, :] - obt_local_vmir[:, np.newaxis, :]
        dist = np.linalg.norm(rv, axis=2)  # (n_mirrors, nseg)
        rv_unit = rv / np.clip(dist[:, :, np.newaxis], 1e-12, None)

        # Cosine of angle between segment-to-mirror vector and mirror pointing vector
        cos_the = np.einsum("mni,mi->mn", rv_unit, self.geometry.mirror_pointing)

        # Angular extent dtheta of each segment as seen by each mirror
        t0 = sv[np.newaxis, :, :] - obt_local_vmir[:, np.newaxis, :]
        t1 = snv[np.newaxis, :, :] - obt_local_vmir[:, np.newaxis, :]
        t0_norm = np.clip(np.linalg.norm(t0, axis=2), 1e-12, None)
        t1_norm = np.clip(np.linalg.norm(t1, axis=2), 1e-12, None)
        cos_dtheta = np.clip(np.sum(t0 * t1, axis=2) / (t0_norm * t1_norm), -1.0, 1.0)
        dtheta = np.arccos(cos_dtheta)

        in_fov = cos_the > COS15 if self.enable_fov_cut else (cos_the > 0.0)

        # Solid angle geometric factor: cq * Area / (4 * pi * r^2)
        safe_dist_sq = np.clip(dist ** 2, 1.0, None)
        geo_factor = np.where(
            in_fov,
            np.clip(cos_the, 0.0, None) * self.mirror_areas[:, np.newaxis] / (4.0 * np.pi * safe_dist_sq),
            0.0,
        )

        # Arriving light fluxes: shape (n_mirrors, nseg, n_wavelengths)
        n_wl = len(self.wavelength_bands)
        if attenuated_fluorescence is not None:
            # Transpose (nseg, n_mirrors, n_wl) -> (n_mirrors, nseg, n_wl) and multiply by geo_factor
            att_fl_mnw = np.swapaxes(attenuated_fluorescence, 0, 1)
            n_fl = att_fl_mnw * geo_factor[:, :, np.newaxis]
        else:
            # Use track.nfl if populated
            n_fl = track.nfl[np.newaxis, :, :] * geo_factor[:, :, np.newaxis]

        zeros_flux = np.zeros((self.n_mirrors, track.nseg, n_wl), dtype=np.float64)
        return ObservedTrack(
            site_id=self.site_id,
            rp=rp,
            psi=psi,
            rp_uv=rp_uv,
            n_pln=n_pln,
            q_view=q_view,
            distance=dist,
            n_fl=n_fl,
            n_cv_dir=zeros_flux.copy(),
            n_cv_mie=zeros_flux.copy(),
            n_cv_ray=zeros_flux.copy(),
            uv=uv_site,
            rimp=rimp,
            local_vsite=obt_local_vsite,
            local_vmir=obt_local_vmir,
            t_gen=t_gen,
            cos_the=cos_the,
            geo_factor=geo_factor,
            dtheta=dtheta,
            in_fov=in_fov,
        )

    def do_ray_trace(
        self,
        cam_id: int,
        origins_site: np.ndarray,
        rp_uv: np.ndarray,
        n_pln: np.ndarray,
        ages: np.ndarray,
        molrads: np.ndarray,
        rng: np.random.Generator,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Vectorized ray-tracing for a batch of N photons aimed at camera `cam_id`.
        Direct vectorized analogue of doRayTrace() (doRayTrace.c) and trace() (trace.c).

        Parameters
        ----------
        cam_id : int, camera/mirror index
        origins_site : (N, 3) photon emission points relative to mirror center in site coords
        rp_uv, n_pln : (3,) shower-detector plane basis vectors in site coords
        ages, molrads : (N,) shower age and Moliere radius for each photon
        rng : numpy.random.Generator

        Returns
        -------
        hit_pmt_ids : (N_valid,) array of hit PMT indices
        travel_times_ns : (N_valid,) array of photon time-of-flight from emission to PMT (ns)
        valid_ray_indices : (N_valid,) indices into the input batch (0..N-1) that hit a PMT
        """
        n_rays = len(origins_site)
        if n_rays == 0:
            return (
                np.empty(0, dtype=np.int64),
                np.empty(0, dtype=np.float64),
                np.empty(0, dtype=np.int64),
            )

        rcurve = float(self.geometry.mirror_radii[cam_id])
        diameter = float(self.geometry.mirror_diameters[cam_id])
        aperture_radius = diameter / 2.0
        vmir_axis = self.geometry.mirror_pointing[cam_id]
        site2cam = self.geometry.camera_rotations[cam_id]
        blur_rad = float(self.geometry.spot_size_blur[cam_id])
        h_chord = float(self.target_h[cam_id])
        sep = float(self.geometry.camera_separations[cam_id])
        mean_time = float(self.mean_time_ns[cam_id])

        # 1. Vector from mirror center of curvature to source in site coordinates (doRayTrace.c:44)
        v_coc = origins_site - rcurve * vmir_axis[np.newaxis, :]

        # 2. Lateral NKG shower displacement (doRayTrace.c:47)
        if self.enable_nkg_spread:
            disp = get_nkg_dr(ages, molrads, rp_uv, n_pln, rng)
            v_coc = v_coc + disp

        # 3. Rotate source vectors into camera coordinates (doRayTrace.c:53)
        # In camera coordinates: center of curvature is (0, 0, 0), mirror vertex is (0, 0, -rcurve)
        rot_site2cam = Rotation.from_matrix(site2cam)
        vsite_cam = rot_site2cam.apply(v_coc)

        # 4. Sample random locations on mirror target chord plane (trace.c:63-68)
        r_disk = aperture_radius * np.sqrt(rng.random(n_rays))
        phi_disk = 2.0 * np.pi * rng.random(n_rays)
        w = np.column_stack([
            r_disk * np.cos(phi_disk),
            r_disk * np.sin(phi_disk),
            np.full(n_rays, -rcurve + h_chord),
        ])

        # Direction unit vector `directions` from source (`vsite_cam`) toward target point (`w`)
        vec_to_mirror = w - vsite_cam
        ray_dirs = vec_to_mirror / np.clip(np.linalg.norm(vec_to_mirror, axis=1, keepdims=True), 1e-12, None)

        # 5. Intersect rays with spherical mirror centered at (0, 0, 0) of radius `rcurve`
        mirror_coc = np.zeros(3, dtype=np.float64)
        mirror_vertex = np.array([0.0, 0.0, -rcurve])
        mirror_axis = np.array([0.0, 0.0, 1.0])

        # Note: mirror_intersection takes origins and finds closer intersection; since vsite_cam is
        # outside the sphere (distance >> rcurve) and shooting toward the back cap (z ~ -rcurve),
        # starting from the chord point `w` along `ray_dirs` directly hits the mirror cap at small t > 0!
        # Let's compute exact sphere intersection of P(s) = w + s * ray_dirs with |P|^2 = rcurve^2:
        b_quad = 2.0 * np.sum(w * ray_dirs, axis=1)
        c_quad = np.sum(w * w, axis=1) - rcurve ** 2
        disc = b_quad ** 2 - 4.0 * c_quad
        sphere_hit = disc >= 0.0
        if not np.any(sphere_hit):
            return (
                np.empty(0, dtype=np.int64),
                np.empty(0, dtype=np.float64),
                np.empty(0, dtype=np.int64),
            )

        active_idx = np.flatnonzero(sphere_hit)
        s_hit = (-b_quad[active_idx] + np.sqrt(disc[active_idx])) / 2.0
        u_hit = w[active_idx] + s_hit[:, np.newaxis] * ray_dirs[active_idx]
        dirs_active = ray_dirs[active_idx]
        vsite_active = vsite_cam[active_idx]

        # 6. Check physical mirror aperture
        in_ap = is_within_aperture(u_hit, mirror_vertex, mirror_axis, aperture_radius)
        if not np.any(in_ap):
            return (
                np.empty(0, dtype=np.int64),
                np.empty(0, dtype=np.float64),
                np.empty(0, dtype=np.int64),
            )

        active_idx = active_idx[in_ap]
        u_hit = u_hit[in_ap]
        dirs_active = dirs_active[in_ap]
        vsite_active = vsite_active[in_ap]

        # Distance from source to mirror intersection and total travel time (trace.c:103-104)
        dist_to_mir = np.linalg.norm(vsite_active - u_hit, axis=1)
        travel_time = dist_to_mir / SPEED_OF_LIGHT + mean_time

        # 7. Optional camera box shadow check (optics.is_shadowed_by_box)
        if self.enable_camera_shadow:
            box_center = np.array([0.0, 0.0, -rcurve + sep + self.camera.depth / 2.0])
            box_half_extents = np.array([
                self.camera.width / 2.0,
                self.camera.height / 2.0,
                self.camera.depth / 2.0,
            ])
            box_axes = np.eye(3, dtype=np.float64)
            shadowed = is_shadowed_by_box(
                vsite_active, dirs_active, dist_to_mir, box_center, box_half_extents, box_axes
            )
            not_shadowed = ~shadowed
            if not np.any(not_shadowed):
                return (
                    np.empty(0, dtype=np.int64),
                    np.empty(0, dtype=np.float64),
                    np.empty(0, dtype=np.int64),
                )
            active_idx = active_idx[not_shadowed]
            u_hit = u_hit[not_shadowed]
            dirs_active = dirs_active[not_shadowed]
            travel_time = travel_time[not_shadowed]

        # 8. Reflect off spherical mirror and apply spot-size Gaussian blur (optics.reflect_and_blur)
        reflected_dirs = reflect_and_blur(
            dirs_active, u_hit, mirror_coc, rcurve, blur_rad, rng
        )

        # 9. Propagate reflected rays to camera focal plane and apply Snell shifts (trace.c:201-235)
        cos_gzen = np.clip(reflected_dirs[:, 2], 1e-8, None)
        total_layer_thickness = sum(layer.thickness for layer in self.camera.optical_layers)
        sep_to_window = sep - total_layer_thickness

        path_to_window = (-u_hit[:, 2] - (rcurve - sep_to_window)) / cos_gzen
        f_window = u_hit + path_to_window[:, np.newaxis] * reflected_dirs

        # Check camera box aperture boundary (trace.c:210-215)
        in_cam_box = (
            (np.abs(f_window[:, 0]) <= (self.camera.width / 2.0))
            & (np.abs(f_window[:, 1]) <= (self.camera.height / 2.0))
            & (reflected_dirs[:, 2] > 0.0)
        )
        if not np.any(in_cam_box):
            return (
                np.empty(0, dtype=np.int64),
                np.empty(0, dtype=np.float64),
                np.empty(0, dtype=np.int64),
            )

        active_idx = active_idx[in_cam_box]
        u_hit = u_hit[in_cam_box]
        reflected_dirs = reflected_dirs[in_cam_box]
        f_window = f_window[in_cam_box]
        cos_gzen = cos_gzen[in_cam_box]
        travel_time = travel_time[in_cam_box]

        # Apply Snell shift for each OpticalLayer (Paraglas, BG3, etc.)
        total_shift = np.zeros_like(f_window)
        for layer in self.camera.optical_layers:
            total_shift += snell_shift(
                reflected_dirs, 1.000293, float(layer.refractive_index), float(layer.thickness)
            )

        f_pmt = (
            f_window
            + (total_layer_thickness / cos_gzen)[:, np.newaxis] * reflected_dirs
            + total_shift
        )

        # 10. Map rays to PMTs (hexagonal focal plane coords if available, else angular map_to_pmts)
        if self.camera.pmt_positions is not None:
            hit_pmts, valid_pmt_mask = pmt_detect_hex(
                f_pmt[:, 0],
                f_pmt[:, 1],
                self.camera.pmt_positions[:, 0],
                self.camera.pmt_positions[:, 1],
                self.camera.pmt_flat2flat,
                self.camera.pmt_point2point,
            )
            hit_pmts = hit_pmts[valid_pmt_mask]
        else:
            hit_pmts, valid_pmt_mask = map_to_pmts(
                reflected_dirs,
                self.geometry.pmt_vectors[cam_id],
                self.camera.fov_limit_rad,
            )

        return (
            hit_pmts,
            travel_time[valid_pmt_mask],
            active_idx[valid_pmt_mask],
        )

    def get_accepted_track(
        self,
        track: Track,
        observed_track: ObservedTrack,
        rng: np.random.Generator,
    ) -> PETimes:
        """
        Vectorized orchestration of ray-tracing across all site cameras and shower segments.
        Direct analogue of getAcceptedTrack() in legacy/src/raytrace/acpttrack.c.
        """
        rot_clf2site = Rotation.from_matrix(self.geometry.site2clf).inv()
        rtop_site = rot_clf2site.apply(track.positions[0] - observed_track.local_vsite)
        tkuv_site = observed_track.uv

        # Total arriving photon flux per (mirror, segment, wavelength)
        total_flux_mnw = (
            observed_track.n_fl
            + observed_track.n_cv_dir
            + observed_track.n_cv_mie
            + observed_track.n_cv_ray
        )
        # Apply wavelength reduction factor (mirror reflectivity * filter trans * QE)
        reduced_flux_mn = np.sum(
            total_flux_mnw * self.reduct_factor[:, np.newaxis, :], axis=2
        )  # shape (n_mirrors, nseg)

        # Segment midpoint vectors wrt site origin for 13-degree camera FOV check (acpttrack.c:135-137)
        seg_pos_1d = track.position
        rmid_site = (
            rtop_site[np.newaxis, :]
            + (seg_pos_1d + track.dlmid)[:, np.newaxis] * tkuv_site[np.newaxis, :]
        )
        rmid_unit = rmid_site / np.clip(np.linalg.norm(rmid_site, axis=1, keepdims=True), 1e-12, None)

        all_cam_t = []
        all_cam_n = []
        participating_mirs = []
        t1_global = np.inf
        t2_global = -np.inf

        for icam in range(self.n_mirrors):
            # Skip segments outside 13-degree cone of mirror axis if FOV cut enabled
            if self.enable_fov_cut:
                dp = rmid_unit @ self.geometry.mirror_pointing[icam]
                seg_mask = dp >= COS13
            else:
                seg_mask = np.ones(track.nseg, dtype=bool)

            cam_flux = np.where(seg_mask, np.clip(reduced_flux_mn[icam], 0.0, None), 0.0)
            nphot_per_seg = rng.poisson(cam_flux)
            total_cam_photons = int(np.sum(nphot_per_seg))

            if total_cam_photons == 0:
                all_cam_t.append([[] for _ in range(self.n_tubes)])
                all_cam_n.append([[] for _ in range(self.n_tubes)])
                continue

            # Expand segment indices for every generated photon in this camera
            seg_indices = np.repeat(np.arange(track.nseg), nphot_per_seg)

            # Sample longitudinal offset along each segment (acpttrack.c:175-178)
            u_rand = rng.random(total_cam_photons)
            lgen = u_rand * track.dlseg[seg_indices]

            rbeg_seg = (
                rtop_site[np.newaxis, :]
                + seg_pos_1d[seg_indices, np.newaxis] * tkuv_site[np.newaxis, :]
            )
            rgen_site = rbeg_seg + lgen[:, np.newaxis] * tkuv_site[np.newaxis, :]
            origins_wrt_cam = rgen_site - observed_track.local_vmir[icam][np.newaxis, :]

            tgen_phot = (
                observed_track.t_gen[seg_indices]
                + (lgen - track.dlmid[seg_indices]) / SPEED_OF_LIGHT
            )

            hit_pmts, tof_ns, valid_idx = self.do_ray_trace(
                cam_id=icam,
                origins_site=origins_wrt_cam,
                rp_uv=observed_track.rp_uv,
                n_pln=observed_track.n_pln,
                ages=track.age[seg_indices],
                molrads=track.molrad[seg_indices],
                rng=rng,
            )

            if len(hit_pmts) == 0:
                all_cam_t.append([[] for _ in range(self.n_tubes)])
                all_cam_n.append([[] for _ in range(self.n_tubes)])
                continue

            arrival_times = tgen_phot[valid_idx] + tof_ns

            # Filter by PMT live flags
            live_flags = self.geometry.pmt_live_flags[icam]
            live_hit_mask = live_flags[hit_pmts]
            hit_pmts = hit_pmts[live_hit_mask]
            arrival_times = arrival_times[live_hit_mask]

            if len(hit_pmts) == 0:
                all_cam_t.append([[] for _ in range(self.n_tubes)])
                all_cam_n.append([[] for _ in range(self.n_tubes)])
                continue

            participating_mirs.append(icam)

            # Group and sort arrival times by tube ID
            order = np.lexsort((arrival_times, hit_pmts))
            sorted_pmts = hit_pmts[order]
            sorted_times = arrival_times[order]

            unique_tubes, split_idx = np.unique(sorted_pmts, return_index=True)
            grouped_times = np.split(sorted_times, split_idx[1:])

            cam_tube_t = [[] for _ in range(self.n_tubes)]
            cam_tube_n = [[] for _ in range(self.n_tubes)]

            for tube_id, t_arr in zip(unique_tubes, grouped_times):
                # Check saturation rate (acpttrack.c:218-239)
                if len(t_arr) > 4096 and (t_arr[-1] > t_arr[0]):
                    pe_rate = 100.0 * len(t_arr) / (t_arr[-1] - t_arr[0])
                    if pe_rate >= SATURATION_RATE:
                        continue

                cam_tube_t[int(tube_id)] = t_arr.tolist()
                cam_tube_n[int(tube_id)] = np.ones(len(t_arr), dtype=np.float64).tolist()
                t1_global = min(t1_global, float(t_arr[0]))
                t2_global = max(t2_global, float(t_arr[-1]))

            all_cam_t.append(cam_tube_t)
            all_cam_n.append(cam_tube_n)

        if not participating_mirs:
            t1_global = 0.0
            t2_global = 0.0

        return PETimes(
            t1=float(t1_global),
            t2=float(t2_global),
            nmir=len(participating_mirs),
            mir=np.array(participating_mirs, dtype=np.int32),
            t=ak.Array(all_cam_t),
            n=ak.Array(all_cam_n),
        )

    def trace_photons(
        self,
        track: Track,
        photons_at_mirrors: np.ndarray,
        rng: np.random.Generator,
    ) -> PETimes:
        """
        Satisfies Experiment.trace_photons() protocol for Simulation.run().
        Takes atmospheric-attenuated photon production (nseg, n_mirrors, n_wavelengths),
        computes ObservedTrack, runs ray-tracing, and returns PETimes.
        """
        observed_track = self.get_observed_track(
            track, attenuated_fluorescence=photons_at_mirrors
        )
        return self.get_accepted_track(track, observed_track, rng)

    def process_electronics(
        self,
        pe_times: Union[PETimes, ak.Array],
        rng: np.random.Generator,
    ) -> tuple[ak.Array, bool]:
        """
        Satisfies Experiment.process_electronics() protocol for Simulation.run().
        Digitizes PE arrival times across cameras and evaluates trigger status.
        """
        if isinstance(pe_times, PETimes):
            if self.n_mirrors == 1:
                readout = self.electronics.digitize(
                    pe_times.for_camera(0), self.pmt_state, rng
                )
                return readout, (pe_times.nmir > 0)

            cam_readouts = []
            for icam in range(self.n_mirrors):
                cam_pe = pe_times.for_camera(icam)
                cam_ro = self.electronics.digitize(cam_pe, self.pmt_state, rng)
                cam_readouts.append(cam_ro)
            return ak.Array(cam_readouts), (pe_times.nmir > 0)

        # Fallback if a 2D (n_tubes, var) jagged ak.Array is passed directly
        readout = self.electronics.digitize(pe_times, self.pmt_state, rng)
        return readout, True
