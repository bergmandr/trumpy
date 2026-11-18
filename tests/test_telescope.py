import numpy as np
import awkward as ak
import pytest

from trumpy.track import Track, ObservedTrack, PETimes
from trumpy.shower import GaisserHillasProfile, get_moliere_radius, get_nkg_dr
from trumpy.camera import CameraAssembly
from trumpy.optics import snell_shift, pmt_detect_hex
from trumpy.electronics import PMTState, TAFADCFrontEnd, generate_ta_impulse_response
from trumpy.atmosphere import HydrostaticProfile, AerosolModel, Atmosphere
from trumpy.fluorescence import KakimotoYield
from trumpy.simulation import Simulation
from trumpy.telescope import (
    TelescopeGeometry,
    FDSite,
    get_sdpn_vector,
    get_rp_origin,
    get_rp_vector,
    get_psi_angle,
    get_zenith_angle,
    get_viewing_angles,
    apply_e_vs_rp_cut,
)


@pytest.fixture
def vertical_track():
    """Creates a vertical shower track at the CLF origin (0, 0, 1400) pointing straight down."""
    nseg = 40
    impact = np.array([0.0, 0.0, 1400.0])
    uv = np.array([0.0, 0.0, -1.0])
    distances = np.linspace(12000.0, 0.0, nseg)
    positions = impact + distances[:, np.newaxis] * (-uv)
    gh = GaisserHillasProfile(x0=-75.8, xmax=773.2, nmax=6.692e9, lambda_inv=59.9)
    slant_depths = np.linspace(100.0, 900.0, nseg)

    return Track(
        species=14,
        log_e=19.0,
        zenith=0.0,
        impact_v=impact,
        track_uv=uv,
        nseg=nseg,
        positions=positions,
        time_gen=np.linspace(0.0, 40000.0, nseg),
        dlseg=np.full(nseg, 300.0),
        dedep=gh.evaluate_dedep(slant_depths),
        age=np.full(nseg, 1.0),
        dlmid=np.full(nseg, 150.0),
        dxseg=np.full(nseg, 20.0),
        nch=gh.evaluate_particles(slant_depths),
        molrad=np.full(nseg, 90.0),
        nfl=np.zeros((nseg, 1)),
        pcv=np.zeros((nseg, 1)),
        ncv=np.zeros((nseg, 1)),
    )


@pytest.fixture
def single_mirror_site():
    """Creates a 1-mirror FDSite located at x = 10 km looking toward the CLF origin (-x direction)."""
    camera = CameraAssembly.default_ta_camera(n_wavelengths=1)
    n_tubes = len(camera.pmt_vectors)

    # Site at (10000, 0, 1400) in CLF coords; mirror pointing along (-1, 0, 0.3) toward shower
    axis_z = np.array([-1.0, 0.0, 0.35])
    axis_z /= np.linalg.norm(axis_z)
    axis_x = np.array([0.0, 1.0, 0.0])
    axis_y = np.cross(axis_z, axis_x)
    site2cam = np.vstack([axis_x, axis_y, axis_z])[np.newaxis, :, :]

    geometry = TelescopeGeometry(
        mirror_centers=np.array([[0.0, 0.0, 0.0]]),
        mirror_radii=np.array([6.067]),
        camera_rotations=site2cam,
        pmt_vectors=camera.pmt_vectors[np.newaxis, :, :],
        pmt_live_flags=np.ones((1, n_tubes), dtype=bool),
        mirror_pointing=axis_z[np.newaxis, :],
        mirror_diameters=np.array([3.3]),
        camera_separations=np.array([3.0]),
        spot_size_blur=np.array([0.005]),
        local_vsite=np.array([10000.0, 0.0, 1400.0]),
        site2clf=np.eye(3),
    )

    pmt_state = PMTState(
        gain=np.ones(n_tubes),
        pedestal=np.full(n_tubes, 100.0),
        pedestal_rms=np.full(n_tubes, 2.0),
        live_flag=np.ones(n_tubes, dtype=bool),
    )

    electronics = TAFADCFrontEnd(
        sample_rate_mhz=10.0,
        trace_length_bins=256,
        impulse_response=generate_ta_impulse_response(),
        nsb_rate_mhz=9.0,
    )

    return FDSite(
        site_id=0,
        geometry=geometry,
        camera=camera,
        pmt_state=pmt_state,
        electronics=electronics,
        wavelength_bands=np.array([350.0]),
    )


def test_fdsite_geometry_vectors(vertical_track, single_mirror_site):
    """Tests vectorized ports of getRpVector, getPsiAngle, getZenithAngle, getViewingAngles, applyEvsRpCut."""
    geom = single_mirror_site.geometry
    rp, rpv_site = get_rp_vector(
        geom.local_vsite, geom.mirror_centers, geom.mirror_pointing, geom.site2clf, vertical_track
    )

    # Vertical shower at (0, 0, 1400) observed from (10000, 0, 1400) has Rp = 10,000 m along -x
    assert rp == pytest.approx(10000.0, rel=1e-6)
    np.testing.assert_allclose(rpv_site, np.array([-10000.0, 0.0, 0.0]), atol=1e-6)

    psi = get_psi_angle(
        geom.local_vsite, geom.mirror_centers, geom.mirror_pointing, geom.site2clf, vertical_track
    )
    assert psi == pytest.approx(np.pi / 2.0, abs=1e-6)

    zen = get_zenith_angle(geom.site2clf, vertical_track)
    assert zen == pytest.approx(0.0, abs=1e-6)

    view_ang, chi_ang = get_viewing_angles(geom.local_vsite, geom.site2clf, vertical_track)
    assert view_ang.shape == (vertical_track.nseg,)
    assert chi_ang.shape == (vertical_track.nseg,)
    # Bottom segment is at same altitude as detector -> 90 deg viewing angle
    assert view_ang[-1] == pytest.approx(np.pi / 2.0, abs=1e-4)

    # E vs Rp cut check
    assert apply_e_vs_rp_cut(19.0, 10000.0) is True
    assert apply_e_vs_rp_cut(15.5, 30000.0) is False


def test_camera_optics_snell_and_hex():
    """Tests vectorized snell_shift and pmt_detect_hex."""
    # Normal incidence ray along +z produces zero lateral shift
    normal_ray = np.array([[0.0, 0.0, 1.0]])
    shift_zero = snell_shift(normal_ray, 1.000293, 1.50, 0.005)
    np.testing.assert_allclose(shift_zero, 0.0, atol=1e-12)

    # Oblique ray produces non-zero radial shift opposite to transverse direction
    oblique_ray = np.array([[0.3, 0.0, np.sqrt(1.0 - 0.09)]])
    shift_oblique = snell_shift(oblique_ray, 1.000293, 1.50, 0.005)
    assert shift_oblique[0, 0] < 0.0
    assert shift_oblique[0, 2] == 0.0

    # Hexagonal PMT detection
    cam = CameraAssembly.default_ta_camera()
    x_tube = cam.pmt_positions[:, 0]
    y_tube = cam.pmt_positions[:, 1]
    # Hitting exact tube centers should match 1-to-1
    hit_ids, valid = pmt_detect_hex(
        x_tube[:10], y_tube[:10], x_tube, y_tube, cam.pmt_flat2flat, cam.pmt_point2point
    )
    assert np.all(valid)
    np.testing.assert_array_equal(hit_ids, np.arange(10))

    # Ray far outside camera misses all tubes
    miss_ids, miss_valid = pmt_detect_hex(
        np.array([5.0]), np.array([5.0]), x_tube, y_tube, cam.pmt_flat2flat, cam.pmt_point2point
    )
    assert not miss_valid[0]
    assert miss_ids[0] == -2


def test_nkg_lateral_displacement():
    """Tests vectorized NKG lateral displacement."""
    rng = np.random.default_rng(123)
    ages = np.ones(200)
    molrads = get_moliere_radius(np.full(200, 1.0))
    rp_uv = np.array([-1.0, 0.0, 0.0])
    n_pln = np.array([0.0, 1.0, 0.0])

    disp = get_nkg_dr(ages, molrads, rp_uv, n_pln, rng)
    assert disp.shape == (200, 3)
    # Displacements lie purely in the (rp_uv, n_pln) plane -> z component is 0
    np.testing.assert_allclose(disp[:, 2], 0.0, atol=1e-12)
    assert np.mean(np.linalg.norm(disp, axis=1)) > 0.0


def test_end_to_end_fdsite_and_petimes(vertical_track, single_mirror_site):
    """Tests ObservedTrack, get_accepted_track -> PETimes, and full Simulation.run()."""
    atmos = Atmosphere(
        profile=HydrostaticProfile(),
        aerosol=AerosolModel(hal=25000.0, scale_height=1000.0, ref_wavelength=355.0),
        fy_model=KakimotoYield(),
    )

    class SingleTrackSource:
        def generate_event(self, rng, timestamp):
            return vertical_track

    class OneTrialSchedule:
        def get_trials(self):
            yield 0, 1000.0

    sim = Simulation(
        experiment=single_mirror_site,
        atmosphere=atmos,
        source=SingleTrackSource(),
        schedule=OneTrialSchedule(),
        seed=42,
    )

    results = sim.run()
    assert len(results) == 1
    assert "readout" in results.fields
    assert len(results[0].readout.pmt_id) == 256
    assert ak.to_numpy(results[0].readout.waveform).shape == (256, 256)
