import numpy as np
import pytest
from trumpy.track import Track
from trumpy.shower import GaisserHillasProfile

@pytest.fixture
def sample_gh_profile():
    """Provides a realistic 10^19 eV proton shower profile."""
    return GaisserHillasProfile(
        x0=-75.8,         # Typical legacy test value for proton
        xmax=773.2,
        nmax=6.692e9,
        lambda_inv=59.9
    )

def test_gaisser_hillas_physical_limits(sample_gh_profile):
    """Verifies the physical boundary conditions of the GH function."""
    gh = sample_gh_profile
    
    # Test 1: Zero particles before the first interaction
    x_before = np.array([-100.0, -80.0, gh.x0])
    n_before = gh.evaluate_particles(x_before)
    np.testing.assert_array_equal(n_before, 0.0)
    
    # Test 2: Peak exactly matches Nmax at Xmax
    x_at_max = np.array([gh.xmax])
    n_at_max = gh.evaluate_particles(x_at_max)
    np.testing.assert_allclose(n_at_max, gh.nmax, rtol=1e-5)
    
    # Test 3: Attenuation tail drops off as expected
    x_tail = np.array([1200.0, 1500.0])
    n_tail = gh.evaluate_particles(x_tail)
    assert n_tail[0] > n_tail[1], "Profile should strictly decrease deep in the tail."


def test_track_vectorized_initialization(sample_gh_profile):
    """Verifies that the Track dataclass correctly handles vectorized SoA inputs."""
    nseg = 500
    
    # Generate mock slant depths mapping to altitude/time
    slant_depths = np.linspace(0, 1200, nseg)
    
    # Populate the array structures
    dedep_array = sample_gh_profile.evaluate_dedep(slant_depths)
    time_array = np.linspace(0, 30000, nseg)  # mock time in ns
    alt_array = np.linspace(30000, 1400, nseg) # mock altitude in m
    dl_array = np.full(nseg, 50.0)            # 50m uniform segments
    
    track = Track(
        species=1, # Proton
        log_e=19.0,
        zenith=np.deg2rad(45.0),
        impact_v=np.array([0.0, 0.0, 1400.0]),
        track_uv=np.array([0.707, 0.0, -0.707]),
        nseg=nseg,
        time_gen=time_array,
        altitude=alt_array,
        dlseg=dl_array,
        dedep=dedep_array
    )
    
    # Assertions to guarantee Structure of Arrays consistency
    assert track.nseg == nseg
    assert track.dedep.shape == (nseg,)
    assert track.time_gen.shape == (nseg,)
    
    # Check that the max energy deposit aligns physically
    max_idx = np.argmax(track.dedep)
    assert slant_depths[max_idx] == pytest.approx(sample_gh_profile.xmax, abs=5.0)