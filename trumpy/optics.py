import numpy as np

def mirror_intersection(origins, directions, center, radius):
    """Calculates intersection distances for N rays hitting a spherical mirror."""
    # V is shape (N, 3)
    V = origins - center 
    
    # Batched dot products: shape (N,)
    V_dot_U = np.sum(V * directions, axis=1)
    V_mag_sq = np.sum(V * V, axis=1)
    
    # Quadratic coefficients
    # a = 1.0 (since directions are unit vectors)
    b = 2.0 * V_dot_U
    c = V_mag_sq - (radius ** 2)
    
    discriminant = (b ** 2) - (4.0 * c)
    
    # Mask out rays that completely miss the mirror sphere
    hit_mask = discriminant >= 0
    
    # Calculate distance t (taking the closer intersection)
    sqrt_disc = np.sqrt(discriminant[hit_mask])
    t = (-b[hit_mask] - sqrt_disc) / 2.0
    
    # Update origins to the intersection point
    hit_points = origins[hit_mask] + t[:, np.newaxis] * directions[hit_mask]
    
    return hit_points, hit_mask, t

def reflect_and_blur(directions, hit_points, center, radius, blur_angle_rad, rng):
    """Reflects N rays off the mirror and applies a Gaussian angular blur."""
    normals = (center - hit_points) / radius
    
    # U_dot_N shape: (N, 1) to allow broadcasting against (N, 3) arrays
    U_dot_N = np.sum(directions * normals, axis=1)[:, np.newaxis]
    
    reflected_dirs = directions - 2.0 * U_dot_N * normals
    
    if blur_angle_rad > 0:
        # Add random 3D Gaussian noise to perturb the direction
        noise = rng.normal(scale=blur_angle_rad, size=reflected_dirs.shape)
        blurred_dirs = reflected_dirs + noise
        
        # Re-normalize
        mags = np.linalg.norm(blurred_dirs, axis=1)[:, np.newaxis]
        reflected_dirs = blurred_dirs / mags
        
    return reflected_dirs

def map_to_pmts(camera_hit_dirs, pmt_dirs, fov_limit_rad):
    """Maps N rays to 256 PMTs using dot products to find the closest tube axis."""
    # camera_hit_dirs is (N, 3)
    # pmt_dirs is (256, 3)
    # Resulting dot product matrix is (N, 256)
    cosine_matrix = np.dot(camera_hit_dirs, pmt_dirs.T)
    
    # Find the PMT with the largest cosine (smallest angle) for each ray
    hit_pmt_indices = np.argmax(cosine_matrix, axis=1)
    max_cosines = np.max(cosine_matrix, axis=1)
    
    # Mask out rays that don't fall within the physical FoV of the closest PMT
    valid_pmt_mask = max_cosines > np.cos(fov_limit_rad)
    
    return hit_pmt_indices[valid_pmt_mask], valid_pmt_mask

def is_within_aperture(hit_points, mirror_vertex, mirror_axis, aperture_radius):
    """
    Vectorized check if rays hitting the mirror sphere fall within the physical rim.
    
    Parameters:
    hit_points      : (N, 3) array of ray intersection points on the mirror sphere.
    mirror_vertex   : (3,) array representing the center point of the mirror face.
    mirror_axis     : (3,) unit vector pointing along the mirror's optical axis.
    aperture_radius : float, the physical transverse radius of the mirror.
    
    Returns:
    valid_aperture_mask : (N,) boolean array where True means the ray hit the physical mirror.
    """
    # Vector from the center of the mirror face to the hit points
    v = hit_points - mirror_vertex
    
    # Projection of v onto the optical axis
    z_proj = np.sum(v * mirror_axis, axis=1)
    
    # Transverse distance squared from the optical axis (Pythagorean theorem)
    r_transverse_sq = np.sum(v**2, axis=1) - (z_proj ** 2)
    
    return r_transverse_sq <= (aperture_radius ** 2)


def is_shadowed(origins, directions, hit_distances, camera_center, shadow_radius):
    """
    Vectorized check if rays pass through the camera cluster before hitting the mirror.
    Models the camera shadow boundary as a sphere to enable fast broadcasted culling.
    
    Parameters:
    origins         : (N, 3) array of emission coordinates.
    directions      : (N, 3) array of unit direction vectors for the rays.
    hit_distances   : (N,) array of distances from origin to the mirror intersection.
    camera_center   : (3,) array representing the geometric center of the PMT cluster.
    shadow_radius   : float, bounding radius of the camera box shadow.
    
    Returns:
    shadow_mask     : (N,) boolean array where True means the ray is blocked.
    """
    # Vector from emission origins to the camera center
    oc = camera_center - origins
    
    # Distance along the ray to the point of closest approach to the camera center
    t_closest = np.sum(oc * directions, axis=1)
    
    # Transverse distance squared from the camera center to the ray line
    d_sq = np.sum(oc**2, axis=1) - (t_closest ** 2)
    
    # A ray is shadowed if:
    # 1. The closest approach is within the camera's bounding radius.
    # 2. The camera is strictly in front of the emission point (t_closest > 0).
    # 3. The camera is reached before the mirror (t_closest < hit_distances).
    shadow_mask = (d_sq <= (shadow_radius ** 2)) & (t_closest > 0) & (t_closest < hit_distances)
    
    return shadow_mask

def is_shadowed_by_box(origins, directions, hit_distances, box_center, box_half_extents, box_axes):
    """
    Vectorized exact shadow check for an Oriented Bounding Box (e.g., Camera Housing).
    Uses the Slab Method.
    
    box_half_extents : (3,) array of the box's half-width, half-height, half-depth
    box_axes         : (3, 3) array where rows are the local u, v, w unit vectors of the box
    """
    N = len(origins)
    t_min = np.full(N, -np.inf)
    t_max = np.full(N, np.inf)
    
    p = box_center - origins  # Vector from ray origin to box center
    
    for i in range(3):
        axis = box_axes[i]
        e = np.sum(axis * p, axis=1)          # Projection of p onto box axis
        f = np.sum(axis * directions, axis=1) # Projection of ray dir onto box axis
        
        # Avoid division by zero for rays exactly parallel to a box face
        f_safe = np.where(np.abs(f) < 1e-8, 1e-8, f)
        
        t1 = (e + box_half_extents[i]) / f_safe
        t2 = (e - box_half_extents[i]) / f_safe
        
        # Slabs min/max
        t_near = np.minimum(t1, t2)
        t_far = np.maximum(t1, t2)
        
        t_min = np.maximum(t_min, t_near)
        t_max = np.minimum(t_max, t_far)
    
    # A ray hits the box if the largest entry time is less than the smallest exit time
    # AND the box is in front of the ray (t_max > 0)
    # AND the box is reached before the mirror (t_min < hit_distances)
    hits_box = (t_min <= t_max) & (t_max > 0) & (t_min < hit_distances)
    
    return hits_box

def is_shadowed_by_cylinder_strut(origins, directions, hit_distances, strut_pt1, strut_pt2, radius):
    """
    Vectorized exact shadow check for a structural strut modeled as a cylinder.
    """
    # Vector defining the cylinder axis
    strut_axis = strut_pt2 - strut_pt1
    strut_length_sq = np.sum(strut_axis**2)
    strut_dir = strut_axis / np.sqrt(strut_length_sq)
    
    # Vector from cylinder start to ray origin
    dp = origins - strut_pt1
    
    # Projections
    v_dot_axis = np.sum(directions * strut_dir, axis=1)
    dp_dot_axis = np.sum(dp * strut_dir, axis=1)
    
    # Quadratic coefficients for infinite cylinder intersection
    a = 1.0 - v_dot_axis**2
    b = 2.0 * (np.sum(directions * dp, axis=1) - v_dot_axis * dp_dot_axis)
    c = np.sum(dp**2, axis=1) - dp_dot_axis**2 - radius**2
    
    discriminant = b**2 - 4*a*c
    hit_inf_cylinder = discriminant >= 0
    
    # Calculate intersection distances for rays that hit the infinite cylinder
    sqrt_disc = np.sqrt(np.maximum(discriminant, 0))
    t1 = (-b - sqrt_disc) / (2*a)
    t2 = (-b + sqrt_disc) / (2*a)
    
    t_first = np.minimum(t1, t2)
    
    # Check if the intersection point lies within the finite length of the strut
    hit_points = origins + t_first[:, np.newaxis] * directions
    hit_proj = np.sum((hit_points - strut_pt1) * strut_dir, axis=1)
    
    valid_finite_hit = (hit_proj >= 0) & (hit_proj**2 <= strut_length_sq)
    
    # Shadowed if it hits the finite cylinder, in front of the origin, before the mirror
    shadowed = hit_inf_cylinder & valid_finite_hit & (t_first > 0) & (t_first < hit_distances)
    
    return shadowed


def snell_shift(ray_dirs: np.ndarray, n_air: float, n_medium: float, thickness: float) -> np.ndarray:
    """
    Vectorized refraction shift through a planar optical layer (Paraglas cover or BG3 filter).
    Direct port of snellShift() in legacy/src/raytrace/camera.c.

    Parameters
    ----------
    ray_dirs : (N, 3) array of unit direction vectors of reflected rays
    n_air : refractive index of air (e.g., 1.000293)
    n_medium : refractive index of optical layer (e.g., 1.50 for Paraglas, 1.53 for BG3)
    thickness : thickness of layer in meters

    Returns
    -------
    shift : (N, 3) lateral displacement vectors in the camera plane (z-component is 0)
    """
    sinthe1 = np.sqrt(np.clip(1.0 - ray_dirs[:, 2] ** 2, 0.0, 1.0))
    sinthe2 = (n_air / n_medium) * sinthe1

    dist = thickness * (sinthe1 - sinthe2)
    phi = np.arctan2(ray_dirs[:, 1], ray_dirs[:, 0])

    shift = np.zeros_like(ray_dirs)
    shift[:, 0] = -dist * np.cos(phi)
    shift[:, 1] = -dist * np.sin(phi)
    return shift


def pmt_detect_hex(
    x_cam: np.ndarray,
    y_cam: np.ndarray,
    x_tube: np.ndarray,
    y_tube: np.ndarray,
    flat2flat: float,
    point2point: float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Vectorized hexagonal PMT detection on the camera focal plane.
    Direct port of pmtDetect() in legacy/src/raytrace/camera.c.

    Parameters
    ----------
    x_cam, y_cam : (N,) arrays of ray impact coordinates on the PMT plane (m)
    x_tube, y_tube : (M,) arrays of PMT center coordinates on the PMT plane (m)
    flat2flat : float, hexagonal PMT flat-to-flat width (m)
    point2point : float, hexagonal PMT point-to-point diameter (m)

    Returns
    -------
    hit_pmt_indices : (N,) array of PMT indices (or -2 for MISSED_TUBE)
    valid_mask : (N,) boolean mask where True indicates a valid PMT hit
    """
    if len(x_cam) == 0:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=bool)

    dx = x_cam[:, np.newaxis] - x_tube[np.newaxis, :]
    dy = y_cam[:, np.newaxis] - y_tube[np.newaxis, :]
    dist = np.sqrt(dx ** 2 + dy ** 2)

    closest_tube = np.argmin(dist, axis=1)
    min_dist = dist[np.arange(len(x_cam)), closest_tube]

    # Fast inner-circle accept and outer-circle reject
    inside_inscribed = min_dist < (flat2flat / 2.0)
    within_circumscribed = min_dist <= (point2point / 2.0)

    # Hexagonal boundary check for rays between inscribed and circumscribed radii
    angles = np.deg2rad(30.0 + np.arange(6) * 60.0)
    hex_x = -np.cos(angles)  # shape (6,)
    hex_y = -np.sin(angles)  # shape (6,)

    rel_x = (x_cam - x_tube[closest_tube])[:, np.newaxis]  # (N, 1)
    rel_y = (y_cam - y_tube[closest_tube])[:, np.newaxis]  # (N, 1)

    seg_x = (point2point / 2.0) * hex_x[np.newaxis, :] + rel_x  # (N, 6)
    seg_y = (point2point / 2.0) * hex_y[np.newaxis, :] + rel_y  # (N, 6)
    seg_norm = np. clip(np.sqrt(seg_x ** 2 + seg_y ** 2), 1e-12, None)

    dots = (hex_x[np.newaxis, :] * seg_x + hex_y[np.newaxis, :] * seg_y) / seg_norm
    inside_hex = np.all(dots >= 0.5, axis=1)

    valid_mask = inside_inscribed | (within_circumscribed & inside_hex)
    hit_pmt_indices = np.where(valid_mask, closest_tube, -2)
    return hit_pmt_indices, valid_mask
