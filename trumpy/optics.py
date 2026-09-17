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