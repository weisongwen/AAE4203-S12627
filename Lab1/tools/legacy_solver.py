# Numerical functions from GAO Yixin's 2025 MIT-licensed example.
# Kept unchanged for historical numerical regression, without plotting imports.
import numpy as np

def compute_elevation_angle(satellite_pos, receiver_pos):
    """Compute elevation angle between satellite and receiver"""
    # Convert to LLA for proper elevation calculation
    lat_ref, lon_ref, h_ref = ecef_to_lla(receiver_pos[0], receiver_pos[1], receiver_pos[2])
    
    # Vector from receiver to satellite
    sat_vector = satellite_pos - receiver_pos
    
    # Convert to ENU frame
    enu_vector = ecef_to_enu_vector(sat_vector, lat_ref, lon_ref)
    
    # Compute elevation angle
    horizontal_dist = np.sqrt(enu_vector[0]**2 + enu_vector[1]**2)
    elevation = np.arctan2(enu_vector[2], horizontal_dist)
    
    return elevation

def ecef_to_enu_vector(vector, lat_ref, lon_ref):
    """Convert ECEF vector to ENU frame"""
    lat_ref_rad = np.radians(lat_ref)
    lon_ref_rad = np.radians(lon_ref)
    
    # Transformation matrix
    t = np.array([
        [-np.sin(lon_ref_rad), np.cos(lon_ref_rad), 0],
        [-np.sin(lat_ref_rad)*np.cos(lon_ref_rad), -np.sin(lat_ref_rad)*np.sin(lon_ref_rad), np.cos(lat_ref_rad)],
        [np.cos(lat_ref_rad)*np.cos(lon_ref_rad), np.cos(lat_ref_rad)*np.sin(lon_ref_rad), np.sin(lat_ref_rad)]
    ])
    
    return t @ vector

def compute_weight_matrix(satellite_positions, receiver_position):
    """Compute weight matrix based on elevation angles"""
    n_sats = len(satellite_positions)
    weights = np.zeros(n_sats)
    
    for i, sat_pos in enumerate(satellite_positions):
        elev_angle = compute_elevation_angle(sat_pos, receiver_position)
        # Weight based on sine of elevation angle (higher elevation = higher weight)
        if elev_angle > np.radians(5):  # Mask low elevation satellites
            weights[i] = np.sin(elev_angle)**2
        else:
            weights[i] = 0.01  # Very low weight for low elevation
    
    return np.diag(weights)

def weighted_least_squares_solution(satellite_positions, receiver_position, pseudoranges_meas, 
                                   satellite_clock_bias, ionospheric_delay, tropospheric_delay):
    """Weighted least squares solution for receiver position"""
    receiver_clock_bias = 0.0
    
    for j in range(10):  # Maximum iterations
        # Compute geometric distances
        estimated_distances = np.linalg.norm(satellite_positions - receiver_position, axis=1)
        
        # Correct pseudorange measurements
        corrected_pseudoranges = pseudoranges_meas + satellite_clock_bias - ionospheric_delay - tropospheric_delay
        
        # Compute residuals
        pseudoranges_diff = corrected_pseudoranges - (estimated_distances + receiver_clock_bias)
        
        # Build design matrix G
        G = np.zeros((len(satellite_positions), 4))
        for i in range(len(satellite_positions)):
            p_i = satellite_positions[i] - receiver_position
            r_i = estimated_distances[i]
            G[i, :3] = -p_i / r_i
            G[i, 3] = 1.0
        
        # Compute weight matrix
        W = compute_weight_matrix(satellite_positions, receiver_position)
        
        # Weighted least squares solution
        try:
            delta_p = np.linalg.inv(G.T @ W @ G) @ G.T @ W @ pseudoranges_diff
        except np.linalg.LinAlgError:
            # Fallback to regular least squares if weights cause singularity
            delta_p, _, _, _ = np.linalg.lstsq(G, pseudoranges_diff, rcond=None)
        
        receiver_position += delta_p[:3]
        receiver_clock_bias += delta_p[3]
        
        # Check convergence
        if np.linalg.norm(delta_p[:3]) < 1e-4:
            break
    
    return receiver_position, receiver_clock_bias

def ecef_to_lla(x, y, z):
    """Convert ECEF coordinates to Latitude, Longitude, Altitude"""
    # WGS84 ellipsoid constants
    a = 6378137.0  # Semi-major axis in meters
    e2 = 6.69437999014e-3  # Square of eccentricity
    
    # Longitude calculation
    lon = np.arctan2(y, x)
    
    # Latitude and altitude iterative calculation
    p = np.sqrt(x**2 + y**2)
    lat = np.arctan2(z, p * (1 - e2))
    lat_prev = 0
    
    while np.abs(lat - lat_prev) > 1e-12:
        lat_prev = lat
        N = a / np.sqrt(1 - e2 * np.sin(lat)**2)
        h = p / np.cos(lat) - N
        lat = np.arctan2(z, p * (1 - e2 * N / (N + h)))
    
    N = a / np.sqrt(1 - e2 * np.sin(lat)**2)
    h = p / np.cos(lat) - N
    
    return np.degrees(lat), np.degrees(lon), h

def ecef_to_enu(x, y, z, x_ref, y_ref, z_ref):
    """Convert ECEF to local ENU coordinates"""
    lat_ref, lon_ref, h_ref = ecef_to_lla(x_ref, y_ref, z_ref)
    lat_ref_rad = np.radians(lat_ref)
    lon_ref_rad = np.radians(lon_ref)
    
    # Difference vector
    dx = x - x_ref
    dy = y - y_ref
    dz = z - z_ref
    
    # Transformation matrix
    t = np.array([
        [-np.sin(lon_ref_rad), np.cos(lon_ref_rad), 0],
        [-np.sin(lat_ref_rad)*np.cos(lon_ref_rad), -np.sin(lat_ref_rad)*np.sin(lon_ref_rad), np.cos(lat_ref_rad)],
        [np.cos(lat_ref_rad)*np.cos(lon_ref_rad), np.cos(lat_ref_rad)*np.sin(lon_ref_rad), np.sin(lat_ref_rad)]
    ])
    
    return t @ np.array([dx, dy, dz])
