import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from matplotlib.patches import Rectangle
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

# Read data from CSV files
satellite_positions = np.loadtxt('satellite_positions.csv', delimiter=',')  # (max_num_sats, num_epochs*3)
pseudoranges_meas = np.loadtxt('pseudoranges_meas.csv', delimiter=',')  # (max_num_sats, num_epochs)
satellite_clock_bias = np.loadtxt('satellite_clock_bias.csv', delimiter=',')  # (max_num_sats, num_epochs)
ionospheric_delay = np.loadtxt('ionospheric_delay.csv', delimiter=',')  # (max_num_sats, num_epochs)
tropospheric_delay = np.loadtxt('tropospheric_delay.csv', delimiter=',')  # (max_num_sats, num_epochs)
ground_truth = np.genfromtxt('NAV-HPPOSECEF.csv', delimiter=',', skip_header=1, usecols=[2,3,4,5])  # Ground truth positions

# Constants
num_epochs = pseudoranges_meas.shape[1]
max_num_sats = pseudoranges_meas.shape[0]
c = 299792458.0  # Speed of light in m/s

# Initialize storage arrays
estimated_positions = []
estimated_clock_biases = []
receiver_position = np.array([0.0, 0.0, 0.0])  # Initial position estimate
ground_truth_position = ground_truth[:,1:4]/100  # Convert to meters

def degrees_to_dms(degrees, is_longitude=False):
    """Convert decimal degrees to degrees-minutes-seconds format with N/S/E/W"""
    abs_degrees = abs(degrees)
    d = int(abs_degrees)
    m = int((abs_degrees - d) * 60)
    s = ((abs_degrees - d) * 60 - m) * 60
    
    if is_longitude:
        direction = 'E' if degrees >= 0 else 'W'
    else:
        direction = 'N' if degrees >= 0 else 'S'
    
    return f"{d}°{m:02d}'{s:04.1f}\"{direction}"

def format_coordinate_axis(ax, axis='lon'):
    """Format axis with degree-minute-second labels"""
    if axis == 'lon':
        formatter = LongitudeFormatter(degree_symbol='°', minute_symbol="'", second_symbol='"')
        ax.xaxis.set_major_formatter(formatter)
    else:
        formatter = LatitudeFormatter(degree_symbol='°', minute_symbol="'", second_symbol='"')
        ax.yaxis.set_major_formatter(formatter)

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

# Main processing loop
print("Processing GNSS positioning...")
for epoch in range(num_epochs):
    # Extract current epoch data
    p_l1_epoch = pseudoranges_meas[:, epoch]
    sat_clock_err_epoch = satellite_clock_bias[:, epoch]
    ion_error_l1_epoch = ionospheric_delay[:, epoch]
    tropo_error_epoch = tropospheric_delay[:, epoch]
    sat_pos_epoch = satellite_positions[:, epoch*3:(epoch+1)*3]
    
    # Filter valid measurements
    valid_idx = (~np.isnan(p_l1_epoch) & 
                 ~np.isnan(sat_clock_err_epoch) & 
                 ~np.isnan(ion_error_l1_epoch) & 
                 ~np.isnan(tropo_error_epoch) & 
                 ~np.isnan(sat_pos_epoch).any(axis=1))
    
    # Check minimum satellite requirement
    if np.sum(valid_idx) < 4:
        print(f"Epoch {epoch+1}: Insufficient satellites ({np.sum(valid_idx)})")
        if epoch > 0:
            estimated_positions.append(estimated_positions[-1])
            estimated_clock_biases.append(estimated_clock_biases[-1])
        else:
            estimated_positions.append(receiver_position.copy())
            estimated_clock_biases.append(0.0)
        continue
    
    # Extract valid data
    p_l1_valid = p_l1_epoch[valid_idx]
    sat_clock_err_valid = sat_clock_err_epoch[valid_idx]
    ion_error_l1_valid = ion_error_l1_epoch[valid_idx]
    tropo_error_valid = tropo_error_epoch[valid_idx]
    sat_pos_valid = sat_pos_epoch[valid_idx, :]
    
    # Use previous position as initial guess
    if epoch > 0:
        receiver_position = estimated_positions[-1].copy()
    else:
        receiver_position = np.array([0.0, 0.0, 0.0])
    
    # Perform weighted least squares estimation
    estimated_position, estimated_receiver_clock_bias = weighted_least_squares_solution(
        sat_pos_valid, receiver_position, p_l1_valid, sat_clock_err_valid, 
        ion_error_l1_valid, tropo_error_valid
    )
    
    # Store results
    estimated_positions.append(estimated_position.copy())
    estimated_clock_biases.append(estimated_receiver_clock_bias)

print(f"Processed {len(estimated_positions)} epochs")

# Convert positions to LLA coordinates
lat_est, lon_est, alt_est = [], [], []
lat_gt, lon_gt, alt_gt = [], [], []

for pos in estimated_positions:
    lat, lon, alt = ecef_to_lla(pos[0], pos[1], pos[2])
    lat_est.append(lat)
    lon_est.append(lon)
    alt_est.append(alt)

for pos_gt in ground_truth_position:
    lat, lon, alt = ecef_to_lla(pos_gt[0], pos_gt[1], pos_gt[2])
    lat_gt.append(lat)
    lon_gt.append(lon)
    alt_gt.append(alt)

# Convert to ENU coordinates for error analysis
x_ref, y_ref, z_ref = ground_truth_position[0]  # Use first ground truth as reference
enu_est, enu_gt = [], []

for pos in estimated_positions:
    enu = ecef_to_enu(pos[0], pos[1], pos[2], x_ref, y_ref, z_ref)
    enu_est.append(enu)

for pos_gt in ground_truth_position:
    enu = ecef_to_enu(pos_gt[0], pos_gt[1], pos_gt[2], x_ref, y_ref, z_ref)
    enu_gt.append(enu)

# Calculate positioning errors
errors_e = [est[0] - gt[0] for est, gt in zip(enu_est, enu_gt)]
errors_n = [est[1] - gt[1] for est, gt in zip(enu_est, enu_gt)]
errors_u = [est[2] - gt[2] for est, gt in zip(enu_est, enu_gt)]
errors_2d = [np.sqrt(e**2 + n**2) for e, n in zip(errors_e, errors_n)]
errors_3d = [np.sqrt(e**2 + n**2 + u**2) for e, n, u in zip(errors_e, errors_n, errors_u)]

# Save estimated positions to CSV
results_df = pd.DataFrame({
    'Epoch': range(1, len(estimated_positions) + 1),
    'Latitude_deg': lat_est,
    'Longitude_deg': lon_est,
    'Altitude_m': alt_est,
    'Clock_Bias_m': estimated_clock_biases,
    'East_m': [enu[0] for enu in enu_est],
    'North_m': [enu[1] for enu in enu_est],
    'Up_m': [enu[2] for enu in enu_est],
    'ECEF_X_m': [pos[0] for pos in estimated_positions],
    'ECEF_Y_m': [pos[1] for pos in estimated_positions],
    'ECEF_Z_m': [pos[2] for pos in estimated_positions]
})
results_df.to_csv('estimated_positions.csv', index=False)

# Save positioning errors to CSV
errors_df = pd.DataFrame({
    'Epoch': range(1, len(errors_e) + 1),
    'Error_East_m': errors_e,
    'Error_North_m': errors_n,
    'Error_Up_m': errors_u,
    'Error_2D_m': errors_2d,
    'Error_3D_m': errors_3d
})
errors_df.to_csv('positioning_errors.csv', index=False)

# Print statistics
print(f"\nPositioning Statistics:")
print(f"Mean 2D Error: {np.mean(errors_2d):.3f} m")
print(f"RMS 2D Error: {np.sqrt(np.mean(np.array(errors_2d)**2)):.3f} m")
print(f"Mean 3D Error: {np.mean(errors_3d):.3f} m")
print(f"RMS 3D Error: {np.sqrt(np.mean(np.array(errors_3d)**2)):.3f} m")

# Enhanced plotting with modern style
plt.style.use('default')
plt.rcParams.update({
    'font.size': 12,
    'font.family': 'Arial',
    'axes.linewidth': 0.8,
    'axes.spines.right': False,
    'axes.spines.top': False,
    'grid.alpha': 0.3,
    'grid.linewidth': 0.5,
})

# Plot 1: Trajectory on map with cartopy
try:
    fig = plt.figure(figsize=(18, 8))
    
    # Map plot with cartopy
    ax1 = fig.add_subplot(121, projection=ccrs.PlateCarree())
    
    # Calculate bounds for the map
    lon_min, lon_max = min(min(lon_est), min(lon_gt)), max(max(lon_est), max(lon_gt))
    lat_min, lat_max = min(min(lat_est), min(lat_gt)), max(max(lat_est), max(lat_gt))
    
    # Add margin
    margin = 0.001  # degrees
    lon_min -= margin
    lon_max += margin
    lat_min -= margin
    lat_max += margin
    
    ax1.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    
    # Add map features
    ax1.add_feature(cfeature.LAND, alpha=0.7, color='lightgray')
    ax1.add_feature(cfeature.OCEAN, alpha=0.7, color='lightblue')
    ax1.add_feature(cfeature.COASTLINE, linewidth=0.8, color='black')
    ax1.add_feature(cfeature.BORDERS, linewidth=0.5, color='gray')
    ax1.add_feature(cfeature.LAKES, alpha=0.7, color='lightblue')
    ax1.add_feature(cfeature.RIVERS, linewidth=0.5, color='blue')
    
    # Plot trajectories
    ax1.plot(lon_est, lat_est, 'r-', linewidth=3, label='Estimated', alpha=0.8, transform=ccrs.PlateCarree())
    #ax1.plot(lon_gt, lat_gt, 'g-', linewidth=3, label='Ground Truth', alpha=0.8, transform=ccrs.PlateCarree())
    ax1.plot(lon_est[0], lat_est[0], 'ro', markersize=10, label='Start (Est)', transform=ccrs.PlateCarree(), zorder=10)
    ax1.plot(lon_est[-1], lat_est[-1], 'rs', markersize=10, label='End (Est)', transform=ccrs.PlateCarree(), zorder=10)
    #ax1.plot(lon_gt[0], lat_gt[0], 'go', markersize=10, label='Start (GT)', transform=ccrs.PlateCarree(), zorder=10)
    #ax1.plot(lon_gt[-1], lat_gt[-1], 'gs', markersize=10, label='End (GT)', transform=ccrs.PlateCarree(), zorder=10)
    
    # Format coordinate labels with DMS
    format_coordinate_axis(ax1, 'lon')
    format_coordinate_axis(ax1, 'lat')
    
    ax1.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, alpha=0.5)
    ax1.set_title('GNSS Receiver Trajectory (Geographic View)', fontsize=14, fontweight='bold', pad=20)
    ax1.legend(loc='upper right', frameon=True, fancybox=True, shadow=True)
    
    # 3D trajectory plot without outer frame labels
    ax2 = fig.add_subplot(122, projection='3d')
    east_est = [enu[0] for enu in enu_est]
    north_est = [enu[1] for enu in enu_est]
    up_est = [enu[2] for enu in enu_est]
    east_gt = [enu[0] for enu in enu_gt]
    north_gt = [enu[1] for enu in enu_gt]
    up_gt = [enu[2] for enu in enu_gt]
    
    ax2.plot(east_est, north_est, up_est, 'r-', linewidth=3, label='Estimated', alpha=0.8)
    ax2.plot(east_gt, north_gt, up_gt, 'g-', linewidth=3, label='Ground Truth', alpha=0.8)
    ax2.scatter(east_est[0], north_est[0], up_est[0], color='red', s=100, marker='o', label='Start', alpha=0.9)
    ax2.scatter(east_est[-1], north_est[-1], up_est[-1], color='red', s=100, marker='s', alpha=0.9)
    ax2.scatter(east_gt[0], north_gt[0], up_gt[0], color='green', s=100, marker='o', alpha=0.9)
    ax2.scatter(east_gt[-1], north_gt[-1], up_gt[-1], color='green', s=100, marker='s', alpha=0.9)
    
    ax2.set_xlabel('East (m)', fontsize=12, labelpad=10)
    ax2.set_ylabel('North (m)', fontsize=12, labelpad=10)
    ax2.set_zlabel('Up (m)', fontsize=12, labelpad=10)
    ax2.set_title('3D Trajectory (ENU Frame)', fontsize=14, fontweight='bold', pad=20)
    ax2.legend(loc='upper right', frameon=True)
    
    # Remove outer frame tick labels (keep only axis labels)
    ax2.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=True)
    ax2.tick_params(axis='y', which='both', left=False, right=False, labelleft=True)
    ax2.tick_params(axis='z', which='both', labelleft=True)
    
    # Set background color
    ax2.xaxis.pane.fill = False
    ax2.yaxis.pane.fill = False
    ax2.zaxis.pane.fill = False
    ax2.xaxis.pane.set_edgecolor('gray')
    ax2.yaxis.pane.set_edgecolor('gray')
    ax2.zaxis.pane.set_edgecolor('gray')
    ax2.xaxis.pane.set_alpha(0.1)
    ax2.yaxis.pane.set_alpha(0.1)
    ax2.zaxis.pane.set_alpha(0.1)
    
    plt.tight_layout()
    plt.show()

except ImportError:
    # Fallback plot without cartopy
    print("Cartopy not available, using fallback plotting...")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8))
    
    # Simple lat/lon plot with enhanced styling
    ax1.plot(lon_est, lat_est, 'r-', linewidth=3, label='Estimated', alpha=0.8)
    ax1.plot(lon_gt, lat_gt, 'g-', linewidth=3, label='Ground Truth', alpha=0.8)
    ax1.plot(lon_est[0], lat_est[0], 'ro', markersize=10, label='Start (Est)')
    ax1.plot(lon_est[-1], lat_est[-1], 'rs', markersize=10, label='End (Est)')
    ax1.plot(lon_gt[0], lat_gt[0], 'go', markersize=10, label='Start (GT)')
    ax1.plot(lon_gt[-1], lat_gt[-1], 'gs', markersize=10, label='End (GT)')
    
    # Format axis labels to DMS
    ax1.set_xlabel('Longitude', fontsize=12)
    ax1.set_ylabel('Latitude', fontsize=12)
    ax1.set_title('GNSS Receiver Trajectory', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend(frameon=True, fancybox=True, shadow=True)
    ax1.set_aspect('equal', adjustable='box')
    
    # Custom tick formatter for DMS
    lon_ticks = ax1.get_xticks()
    lat_ticks = ax1.get_yticks()
    ax1.set_xticklabels([degrees_to_dms(lon, True) for lon in lon_ticks])
    ax1.set_yticklabels([degrees_to_dms(lat, False) for lat in lat_ticks])
    
    # 3D trajectory plot
    ax2 = fig.add_subplot(122, projection='3d')
    east_est = [enu[0] for enu in enu_est]
    north_est = [enu[1] for enu in enu_est]
    up_est = [enu[2] for enu in enu_est]
    east_gt = [enu[0] for enu in enu_gt]
    north_gt = [enu[1] for enu in enu_gt]
    up_gt = [enu[2] for enu in enu_gt]
    
    ax2.plot(east_est, north_est, up_est, 'r-', linewidth=3, label='Estimated', alpha=0.8)
    ax2.plot(east_gt, north_gt, up_gt, 'g-', linewidth=3, label='Ground Truth', alpha=0.8)
    ax2.set_xlabel('East (m)', fontsize=12)
    ax2.set_ylabel('North (m)', fontsize=12)
    ax2.set_zlabel('Up (m)', fontsize=12)
    ax2.set_title('3D Trajectory (ENU Frame)', fontsize=14, fontweight='bold')
    ax2.legend()
    
    plt.tight_layout()
    plt.show()

# Plot 2: Enhanced positioning errors with modern styling
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
epochs = range(1, len(errors_e) + 1)

# Color palette
colors = {'east': '#2E86C1', 'north': '#E74C3C', 'up': '#58D68D', '2d': '#AF7AC5', '3d': '#F39C12'}

# East error
ax1.plot(epochs, errors_e, color=colors['east'], linewidth=2, alpha=0.8)
ax1.fill_between(epochs, errors_e, alpha=0.3, color=colors['east'])
ax1.set_xlabel('Epoch', fontsize=12)
ax1.set_ylabel('East Error (m)', fontsize=12)
ax1.set_title('East Positioning Error', fontsize=13, fontweight='bold')
ax1.grid(True, alpha=0.3)
ax1.axhline(y=0, color='k', linestyle='--', alpha=0.7, linewidth=1)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)

# North error
ax2.plot(epochs, errors_n, color=colors['north'], linewidth=2, alpha=0.8)
ax2.fill_between(epochs, errors_n, alpha=0.3, color=colors['north'])
ax2.set_xlabel('Epoch', fontsize=12)
ax2.set_ylabel('North Error (m)', fontsize=12)
ax2.set_title('North Positioning Error', fontsize=13, fontweight='bold')
ax2.grid(True, alpha=0.3)
ax2.axhline(y=0, color='k', linestyle='--', alpha=0.7, linewidth=1)
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)

# Up error
ax3.plot(epochs, errors_u, color=colors['up'], linewidth=2, alpha=0.8)
ax3.fill_between(epochs, errors_u, alpha=0.3, color=colors['up'])
ax3.set_xlabel('Epoch', fontsize=12)
ax3.set_ylabel('Up Error (m)', fontsize=12)
ax3.set_title('Up Positioning Error', fontsize=13, fontweight='bold')
ax3.grid(True, alpha=0.3)
ax3.axhline(y=0, color='k', linestyle='--', alpha=0.7, linewidth=1)
ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)

# 2D and 3D errors
ax4.plot(epochs, errors_2d, color=colors['2d'], linewidth=2, label='2D Error', alpha=0.8)
ax4.plot(epochs, errors_3d, color=colors['3d'], linewidth=2, label='3D Error', alpha=0.8)
ax4.fill_between(epochs, errors_2d, alpha=0.2, color=colors['2d'])
ax4.fill_between(epochs, errors_3d, alpha=0.2, color=colors['3d'])
ax4.set_xlabel('Epoch', fontsize=12)
ax4.set_ylabel('Error (m)', fontsize=12)
ax4.set_title('2D and 3D Positioning Errors', fontsize=13, fontweight='bold')
ax4.grid(True, alpha=0.3)
ax4.legend(frameon=True, fancybox=True, shadow=False)
ax4.spines['right'].set_visible(False)
ax4.spines['top'].set_visible(False)

plt.tight_layout()
plt.show()

# Plot 3: Enhanced error statistics histogram
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# 2D error histogram
n1, bins1, patches1 = ax1.hist(errors_2d, bins=30, alpha=0.8, color=colors['2d'], edgecolor='white', linewidth=0.8)
ax1.set_xlabel('2D Error (m)', fontsize=12)
ax1.set_ylabel('Frequency', fontsize=12)
ax1.set_title('2D Error Distribution', fontsize=13, fontweight='bold')
ax1.axvline(np.mean(errors_2d), color='red', linestyle='--', linewidth=2, alpha=0.8, 
           label=f'Mean: {np.mean(errors_2d):.3f} m')
ax1.axvline(np.sqrt(np.mean(np.array(errors_2d)**2)), color='orange', linestyle='--', linewidth=2, alpha=0.8,
           label=f'RMS: {np.sqrt(np.mean(np.array(errors_2d)**2)):.3f} m')
ax1.legend(frameon=True, fancybox=True)
ax1.grid(True, alpha=0.3)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)

# Add gradient coloring to histogram bars
for i, p in enumerate(patches1):
    p.set_facecolor(plt.cm.viridis(i / len(patches1)))

# 3D error histogram  
n2, bins2, patches2 = ax2.hist(errors_3d, bins=30, alpha=0.8, color=colors['3d'], edgecolor='white', linewidth=0.8)
ax2.set_xlabel('3D Error (m)', fontsize=12)
ax2.set_ylabel('Frequency', fontsize=12)
ax2.set_title('3D Error Distribution', fontsize=13, fontweight='bold')
ax2.axvline(np.mean(errors_3d), color='red', linestyle='--', linewidth=2, alpha=0.8,
           label=f'Mean: {np.mean(errors_3d):.3f} m')
ax2.axvline(np.sqrt(np.mean(np.array(errors_3d)**2)), color='orange', linestyle='--', linewidth=2, alpha=0.8,
           label=f'RMS: {np.sqrt(np.mean(np.array(errors_3d)**2)):.3f} m')
ax2.legend(frameon=True, fancybox=True)
ax2.grid(True, alpha=0.3)
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)

# Add gradient coloring to histogram bars
for i, p in enumerate(patches2):
    p.set_facecolor(plt.cm.plasma(i / len(patches2)))

plt.tight_layout()
plt.show()

# Plot 4: Additional analysis - Error vs Time and Scatter Plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

# Error vs Time with rolling average
window_size = max(1, len(errors_2d) // 20)  # Adaptive window size
if len(errors_2d) > window_size:
    errors_2d_smooth = pd.Series(errors_2d).rolling(window=window_size, center=True).mean()
    errors_3d_smooth = pd.Series(errors_3d).rolling(window=window_size, center=True).mean()
    
    ax1.plot(epochs, errors_2d, color=colors['2d'], alpha=0.4, linewidth=1, label='2D Raw')
    ax1.plot(epochs, errors_2d_smooth, color=colors['2d'], linewidth=3, label='2D Smoothed')
    ax1.plot(epochs, errors_3d, color=colors['3d'], alpha=0.4, linewidth=1, label='3D Raw')  
    ax1.plot(epochs, errors_3d_smooth, color=colors['3d'], linewidth=3, label='3D Smoothed')
else:
    ax1.plot(epochs, errors_2d, color=colors['2d'], linewidth=2, label='2D Error')
    ax1.plot(epochs, errors_3d, color=colors['3d'], linewidth=2, label='3D Error')

ax1.set_xlabel('Epoch', fontsize=12)
ax1.set_ylabel('Positioning Error (m)', fontsize=12)
ax1.set_title('Positioning Error Evolution', fontsize=13, fontweight='bold')
ax1.legend(frameon=True, fancybox=True)
ax1.grid(True, alpha=0.3)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)

# East vs North Error Scatter Plot
scatter = ax2.scatter(errors_e, errors_n, c=epochs, cmap='viridis', alpha=0.7, s=30, edgecolors='white', linewidth=0.5)
ax2.set_xlabel('East Error (m)', fontsize=12)
ax2.set_ylabel('North Error (m)', fontsize=12)
ax2.set_title('East vs North Error Distribution', fontsize=13, fontweight='bold')
ax2.axhline(y=0, color='k', linestyle='--', alpha=0.5, linewidth=1)
ax2.axvline(x=0, color='k', linestyle='--', alpha=0.5, linewidth=1)
ax2.grid(True, alpha=0.3)
ax2.set_aspect('equal', adjustable='box')
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)

# Add colorbar
cbar = plt.colorbar(scatter, ax=ax2)
cbar.set_label('Epoch', fontsize=12)

plt.tight_layout()
plt.show()

# Summary table
print(f"\n{'='*60}")
print(f"{'GNSS SPP POSITIONING RESULTS SUMMARY':^60}")
print(f"{'='*60}")
print(f"Total Epochs Processed: {len(estimated_positions)}")
print(f"{'='*60}")
print(f"{'POSITIONING ACCURACY STATISTICS':^60}")
print(f"{'='*60}")
print(f"{'Metric':<20} {'2D Error (m)':<15} {'3D Error (m)':<15}")
print(f"{'-'*50}")
print(f"{'Mean':<20} {np.mean(errors_2d):<15.3f} {np.mean(errors_3d):<15.3f}")
print(f"{'RMS':<20} {np.sqrt(np.mean(np.array(errors_2d)**2)):<15.3f} {np.sqrt(np.mean(np.array(errors_3d)**2)):<15.3f}")
print(f"{'Standard Deviation':<20} {np.std(errors_2d):<15.3f} {np.std(errors_3d):<15.3f}")
print(f"{'Maximum':<20} {np.max(errors_2d):<15.3f} {np.max(errors_3d):<15.3f}")
print(f"{'Minimum':<20} {np.min(errors_2d):<15.3f} {np.min(errors_3d):<15.3f}")
print(f"{'95th Percentile':<20} {np.percentile(errors_2d, 95):<15.3f} {np.percentile(errors_3d, 95):<15.3f}")
print(f"{'='*60}")

# Print coordinate range information
print(f"{'COORDINATE RANGES':^60}")
print(f"{'='*60}")
print(f"Latitude Range:  {min(lat_est):.6f}° to {max(lat_est):.6f}°")
print(f"Longitude Range: {min(lon_est):.6f}° to {max(lon_est):.6f}°") 
print(f"Altitude Range:  {min(alt_est):.3f} m to {max(alt_est):.3f} m")
print(f"{'='*60}")

print(f"\nResults saved to:")
print(f"- estimated_positions.csv: Complete estimated positions with ECEF, LLA, and ENU coordinates")
print(f"- positioning_errors.csv: Detailed error analysis for each epoch")
print(f"\nEstimated positions CSV contains the following columns:")
print(f"- Epoch, Latitude_deg, Longitude_deg, Altitude_m, Clock_Bias_m")
print(f"- East_m, North_m, Up_m (ENU coordinates)")  
print(f"- ECEF_X_m, ECEF_Y_m, ECEF_Z_m (ECEF coordinates)")
print(f"{'='*60}")

# Print first few rows as example
print(f"\nFirst 5 estimated positions:")
print(results_df.head().to_string(index=False))