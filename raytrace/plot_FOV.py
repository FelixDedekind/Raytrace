import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

# Read data from file
data = np.loadtxt('evaluation.txt')

# Extract columns
arcsin_r = data[:, 0]
phi = data[:, 1]
rgb = data[:, 2:5]

# Convert RGB to normalized 0-1 range
rgb_norm = rgb / 255.0

# Convert arcsin_r back to r
r = np.sin(arcsin_r)

# Downsample the data
downsample_factor = 1  # Adjust this value to fit your needs - IF IT'S TAKING TOO LONG, SET TO LIKE 1000
r = r[::downsample_factor]
phi = phi[::downsample_factor]
rgb_norm = rgb_norm[::downsample_factor]

# Create a 2D grid
r_grid, phi_grid = np.meshgrid(np.unique(r), np.unique(phi))

# Interpolate RGB values onto the grid
color_grid = np.empty(r_grid.shape + (3,))
for i in range(3):
    color_grid[..., i] = griddata((r, phi), rgb_norm[:, i], (r_grid, phi_grid), method='linear')

# Create a polar plot
ax = plt.subplot(projection="polar")

# Plot the data
plt.pcolormesh(phi_grid, r_grid, color_grid, shading='auto')

# Remove ticks and labels
ax.set_yticklabels([])
ax.set_xticklabels([])

# Show the plot
plt.savefig("./plot.png", dpi=400)
