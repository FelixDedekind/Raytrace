import numpy as np
import matplotlib.pyplot as plt

# Load the data from the file
data = np.loadtxt('intensity_grid.txt')

# Separate the data into radius and intensity
radius = data[:, 0]
intensity = data[:, 1]

# Normalize intensity to [0, 1]
intensity = intensity / np.max(intensity)

# Find maximum radius for plot
max_radius = max(radius)

# Create a grid of x, y coordinates
x = np.linspace(-max_radius, max_radius, 1000)
y = np.linspace(-max_radius, max_radius, 1000)
x, y = np.meshgrid(x, y)

# Calculate the radius of each point in the grid
grid_r = np.sqrt(x**2 + y**2)

# Create an empty intensity image
image = np.zeros((1000, 1000))

# For each point in the grid, find the closest radius in the data and use its intensity
for i in range(1000):
    for j in range(1000):
        idx = np.abs(radius - grid_r[i, j]).argmin()
        image[i, j] = intensity[idx]

# Display the image as a heatmap
plt.imshow(image, extent=(-max_radius, max_radius, -max_radius, max_radius), cmap='hot', interpolation='nearest')
plt.colorbar()
plt.show()
