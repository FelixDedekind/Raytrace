import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata



# Read data from file
data = np.loadtxt('evaluation.txt')
print("Data loaded")

cutoff = -1 #int(np.shape(data)[0]/100)
resolution = 1000

# Extract angles and convert to radius
theta = data[0:cutoff, 0]
r = np.sin(theta)

# Extract RGB values and normalize to [0, 1]
rgb = data[0:cutoff, 2:5] / 255.0

# Find maximum radius for plot
maxradius = max(r)
print("Data extracted")


# Create a grid of x, y coordinates
x = np.linspace(-maxradius, maxradius, resolution)
y = np.linspace(-maxradius, maxradius, resolution)
x, y = np.meshgrid(x, y)


# Calculate the radius of each point in the grid
grid_r = np.sqrt(x**2 + y**2)

# Create an empty RGB image
image = np.zeros((resolution, resolution, 3))

print("Arrays created")
print("interpolating...")
# For each point in the grid, find the closest radius in the data and use its color
for i in range(resolution):
    print(i/resolution)
    for j in range(resolution):
        idx = np.abs(r - grid_r[i, j]).argmin()
        image[i, j] = rgb[idx]
print("interpolated")

# Display the image
plt.imshow(image)
plt.show()
