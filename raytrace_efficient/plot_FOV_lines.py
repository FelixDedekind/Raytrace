import numpy as np
import matplotlib.pyplot as plt

# Read data from file
data = np.loadtxt('evaluation.txt')

# Ignore the second column
theta = data[:, 0]
r = np.sin(theta)
rgb = data[:, 2:5] / 255.0

# Set the resolution and create an x array from -r to r
resolution = 1000
x = np.linspace(-max(r), max(r), resolution)

# Create an array for RGB values mirrored around zero
rgb_mirrored = np.zeros((resolution, 3))

# Fill the RGB array
half_res = resolution // 2
for i in range(half_res):
    idx = np.abs(r - x[half_res + i]).argmin()
    rgb_mirrored[half_res + i] = rgb[idx]
    rgb_mirrored[half_res - i - 1] = rgb[idx]

# Extend the RGB values to a second dimension
y_resolution = 1000
image = np.zeros((y_resolution, resolution, 3))
for i in range(y_resolution):
    image[i] = rgb_mirrored

# Display the image
plt.imshow(image, extent=[-max(r), max(r), -max(r), max(r)], aspect='auto')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('RGB Image')
plt.show()

# Overlay a plot of the red value over the x coordinate
plt.plot(x, rgb_mirrored[:, 0], color='red')
plt.xlabel('X')
plt.ylabel('Red Value')
plt.title('Red Value Over X Coordinate')
plt.show()
