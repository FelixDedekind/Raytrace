import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Read the points from the file
points = np.loadtxt('../directions.txt') 

# Create a new figure and add a 3D subplot to it
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Set the aspect ratio to 1:1:1
ax.set_box_aspect([1,1,1])

# Plot the points
for x, y, z in points:
    ax.scatter(x, y, z, color='b')

# Set labels for the axes
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# Display the plot
plt.show()
