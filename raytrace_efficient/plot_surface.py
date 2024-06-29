import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# Read the file
df = pd.read_csv('surfaces.txt', header=None)

# Parse the vertices
vertices = []
with open('surfaces.txt', 'r') as f:
    for line in f:
        triangle = [list(map(float, vertex.split())) for vertex in line.split(' , ')]
        if len(triangle) == 3:  # Only add complete triangles
            vertices.append(triangle)



# Create a 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Create a Poly3DCollection object
poly = Poly3DCollection(vertices, alpha=.25, linewidths=1, edgecolors='r', facecolors='blue')  # Set facecolors to 'white'
poly.set_alpha(0.1)

# Add the collection to the plot
ax.add_collection3d(poly)

# # Auto scale to the mesh size
# scale = np.concatenate([np.ravel(v) for v in vertices])  # concatenate all vertices
# ax.auto_scale_xyz(scale, scale, scale)

# Calculate the min and max for each axis
# x_min = np.min([v[0] for v in np.concatenate(vertices)])
# x_max = np.max([v[0] for v in np.concatenate(vertices)])
# y_min = np.min([v[1] for v in np.concatenate(vertices)])
# y_max = np.max([v[1] for v in np.concatenate(vertices)])
# z_min = np.min([v[2] for v in np.concatenate(vertices)])
# z_max = np.max([v[2] for v in np.concatenate(vertices)])

# Set the limits
# xylim = 0.001
# zlim = 0.005 # around zavg
# zavg = 0.02
# ax.set_xlim([-xylim, xylim])
# ax.set_ylim([-xylim, xylim])
# ax.set_zlim([zavg - zlim, zavg + zlim])

# Read the rays file
plot_every_nth_ray = 100
rays = []
with open('rays.txt', 'r') as f:
    for i, line in enumerate(f):
        if i % plot_every_nth_ray == 0:  # Only add every nth point
            ray = [list(map(float, point.split())) for point in line.split(',') if point.strip()]
            if len(ray) > 1 and all(len(point) == 3 for point in ray):  # Only add complete rays
                rays.append(ray)

# Plot lines for each ray
for ray in rays:
    X = [point[0] for point in ray]
    Y = [point[1] for point in ray]
    Z = [point[2] for point in ray]
    ax.plot(X, Y, Z, color='r')

plt.show()



# plot points at rays
# X = [ray[0] for ray in rays]
# Y = [ray[1] for ray in rays]
# Z = [ray[2] for ray in rays]
# ax.scatter(X, Y, Z, color='green', s=20)


