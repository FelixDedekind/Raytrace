import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Read the file
with open('surface_normals.txt', 'r') as f:
    lines = f.readlines()

centers = []
normals = []
for line in lines:
    center_str, normal_str = line.split(',')
    centers.append(list(map(float, center_str.split())))
    normals.append(list(map(float, normal_str.split())))

centers = np.array(centers)
normals = np.array(normals)

# Create a 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Read the file
with open('surface_normals.txt', 'r') as f:
    lines = f.readlines()

centers = []
normals = []
for line in lines:
    center_str, normal_str = line.split(',')
    centers.append(list(map(float, center_str.split())))
    normals.append(list(map(float, normal_str.split())))

centers = np.array(centers)
normals = np.array(normals)

# Create a 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot the normals as vectors at the surface centers
ax.quiver(centers[:, 0], centers[:, 1], centers[:, 2], normals[:, 0], normals[:, 1], normals[:, 2], length=0.0001, normalize=True)

plt.show()
