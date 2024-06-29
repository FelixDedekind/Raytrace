import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Read the data from the file
data = []
with open('beams.txt', 'r') as file:
    for line in file:
        x, y, z, vx, vy, vz = map(float, line.split())
        data.append((x, y, z, vx, vy, vz))

# Separate the data into different lists for plotting
x, y, z, vx, vy, vz = zip(*data)

# Create a 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot the velocities as vectors at the corresponding points
ax.quiver(x, y, z, vx, vy, vz, color='b')

xylim = 0.001
zlim = 1 # around zavg
zavg = 0
ax.set_xlim([-xylim, xylim])
ax.set_ylim([-xylim, xylim])
ax.set_zlim([zavg - zlim, zavg + zlim])

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()