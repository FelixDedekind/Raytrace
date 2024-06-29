import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np


# Create a 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

def plot_surface():

    # Parse the vertices
    vertices = []
    with open('surfaces.txt', 'r') as f:
        for line in f:
            triangle = [list(map(float, vertex.split())) for vertex in line.split(' , ')]
            if len(triangle) == 3:  # Only add complete triangles
                vertices.append(triangle)

    # Create a Poly3DCollection object
    poly = Poly3DCollection(vertices, alpha=.25, linewidths=1, edgecolors='r', facecolors='blue')  # Set facecolors to 'white'
    poly.set_alpha(0.1)

    # Add the collection to the plot
    ax.add_collection3d(poly)


def plot_every_nth_beam(n):
    lines = []
    counter = 0
    with open('rays.txt', 'r') as file:
        for line in file:
            if counter%n == 0:
                vertices = line.split(',')
                # Ignore the last element if it's an empty string
                vertices = vertices[:-1]
                line_data = [tuple(map(float, vertex.split())) for vertex in vertices]
                lines.append(line_data)
            counter+=1

    # Plot the lines between subsequent vertices in each line
    for line_data in lines:
        x, y, z = zip(*line_data)
        ax.plot(x, y, z)




plot_every_nth_beam(11711)
plot_surface()



xylim = 0.005
zlim = 0.005 # around zavg
zavg = 0.05
ax.set_xlim([-xylim, xylim])
ax.set_ylim([-xylim, xylim])
ax.set_zlim([zavg - zlim, zavg + zlim])


ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()
