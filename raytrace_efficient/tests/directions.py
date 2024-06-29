import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

maxt = np.pi/4
maxp = 2*np.pi

Nt = 4
Np = int(maxp / maxt * Nt)

dt = maxt/Nt
dp = maxp/Np

# Create a new figure and add a 3D subplot to it
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')



for ii in range(Nt):
    for jj in range(Np):
        if ii == 0 and jj > 0: 
            break

        t = ii*dt
        p = jj*dp

        x = np.sin(t)*np.cos(p)
        y = np.sin(t)*np.sin(p)
        z = np.cos(t)

        # Plot the point (x, y, z)
        ax.scatter(x, y, z, color='b')

# Set labels for the axes
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# Display the plot
plt.show()
