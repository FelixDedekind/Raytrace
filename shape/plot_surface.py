import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

# Read data from file
with open('surfaces.txt', 'r') as f:
    lines = f.readlines()

# Initialize a list to hold all configurations
configs = []
x_coords = []
y_coords = []
z_coords = []

# Parse the data
for line in lines:
    if line.strip():
        v1, v2, v3 = line.split(',')
        x1, y1, z1 = map(float, v1.split())
        x2, y2, z2 = map(float, v2.split())
        x3, y3, z3 = map(float, v3.split())
        x_coords.extend([x1, x2, x3])
        y_coords.extend([y1, y2, y3])
        z_coords.extend([z1, z2, z3])
    else:
        configs.append((x_coords, y_coords, z_coords))
        x_coords = []
        y_coords = []
        z_coords = []

# Calculate the limits
x_all = [x for config in configs for x in config[0]]
y_all = [y for config in configs for y in config[1]]
z_all = [z for config in configs for z in config[2]]
max_range = max(max(x_all) - min(x_all), max(y_all) - min(y_all), max(z_all) - min(z_all))
xlim = [min(x_all) - 0.1 * max_range, max(x_all) + 0.1 * max_range]
ylim = [min(y_all) - 0.1 * max_range, max(y_all) + 0.1 * max_range]
zlim = [min(z_all) - 0.1 * max_range, max(z_all) + 0.1 * max_range]

# Create a 3D scatter plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Set labels for the axes
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# Animation update function
def update(num, configs, scatter):
    ax.clear()
    x_coords, y_coords, z_coords = configs[num]
    scatter = ax.scatter(x_coords, y_coords, z_coords)

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_zlim(zlim)

    return scatter,

# Create the animation
ani = animation.FuncAnimation(fig, update, frames=len(configs), fargs=(configs, ax), interval=20, repeat = True)

plt.show()
