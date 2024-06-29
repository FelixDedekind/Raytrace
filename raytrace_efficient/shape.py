# (x-0.0308269482009)^(2)+(y+3.201303162077)^(2)=20.6242818446884
import numpy as np
import matplotlib.pyplot as plt

centre = [0.0,0.0]
radius = 0.005

real_radius = 0.005
nx = 150
effective_radius = radius

scale_factor = real_radius/radius
radius *= scale_factor
centre = [cc*scale_factor for cc in centre]



def circ(centre, radius, x):
    y = np.sqrt(radius**2-(x-centre[0])**2)+centre[1]
    return y

print(effective_radius)

xspace = np.linspace(0,0.0005,nx)
yspace = circ(centre, effective_radius, xspace)

plt.plot(xspace,yspace)
plt.show()

yspace_zero = np.zeros_like(yspace)  # y is always zero
data = np.column_stack((xspace, yspace_zero, yspace))

# Write the data to the file
with open("vertices.txt", "w") as f:
    f.write(f"{nx} {nx}\n")  # Write nx repeated twice with a space in between
    np.savetxt(f, data, fmt='%.20f')  # Write the data with 20 decimal places