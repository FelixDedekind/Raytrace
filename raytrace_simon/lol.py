import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# Define lens surface function
def lens_surface(x, R, k, aspheric_coefficients):
    # Example: R is the vertex radius, k is the conic constant, aspheric_coefficients are higher-order terms
    z = (x*2 / (R * (1 + np.sqrt(1 - (1 + k) * (x / R*2)))))
    for i, a in enumerate(aspheric_coefficients):
        z += a * (x ** (2 * (i + 1)))
    return z

# Snell's law to calculate refraction
def snells_law(n1, n2, theta1):
    return np.arcsin((n1 / n2) * np.sin(theta1))

# Function to compute the intersection of a ray with the lens surface
def find_intersection(ray, surface_func, guess):
    def func(x):
        ray_x, ray_y = ray
        return ray_y - surface_func(x) - (x - ray_x) * (np.tan(np.arctan2(ray_y, x)))
    return fsolve(func, guess)

# Define parameters for lens
R = 10
k = -1
aspheric_coefficients = [0.01, -0.0005]  # Example coefficients
n1 = 1.0  # Refractive index of air
n2 = 1.5  # Refractive index of lens material

# Generate incident rays
num_rays = 10
rays = np.linspace(-5, 5, num_rays)
incident_rays = [(0, y) for y in rays]

# Plot lens surface
x_vals = np.linspace(-5, 5, 400)
y_vals = lens_surface(x_vals, R, k, aspheric_coefficients)

plt.plot(x_vals, y_vals, label='Lens Surface')

# Trace rays through the lens
focal_points = []
for ray in incident_rays:
    x0, y0 = ray
    intersection = find_intersection(ray, lambda x: lens_surface(x, R, k, aspheric_coefficients), x0 + 1)
    theta1 = np.arctan2(y0, intersection - x0)
    theta2 = snells_law(n1, n2, theta1)
    focal_x = intersection + 5 * np.tan(theta2)
    focal_points.append(focal_x)
    plt.plot([x0, intersection], [y0, lens_surface(intersection, R, k, aspheric_coefficients)], 'r-')

# Plot focal points
for fp in focal_points:
    plt.plot(fp, 0, 'bo')

# Calculate magnification
object_distance = 20  # Example object distance
image_distance = np.mean(focal_points)  # Average focal point distance
magnification = -image_distance / object_distance

print("Magnification:", magnification)

plt.xlabel('x')
plt.ylabel('y')
plt.title('Ray Tracing through Lens')
plt.legend()
plt.grid(True)
plt.axis('equal')
plt.show()