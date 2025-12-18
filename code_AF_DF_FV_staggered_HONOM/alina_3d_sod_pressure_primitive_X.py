import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
from optparse import OptionParser

# parser = OptionParser()
# parser.add_option("-d", "--dir", action="store", type="string", dest="d", default="tsunami")
# options, args = parser.parse_args()
# folder = options.d
    
folder="simulations_paper/sod2d/space_reconstruction-27/time_scheme_3/PP1/CFL0.25/"
folder_name = folder

# Define file paths
mesh_x_file = os.path.join(folder, "W_X_MESH_CoordinateX.dat")
mesh_y_file = os.path.join(folder, "W_X_MESH_CoordinateY.dat")
solution_file = os.path.join(folder, "W_X_SOLUTION_OCT_0000000.2500000.dat")  # Specific time step file

# Check if files exist
assert os.path.exists(mesh_x_file), f"File not found: {mesh_x_file}"
assert os.path.exists(mesh_y_file), f"File not found: {mesh_y_file}"
assert os.path.exists(solution_file), f"File not found: {solution_file}"

# Load mesh data
x = np.loadtxt(mesh_x_file, delimiter=' ', skiprows=1)
y = np.loadtxt(mesh_y_file, delimiter=' ', skiprows=1)
X, Y = np.meshgrid(x, y)

# Load solution data
solution_data = np.loadtxt(solution_file, skiprows=1)
ro = solution_data[:, 0]
u = solution_data[:, 1]
v = solution_data[:, 2]
p = solution_data[:, 3]

# Compute derived quantities
Gmm = 1.4

RO  = np.reshape(ro , (X.shape[0], Y.shape[1]))
U   = np.reshape(u  , (X.shape[0], Y.shape[1]))
V   = np.reshape(v  , (X.shape[0], Y.shape[1]))
P   = np.reshape(p  , (X.shape[0], Y.shape[1]))

ROU=RO*U
ROV=RO*V

ENERGY=P/(Gmm-1)+0.5*RO*(U**2+V**2)

aspect_ratio_dimensions=[4, 4, 4]

# Plotting
fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, P, color='red', alpha=1.0, edgecolor='none', linewidth=0.5)
ax.set_title(r"$p$")
ax.set_xlabel('X')
ax.set_ylabel('Y')
# ax.set_zlabel(r"Density")
ax.set_box_aspect(aspect_ratio_dimensions)  # Set aspect ratio
ax.view_init(elev=30, azim=225)
ax.grid(False)  # Disable grid lines
# Reduce the number of ticks on x and y axes
ax.set_xticks([-1, -0.5, 0, 0.5, 1])  # Manually specify x-ticks
ax.set_yticks([-1, -0.5, 0, 0.5, 1])  # Manually specify y-ticks

plt.savefig("sod2d_V_X_3d_pressure_without_grids.pdf", format="pdf", bbox_inches="tight", transparent=True)
# plt.show()
