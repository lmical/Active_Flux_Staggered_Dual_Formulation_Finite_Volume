import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser

# =========================
# User input / options
# =========================
parser = OptionParser()
parser.add_option("-d", "--dir", action="store", type="string", dest="folder", default="tsunami")
parser.add_option("-o", "--outdir", action="store", type="string", dest="outfolder",
                  default="final_results_multifld_CFL0.475")
options, args = parser.parse_args()
folder = options.folder
outfolder = options.outfolder

aspect_ratio = "quad"
if aspect_ratio == "quad":
    aspect_ratio_dims = [1, 1]
elif aspect_ratio == "rect_X":
    aspect_ratio_dims = [14, 4]
else:
    print("Aspect ratio undefined")
    sys.exit(1)

gmm = 1.4
K = 80       # exponential scaling
levels = 30  # number of contour levels
cmap = 'Greys_r'
colorbar_shrink = 0.5

# =========================
# Check folders
# =========================
if not os.path.isdir(folder):
    raise ValueError(f"Input folder '{folder}' not valid.")
if not os.path.isdir(outfolder):
    os.makedirs(outfolder)

# =========================
# Load mesh
# =========================
mesh_files = sorted([f for f in os.listdir(folder) if f.startswith('W_X_MESH_')])
if len(mesh_files) < 2:
    raise ValueError("Mesh files missing.")

x = np.loadtxt(os.path.join(folder, mesh_files[0]), skiprows=1)
y = np.loadtxt(os.path.join(folder, mesh_files[1]), skiprows=1)
X, Y = np.meshgrid(x, y)
dx = X[0, 1] - X[0, 0]
dy = Y[1, 0] - Y[0, 0]

# Domain limits exactly at the mesh boundaries
x_lim = [x.min(), x.max()]
y_lim = [y.min(), y.max()]

# =========================
# Load solution files
# =========================
solution_files = sorted([f for f in os.listdir(folder) if f.startswith('W_X_SOLUTION_')])
if len(solution_files) == 0:
    raise ValueError("No solution files found.")

# =========================
# Loop over solution files
# =========================
for i, filename in enumerate(solution_files):
    print(f"Processing {filename} ({i+1}/{len(solution_files)})")

    data = np.loadtxt(os.path.join(folder, filename), skiprows=1)
    ro, rou, rov, energy = data[:, 0], data[:, 1], data[:, 2], data[:, 3]

    # Reshape to mesh
    RO = np.reshape(ro, X.shape)
    ROU = np.reshape(rou, X.shape)
    ROV = np.reshape(rov, X.shape)
    ENERGY = np.reshape(energy, X.shape)

    # Vectorized gradient magnitude (VORT)
    dRO_dy, dRO_dx = np.gradient(RO, dy, dx)
    VORT = np.sqrt(dRO_dx**2 + dRO_dy**2)
    VORT = np.exp(-K * VORT / np.max(VORT))

    # Plot
    fig, ax = plt.subplots(figsize=(12, 6))
    cont = ax.contourf(X, Y, VORT, cmap=cmap, levels=levels)

    # Colorbar next to plot
    cbar = fig.colorbar(cont, ax=ax, shrink=colorbar_shrink, pad=0.02)
    cbar.set_label('')  # remove label

    # Set axis limits exactly to mesh
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)
    ax.set_aspect(aspect_ratio_dims[0] / aspect_ratio_dims[1])

    # Show axis ticks
    ax.tick_params(direction='in', top=True, right=True)
    # ax.set_xlabel('X')
    # ax.set_ylabel('Y')

    # Save PNG
    save_path = os.path.join(outfolder, f"{filename}.png")
    plt.savefig(save_path, format='png', dpi=600, bbox_inches='tight')
    plt.close(fig)

print(f"All images saved in folder '{outfolder}'.")
