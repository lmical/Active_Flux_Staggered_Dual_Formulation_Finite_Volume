import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser

which_tone = "gist_rainbow"

# parser = OptionParser()
# parser.add_option("-d", "--dir", action="store", type="string", dest="d", default="tsunami")
# options, args = parser.parse_args()
# folder = options.d
    
folder="simulations_paper/sod2d/space_reconstruction-27/time_scheme_3/PP1/CFL0.25/"
filename = "SOLUTION_OCT_0000000.2500000.dat"



folder_name = folder
aspect_ratio = "quad"

if aspect_ratio == "quad":
    aspect_ratio_dimensions = [1, 1]
elif aspect_ratio == "rect_X":
    aspect_ratio_dimensions = [14/4, 1]
else:
    print("Aspect ratio undefined")
    quit()

Gmm = 1.4

print(folder_name)
if not os.path.isdir(folder_name):
    raise ValueError("Folder not valid. Run as 'python plotting.py --dir folder_name'")

def plotting(folder_name, filename, what_plot='E', n_levels=20):
    global contours, colorbars
    contours = {}
    colorbars = {}

    delimiter_in = ' '
    headerlines_in = 1

    import_mesh = [file for file in os.listdir(folder_name) if file.startswith('MESH_')]
    import_mesh.sort()
    nfiles_mesh = len(import_mesh)
    for i in range(nfiles_mesh):
        if i == 0:
            mesh_filename = import_mesh[i]
            mydata_mesh = np.loadtxt(folder_name + mesh_filename, delimiter=delimiter_in, skiprows=headerlines_in)
            x = mydata_mesh
        if i == 1:
            mesh_filename = import_mesh[i]
            mydata_mesh = np.loadtxt(folder_name + mesh_filename, delimiter=delimiter_in, skiprows=headerlines_in)
            y = mydata_mesh

    X, Y = np.meshgrid(x, y)

    mydata_solution = np.loadtxt(folder_name + filename, skiprows=headerlines_in)
    ro = mydata_solution[:, 0]
    rou = mydata_solution[:, 1]
    rov = mydata_solution[:, 2]
    Energy = mydata_solution[:, 3]
    phi = mydata_solution[:, 4]

    RO = np.reshape(ro, (X.shape[0], Y.shape[1]))
    ROU = np.reshape(rou, (X.shape[0], Y.shape[1]))
    ROV = np.reshape(rov, (X.shape[0], Y.shape[1]))
    ENERGY = np.reshape(Energy, (X.shape[0], Y.shape[1]))
    U=ROU/RO
    V=ROV/RO
    P = (Gmm-1)*(ENERGY-0.5*RO*(U**2+V**2))

    PHI = np.reshape(phi, (X.shape[0], Y.shape[1]))

    plot_data = {}
    plot_data['E'] = P / (Gmm - 1.) + 0.5 * RO * (U**2 + V**2)

    fig = plt.figure(1, figsize=(8, 6))
    axes = {}
    axes['E'] = fig.add_subplot(111)

    def plot_data_step():
        global contours, colorbars

        for key in axes:
            if key in contours:
                for c in contours[key].collections:
                    c.remove()
                contours[key] = None
            if key in colorbars:
                colorbars[key].remove()
                colorbars[key] = None

        for key, ax in axes.items():
            if key in plot_data:
                data = plot_data[key]
                vmin, vmax = np.min(data), np.max(data)
                levels = np.linspace(vmin, vmax, n_levels)
                contours[key] = ax.contourf(X, Y, data, levels=levels, cmap=which_tone, vmin=vmin, vmax=vmax)
                ax.contour(X, Y, data, levels=levels, colors='black', linewidths=0.5)
                
                if key in colorbars:
                    colorbars[key].update_normal(contours[key])
                else:
                    colorbars[key] = fig.colorbar(contours[key], ax=ax)
                
                colorbar = colorbars[key]
                ticks = np.linspace(vmin, vmax, n_levels)
                colorbar.set_ticks(np.append(ticks, vmax))
                colorbar.set_ticklabels([f'{tick:.3f}' for tick in np.append(ticks, vmax)])
                
                ax.set_title(r'$E$')
                ax.set_xlabel('x')
                ax.set_ylabel('y')
                ax.set_aspect(aspect_ratio_dimensions[0] / aspect_ratio_dimensions[1])

        plt.draw()

    plot_data_step()

    plt.tight_layout()
    plt.savefig("sod2d_U_contour_energy.pdf", format="pdf", bbox_inches="tight")
    # plt.show()

if __name__ == '__main__':
    plotting(folder_name, filename, 'E', n_levels=16) #52 #70
