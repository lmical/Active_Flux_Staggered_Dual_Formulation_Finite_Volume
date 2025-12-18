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
# folder="simulations_paper/sod2d_200_200/space_reconstruction-27/time_scheme_3/PP1/CFL0.25/"
filename = "SOLUTION_OCT_0000000.2500000.dat"

folder_name = folder
aspect_ratio = "quad"

if aspect_ratio == "quad":
    aspect_ratio_dimensions = [1, 1]  # Adjusted for 2D plots
elif aspect_ratio == "rect_X":
    aspect_ratio_dimensions = [14/4, 1]  # Maintain aspect ratio in 2D
else:
    print("Aspect ratio undefined")
    quit()

Gmm = 1.4

print(folder_name)
if not os.path.isdir(folder_name):
    raise ValueError("Folder not valid. Run as 'python plotting.py --dir folder_name --file filename.dat'")

def plotting(folder_name, filename, what_plot='all', n_levels=20):
    """
    Generate 2D contour plots from a single simulation data file.

    Parameters:
    - folder_name: str, directory containing the data files.
    - filename: str, name of the data file to plot.
    - what_plot: str, type of plots to generate ('density', 'primitive', 'conservative', 'all').
    - n_levels: int or list, number of contour levels or specific level boundaries.
    """
    global contours, colorbars
    contours = {}  # Define contours as a global variable
    colorbars = {}  # Define colorbars as a global variable

    delimiter_in = ' '
    headerlines_in = 1

    # Importing the MESH DATA
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

    # Importing the SOLUTION DATA
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

    print(np.max(RO), np.min(RO))

    # Define data for plots based on 'what_plot'
    plot_data = {}
    plot_data['ro'] = RO
    if what_plot in ['primitive', 'all']:
        plot_data['u'] = U
        plot_data['v'] = V
        plot_data['p'] = P
    if what_plot in ['conservative', 'all']:
        plot_data['qx'] = RO * U
        plot_data['qy'] = RO * V
        plot_data['E'] = P / (Gmm - 1.) + 0.5 * RO * (U**2 + V**2)

    # Create the figure
    fig = plt.figure(1, figsize=(20, 8))

    # Define subplots based on the 'what_plot' parameter without 3D projection
    axes = {}
    if what_plot == 'density':
        axes['ro'] = fig.add_subplot(111)
    elif what_plot == 'primitive':
        axes['ro'] = fig.add_subplot(141)
        axes['u'] = fig.add_subplot(142)
        axes['v'] = fig.add_subplot(143)
        axes['p'] = fig.add_subplot(144)
    elif what_plot == 'conservative':
        axes['ro'] = fig.add_subplot(141)
        axes['qx'] = fig.add_subplot(142)
        axes['qy'] = fig.add_subplot(143)
        axes['E'] = fig.add_subplot(144)
    elif what_plot == 'all':
        axes['ro'] = fig.add_subplot(241)
        axes['u'] = fig.add_subplot(242)
        axes['v'] = fig.add_subplot(243)
        axes['p'] = fig.add_subplot(244)
        axes['qx'] = fig.add_subplot(245)
        axes['qy'] = fig.add_subplot(246)
        axes['E'] = fig.add_subplot(247)
    else:
        print("Plotting option not available")
        print("It was")
        print(what_plot)
        quit()

    def plot_data_step():
        global contours, colorbars

        # Clear old contours and colorbars
        for key in axes:
            if key in contours:
                for c in contours[key].collections:
                    c.remove()
                contours[key] = None
            if key in colorbars:
                colorbars[key].remove()
                colorbars[key] = None

        # Plot new contours and colorbars
        for key, ax in axes.items():
            if key in plot_data:
                data = plot_data[key]
                vmin, vmax = np.min(data), np.max(data)
                levels = np.linspace(vmin, vmax, n_levels)  # Create intermediate levels
                contours[key] = ax.contourf(X, Y, data, levels=levels, cmap=which_tone, vmin=vmin, vmax=vmax)
                ax.contour(X, Y, data, levels=levels, colors='black', linewidths=0.5)
                
                if key in colorbars:
                    colorbars[key].update_normal(contours[key])
                else:
                    colorbars[key] = fig.colorbar(contours[key], ax=ax)
                
                # Set colorbar ticks and labels
                colorbar = colorbars[key]
                ticks = np.linspace(vmin, vmax, n_levels)  # Create intermediate levels
                colorbar.set_ticks(np.append(ticks, vmax))  # Include the maximum value
                colorbar.set_ticklabels([f'{tick:.3f}' for tick in np.append(ticks, vmax)])  # Format tick labels
                # colorbar.set_label(f'{key} Value')                

                # Set plot titles and labels
                title_dict = {
                    'ro': r"$\rho$",
                    'u': 'u',
                    'v': 'v',
                    'p': 'p',
                    'qx': r'$q_x$',
                    'qy': r'$q_y$',
                    'E': r'$E$'
                }
                ax.set_title(title_dict.get(key, key))
                ax.set_xlabel('x')
                ax.set_ylabel('y')
                ax.set_aspect(aspect_ratio_dimensions[0] / aspect_ratio_dimensions[1])

        plt.draw()

    plot_data_step()

    plt.tight_layout()
    plt.savefig("sod2d_U_contour_density.pdf", format="pdf", bbox_inches="tight")
    # plt.show()

if __name__ == '__main__':
    # You can change the second argument to 'density', 'primitive', 'conservative', or 'all'
    # You can also specify the number of levels, e.g., n_levels=30
    plotting(folder_name, filename, 'density', n_levels=16) #48 #70
    # plotting(folder_name, filename, 'primitive', n_levels=20)
    # plotting(folder_name, filename, 'conservative', n_levels=20)
    # plotting(folder_name, filename, 'all', n_levels=20)
