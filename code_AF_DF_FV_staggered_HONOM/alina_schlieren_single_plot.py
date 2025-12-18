#####################
#CFL 0.95
#200 elements
#####################
nametest="shock_vortex_1200_601"
filename="SOLUTION_OCT_0000000.7000000.dat"
# filename="SOLUTION_OCT_0000000.6900000.dat"
CFL_test=0.475 #OK
word=None #"twilight_shifted" #"cubehelix" #"viridis" #"magma" #"inferno" #"plasma" #"rainbow" #"hsv" #"nipy_spectral" #"turbo" #"jet" #None

gmm=1.4

test=[-27,3, 12, CFL_test]

name_base_folder="simulations_paper"

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser
import matplotlib.ticker as ticker
from matplotlib.ticker import FuncFormatter


space_reconstruction=test[0]
time_scheme         =test[1]
post_processing     =test[2]
CFL                 =test[3]
foldName=name_base_folder+"/"+nametest+"/space_reconstruction"+str(space_reconstruction)+"/time_scheme_"+str(time_scheme)+"/PP"+str(post_processing)+"/CFL"+str(CFL)

gmm = 1.4

if not os.path.isdir(foldName):
    raise ValueError("Folder not valid. Run as 'python plotting.py --dir folder_name --file file_name'")

def plotting(folder_name, file_name, levels):
    folder_name = folder_name + "/"
    delimiter_in = ' '
    headerlines_in = 1

    # Importing the MESH DATA
    import_mesh = sorted([file for file in os.listdir(folder_name) if file.startswith('MESH_')])
    if len(import_mesh) < 2:
        raise ValueError("Mesh files missing")
    
    x = np.loadtxt(folder_name + import_mesh[0], delimiter=delimiter_in, skiprows=headerlines_in)
    y = np.loadtxt(folder_name + import_mesh[1], delimiter=delimiter_in, skiprows=headerlines_in)
    X, Y = np.meshgrid(x, y)

    # Importing the specified SOLUTION DATA
    if not file_name:
        raise ValueError("No solution file specified")
    if not os.path.isfile(folder_name + file_name):
        raise ValueError("Specified solution file not found")
    
    mydata_solution = np.loadtxt(folder_name + file_name, skiprows=headerlines_in)
    ro = mydata_solution[:, 0]
    rou = mydata_solution[:, 1]
    rov = mydata_solution[:, 2]
    energy = mydata_solution[:, 3]
    
    RO = np.reshape(ro, X.shape)
    ROU = np.reshape(rou, X.shape)
    ROV = np.reshape(rov, X.shape)
    ENERGY = np.reshape(energy, X.shape)

    U = ROU / RO
    V = ROV / RO
    P = (gmm - 1.0) * (ENERGY - 0.5 * RO * (U**2 + V**2))

    VORT = np.zeros_like(RO)
    dx = X[1, 1] - X[1, 0]
    dy = Y[1, 1] - Y[0, 1]
    
    #CHECKED
    print(X[1, 1],X[1, 0])
    print(Y[1, 1],Y[0, 1])
    print(X.shape[0],Y.shape[0])
    print(X.shape[1],Y.shape[1])

    for indy in range(1, X.shape[0] - 1):
        for indx in range(1, X.shape[1] - 1):
            VORT[indy, indx] = np.sqrt(((RO[indy + 1, indx] - RO[indy - 1, indx]) / (2 * dy))**2 + 
                                        ((RO[indy, indx + 1] - RO[indy, indx - 1]) / (2 * dx))**2)
    

    maxgrad = np.max(VORT)
    K = 80
    VORT = np.exp(-K * VORT / maxgrad)

    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))
    if word != None:
        cont_vort = ax.contourf(X, Y, VORT, cmap=word+'_r', levels=levels)    
    else:
        cont_vort = ax.contourf(X, Y, VORT, cmap='Greys_r', levels=levels)
    cbar = fig.colorbar(cont_vort, ax=ax, shrink=0.5)  # Reduce colorbar size
    # cbar.set_label('vort')
    # ax.set_title(r"$schilieren$")
    # ax.set_xlabel('x')
    # ax.set_ylabel('y')
    ax.set_aspect('equal')
    def smart_format(x, _):
        return str(float(x))  # converte tipo 0.50 → '0.5', 1.00 → '1.0' → '1'
    ax.xaxis.set_major_formatter(FuncFormatter(smart_format))
    cbar.formatter = FuncFormatter(lambda x, _: str(float(x)))
    cbar.update_ticks()

    print("VORT min, max =", VORT.min(), VORT.max())
    print("Colorbar limits =", cont_vort.get_clim())

    if word != None:
        plt.savefig(nametest+"_schlieren_K"+str(K)+"_PP"+str(post_processing)+"_timescheme"+str(time_scheme)+"_CFL"+str(CFL)+"_"+filename+"_"+word+".pdf", format="pdf", bbox_inches="tight")
    else:
        plt.savefig(nametest+"_schlieren_K"+str(K)+"_PP"+str(post_processing)+"_timescheme"+str(time_scheme)+"_CFL"+str(CFL)+"_"+filename+".pdf", format="pdf", bbox_inches="tight")
    # plt.show()

if __name__ == '__main__':
    plotting(foldName, filename, 50)