coordinate="X" #"X" #"Y" #NB:Actually y is not implemented
lw=2. #2. for shock interaction
ms=3.0 #1. 
fs=10
#overwritten later if needed


# nametest="sod"
# file="WS_SOLUTION_OCT_0000000.2000000.dat"
# CFL_test=0.475 #0.8 #original #0.95 #numerical fluxes
# solution_type="exact solution"
# xspan=[0,1]



# nametest="shock_turbulence_interaction_shu_osher"
# file="WS_SOLUTION_OCT_0000001.8000000.dat"
# CFL_test=0.475 #OK
# lw=2. #for shock interaction
# ms=2. 
# solution_type="reference solution"
# xspan=[-5,5]

# nametest="woodward_colella"
# file="WS_SOLUTION_OCT_0000000.0380000.dat"
# CFL_test=0.475
# solution_type="reference solution"
# xspan=[0,1]


# nametest="woodward_colella_theta1.1"
# file="WS_SOLUTION_OCT_0000000.0380000.dat"
# CFL_test=0.25
# solution_type="reference solution"
# xspan=[0,1]


nametest="RP2"
file="WS_SOLUTION_OCT_0000000.1500000.dat"
CFL_test=0.475
solution_type="exact solution"
xspan=[0,1]

# nametest="RP2_200"
# file="WS_SOLUTION_OCT_0000000.1500000.dat"
# CFL_test=0.25
# solution_type="exact solution"
# xspan=[0,1]

# nametest="lax"
# file="WS_SOLUTION_OCT_0000001.3000000.dat"
# CFL_test=0.25 #OK
# solution_type="exact solution"
# xspan=[0,10]




if not(nametest.startswith("shock_turbulence_interaction")):
    lwex=lw+1.5
else:
    lwex=lw


gmm=1.4


teststocompare=[] 
#Reconstructed variable, riemann_solver, speed_estimate, order, CFL, RelaxedCFL, NRelaxedTimeSteps


teststocompare=[] 
#Reconstructed variable, riemann_solver, speed_estimate, order, CFL, RelaxedCFL, NRelaxedTimeSteps
# teststocompare.append([-25,12, 1, CFL_test]) 
# teststocompare.append([-24,12, 1, CFL_test]) 
# teststocompare.append([-23,12, 1, CFL_test]) 
# teststocompare.append([-22,12, 1, CFL_test]) 
# teststocompare.append([-21,12, 1, CFL_test]) 
# teststocompare.append([-20,12, 1, CFL_test]) 
# teststocompare.append([-1,  12, 1, CFL_test]) 
# teststocompare.append([25,  12, 1, CFL_test]) 
# teststocompare.append([24,  12, 1, CFL_test]) 
# teststocompare.append([23,  12, 1, CFL_test]) 
# teststocompare.append([22,  12, 1, CFL_test]) 
# teststocompare.append([21,  12, 1, CFL_test]) 
# teststocompare.append([20,  12, 1, CFL_test]) 
# teststocompare.append([1,   12, 1, CFL_test]) 

teststocompare.append([-27,   3, -1, CFL_test]) 


import os
from glob import glob
import numpy as np
import matplotlib.pyplot as pl

pl.rcParams['xtick.labelsize'] = 4
pl.rcParams['ytick.labelsize'] = 4


linestyles_space_reconstruction = {1: "--",
                                   10:"--",
                                   2: "--",
                                   20:"--",
                                   21:"--",
                                   22:"--",
                                   23:"--",
                                   24:"--",
                                   25:"--",
                                  -1: "-",
                                  -10:"-",
                                  -2: "-",
                                  -20:"-",
                                  -21:"-",
                                  -22:"-",
                                  -23:"-",
                                  -24:"-",
                                  -25:"-",
                                  -26:"-",
                                  -27:"-"}

markers_space_reconstruction   = { 1: "x",
                                   10:"x",
                                   2: "x",
                                   20:"x",
                                   21:"x",
                                   22:"x",
                                   23:"x",
                                   24:"x",
                                   25:"x",
                                  -1: "o",
                                  -10:"o",
                                  -2: "o",
                                  -20:"o",
                                  -21:"o",
                                  -22:"o",
                                  -23:"o",
                                  -24:"o",
                                  -25:"o",
                                  -26:"o",
                                  -27:"o"}


#Updated col
colors = {
     1:"#17becf",  #PWC
     2:"#d62728",  #MINMOD
     20:"#d62728", #MINMOD
     21:"#1f77b4", #
     22:"#ff7f0e", #
     23:"#2ca02c", #
     24:"#9467bd",  #

    -1:"#17becf",  #PWC
    -2:"#d62728",  #MINMOD
    -20:"#d62728", #MINMOD
    -21:"#1f77b4", #
    -22:"#ff7f0e", #
    -23:"#2ca02c", #
    -24:"#9467bd",  #
    -26:"red", #AAAAAAAAAAAAAAAAAAAAAA
    -27:"red", #AAAAAAAAAAAAAAAAAAAAAA

     93:"#e377c2", #

     6:"#8c564b",  #
    -2:"#bcbd22",  #

    -6:"#ff9896",  #
     7:"#7f7f7f"   #
}

# colors = {
#     3:  "#1f77b4",   # 
#     5:  "#ff7f0e",   # 
#     7:  "#2ca02c",   # 
#     9:  "#d62728",   # 
#     11: "#9467bd",   # 
#     13: "#8c564b",   # 
#     15: "#17becf",   # 
#     17: "#bcbd22",   # 
#     19: "#e377c2",   # 
#     21: "#7f7f7f",   # 
#     23: "#ff9896",   # 
#     25: "#98df8a",   # New
#     27: "#c5b0d5",   # New
#     29: "#c49c94",   # New
#     31: "#f7b6d2",   # New
#     33: "#aec7e8",   # New
#     35: "#ffbb78",   # New
#     37: "#9edae5",   # New
#     39: "#dbdb8d",   # New
#     41: "#c7c7c7",   # New
#     43: "#bc80bd",   # New
#     45: "#6b6ecf",   # New
#     47: "#b5cf6b",   # New
#     49: "#e7969c",   # New
#     51: "#d6616b",   # New
#     53: "#843c39",   # New
#     55: "#8c6d31",   # New
#     57: "#d9d9d9",   # New
#     59: "#7b4173",   # New
#     61: "#b15928",   # New
#     63: "#6baed6"    # New
# }

####################################################
name_reconstructed_variable = {0:"cons",    1:"char"}
name_riemann_solver         = {-2:"central", -1:"LF", 0:"rusanov", 1:"exact", 2:"HLL", 3:"CU", 4:"LDCU", 5:"HLLC", 6:"FORCE",7:"FORCEalpha"}
name_speed_estimate         = {0:"standard",1:"riemann"}
order_of_plotting           = {1:9,2:8,3:7,5:6,7:5,9:4,11:3,13:2}
####################################################

####################################################
label_reconstructed_variable = {0:"cons.",    1:"char."}
label_riemann_solver         = {-2:"Central", -1:"LxF", 0:"Rus", 1:"Ex.RS", 2:"HLL", 3:"CU", 4:"LDCU", 5:"HLLC", 6:"FORCE",7:r"FORCE-$\alpha$"}
label_space_reconstruction   =   { 1: "FV-o1"          ,
                                   10:"FV-o1-PWL in PC",
                                   2: "FV-o2-MUSCL"    ,
                                   20:"FV-o2-MUSCL"    ,
                                   21:"FV-o2-k-MUSCL"  ,
                                   22:"FV-o2-CO-MUSCL" ,
                                   23:"FV-o2-VL-MUSCL" ,
                                   24:"FV-o2-M-MUSCL"  ,
                                   25:"FV-o2-VA-MUSCL" ,
                                  -1: "AF-o1"          ,
                                  -10:"AF-o1-PWL in PC",
                                  -2: "AF-o2-MUSCL"    ,
                                  -20:"AF-o2-MUSCL"    ,
                                  -21:"AF-o2-k-MUSCL"  ,
                                  -22:"AF-o2-CO-MUSCL" ,
                                  -23:"AF-o2-VL-MUSCL" ,
                                  -24:"AF-o2-M-MUSCL"  ,
                                  -25:"AF-o2-VA-MUSCL" ,
                                  -26:"AF-o2-SBM",
                                  -27:"DF-FV" }
label_time_scheme            =   { 1: "Euler"          ,
                                   2: "SSPRK2"         ,
                                   3: "SSPRK3"         ,
                                   4: "SSPRK4"         ,
                                   5: "RK65"           ,
                                   12:"SSPRK2"         , #NB:DeC2 is SSPRK2
                                   13:"DeC3"           ,
                                   14:"DeC4"           ,
                                   15:"DeC5" }
####################################################

nVar=3


params = {'mathtext.default': 'regular' } #parameters plot




fig = pl.figure(1, figsize=(20,5))
ax_ro = fig.add_subplot(131)
ax_ro.set_title( r'$\rho$',fontsize=fs+5)
# ax_ro.set_xlabel('x',fontsize=fs)
# ax_ro.set_ylabel('density',fontsize=fs)
ax_v  = fig.add_subplot(132)
ax_v.set_title( 'u',fontsize=fs+5)
# ax_v.set_xlabel('x',fontsize=fs)
# ax_v.set_ylabel('velocity',fontsize=fs)
ax_p  = fig.add_subplot(133)
ax_p.set_title( 'p',fontsize=fs+5)
# ax_p.set_xlabel('x',fontsize=fs)
# ax_p.set_ylabel('pressure',fontsize=fs)
# ax_q  = fig.add_subplot(234)
# ax_q.set_title( 'momentum')
# ax_q.set_xlabel('x')
# ax_q.set_ylabel('momentum')
# ax_E  = fig.add_subplot(235)
# ax_E.set_title( 'total energy')
# ax_E.set_xlabel('x')
# ax_E.set_ylabel('total energy')

firsttime=True
index_plot=0

ax_ro.tick_params(axis='both', which='major', labelsize=15)
ax_v.tick_params(axis='both', which='major',  labelsize=15)
ax_p.tick_params(axis='both', which='major',  labelsize=15)



##############################
# READ AND PLOT EXACT SOLUTION
##############################
# FILLING x_ref AND ref_sol
##############################
if (nametest.startswith("composite_wave")):
    solfoldName="reference_solutions/composite_wave"
elif (nametest.startswith("woodward_colella")):
    solfoldName="reference_solutions/woodward_colella_200000_theta1.3_char_LDCU_CFL0.25"
elif (nametest.startswith("RP2")):
    solfoldName="reference_solutions/RP2"
elif (nametest.startswith("shock_turbulence_interaction_shu_osher")):
    solfoldName="reference_solutions/shock_turbulence_interaction_shu_osher_200000_theta1.3_char_LDCU_CFL0.25"
else:
    solfoldName="reference_solutions/"+nametest

    
if (1==1): #IN CASE YOU DO NOT WANT THE REFERENCE
    if os.path.isdir(solfoldName):  #CONDITION: Is it a folder? If yes go on
        namefile=solfoldName+'/'+"reference.out"
        if os.path.isfile(namefile): #If the file exists, read it
            print("You are in the folder "+solfoldName+" dealing with "+namefile)
            lines=[]
            with open(namefile, 'r') as a:
                for idx, line in enumerate(a.readlines()):
                    lines.append(line)

            x_ref=np.array([])
            ro_ref=np.array([])
            q_ref=np.array([])
            E_ref=np.array([])
            v_ref=np.array([])
            p_ref=np.array([])

            #Loop over the lines that we have read
            for idx, line in enumerate(lines): #rmk: 0-based numeration
                data=line.split()
                x_ref=np.append(x_ref,float(data[0]))
                ro_ref=np.append(ro_ref,float(data[1]))
                v_ref=np.append(v_ref,float(data[2]))
                p_ref=np.append(p_ref,float(data[3]))
                # E_ref=np.append(E_ref,float(data[4]))

            q_ref=ro_ref*v_ref
            E_ref=p_ref/(gmm-1)+0.5*ro_ref*v_ref**2

            ax_ro.plot(x_ref,ro_ref,linestyle="-",linewidth=lwex,label=solution_type,color="k")
            ax_v.plot( x_ref, v_ref,linestyle="-",linewidth=lwex,label=solution_type,color="k")
            ax_p.plot( x_ref, p_ref,linestyle="-",linewidth=lwex,label=solution_type,color="k")
            # ax_q.plot( x_ref, q_ref,linestyle="-",linewidth=lwex,label=solution_type,color="k")
            # ax_E.plot( x_ref, E_ref,linestyle="-",linewidth=lwex,label=solution_type,color="k")
            index_plot=index_plot+1


name_base_folder="simulations_paper"
# name_base_folder="simulations_splitting"


for indt, test in enumerate(teststocompare): #Loop on the schemes
    space_reconstruction=test[0]
    time_scheme         =test[1]
    post_processing     =test[2]
    CFL                 =test[3]
    foldName=name_base_folder+"/"+nametest+"/space_reconstruction"+str(space_reconstruction)+"/time_scheme_"+str(time_scheme)+"/PP"+str(post_processing)+"/CFL"+str(CFL)
    print(foldName)
    if os.path.isdir(foldName):  #CONDITION: Is it a folder? If yes go on
        filename = "WS_MESH_Coordinate"+coordinate+".dat"
        if os.path.isfile(foldName+"/"+ filename):
            delimiter_in = ' '
            headerlines_in = 1
            mydata_mesh = np.loadtxt(foldName+"/"+ filename, delimiter=delimiter_in, skiprows=headerlines_in)
            z = mydata_mesh
            z_ini = min(z)
            z_end = max(z)
            Nz    = len(z)

            # Importing the SOLUTION DATA
            if os.path.isfile(foldName+"/"+ file):

                mydata_solution = np.loadtxt(foldName+"/"+ file, skiprows=headerlines_in)

                ro   = mydata_solution[:, 0]
                u    = mydata_solution[:, 1]
                v    = mydata_solution[:, 2]
                p    = mydata_solution[:, 3]
                phi  = mydata_solution[:, 4]

                RO     = ro[:Nz]
                U      = u[:Nz]
                V      = v[:Nz]
                P      = p[:Nz]
                PHI    = phi[:Nz]


                ax_ro.plot(z,RO,marker=markers_space_reconstruction[space_reconstruction],linestyle=linestyles_space_reconstruction[space_reconstruction],linewidth=lw, markersize=ms,label=label_space_reconstruction[space_reconstruction],color=colors[space_reconstruction],alpha=0.7,zorder=3)
                ax_v.plot( z, U,marker=markers_space_reconstruction[space_reconstruction],linestyle=linestyles_space_reconstruction[space_reconstruction],linewidth=lw, markersize=ms,label=label_space_reconstruction[space_reconstruction],color=colors[space_reconstruction],alpha=0.7,zorder=3)
                ax_p.plot( z, P,marker=markers_space_reconstruction[space_reconstruction],linestyle=linestyles_space_reconstruction[space_reconstruction],linewidth=lw, markersize=ms,label=label_space_reconstruction[space_reconstruction],color=colors[space_reconstruction],alpha=0.7,zorder=3)
                ax_ro.set_xlim(xspan)
                ax_v.set_xlim(xspan)
                ax_p.set_xlim(xspan)
                index_plot=index_plot+1
            



import numpy as np
import matplotlib.pyplot as plt
import os

# Gas constant
gamma = 1.4

# Use the same folder as the script for the data file
file_name = 'RP2_HLL_200_-230_SSPRK3_CFL0.475/SOLUTION_0000000.1500000.dat'
file_path = os.path.join(os.path.dirname(__file__), file_name)

# Load data: columns are x, rho, rho*u, E
data = np.loadtxt(file_path)

x = data[:, 0]
rho = data[:, 1]
rho_u = data[:, 2]
E = data[:, 3]

# Compute primitive variables
u = rho_u / rho
p = (gamma - 1) * (E - 0.5 * rho * u**2)

ax_ro.plot(x,rho,marker=markers_space_reconstruction[space_reconstruction],linestyle=linestyles_space_reconstruction[space_reconstruction],linewidth=lw, markersize=ms,label="HLL",color="blue",alpha=0.7,zorder=3)
ax_v.plot( x, u,marker=markers_space_reconstruction[space_reconstruction],linestyle=linestyles_space_reconstruction[space_reconstruction],linewidth=lw, markersize=ms,label="HLL",color="blue",alpha=0.7,zorder=3)
ax_p.plot( x, p,marker=markers_space_reconstruction[space_reconstruction],linestyle=linestyles_space_reconstruction[space_reconstruction],linewidth=lw, markersize=ms,label="HLL",color="blue",alpha=0.7,zorder=3)

index_plot=index_plot+1


plt.tight_layout()
# plt.show()

# ax_ro.grid()
# ax_v.grid()
# ax_p.grid()
# ax_q.grid()
# ax_E.grid()
# pl.grid()

# Capturing handles and labels
handles, labels = pl.gca().get_legend_handles_labels()
# Define the desired order
desired_order=np.arange(index_plot)
# print(np.flip(desired_order))
desired_order=np.flip(desired_order)

# Reorder handles and labels
handles[:] = [handles[i] for i in desired_order[:]]
labels[:] = [labels[i] for i in desired_order[:]]



# Create the legend with reordered handles and labels
pl.legend(handles, labels,fontsize=fs+2)#,loc='center left', bbox_to_anchor=(1, 0.5),fontsize=fs)

pl.savefig(nametest+"_compare_HLL_V_SSPRK"+str(time_scheme)+"_PP"+str(post_processing)+"_CFL"+str(CFL)+".pdf", format="pdf", bbox_inches="tight")
# pl.savefig(nametest+"_CFL"+str(CFL)+"_compare_char_plot.pdf", format="pdf", bbox_inches="tight")
# pl.show()

fig = pl.gcf()
fig.canvas.draw()  # Necessary for accurate bbox

def save_axis(ax, filename):
    extent = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(filename, bbox_inches=extent)


ax_ro.legend(handles, labels,fontsize=fs+2)#,loc='center left', bbox_to_anchor=(1, 0.5),fontsize=fs)
ax_v.legend(handles, labels,fontsize=fs+2)#,loc='center left', bbox_to_anchor=(1, 0.5),fontsize=fs)
ax_p.legend(handles, labels,fontsize=fs+2)#,loc='center left', bbox_to_anchor=(1, 0.5),fontsize=fs)


save_axis(ax_ro, nametest+"_compare_HLL_V_SSPRK"+str(time_scheme)+"_PP"+str(post_processing)+"_CFL"+str(CFL)+"_rho.pdf")
save_axis(ax_v, nametest+"_compare_HLL_V_SSPRK"+str(time_scheme)+"_PP"+str(post_processing)+"_CFL"+str(CFL)+"_u.pdf")
save_axis(ax_p, nametest+"_compare_HLL_V_SSPRK"+str(time_scheme)+"_PP"+str(post_processing)+"_CFL"+str(CFL)+"_p.pdf")