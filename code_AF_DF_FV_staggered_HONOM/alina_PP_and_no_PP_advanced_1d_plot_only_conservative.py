coordinate="X" #"X" #"Y" #NB:Actually y is not implemented
lw=2. #2. for shock interaction
ms=3.0 #1. 
fs=10
#overwritten later if needed



#####################
#CFL 0.8
#100 elements
#####################
#Let us assume standard speed estimate.
#The only simulation crashing is cons with exact riemann solver for order 13.
#All the others manage to arrive to the final time.
#It is possible to make also this setting work if using the exact riemann solver speed estimate.
#In order to run with standard speed estimate, for order 13, with such setting, the CFL must be decreased to 0.8.
#We perform the comparison for CFL=0.8
#The results with reconstruction of conserved variables present some slight oscillations in the central plateau of velocity and pressure.
#Results with rusanov and exact rieman solver are comparable. 
#Rusanov is a bit more diffusive and prevents over- and undershoots when reconstructing characteristic variables 
#Slight over and undershoots are present for the exact riemann solver even when characteristic variables are reconstructed
#Overall the best setting seems to be characteristic variables with Rusanov
#####################
nametest="sod"
file="SOLUTION_OCT_0000000.2000000.dat"
CFL_test=0.475 #0.8 #original #0.95 #numerical fluxes
solution_type="exact solution"
xspan=[0,1]





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

teststocompare.append([-27,   3, 12, CFL_test]) 
teststocompare.append([-27,   3, 0, CFL_test]) 


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

colors_PP = {
    -1:"red", #AAAAAAAAAAAAAAAAAAAAAA
    1:"red", #AAAAAAAAAAAAAAAAAAAAAA
    12:"red", #AAAAAAAAAAAAAAAAAAAAAA
    0:"blue", #AAAAAAAAAAAAAAAAAAAAAA

}

linestyles_PP = {
    -1:"-",
    1:"-", 
    12:"-", 
    0:"--", 
}

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

label_PP   =   { 0: "DF-FV without PP"          ,
                 1:"DF-FV",
                 -1:"DF-FV",
                 12:"DF-FV",
                }

nVar=3


params = {'mathtext.default': 'regular' } #parameters plot




fig = pl.figure(1, figsize=(20,5))
ax_ro = fig.add_subplot(131)
ax_ro.set_title( r'$\rho$',fontsize=fs+5)
# ax_ro.set_xlabel('x',fontsize=fs)
# ax_ro.set_ylabel('density',fontsize=fs)
ax_q  = fig.add_subplot(132)
ax_q.set_title( r'$\rho u$',fontsize=fs+5)
# ax_q.set_xlabel('x',fontsize=fs)
# ax_v.set_ylabel('velocity',fontsize=fs)
ax_E  = fig.add_subplot(133)
ax_E.set_title( r'$E$',fontsize=fs+5)
# ax_E.set_xlabel('x',fontsize=fs)
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
ax_q.tick_params(axis='both', which='major',  labelsize=15)
ax_E.tick_params(axis='both', which='major',  labelsize=15)



##############################
# READ AND PLOT EXACT SOLUTION
##############################
# FILLING x_ref AND ref_sol
##############################
if (nametest.startswith("composite_wave")):
    solfoldName="reference_solutions/composite_wave"
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
            ax_q.plot( x_ref, q_ref,linestyle="-",linewidth=lwex,label=solution_type,color="k")
            ax_E.plot( x_ref, E_ref,linestyle="-",linewidth=lwex,label=solution_type,color="k")
            # ax_q.plot( x_ref, q_ref,linestyle="-",linewidth=lwex,label=solution_type,color="k")
            # ax_E.plot( x_ref, E_ref,linestyle="-",linewidth=lwex,label=solution_type,color="k")
            index_plot=index_plot+1


name_base_folder="simulations_paper"


for indt, test in enumerate(teststocompare): #Loop on the schemes
    space_reconstruction=test[0]
    time_scheme         =test[1]
    post_processing     =test[2]
    CFL                 =test[3]
    foldName=name_base_folder+"/"+nametest+"/space_reconstruction"+str(space_reconstruction)+"/time_scheme_"+str(time_scheme)+"/PP"+str(post_processing)+"/CFL"+str(CFL)
    print(foldName)
    if os.path.isdir(foldName):  #CONDITION: Is it a folder? If yes go on
        filename = "MESH_Coordinate"+coordinate+".dat"
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
                rou    = mydata_solution[:, 1]
                rov    = mydata_solution[:, 2]
                energy    = mydata_solution[:, 3]
                phi  = mydata_solution[:, 4]

                RO     = ro[:Nz]
                ROU    = rou[:Nz]
                ROV    = rov[:Nz]
                ENERGY = energy[:Nz]
                PHI    = phi[:Nz]

                U=ROU/RO
                V=ROV/RO

                P=(gmm-1.0)*(ENERGY-0.5*RO*(U**2+V**2))

                ax_ro.plot(z,RO,marker=markers_space_reconstruction[space_reconstruction],linestyle=linestyles_PP[post_processing],linewidth=lw, markersize=ms,label=label_PP[post_processing],color=colors_PP[post_processing],alpha=0.7,zorder=3)
                ax_q.plot( z,ROU,marker=markers_space_reconstruction[space_reconstruction],linestyle=linestyles_PP[post_processing],linewidth=lw, markersize=ms,label=label_PP[post_processing],color=colors_PP[post_processing],alpha=0.7,zorder=3)
                ax_E.plot( z,ENERGY,marker=markers_space_reconstruction[space_reconstruction],linestyle=linestyles_PP[post_processing],linewidth=lw, markersize=ms,label=label_PP[post_processing],color=colors_PP[post_processing],alpha=0.7,zorder=3)
                ax_ro.set_xlim(xspan)
                ax_q.set_xlim(xspan)
                ax_q.set_ylim([-0.05, np.max(ROU)*1.2])
                ax_E.set_xlim(xspan)
                index_plot=index_plot+1
            




# ax_ro.grid()
# ax_q.grid()
# ax_E.grid()
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
# handles[:] = [handles[i] for i in desired_order[:]]
# labels[:] = [labels[i] for i in desired_order[:]]



# Create the legend with reordered handles and labels
pl.legend(handles, labels,fontsize=fs+2)#,loc='center left', bbox_to_anchor=(1, 0.5),fontsize=fs)

pl.savefig(nametest+"_U_SSPRK"+str(time_scheme)+"_CFL"+str(CFL)+"_compare_PP_and_not_PP.pdf", format="pdf", bbox_inches="tight")
# pl.savefig(nametest+"_CFL"+str(CFL)+"_compare_char_plot.pdf", format="pdf", bbox_inches="tight")
# pl.show()

#SAVE FIGURES INDIVIDUALLY
from matplotlib.transforms import Bbox

fig = pl.gcf()
fig.canvas.draw()  # Necessary for accurate bbox

def save_axis(ax, filename):
    extent = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(filename, bbox_inches=extent)


ax_ro.legend(handles, labels,fontsize=fs+2)#,loc='center left', bbox_to_anchor=(1, 0.5),fontsize=fs)
ax_q.legend(handles, labels,fontsize=fs+2,loc='upper left')#,loc='center left', bbox_to_anchor=(1, 0.5),fontsize=fs)
ax_E.legend(handles, labels,fontsize=fs+2)#,loc='center left', bbox_to_anchor=(1, 0.5),fontsize=fs)



save_axis(ax_ro, nametest+"_U_SSPRK"+str(time_scheme)+"_CFL"+str(CFL)+"_compare_PP_and_not_PP_rho.pdf")
save_axis(ax_q, nametest+"_U_SSPRK"+str(time_scheme)+"_CFL"+str(CFL)+"_compare_PP_and_not_PP_rhou.pdf")
save_axis(ax_E, nametest+"_U_SSPRK"+str(time_scheme)+"_CFL"+str(CFL)+"_compare_PP_and_not_PP_E.pdf")