import numpy as np
import shutil, os
import re
from joblib import Parallel, delayed
import multiprocessing
import sys




numbertest=40
starting_elements_X=2000
starting_elements_Y=500
nametest="bubble_helium_"+str(starting_elements_X)+"_"+str(starting_elements_Y)
space_reconstructions = [-27] #-np.array([-1, -20, -21, -22, -23, -24, -26]) #put sign -
time_schemes           =[3]
post_processings       =[12]
reconstructed_variables=[1]
N_refinements = {1:1,2:1,20:1,21:1,22:1,23:1,24:1,25:1,26:1,27:1,-1:1,-2:1,-20:1,-21:1,-22:1,-23:1,-24:1,-25:1,-26:1,-27:1}
starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
test_type="2D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
CFL=0.25 #0.95

# numbertest=41
# starting_elements_X=2000
# starting_elements_Y=500
# nametest="bubble_r22_"+str(starting_elements_X)+"_"+str(starting_elements_Y)
# space_reconstructions = [-27] #-np.array([-1, -20, -21, -22, -23, -24, -26]) #put sign -
# time_schemes           =[3]
# post_processings       =[12]
# reconstructed_variables=[1]
# N_refinements = {1:1,2:1,20:1,21:1,22:1,23:1,24:1,25:1,26:1,27:1,-1:1,-2:1,-20:1,-21:1,-22:1,-23:1,-24:1,-25:1,-26:1,-27:1}
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="2D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# CFL=0.25 #0.95





# numbertest=28
# nametest="multifluid_sod"
# space_reconstructions = [-27] #-np.array([-1, -20, -21, -22, -23, -24, -26]) #put sign -
# time_schemes         =[3]
# post_processings     =[12]
# N_refinements = {1:1,2:1,20:1,21:1,22:1,23:1,24:1,25:1,26:1,27:1,-1:1,-2:1,-20:1,-21:1,-22:1,-23:1,-24:1,-25:1,-26:1,-27:1}
# starting_elements_X=200
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="1D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=5
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.475 #0.95


# numbertest=30
# nametest="multifluid_stiff_sod"
# space_reconstructions = [-27] #-np.array([-1, -20, -21, -22, -23, -24, -26]) #put sign -
# time_schemes         =[3]
# post_processings     =[12]
# N_refinements = {1:1,2:1,20:1,21:1,22:1,23:1,24:1,25:1,26:1,27:1,-1:1,-2:1,-20:1,-21:1,-22:1,-23:1,-24:1,-25:1,-26:1,-27:1}
# starting_elements_X=400
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="1D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=5
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.475 #0.95

# numbertest=31
# nametest="multifluid_water_air"
# space_reconstructions = [-27] #-np.array([-1, -20, -21, -22, -23, -24, -26]) #put sign -
# time_schemes         =[3]
# post_processings     =[12]
# N_refinements = {1:1,2:1,20:1,21:1,22:1,23:1,24:1,25:1,26:1,27:1,-1:1,-2:1,-20:1,-21:1,-22:1,-23:1,-24:1,-25:1,-26:1,-27:1}
# starting_elements_X=400
# starting_elements_X_dict = {key: starting_elements_X for key in space_reconstructions}  # use elements of 'orders' as keys
# test_type="1D" #"1D" #"2D" #Depending on this I will modify or not the elements in the Y direction
# starting_elements_Y=5
# starting_elements_Y_dict = {key: starting_elements_Y for key in space_reconstructions}  # use elements of 'orders' as keys
# speed_estimate=0
# CFL=0.475 #0.95

####################################################
name_reconstructed_variable = {0:"cons",    1:"char",    2:"prim",    -1:"char_interface"}
name_riemann_solver         = {-2:"central", -1:"LF", 0:"rusanov", 1:"exact", 2:"HLL", 3:"CU", 4:"LDCU", 5:"HLLC", 6:"FORCE",7:"FORCEalpha"}
name_speed_estimate         = {0:"standard",1:"riemann"}
####################################################

name_base_folder="simulations_multifluid"

# #Clean, Compile
# instruction="make clean; make"
# os.system(instruction)

#Create simulations folders and copy main there
for space_reconstruction in space_reconstructions:
    for time_scheme in time_schemes:
        for post_processing in post_processings:
            foldName=name_base_folder+"/"+nametest+"/space_reconstruction"+str(space_reconstruction)+"/time_scheme_"+str(time_scheme)+"/PP"+str(post_processing)+"/CFL"+str(CFL)
            instruction =""
            instruction+="mkdir -p "+foldName+" \n" #creation of the folder
            os.system(instruction)
            instruction =""
            instruction+="cp bin/main "+foldName+"/main \n" #copying DATA/don1d
            os.system(instruction)

#Run simulations
for space_reconstruction in space_reconstructions:
    for time_scheme in time_schemes:
        for post_processing in post_processings:
            foldName=name_base_folder+"/"+nametest+"/space_reconstruction"+str(space_reconstruction)+"/time_scheme_"+str(time_scheme)+"/PP"+str(post_processing)+"/CFL"+str(CFL)

            elements_X=np.zeros(N_refinements[space_reconstruction])
            elements_X[0]=starting_elements_X_dict[space_reconstruction]
            for ind in range(1,N_refinements[space_reconstruction]):
                elements_X[ind]=elements_X[ind-1]*2

            elements_Y=np.zeros(N_refinements[space_reconstruction])
            elements_Y[0]=starting_elements_Y_dict[space_reconstruction]
            for ind in range(1,N_refinements[space_reconstruction]):
                elements_Y[ind]=elements_Y[ind-1]*2

            instruction =""
            instruction+="cd "+foldName+" \n" #move in the folder
            for indi,element_X in enumerate(elements_X):        
                if test_type=="2D":
                    instruction+="./main "+str(numbertest)+" "+str(int(element_X))+" "+str(int(elements_Y[indi]))+" "+str(space_reconstruction)+" "+str(space_reconstruction)+" "+str(time_scheme)+" "+str(CFL)+" "+str(post_processing)+" "+" \n" 
                elif test_type=="1D":
                    instruction+="./main "+str(numbertest)+" "+str(int(element_X))+" "+str(int(elements_Y[0]))+" "+str(space_reconstruction)+" "+str(space_reconstruction)+" "+str(time_scheme)+" "+str(CFL)+" "+str(post_processing)+" "+" \n" 
                else:
                    print("Test type not recognized")
                    quit()


            print(instruction)
            os.system(instruction)
