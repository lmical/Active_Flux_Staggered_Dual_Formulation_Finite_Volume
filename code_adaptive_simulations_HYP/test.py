import numpy as np
import shutil, os
import re
from joblib import Parallel, delayed
import multiprocessing
import sys

CodeFolder="src"
ReferenceTestFolder="tests/test"

#Change folder
instruction="cd "+CodeFolder+" \n"
#Clean
instruction+="make clean \n"
#Compile
instruction+="make"
os.system(instruction)

#Change folder
instruction="cd "+ReferenceTestFolder+" \n"

#Remove existing .dat
instruction+="rm *.dat\n"

#Run





# instruction+="../../bin/main 8 200 200 0.5 2 -2\n"

# instruction+="../../bin/main 6 200 200 0.5 2 -2\n"

#SAFETY CHECK
# PP1
# instruction+="../../bin/main 3 25 25 2 2 0.5 1 1.0\n"
# instruction+="../../bin/main 3 25 25 2 -5 0.5 1 1.0\n"

# PP0
# instruction+="../../bin/main 3 25 25 2 2 0.5 0 1.0\n"
# instruction+="../../bin/main 3 25 25 2 -5 0.5 0 1.0\n"

# instruction+="../../bin/main 3 400 400 2 2 0.5 0 1.0\n"
# instruction+="../../bin/main 3 400 400 2 -5 0.5 0 1.0\n"

#SOD
# instruction+="../../bin/main 7 200 200 2 2 0.5 1 1.0\n"
# instruction+="../../bin/main 7 200 200 2 -2 0.5 1 1.0\n"

#SOD 1D
# instruction+="../../bin/main 4 50 50 -1 1 0.3 0 1.0\n"
# instruction+="../../bin/main 8 200 3 26 3 0.3 1 1.0\n"


# instruction+="../../bin/main 8 200 3 26 26 3 0.3 1 1.0\n"
# instruction+="../../bin/main 8 200 3 28 28 3 0.3 1 1.0\n"
# instruction+="../../bin/main 8 200 3 28 29 3 0.3 1 1.0\n"

# instruction+="../../bin/main 350 400 3 28 28 3 0.3 1 1.0\n" #Woodward colella
# instruction+="../../bin/main 350 400 3 28 29 3 0.3 1 1.0\n" #Woodward colella


# instruction+="../../bin/main -302 800 3 28 28 3 0.3 1 1.0\n" #Woodward colella


# instruction+="../../bin/main 8 200 3 28 29 3 0.45 0.5 1\n" 
# instruction+="../../bin/main 350 400 3 28 29 3 0.45 0.05 1\n" 
# instruction+="../../bin/main 350 1600 3 28 29 3 0.4 0.05 1\n" 



# instruction+="../../bin/main 3 200 200 28 28 3 0.4 0.025 1\n" 
# instruction+="../../bin/main 10 250 250 28 29 3 0.4 0.025 1\n" 

# instruction+="../../bin/main 7 100 100 28 28 2 0.25 1\n"
# Error  0.45349579196493206       0.31955188388933026       0.31955188388933103        1.2962410225134924     


# instruction+="../../bin/main 7 100 100 28 29 2 0.25 1\n"
# FCFLAGS += -I. -DCENTEREDPRIMITIVE
# FCFLAGS += -I. -DRECONSTRUCTFROMCONSERVED
# FCFLAGS += -I. -DPATHCONSERVATIVESHOCKDETECTION
# FCFLAGS += -I. -DSMOOTHINGLIMITER
# Error  0.45451303561057155       0.32005349359397323       0.32005349359397256        1.2942214233189244     

# instruction+="../../bin/main 350 400 3 28 29 3 0.3 1 1.0\n" #Woodward colella


# instruction+="../../bin/main 7 200 200 27 3 0.3 1 1.0\n"

# instruction+="../../bin/main 8 200 3 5 4 0.3 1 1.0\n"
# instruction+="../../bin/main 7 200 200 5 4 0.3 1 1.0\n"


# instruction+="../../bin/main 101 50 50 23 -22 0.25 5 0.01\n"

# instruction+="../../bin/main 7 200 200 27 3 0.3 1 1.0\n"

# instruction+="../../bin/main 101 100 100 -23 -22 0.25 5 0.1\n"
# instruction+="../../bin/main 3 100 100 -23 12 0.25 0 1.0\n"
# instruction+="../../bin/main 3 100 100 -23 12 0.25 5 1.0\n"

# instruction+="../../bin/main 3 200 200 -26 -102 0.3 5 1.0\n"

# instruction+="../../bin/main 3 100 100 23 12 0.25 5 1.0\n"
# instruction+="../../bin/main 3 100 100 -23 12 0.25 5 1.0\n"

# instruction+="../../bin/main 18 200 10 23 12 0.25 0 1.0\n"

#X
# instruction+="../../bin/main 28 200 3 -26 3 0.3 1 1.0\n"
# instruction+="../../bin/main 30 400 3 -26 3 0.3 1 1.0\n"
# instruction+="../../bin/main 31 400 3 -26 3 0.3 1 1.0\n"

# instruction+="../../bin/main 30 800 2 -26 3 0.3 1 1.0\n"

# instruction+="../../bin/main 101 100 100 -1 -101 0.25 0 0.1\n"

#8 Sod
# instruction+="../../bin/main 8 200 5 -5 2 0.25 1 1.0\n"
#28 Sod multifluid
# instruction+="../../bin/main 28 200 2 -26 3 0.3 0 1.0\n" #######################
#30 Stiff Sod multifluid
# instruction+="../../bin/main 30 400 5 -23 12 0.25 1 1.0\n"
# instruction+="../../bin/main 30 400 5 -23 12 0.25 1 1.0\n"
# instruction+="../../bin/main 30 400 5 -1 1 0.25 1 1.0\n"

#31 Water air multifluid
# instruction+="../../bin/main 31 400 3 -1 12 0.25 0 1.0\n"

# instruction+="../../bin/main 31 20 3 -26 3 0.5 0 1.0\n"

# instruction+="../../bin/main 28 400 3 -23 3 0.5 1 1.0\n"

# instruction+="../../bin/main 28 200 2 -26 -26 3 0.3 1 1.0\n"

# instruction+="../../bin/main 28 200 2 -26 3 0.3 1 1.0\n"

# instruction+="../../bin/main 30 400 2 -26 3 0.3 1 1.0\n"

# instruction+="../../bin/main 31 400 2 -26 3 0.3 1 1.0\n"

# instruction+="../../bin/main 31 800 2 -26 3 0.3 1 1.0\n"

# instruction+="../../bin/main 30 1600 2 -23 3 0.5 1 1.0\n"

# instruction+="../../bin/main 7 200 200 -26 3 0.3 1 1.0\n"

# instruction+="../../bin/main 3 100 100 -26 3 0.3 1 1.0\n"


# instruction+="../../bin/main 28 200 3 -26 3 0.3 1 1.0\n"

# instruction+="../../bin/main 29 3 200 -26 3 0.3 1 1.0\n"
# instruction+="../../bin/main 8 200 2 -26 3 0.3 1 1.0\n"

# instruction+="../../bin/main 28 200 3 -26 3 0.3 1 1.0\n"


# instruction+="../../bin/main 30 1600 2 -26 3 0.3 1 1.0\n"

# instruction+="../../bin/main 30 400 2 26 3 0.3 1 1.0\n"

# instruction+="../../bin/main 30 400 2 -26 3 0.3 1 1.0\n"
# instruction+="../../bin/main 31 400 2 -26 3 0.3 1 1.0\n"

# instruction+="../../bin/main 28 400 3 -23 3 0.3 1 1.0\n"
# instruction+="../../bin/main 30 400 3 -26 3 0.3 1 1.0\n"



# instruction+="../../bin/main 31 100 3 -26 3 0.5 0 1.0\n"


# instruction+="../../bin/main 31 400 3 -26 3 0.5 0 1.0\n" #play


# instruction+="../../bin/main 31 50 3 -23 12 0.25 0 1.0\n"
#40 helium bubble multifluid
# instruction+="../../bin/main 40 500 200 23 12 0.25 1 1.0\n"

# instruction+="../../bin/main 31 200 5 21 12 0.05 0 1.0\n"


#Y
# instruction+="../../bin/main 29 10 400 1 1 0.25 0 1.0\n"

#EXPLICIT
#SOD 1D
# instruction+="../../bin/main 8 200 3 -27 -27 3 0.3 1 1.0\n" #ACTIVEFLUX
# instruction+="../../bin/main 8 200 3 26 3 0.3 1 1.0\n"

# VORTEX
# instruction+="../../bin/main 3 200 200 -26 3 0.3 -1 1.0\n"

# instruction+="../../bin/main 3 200 200 -26 3 0.3 -1 1.0\n"


# instruction+="../../bin/main 350 400 3 -27 -27 3 0.3 1 1.0\n" #Woodward colella
# instruction+="../../bin/main 401 200 3 -27 -27 3 0.3 1 1.0\n" #RP1
# instruction+="../../bin/main 402 200 3 -27 -27 3 0.3 1 1.0\n" #RP2
# instruction+="../../bin/main 403 200 3 -27 -27 3 0.3 1 1.0\n" #RP3
# instruction+="../../bin/main 404 200 3 -27 -27 3 0.3 1 1.0\n" #RP4
# instruction+="../../bin/main 405 200 3 -27 -27 3 0.3 1 1.0\n" #RP5
# instruction+="../../bin/main 406 200 3 -27 -27 3 0.3 1 1.0\n" #RP6
# instruction+="../../bin/main 500 200 3 -27 -27 3 0.3 1 1.0\n" #RP6
# instruction+="../../bin/main 501 200 3 -27 -27 3 0.3 1 1.0\n" #RP6

# instruction+="../../bin/main 3 100 100 -5 12 0.25 5 1.0\n"


# instruction+="../../bin/main 3 100 100 -26 3 0.3 1 1.0\n" #DEEPACK

# instruction+="../../bin/main 3 400 400 -26 3 0.3 1 1.0\n"

# instruction+="../../bin/main 350 400 3 27 3 0.3 1 1.0\n" #Woodward colella

# instruction+="../../bin/main 10 50 50 28 28 3 0.45 1\n"

# instruction+="../../bin/main 892 500 500 28 28 3 0.45 1\n"


# instruction+="../../bin/main 892 50 50 -26 -26 3 0.3 1 1.0\n"

# instruction+="../../bin/main 8 200 2 -26 -26 3 0.3 1 1.0\n"

# instruction+="../../bin/main 892 500 500 -28 -28 3 0.45 1 1.0\n"


# instruction+="../../bin/main 10 250 250 28 28 3 0.45 1.0\n"

#T=0.01
#ALEX TIMESTEP
#PRIMTIVE CENTERED
# instruction+="../../bin/main 10 250 250 28 28 3 0.45 1.0\n"

        #    1   1.0000125258931736       0.12499996683328857     
        #    2   4.4204584463481637E-003 -0.45753052141369521     
        #    3   4.4204584463482296E-003 -0.45753052141369516     
        #    4   2.5000000762210508       0.34999986570604541    

#NO FLAGS
# instruction+="../../bin/main 10 250 250 28 28 3 0.45\n"

        #SAME

#T=0.01
#MY TIMESTEP
#PRIMTIVE CENTERED
# instruction+="../../bin/main 10 250 250 28 28 3 0.45 1.0 3\n"
        #    1   1.0000233001944101       0.12499997154347926     
        #    2   4.5249735942836621E-003 -0.45768110077329782     
        #    3   4.5249735942837020E-003 -0.45768110077329777     
        #    4   2.5000000551828574       0.34999988671351362     

#NO FLAGS
# instruction+="../../bin/main 10 250 250 28 28 3 0.45\n"
        #SAME


#SOD
# CHAR
# instruction+="../../bin/main 7 100 100 -28 -28 3 0.45 1 1.0\n" #Sod
# FCFLAGS += -I. -DACTIVEFLUX
# FCFLAGS += -I. -DALTERNATIVEFORMULATIONPRESSURE

#  Error  0.45326470896227317       0.31880091886883571       0.31892626202210961        1.3002943694904594     
#  Error W_X  0.45089691502546814        1.0475524083881482        1.0491448667914711       0.47397752608561033     
#  Error W_Y  0.45089126831853221        1.0484827158723624        1.0489755808593624       0.47399565571055585     


#PRIM
#MULTI
# instruction+="../../bin/main 7 100 100 -28 -28 3 0.45 1 1.0\n" #Sod
#  Error  0.45391085609879256       0.31878155132184244       0.31896263802376829        1.3013752392510036        0.0000000000000000     
#  Error W_X  0.45152006480347096        1.0483546856574315        1.0499708503814407       0.47456209698625801        0.0000000000000000     
#  Error W_Y  0.45151485244266204        1.0491680436692283        1.0499387095890862       0.47458559733974415        0.0000000000000000     

#NOMULTI
#  Error  0.45391085609879256       0.31878155132184249       0.31896263802376823        1.3013752392510034     
#  Error W_X  0.45152006480347096        1.0483546856574320        1.0499708503814409       0.47456209698625801     
#  Error W_Y  0.45151485244266198        1.0491680436692279        1.0499387095890860       0.47458559733974398     


# instruction+="../../bin/main 7 100 100 -27 -27 3 0.25 1 1.0\n" #Sod YES
# instruction+="../../bin/main 8 200 3 -27 -27 3 0.25 1 1.0\n" #Sod 1D YES
# instruction+="../../bin/main 401 200 3 -27 -27 3 0.25 1 1.0\n" #Modified Sod 1D YES
# instruction+="../../bin/main 402 201 3 -27 -27 3 0.25 1 1.0\n" #Double rarefaction YES
# instruction+="../../bin/main 403 400 3 -27 -27 3 0.25 1 1.0\n" #Spurious spike in u YES
# instruction+="../../bin/main 404 200 3 -27 -27 3 0.25 1 1.0\n" #NO
# instruction+="../../bin/main 405 200 3 -27 -27 3 0.25 1 1.0\n" #NOOOOOOOOOOO NEGATIVE PRESSURE
# instruction+="../../bin/main 500 200 3 -27 -27 3 0.3 1 1.0\n" #OK
# instruction+="../../bin/main 350 1000 3 -27 -27 3 0.25 1 1.0\n" #WC OK
# instruction+="../../bin/main -300 1000 3 -27 -27 3 0.3 1 1.0\n" 
# instruction+="../../bin/main 412 233 3 -27 -27 3 0.25 1 1.0\n" #Double rarefaction YES


# instruction+="../../bin/main 450 400 3 -27 -27 3 0.1 1 1.0\n" #Double rarefaction NO
# instruction+="../../bin/main 451 800 3 -27 -27 3 0.1 1 1.0\n" #Leblanc



# instruction+="../../bin/main 8 200 3 27 27 3 0.3 1 1.0\n" #Sod 1D YES
# instruction+="../../bin/main 405 200 3 27 27 3 0.3 1 1.0\n" #NOOOOOOOOOOO

# instruction+="../../bin/main 8 200 3 -27 -27 3 0.25 1 1.0\n" #Sod 1D YES
# instruction+="../../bin/main 401 200 3 -27 -27 3 0.25 1 1.0\n" #Modified Sod 1D YES
# instruction+="../../bin/main 402 201 3 -27 -27 3 0.25 1 1.0\n" #Double rarefaction YES
# instruction+="../../bin/main 350 1000 3 -27 -27 3 0.25 1 1.0\n" #WC OK #NB: Positivity checks
# instruction+="../../bin/main -300 800 3 -27 -27 3 0.25 1 1.0\n" #WC OK #NB: Positivity checks

# instruction+="../../bin/main 3 50 50 -27 -27 3 0.25 1 1.0\n" #WC OK #NB: Positivity checks

# instruction+="../../bin/main 7 200 200 -27 -27 3 0.25 1 1.0\n" #Sod 2D YES

# instruction+="../../bin/main 350 1000 3 -27 -27 3 0.25 1 1.0\n" #Sod 2D YES


# instruction+="../../bin/main 7 200 200 -27 -27 3 0.25 1 1.0\n" #Sod 2D YES

# instruction+="../../bin/main 8 200 3 -27 -27 3 0.25 1 1.0\n" #Sod 2D YES

# instruction+="../../bin/main 8 200 3 -27 -27 3 0.25 1 1.0\n" #Sod 2D YES

#  Error  0.15725985556491201       0.17999999999999988        0.0000000000000000       0.46144287067397288     
#  Error W_X  0.15792795984835936       0.44164503931033683        0.0000000000000000       0.17256053325018000     
#  Error W_Y  0.20967980741988271       0.58916933512794534        0.0000000000000000       0.22846314570626702 

# instruction+="../../bin/main 350 400 3 -27 -27 3 0.25 1 1.0\n" #Sod 2D YES


instruction+="../../bin/main 403 200 3 -27 -27 3 0.25 1 1.0\n" #Sod 2D YES

#  Error  0.45447643958372086        11.999879999974214        0.0000000000000000        814.00863319631662     
#  Error W_X  0.45511669192626064        11.676368281929967        0.0000000000000000        306.43842320599805     
#  Error W_Y  0.60596858611162718        15.594770124777629        0.0000000000000000        408.45716657506830     

print(instruction)
os.system(instruction)

#Plot
instruction="python3 tools/plotting_gif.py --dir "+ReferenceTestFolder


os.system(instruction)

sys.exit()
