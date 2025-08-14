


############### VERSION OF ELMER/ICE ##########################

MAde to run on 8 partitions


###### SIF #######

SIF is the directory where all the SIF are stored, 
 * TOPOGRAPHY* to make the TOPOGRAPHY restart file that you will use to lauch the simulation
 * OPTI* to do the optimization on the variable eta
 * STRESS* to calculate the Cauchy stress tensor from the calculated velocities, need to be calculated from an optimized beta, and set a Weertman law !

 
###### TEST ####

You can launch an optimization on eta, the effective viscosity from RUN_TEST.sh 
