#!/bin/bash

#OAR -n PIG_IS_Optim_visco
#OAR -l /nodes=1/core=10,walltime=47:00:00
#OAR --stdout PIG_run_stress2.out
#OAR --stderr PIG_run_stress2.err
#OAR --project sno-elmerice

source ~/.nix-profile/setvars.sh > /dev/null
module purge
module load  netcdf/netcdf-4.7.2_intel21_hdf5_MPIO xios/xios-2.0_rev2320_intel21 elmerfem/elmerfem-6aa7bc75c_intel21



## compile required USFs
make USFs_PIG

years='2016 '
lambda='0.0'

#echo SIF/TOPOGRAPHY.sif > ELMERSOLVER_STARTINFO
#mpirun -np 8 ElmerSolver_mpi


for year in $years
do 

	#mkdir Simu_$year
	rm -rf LCurve_$year.dat
	DATAFILE="..\/..\/..\/..\/DATA\/VELOCITY\/velocity_lucille_$year.nc"
	
	echo "open"
	
	
	for i in $lambda
 	 do

  	 c=$((c+1))

  	 echo $i
	 NAME=OPT_"$year-$c"
  	 #sed  "s/<Lambda>/"$i"/g;s/<NAME>/$NAME/g;s/<OBS_FILE>/$DATAFILE/g;s/<YEAR>/$year/g" SIF/OPTIM_ETA_STRESS.sif > Simu_$year/OPTIM_STRESS_$year-$c.sif
	 sed  "s/<Lambda>/"$i"/g;s/<NAME>/$NAME/g;s/<OBS_FILE>/$DATAFILE/g;s/<YEAR>/$year/g" SIF/OPTIM_ETA_L.sif > Simu_$year/OPTIM_L_$year-$c.sif


  	 echo Simu_$year/OPTIM_L_$year-$c.sif > ELMERSOLVER_STARTINFO
  	 # Has to be parallel on 2 partition to restart initial file
  	 mpirun -np 8 ElmerSolver_mpi

		
 	done
done
