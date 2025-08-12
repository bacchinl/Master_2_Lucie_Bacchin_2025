#!/bin/bash

#OAR -n PIG_IS_ComputeStress
#OAR -l /nodes=1/core=10,walltime=20:00:00
#OAR --stdout PIG_run_only_stress.out
#OAR --stderr PIG_run_only_stress.err
#OAR --project sno-elmerice

source ~/.nix-profile/setvars.sh > /dev/null
module purge
module load  netcdf/netcdf-4.7.2_intel21_hdf5_MPIO xios/xios-2.0_rev2320_intel21 elmerfem/elmerfem-6aa7bc75c_intel21



## compile required USFs
make USFs_PIG

## Data file
#DATAFILE="..\/..\/DATA\/VELOCITY\/ASE_TimeSeries_1973-2018.nc"
  

## regularisation parameters

#lambda='5.0e05 1.0e06 2.0e06 5.0e06 1.0e07 2.0e07 5.0e07 1.0e08 1.0e09 1.0e10 1.0e11'

#lambda='0.0e00 1.0e03 1.0e04 2.0e04 5.0e04 1.0e05 2.0e05 5.0e05 1.0e06 2.0e06 5.0e06 1.0e07 1.0e08 1.0e09 1.0e10 1.0e11'
years=' 2007 2009 2012 2016 2018 2020 2022'
#years='2012'
lambda='0.0e00'

#echo SIF/TOPOGRAPHY.sif > ELMERSOLVER_STARTINFO
#mpirun -np 8 ElmerSolver_mpi


for year in $years
do 

	for i in $lambda
 	 do
	

  	 
	NAME=OPT_"$year"

	sed  "s/<Lambda>/"$i"/g;s/<NAME>/$NAME/g;s/<OBS_FILE>/$DATAFILE/g;s/<YEAR>/$year/g" SIF/STRESS.sif > Simu_$year/STRESS_$year.sif
	echo Simu_$year/STRESS_$year.sif > ELMERSOLVER_STARTINFO
  	## Has to be parallel on 2 partition to restart initial file
  	mpirun -np 8 ElmerSolver_mpi


  	 #python ../SCRIPTS/MakeReport.py $NAME
  	 #echo $(tail -n 1 Cost_"$NAME".dat | awk '{print $3}') $(tail -n 1 CostReg_"$NAME".dat | awk '{print $2}') $i $c >> Simu_$year/LCurve_$year.dat
 	done
done
