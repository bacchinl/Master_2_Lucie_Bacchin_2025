#!/bin/bash

#OAR -n PIG_IS_viscosity
#OAR -l /nodes=1/core=10,walltime=20:00:00
#OAR --stdout PIG_run.out
#OAR --stderr PIG_run.err
#OAR --project sno-elmerice




## compile required USFs
make USFs_PIG

## Data file
#DATAFILE="..\/..\/DATA\/VELOCITY\/ASE_TimeSeries_1973-2018.nc"
  

## regularisation parameters

#lambda='5.0e05 1.0e06 2.0e06 5.0e06 1.0e07 2.0e07 5.0e07 1.0e08 1.0e09 1.0e10 1.0e11'

#lambda='0.0e00 1.0e03 1.0e04 2.0e04 5.0e04 1.0e05 2.0e05 5.0e05 1.0e06 2.0e06 5.0e06 1.0e07 1.0e08 1.0e09 1.0e10 1.0e11'
years='2005 2011 2014 2017'
lambda='0.0e00'

#echo SIF/TOPOGRAPHY.sif > ELMERSOLVER_STARTINFO
#mpirun -np 8 ElmerSolver_mpi

c=0
for year in $years
do 

	#mkdir Simu_$year
	rm -rf LCurve_$year.dat
	DATAFILE="..\/..\/..\/..\/DATA\/VELOCITY\/velocity_lucille_$year.nc"
	
	echo "open"
	sed  "s/<YEAR>/$year/g" SIF/TOPOGRAPHY_Paolo.sif > SIF/TOPOGRAPHY_$year.sif
	echo SIF/TOPOGRAPHY_$year.sif > ELMERSOLVER_STARTINFO
	mpirun -np 8 ElmerSolver_mpi
	c=0
	

	for i in $lambda
 	 do

  	 c=$((c+1))

  	 echo $i
	 NAME=OPT_"$year-$c"
  	 #sed  "s/<Lambda>/"$i"/g;s/<NAME>/$NAME/g;s/<OBS_FILE>/$DATAFILE/g;s/<YEAR>/$year/g" SIF/OPTIM_ETA_STRESS.sif > Simu_$year/OPTIM_STRESS_$year-$c.sif
	 sed  "s/<Lambda>/"$i"/g;s/<NAME>/$NAME/g;s/<OBS_FILE>/$DATAFILE/g;s/<YEAR>/$year/g" SIF/OPTIM_ETA.sif > Simu_$year/OPTIM_$year-$c.sif


  	 echo Simu_$year/OPTIM_$year-$c.sif > ELMERSOLVER_STARTINFO
  	 # Has to be parallel on 2 partition to restart initial file
  	 mpirun -np 8 ElmerSolver_mpi

		

	#sed  "s/<Lambda>/"$i"/g;s/<NAME>/$NAME/g;s/<OBS_FILE>/$DATAFILE/g;s/<YEAR>/$year/g" SIF/STRESS.sif > Simu_$year/STRESS_$year-$c.sif
	#echo Simu_$year/STRESS_$year-$c.sif > ELMERSOLVER_STARTINFO
  	## Has to be parallel on 2 partition to restart initial file
  	#mpirun -np 8 ElmerSolver_mpi


  	 python ../SCRIPTS/MakeReport.py $NAME
  	 echo $(tail -n 1 Cost_"$NAME".dat | awk '{print $3}') $(tail -n 1 CostReg_"$NAME".dat | awk '{print $2}') $i $c >> Simu_$year/LCurve_$year.dat
 	done
done
