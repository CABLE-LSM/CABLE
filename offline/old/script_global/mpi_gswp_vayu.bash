#!/bin/bash
###!!!#PBS -l nodes=1:ppn=8  #  8 CPUs, 1.5 GB vmem, ~22 minutes
###!!!#PBS -l nodes=2:ppn=8  # 16 CPUs, 3.0 GB vmem, ~13 minutes  
###!!!#PBS -l nodes=3:ppn=8  # 24 CPUs, 5.1 GB vmem, ~12 minutes
#PBS -l nodes=2:ppn=8
#PBS -l vmem=6gb
#PBS -l walltime=40:00
module add netcdf openmpi
cd /short/p66/YOUR_IDENT/runCABLE2.0/
ulimit -s 65536
yr=1986
while [ $yr -le 1995 ]
do
  cp namelistDir/icycle3/cable.nml_gswp${yr} cable.nml
  ## cp namelistDir/icycle1/cable.nml_gswp${yr} cable.nml
  mpirun -np 16 ./cable-mpi
  ###!!! mpirun -np  8 ./cable-mpi   #  8 CPUs
  ###!!! mpirun -np 16 ./cable-mpi   # 16 CPUs
  ###!!! mpirun -np 24 ./cable-mpi   # 24 CPUs
  mv out_cable.nc      out_gswp/out_gswp${yr}.nc
  mv log_cable.txt     out_gswp/log_gswp${yr}.txt
  cp -p restart_out.nc out_gswp/restart_gswp${yr}.nc
  mv restart_out.nc    restart_in_gswp.nc
  yr=`expr $yr + 1`
done

