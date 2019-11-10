#!/bin/bash
###!!! #PBS -l nodes=1:ppn=12   # 12 CPUs, 7.2 GB vmem, ~17 minutes
###!!! #PBS -l nodes=2:ppn=12   # 24 CPUs, 11  GB vmem, ~11 minutes
#PBS -l nodes=2:ppn=12
#PBS -l vmem=12gb
#PBS -l walltime=00:25:00
module add netcdf openmpi
cd /data/flush/YOUR_IDENT/CABLE-run
ulimit -s 65536
ulimit -c unlimited
ulimit -a
yr=1986
while [ $yr -le 1995 ]
do
  cp namelist_r491/350ppm/cable.nml_gswp${yr} cable.nml
  mpirun -np 24 ./cable-mpi
  #tar -cvf out_gswp/${yr}fort.tar fort.*
  #rm fort.*
  mv out_cable.nc      out_gswp/out_gswp${yr}.nc
  mv log_cable.txt     out_gswp/log_gswp${yr}.txt
  cp -p restart_out.nc out_gswp/restart_gswp${yr}.nc
  mv restart_out.nc    restart_in_gswp.nc
  yr=`expr $yr + 1`
done

