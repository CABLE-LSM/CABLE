#!/bin/bash
#PBS -l vmem=4gb
#PBS -l walltime=03:20:00   # 1 CPU, ~2:20:00
cd /short/p66/YOUR_IDENT/runCABLE2.0
#ulimit -c unlimited
ulimit -s 65536
yr=1986
while [ $yr -le 1995 ]
do
  cp namelistDir/icycle1/cable.nml_gswp${yr} cable.nml
  ## cp namelistDir/icycle3/cable.nml_gswp${yr} cable.nml
  ./cable
  mv out_cable.nc      out_gswp/out_gswp${yr}.nc
  mv log_cable.txt     out_gswp/log_gswp${yr}.txt
  cp -p restart_out.nc out_gswp/restart_gswp${yr}.nc
  mv restart_out.nc    restart_in_gswp.nc
  yr=`expr $yr + 1`
done

