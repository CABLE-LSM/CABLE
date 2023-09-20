#!/usr/bin/env bash

# Gadi
# https://opus.nci.org.au/display/Help/How+to+submit+a+job
#PBS -N CABLE_merge
#PBS -P pr09
#PBS -q express
#PBS -l walltime=09:30:00
#PBS -l mem=16GB
#PBS -l ncpus=1
#PBS -l storage=gdata/x45+scratch/pt17+gdata/pr09
#PBS -l software=netCDF:MPI:Intel:GNU:scorep
#PBS -r y
#PBS -l wd
#PBS -j oe
#PBS -S /bin/bash
#PBS -v PYTHONPATH

python3 /home/599/jk8585/CABLE_run/TRENDY_v12/scripts/aux/merge_to_output2d.py -o /g/data/pr09/TRENDY_v12/S0/output/cru_out_casa_1901_2022.nc /g/data/pr09/TRENDY_v12/S0/run*/outputs/cru_out_casa_1901_2022.nc