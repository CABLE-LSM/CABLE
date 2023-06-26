#!/usr/bin/env bash

# Gadi
# https://opus.nci.org.au/display/Help/How+to+submit+a+job
#PBS -N CABLE_merge
#PBS -P x45
#PBS -q express
#PBS -l walltime=06:30:00
#PBS -l mem=16GB
#PBS -l ncpus=1
#PBS -l storage=gdata/x45+scratch/pt17+gdata/vl59
#PBS -l software=netCDF:MPI:Intel:GNU:scorep
#PBS -r y
#PBS -l wd
#PBS -j oe
#PBS -S /bin/bash

python3 merge_to_output2d.py -o /scratch/pt17/jk8585/CABLE_tests/S3_test8_serial_global/output/cru_out_LUC_1901_2021.nc /scratch/pt17/jk8585/CABLE_tests/S3_test8_serial_global/run*/outputs/cru_out_LUC_1901_2021.nc