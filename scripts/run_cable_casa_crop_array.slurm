#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --mem=5GB
#SBATCH --job-name=CABLE_site

module del intel-cc intel-fc
module add intel-cc/16.0.1.150 intel-fc/16.0.1.150
module unload intel-mpi/5.0.1.035
module add netcdf/4.3.3.1 openmpi/1.8.8

read -a settings < run_settings_crops.txt  
(( line = SLURM_ARRAY_TASK_ID + 1 ))               # to avoid header
read -a site < <(sed -n ${line}\p ${settings[6]})  # reading site list (=settings[5] (counting starts at 0))

python ./run_cable_site_CNP_meta.py ${site[*]:0:3} ${settings[*]:0:6}
