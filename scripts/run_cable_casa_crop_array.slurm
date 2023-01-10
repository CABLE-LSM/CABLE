#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --mem=5GB
#SBATCH --job-name=CABLE_site

## command line arguments to 'run_cable_site_CNP_meta.py'
# 1: site name
# 2: start year
# 3: end year
# 4: lai_feedback? (TRUE or FALSE)
# 5: site directory
# 6: observations (EC data) directory
# 7: plot directory (location of plotting scripts)
# 8: finite gm? (TRUE or FALSE)
# 9: experiment name

module unload intel-cc
module unload intel-fc
module unload intel-mpi
module load netcdf
module load openmpi
module load hdf5
module load python
# If your system does not already have pandas, xarray, and matplotlib then activate your own env.
source ~/cable-env/bin/activate

read -a settings < run_settings_crops.txt
(( line = SLURM_ARRAY_TASK_ID + 1 ))               # to avoid header
read -a site < <(sed -n ${line}\p ${settings[6]})  # reading site list (=settings[5] (counting starts at 0))

python ./scripts/run_cable_site_CNP_meta.py ${site[*]:0:3} ${settings[*]:0:6}
