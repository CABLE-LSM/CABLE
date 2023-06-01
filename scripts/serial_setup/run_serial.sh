#!/usr/bin/env bash

# This script starts multiple sessions of the CABLE run script, each time using different land masks and output folders
# Before running this script, the following needs to be checked in the CABLE run script:
# PBS settings, Run sequence, CABLE settings 

# Requires split_landmask.R and assumes it is the same directory



#-------------------------------------------------------
# Modules
#-------------------------------------------------------
module purge
module load R/4.2.2
module load proj/6.2.1
module load gdal/3.0.2
module load geos/3.8.0
module load intel-compiler/2021.8.0
export R_LIBS=/g/data/x45/R/libs


#-------------------------------------------------------
# Settings
#-------------------------------------------------------
experiment="S3"
experiment_name="${experiment}_test_serial"
merge_results=1   # after runs are finished, merge results into one folder and backup 
                  # restart, logs, landmasks etc. (1) or keep folder structure as it is (0).
                  # The latter is useful if runs are to be resumed from restart files. 

### Spatial subruns ###
create_landmasks=1           # create new landmask files (1) or use existing ones (0)?
nruns=2                      # number of runs in parallel
extent="0.0,9.0,44.0,50.0"  # "global" or "xmin,xmax,ymin,ymax"
#extent="0.0,2.0,44.0,46.0"
#extent="global"

### Directories ###
# Output directory
outpath="/scratch/pt17/jk8585/CABLE_tests/${experiment_name}"
# Code directory
cablecode="/home/599/jk8585/CABLE_code/SHARE/CABLE-POP_TRENDY"
# Script directory
cablehome="/home/599/jk8585/CABLE_run/TRENDY_v11"
# Cable run script
run_script="${cablehome}/scripts/serial/run_cable.sh"
# Cable executable
#exe="${cablecode}/offline/cable"
exe="${cablecode}/offline/cable"
# CABLE-AUX directory (uses offline/gridinfo_CSIRO_1x1.nc and offline/modis_phenology_csiro.txt)
aux="/g/data/vl59/TRENDY_v11/aux"
# Global Meteorology
GlobalMetPath="/g/data/x45/CRUJRA2022/daily_1deg"
# Global LUC
GlobalTransitionFilePath="/g/data/x45/LUH2/GCB_2022/1deg/EXTRACT"
# Global Surface file 
SurfaceFile="${aux}/gridinfo_CSIRO_1x1.nc"   
# Global Land Mask
GlobalLandMaskFile="${aux}/landmasks/glob_ipsl_1x1.nc"


# -----------------------------------------------------------------------
# 1) Create landmasks (created in folders ${outpath}/runX/landmask)
# -----------------------------------------------------------------------
if [[ ${create_landmasks} -eq 1 ]] ; then
   ./split_landmask.R $GlobalLandMaskFile $nruns $outpath $extent
fi

# -----------------------------------------------------------------------
# 2) Start runs
# -----------------------------------------------------------------------
# 2.1) Write general settings into run script
sed -i "s!^#PBS -N.*!#PBS -N ${experiment_name}!" $run_script
sed -i "s!^experiment=.*!experiment='${experiment}'!" $run_script
sed -i "s!^experiment_name=.*!experiment_name='${experiment_name}'!" $run_script
sed -i "s!^cablecode=.*!cablecode='${cablecode}'!" $run_script
sed -i "s!^cablehome=.*!cablehome='${cablehome}'!" $run_script
sed -i "s!^exe=.*!exe='${exe}'!" $run_script
sed -i "s!^aux=.*!aux='${aux}'!" $run_script
sed -i "s!^MetPath=.*!MetPath='${GlobalMetPath}'!" $run_script
sed -i "s!^TransitionFilePath=.*!TransitionFilePath='${GlobalTransitionFilePath}'!" $run_script
sed -i "s!^SurfaceFile=.*!SurfaceFile='${SurfaceFile}'!" $run_script

# 2.2) Loop over landmasks and start runs
for ((irun=1; irun<=${nruns}; irun++)) ; do
    runpath="${outpath}/run${irun}"
    sed -i "s!^runpath=.*!runpath='${runpath}'!" $run_script
    sed -i "s!^LandMaskFile=.*!LandMaskFile='${runpath}/landmask/landmask${irun}.nc'!" $run_script

    JOB_IDS="${JOB_IDS}:$(qsub $run_script)"
done

echo $JOB_IDS
# -----------------------------------------------------------------------
# 3) Merge outputs (only if all previous runs were OK)
# -----------------------------------------------------------------------
#if [[ ${merge_results} -eq 1 ]] ; then
#    qsub -W "depend=afterok:${JOB_IDS}" --merge_outputs.sh $outpath
#fi 





