#!/usr/bin/env bash

# This script starts multiple sessions of the CABLE run script, each time using different land masks and output folders
# Before running this script, the following needs to be checked in the CABLE run script:
# PBS settings, Run sequence, CABLE settings 

# Script requires the following files
# landmask_script --> creates landmasks
# run_script      --> runs CABLE instances
# merge_script    --> merges CABLE outputs to latlon grid
# cleanup script  --> cleans up folder structure

#-------------------------------------------------------
# Modules
#-------------------------------------------------------
module purge
module load R/4.2.2
module load python3/3.10.4
module load proj/6.2.1
module load gdal/3.0.2
module load geos/3.8.0
module load intel-compiler/2021.8.0
export R_LIBS=/g/data/x45/R/libs


#-------------------------------------------------------
# Settings
#-------------------------------------------------------
experiment="S3"
experiment_name="${experiment}_test8_serial_global"
merge_results=1   # after runs are finished, merge results into one folder and backup 
                  # restart, logs, landmasks etc. (1) or keep folder structure as it is (0).
                  # The latter is useful if runs are to be resumed from restart files. 
mergesteps="zero_biomass spinup_nutrient_limited_1 spinup_nutrient_limited_2 1700_1900 1901_2021"   # sub-steps to be merged

### Spatial subruns ###
create_landmasks=1           # create new landmask files (1) or use existing ones (0)?
nruns=100                    # number of runs in parallel
#extent="0.0,4.0,47.0,50.0"  # "global" or "lon_min,lon_max,lat_min,lat_max"
#extent="0.0,2.0,44.0,46.0"
extent="global"
climate_restart="cru_climate_rst"       # name of climate restart file (without file extension)
keep_dump=1                             # keep dump files (1) or discard (0)? they are always kept for LUC runs


### Directories and files###
# Output directory
outpath="/scratch/pt17/jk8585/CABLE_tests/${experiment_name}"
# Code directory
cablecode="/home/599/jk8585/CABLE_code/SHARE/CABLE-POP_TRENDY"
# Script directory
cablehome="/home/599/jk8585/CABLE_run/TRENDY_v11"
# Scripts
landmask_script="${cablehome}/scripts/serial/split_landmask.R"
run_script="${cablehome}/scripts/serial/run_cable.sh"
merge_script="${cablehome}/scripts/serial/merge_outputs.sh"
cleanup_script="${cablehome}/scripts/serial/cleanup.sh"
# Cable executable
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


## ---------------------------- End Settings ---------------------------------- ## 

# -----------------------------------------------------------------------
# 1) Create landmasks (created in folders ${outpath}/runX/landmask)
# -----------------------------------------------------------------------
if [[ ${create_landmasks} -eq 1 ]] ; then
   $landmask_script $GlobalLandMaskFile $nruns $outpath $extent
fi

# -----------------------------------------------------------------------
# 2) Run CABLE
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

    RUN_IDS="${RUN_IDS}:$(qsub $run_script)"
done

# -----------------------------------------------------------------------
# 3) Merge outputs (only if all previous runs were OK)
# -----------------------------------------------------------------------
if [[ ${merge_results} -eq 1 ]] ; then
    
    ftypes="cable casa LUC"
    outfinal="${outpath}/output"
    if [[ -d ${outfinal} ]] ; then
        rm -r ${outfinal}
    fi
    mkdir -p ${outfinal}

    for mergestep in ${mergesteps} ; do
        for ftype in ${ftypes} ; do
            if [[ ("${ftype}" != "LUC") || ("${experiment}" == "S3" && ("${mergestep}" == "1700_1900" || "${mergestep}" == "1901_"* )) ]] ; then
                sed -i "s!^python3.*!python3 merge_to_output2d.py -o ${outfinal}/cru_out_${ftype}_${mergestep}.nc ${outpath}/run*/outputs/cru_out_${ftype}_${mergestep}.nc!" $merge_script
                MERGE_IDS="${MERGE_IDS}:$(qsub -W "depend=afterok${RUN_IDS}" $merge_script)"
            fi
        done
    done
fi

# -----------------------------------------------------
# 4) Backup and cleanup
# -----------------------------------------------------
sed -i "s!^exp_name=.*!exp_name='${experiment_name}'!" $cleanup_script
sed -i "s!^outpath=.*!outpath='${outpath}'!" $cleanup_script
sed -i "s!^nruns=.*!nruns=${nruns}!" $cleanup_script
sed -i "s!^climate_restart=.*!climate_restart='${climate_restart}'!" $cleanup_script
sed -i "s!^keep_dump=.*!keep_dump=${keep_dump}!" $cleanup_script
sed -i "s!^mergesteps=.*!mergesteps='${mergesteps}'!" $cleanup_script

qsub -W "depend=afterok${MERGE_IDS}" $cleanup_script
