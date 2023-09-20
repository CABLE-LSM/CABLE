#!/usr/bin/env bash

# This script starts multiple sessions of the CABLE run script, each time using different land masks and output folders
# Before running this script, the following needs to be checked in the CABLE run script:
# PBS settings, Run sequence, CABLE settings 

# Script requires the following files
# landmask_script --> creates landmasks
# run_script      --> runs CABLE instances
#     - requires run_cable-pop_lib.sh
# merge_script    --> merges CABLE outputs to latlon grid
#     - requires merge_to_output2d.py 
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
export PYTHONPATH=/g/data/x45/python3.10.4/lib/python3.10/site-packages


#-------------------------------------------------------
# Settings
#-------------------------------------------------------
experiment="S0"
experiment_name="${experiment}"
run_model=1       # run the model or just do other steps (e.g. merging)?
merge_results=1   # after runs are finished, merge results into one folder and backup 
                  # restart, logs, landmasks etc. (1) or keep folder structure as it is (0).
                  # The latter is useful if runs are to be resumed from restart files. 
#mergesteps="zero_biomass spinup_nutrient_limited_1 spinup_nutrient_limited_2 1700_1900 1901_2022"   # sub-steps to be merged
mergesteps="1700_1900 1901_2022"

### Spatial subruns ###
create_landmasks=1               # create new landmask files (1) or use existing ones (0)?
nruns=100                        # number of runs in parallel
#extent="64.0,66.0,60.0,62.0"    # "global" or "lon_min,lon_max,lat_min,lat_max"
extent="global"
climate_restart="cru_climate_rst"       # name of climate restart file (without file extension)
keep_dump=1                             # keep dump files (1) or discard (0)? They are always kept for LUC runs


### Directories and files###
# Output directory
outpath="/g/data/pr09/TRENDY_v12/${experiment_name}"
# Code directory
cablecode="/home/599/jk8585/CABLE_code/SHARE/CABLE-POP_TRENDY"
# Run directory
rundir="/home/599/jk8585/CABLE_run/TRENDY_v12"
# Scripts
landmask_script="${rundir}/scripts/split_landmask.R"
run_script="${rundir}/scripts/run_cable.sh"
merge_script="${rundir}/scripts/merge_outputs.sh"
cleanup_script="${rundir}/scripts/cleanup.sh"
# Cable executable
exe="${cablecode}/offline/cable"
# CABLE-AUX directory (uses offline/gridinfo_CSIRO_1x1.nc and offline/modis_phenology_csiro.txt)
aux="/g/data/pr09/TRENDY_v12/aux"
# Global Meteorology
GlobalMetPath="/g/data/x45/CRUJRA2023/daily_1deg"
MetVersion="CRUJRA_2023"
# Global LUC
GlobalTransitionFilePath="/g/data/x45/LUH2/GCB_2023/1deg/EXTRACT"
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
if [[ ${run_model} -eq 1 ]] ; then
    # 2.1) Write general settings into run script
    sed -i "s!^#PBS -N.*!#PBS -N TRENDY_${experiment_name}!" $run_script
    sed -i "s!^experiment=.*!experiment='${experiment}'!" $run_script
    sed -i "s!^experiment_name=.*!experiment_name='${experiment_name}'!" $run_script
    sed -i "s!^cablecode=.*!cablecode='${cablecode}'!" $run_script
    sed -i "s!^rundir=.*!rundir='${rundir}'!" $run_script
    sed -i "s!^exe=.*!exe='${exe}'!" $run_script
    sed -i "s!^aux=.*!aux='${aux}'!" $run_script
    sed -i "s!^MetPath=.*!MetPath='${GlobalMetPath}'!" $run_script
    sed -i "s!^MetVersion=.*!MetVersion='${MetVersion}'!" $run_script
    sed -i "s!^TransitionFilePath=.*!TransitionFilePath='${GlobalTransitionFilePath}'!" $run_script
    sed -i "s!^SurfaceFile=.*!SurfaceFile='${SurfaceFile}'!" $run_script

    # 2.2) Loop over landmasks and start runs
    for ((irun=1; irun<=${nruns}; irun++)) ; do
        runpath="${outpath}/run${irun}"
        sed -i "s!^runpath=.*!runpath='${runpath}'!" $run_script
        sed -i "s!^LandMaskFile=.*!LandMaskFile='${runpath}/landmask/landmask${irun}.nc'!" $run_script

        RUN_IDS="${RUN_IDS}:$(qsub $run_script)"
    done
fi

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
                sed -i "s!^python3.*!python3 ${rundir}/scripts/aux/merge_to_output2d.py -o ${outfinal}/cru_out_${ftype}_${mergestep}.nc ${outpath}/run*/outputs/cru_out_${ftype}_${mergestep}.nc!" $merge_script
                if [[ ${run_model} -eq 1 ]] ; then
                    MERGE_IDS="${MERGE_IDS}:$(qsub -W "depend=afterok${RUN_IDS}" $merge_script)"
                else
                    MERGE_IDS="${MERGE_IDS}:$(qsub $merge_script)"
                fi
            fi
        done
    done

    # -----------------------------------------------------
    # 4) Backup and cleanup
    # -----------------------------------------------------

    sed -i "s!^exp_name=.*!exp_name='TRENDY_${experiment_name}'!" $cleanup_script
    sed -i "s!^outpath=.*!outpath='${outpath}'!" $cleanup_script
    sed -i "s!^nruns=.*!nruns=${nruns}!" $cleanup_script
    sed -i "s!^climate_restart=.*!climate_restart='${climate_restart}'!" $cleanup_script
    sed -i "s!^keep_dump=.*!keep_dump=${keep_dump}!" $cleanup_script
    sed -i "s!^mergesteps=.*!mergesteps='${mergesteps}'!" $cleanup_script

    qsub -W "depend=afterok${MERGE_IDS}" $cleanup_script

fi