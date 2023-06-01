#!/usr/bin/env bash

# Gadi
# https://opus.nci.org.au/display/Help/How+to+submit+a+job
#PBS -N CABLE_test_merge
#PBS -P x45
#PBS -q express
#PBS -l walltime=04:30:00
#PBS -l mem=128GB
#PBS -l ncpus=1
#PBS -l storage=gdata/x45+scratch/pt17+gdata/vl59
#PBS -l software=netCDF:MPI:Intel:GNU:scorep
#PBS -r y
#PBS -l wd
#PBS -j oe
#PBS -S /bin/bash

## TODO: add script description here!

# -------------------------------------------------------------------------------
# Settings
# -------------------------------------------------------------------------------
# Location of CABLE output files
#outpath=$1
outpath="/scratch/pt17/jk8585/CABLE_tests/S3_test_serial_global"
exp_name=$(basename ${outpath})

basename_cable="cru_out_cable"
basename_casa="cru_out_casa"
basename_LUC="cru_out_LUC"

# run steps to process
#runsteps="zero_biomass spinup_limit_labile_1 spinup_limit_labile_2 spinup_limit_labile_3
#          spinup_nutrient_limited_1 spinup_nutrient_limited_2 spinup_nutrient_limited_3
#          1700_1900 1901_2021"
runsteps="zero_biomass spinup_limit_labile_1 spinup_nutrient_limited_1 1700_1900"

#runsteps="zero_biomass"

# check if LUC output needs to be processed
if [[ "${outpath}" == "S3*" && ("${runsteps}" == "*1700_1900" || "${runsteps}" == "1901_*")  ]] ; then
    process_LUC=1
else
    process_LUC=0
fi

climate_restart="cru_climate_rst"       # name of climate restart file (without file extension)
keep_dump=1                             # keep dump files (1) or discard (0)? they are always kept for LUC runs


# -------------------------------------------------------------------------------
# prepare folder structure and determine number of runs
# -------------------------------------------------------------------------------
outfinal="${outpath}/output"
if [[ -d ${outfinal} ]] ; then
    rm -r ${outfinal}
fi
nruns="$(find ${outpath} -maxdepth 1 -mindepth 1 -name "run*" -type d | wc -l)"
mkdir -p ${outfinal}


# -------------------------------------------------------------------------------
# Merge outputs
# -------------------------------------------------------------------------------

cd ${outfinal}

for runstep in ${runsteps} ; do

    ## merge casa output files
    files="$(find ${outpath} -name ${basename_casa}_${runstep}.nc)"
    for file in ${files} ; do
        echo $file
        if [[ ! -f tmp.nc ]] ; then 
            echo FIRST $file
            cp $file tmp.nc
            continue
        fi
        cdo -O -z zip_4 mergegrid $file tmp.nc tmp1.nc 
        mv tmp1.nc tmp.nc
    done
    mv tmp.nc ${outfinal}/${basename_casa}_${runstep}.nc
    ## zip end product, not intermediates
    
exit 1
    ## merge cable output files
    files="$(find ${outpath} -name ${basename_cable}_${runstep}.nc)"
    cdo -O -z zip_4 mergegrid $files ${outfinal}/${basename_cable}_${runstep}.nc 

    ## if LUC is simulated, merge those files as well
    if [[ $process_LUC -eq 1 ]] ; then
        files="$(find ${outpath} -name ${basename_LUC}_${runstep}.nc)"
        cdo -O -z zip_4 mergegrid $files ${outfinal}/${basename_LUC}_${runstep}.nc
    fi

done


# -------------------------------------------------------------------------------
# Back up the rest and clean up
# -------------------------------------------------------------------------------

# executable
mkdir -p ${outpath}/exe
mv ${outpath}/run1/cable ${outpath}/exe/.

for ((irun=1; irun<=${nruns}; irun++)) ; do

    # PBS log files
    mkdir -p ${outpath}/PBS
    cp ${PWD}/${exp_name}.o* ${outpath}/PBS/.

    # Restart files (includes namelists)
    mkdir -p ${outpath}/restart/run${irun}
    mv ${outpath}/run${irun}/restart/* ${outpath}/restart/run${irun}/.

    # Climate restart
    mkdir -p ${outpath}/climate_restart
    mv ${outpath}/run${irun}/${climate_restart}.nc ${outpath}/climate_restart/${climate_restart}${irun}.nc

    # log files
    mkdir -p ${outpath}/logs/run${irun}
    mv ${outpath}/run${irun}/logs/* ${outpath}/log/run${irun}/.

    # landmasks
    mkdir -p ${outpath}/landmasks
    mv ${outpath}/run${irun}/landmask/landmask${irun}.nc ${outpath}/landmasks/.

    # dump files
    if [[ ${keep_dump} -eq 1 || ${process_LUC} -eq 1 ]] ; then
        mkdir -p ${outpath}/dump_files/run${irun}
        mv ${outpath}/run${irun}/*_dump.nc ${outpath}/dump_files/run${irun}/.
    fi

    # delete all the rest
    #rm -r ${outpath}/run${irun}/
done

# finished script!