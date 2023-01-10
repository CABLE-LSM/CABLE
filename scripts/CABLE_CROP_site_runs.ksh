#!/bin/ksh

# This script runs CABLE at multiple OzFLUX/FLUXNET sites
# Juergen Knauer, February/March 2019
#

#set -x # enable ksh debug output

# Global settings:
PROJECT_ID=OD-202890 # CSIROproject ID for slurm submission on petrichor
SITE_LIST=/datasets/work/oa-globalcable/work/Juergen/CABLE_run/CROPS/crop_sites_v2.txt
SITE_DIR=/home/lou059/sources/NESP2pt9_CROP
AUX_DIR=/datasets/work/oa-globalcable/work/Vanessa/CABLE-AUX
LOG_DIR=${SITE_DIR}/logs
CODE_DIR=/home/lou059/sources/NESP2pt9_CROP
FORCING_DIR=/datasets/work/oa-globalcable/work/BIOS3_forcing/site_met/FLUXNET2015/Crops
OBS_DIR=/datasets/work/oa-globalcable/work/Data_EC/OzFlux
PLOT_DIR=/datasets/work/oa-globalcable/work/CABLE_files/plots_R

# Settings
call_crop=TRUE
LAI_feedback=TRUE
RUN_ONLY=FALSE       # no spinup and transient (restart files must exist)
RESET_RESTART=TRUE   # reset restart files? (making sure the same files are used)?
                     # only used if RUN_ONLY=TRUE

# Experiment naming
if [ $call_crop = 'TRUE' ]; then
    EXP_NAME=crop
else
    EXP_NAME=default
fi

### no changes needed beyond that point ###

echo ${call_crop} ${LAI_feedback} ${SITE_DIR} ${OBS_DIR} ${PLOT_DIR} ${EXP_NAME} ${SITE_LIST} \
        > run_settings_crops.txt

# site names and number of sites
site=$(cut --fields=1 $SITE_LIST)
let nr_sites=$(wc --lines $SITE_LIST | awk '{print $1}')-1

# set current directory as base directory
BASE_DIR=$PWD

# define (and create) log folder
if [ ! -d $LOG_DIR ]; then
    mkdir $LOG_DIR
fi

for s in $site; do
    echo starting site $s
    # create working directory if not existing (copy from template directory)
    WD=./${s}
    if [ ! -d $WD ]; then
        mkdir $WD
        cd $WD

        # create links
        ln --symbolic $FORCING_DIR met  # maybe better below (met location could change)
        ln --symbolic $AUX_DIR ./CABLE-AUX
        #chmod -R 755 *  # permissions
    fi

    cd $WD
    # take crop params from code dir, the rest from basedir
    #cp ${CODE_DIR}/params/def_crop_params.txt ${BASE_DIR}/params/. # This is redundant
    cp -r ${BASE_DIR}/params .
    cp -r ${BASE_DIR}/namelists .
    cp namelists/cable.nml .
    cp namelists/site.nml .
    rm cable
    ln -s ${CODE_DIR}/offline/cable cable
    chmod +x cable
    cp ${BASE_DIR}/scripts/reset_restart.bash .
    cp ${BASE_DIR}/scripts/run_cable_site_CNP.py .
    cp ${BASE_DIR}/scripts/run_cable_site_CNP_meta.py .

    if [ $RUN_ONLY = "TRUE" ]; then
        # check for restart files
        if [ ! -e restart_files/${s}_CNP_${EXP_NAME}_cable_rst_transient.nc ]; then
            print "no existing restart files! Exiting..."
            exit 1
        else
            if [ $RESET_RESTART = "TRUE" ]; then
                sed -i "s!^file_ext=.*!file_ext=${s}_CNP_${EXP_NAME}!" reset_restart.bash
                ./reset_restart.bash
            fi
            sed -i "$ s!C.main(.*!C.main(SPIN_UP=False, TRANSIENT=False, SIMULATION=True)!" \
                    ${BASE_DIR}/scripts/run_cable_site_CNP_meta.py
        fi
    else
        sed -i "$ s!C.main(.*!C.main(SPIN_UP=True, TRANSIENT=True, SIMULATION=True)!" \
                ${BASE_DIR}/scripts/run_cable_site_CNP_meta.py
    fi

    cd $BASE_DIR
done

## extract start- and endyear from sitelist (alternative: from metfile)
#startyear=`awk -v site=${site} '$0~site {print $2}' $SITE_LIST`
#endyear=`awk -v site=${site} '$0~site {print $3}' $SITE_LIST`

#sed -i "s!^#SBATCH --ntasks-per-node.*!#SBATCH --ntasks-per-node=${nr_sites}!" run_cable_casa.slurm

## send to cluster
sbatch \
        --account=$PROJECT_ID \
        --array=1-${nr_sites} \
        --output=${LOG_DIR}/%x_%a.out \
        --error=${LOG_DIR}/%x_%a.err \
        scripts/run_cable_casa_crop_array.slurm

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

