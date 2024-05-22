#!/usr/bin/env bash

# biocomp
#SBATCH --ignore-pbs
#SBATCH --job-name=llmeteo
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.out
# #SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --mail-type=FAIL,STAGE_OUT,TIME_LIMIT,INVALID_DEPEND,END
#SBATCH --mail-user=matthias.cuntz@inrae.fr

# This script extracts one lat,lon point from meteo, land use, and mask files.

set -e

# --------------------------------------------------------------------
# Load Modules
#
eval "$(${HOME}/miniconda3/bin/conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate pystd


# --------------------------------------------------------------------
# lat,lon and directories
#

# lat,lon  or  latmin,latmax,lonmin,lonmax
# must have . in numbers otherwise indexes are taken
# FR-Hes
latlon=48.6742166667,7.06461666667
# latlon=-34.5,-33.5,149.5,156.5
# latlon=42.5,43.5,109.5,110.5
# latlon=-44.0,-10.0,110.0,155.0  # Australia

# Global Meteorology
GlobalMetPath='/home/mcuntz/data/met_forcing/CRUJRA2023/daily_1deg'
# Global LUC
GlobalTransitionFilePath='/home/mcuntz/data/cable/LUH2/GCB_2023/1deg/EXTRACT'
# Global Landmask file
GlobalLandMaskFile="/home/mcuntz/data/cable/ipbes/masks/glob_ipsl_1x1.nc"

# Experiment directory
expdir='/home/mcuntz/projects/cable/sites/FR-Hes'
# Local Meteorology
MetPath="${expdir}/met"
# Local LUC
TransitionFilePath="${expdir}/luh"
# Local Mask
LandMaskFile="${expdir}/mask/landmask_FR-Hes.nc"


# --------------------------------------------------------------------
# Setup
#

set -e

trap cleanup 1 2 3 6

pid=$$
isdir="${PWD}"
prog=$0
pprog=$(basename ${prog})
pdir=$(dirname ${prog})
tmp=${TMPDIR:-"/tmp"}

cd ${isdir}

source ${isdir}/run_cable-pop_lib.sh

# usage of script
function usage()
{
    printf "${pprog} [-h]\n"
    printf "Extract lat,lon from meteo, land use, and mask files.\n"
    printf "\n"
    printf "Options\n"
    printf "    -h    Prints this help screen.\n"
}

# cleanup at end or at trap
function cleanup()
{
    \rm -f ${tmp}/*.${pid}*
    exit 1
}

# Get options
while getopts "h" option ; do
    case ${option} in
        h) usage
	   exit
	   ;;
        *) printf "Error ${pprog}: unimplemented option.\n\n" 1>&2
	   usage 1>&2
	   exit 1
	   ;;
    esac
done
shift $((${OPTIND} - 1))


# --------------------------------------------------------------------
# Extract
#
t1=$(date +%s)
printf "Started at %s\n" "$(date)"

# Extract meteo, land use and mask for one specific site from global files
echo "Extract local meteo and mask"

# meteorology
met_list="dlwrf fd ndep pre pres spfh tmax tmin tswrf ugrd vgrd"
mkdir -p ${MetPath}
MetPath=$(abspath ${MetPath})
for mm in ${met_list} ; do
    echo ${GlobalMetPath}/${mm}
    mkdir -p ${MetPath}/${mm}
    for nc in ${GlobalMetPath}/${mm}/*.nc ; do
        ff=$(basename ${nc})
        echo "    ${ff}"
        ncks -O $(nckslatlon ${nc} ${latlon}) ${nc} ${MetPath}/${mm}/${ff}
    done
done
echo ${GlobalMetPath}/co2/*.txt
mkdir -p ${MetPath}/co2
cp ${GlobalMetPath}/co2/*.txt ${MetPath}/co2/

# land use
mkdir -p ${TransitionFilePath}
TransitionFilePath=$(abspath ${TransitionFilePath})
echo ${GlobalTransitionFilePath}
for nc in ${GlobalTransitionFilePath}/*.nc ; do
    ff=$(basename ${nc})
    echo "    ${ff}"
    ncks -O $(nckslatlon ${nc} ${latlon}) ${nc} ${TransitionFilePath}/${ff}
done

# mask
LandMaskFilePath=$(dirname ${LandMaskFile})
mkdir -p ${LandMaskFilePath}
LandMaskFile=$(absfile ${LandMaskFile})
echo $(basename ${LandMaskFile})
ncks -O $(nckslatlon ${GlobalLandMaskFile} ${latlon}) ${GlobalLandMaskFile} ${LandMaskFile}


# --------------------------------------------------------------------
# Finish
#
cd ${isdir}

t2=$(date +%s)
dt=$((t2-t1))
printf "\n"
if [[ ${dt} -lt 60 ]] ; then
    printf "Finished at %s   in %i seconds.\n" "$(date)" ${dt}
else
    dm=$(echo "(${t2}-${t1})/60." | bc -l)
    printf "Finished at %s   in %.2f minutes.\n" "$(date)" ${dm}
fi

exit
