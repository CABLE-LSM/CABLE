#!/bin/bash
#------------------------------------------------------------------------------
# Purpose:
# Retrieves input files required for CABLE offline test runs, and fills in the
# namefile template.
#
# This script is used for testing, and should be modified to point to updated
# input files as required in branches. Once merged, files in CABLE-AUX should
# be updated, and this file reverted to the default.
#
# Contact: ned haughton <ned@nedhaughton.com>
#------------------------------------------------------------------------------

BASENAME=$0

function usage {
    echo "usage: $BASENAME [offline|offline-mpi|UM] <out-dir>"
exit 1
}

[ 2 -gt "$#" ] && { usage; }

MODE=$1
OUTDIR=$2

mkdir -p ${OUTDIR}

CABLE_AUX="https://trac.nci.org.au/svn/cable/branches/Share/CABLE-AUX"

# Modify these lines as necessary. Use -r REV if you need a specific version.
# See `svn help export` for more detail.

# Required for all modes
svn export --force ${CABLE_AUX}/core/biogeophys/def_veg_params_zr_clitt_albedo_fix.txt ${OUTDIR}/def_veg_params.txt
svn export --force ${CABLE_AUX}/core/biogeophys/def_soil_params.txt ${OUTDIR}/def_soil_params.txt

if [[ $MODE == offline* ]]; then
    # Files required for offline modes
    svn export --force ${CABLE_AUX}/offline/cable.nml ${OUTDIR}/cable.nml
    svn export --force ${CABLE_AUX}/offline/TumbaFluxnet.1.3_met.nc ${OUTDIR}/TumbaFluxnet.1.3_met.nc
    svn export --force ${CABLE_AUX}/offline/gridinfo_CSIRO_1x1.nc ${OUTDIR}/gridinfo_CSIRO_1x1.nc

    # Change the namelist to match the downloaded files - modify as required.
    sed -e "s/\(filename%met *= '\).*'/\1TumbaFluxnet.1.3_met.nc'/" \
        -e "s/\(filename%type *= '\).*'/\1gridinfo_CSIRO_1x1.nc'/" \
        -e "s/\(filename%veg *= '\).*'/\1def_veg_params.txt'/" \
        -e "s/\(filename%soil *= '\).*'/\1def_soil_params.txt'/" \
        -i ${OUTDIR}/cable.nml

fi  # offline*
