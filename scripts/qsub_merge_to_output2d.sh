#!/usr/bin/env bash

# Gadi
# https://opus.nci.org.au/display/Help/How+to+submit+a+job
#PBS -N merge2d
#PBS -P x45
# express / normal / copyq (2x24, cascadelake)
# expressbw / normalbw (2x14, broadwell) / normalsl (2x16, skylake)- ex-Raijin nodes
#PBS -q normal
# casa < 1 h, cable < 2 h, sum_patchfrac < 25 min
#PBS -l walltime=04:00:00
#PBS -l mem=4GB
#PBS -l ncpus=1
# #PBS -l jobfs=1GB
#PBS -l storage=scratch/x45+gdata/x45
#PBS -l software=netCDF:python
#PBS -r y
#PBS -l wd
#PBS -j oe
#PBS -S /bin/bash
#PBS -M matthias.cuntz@inrae.fr

set -e

# setup
rundir=/scratch/x45/jk8585/S3_test_serial2
outdir=/scratch/x45/moc801
# gives files  cru_out_${otype}_${oend}.nc  such as  cru_out_cable_1901_2021.nc
otypes="casa cable"
oend="1901_2021"

# Python
module load python3/3.10.4
module load python3-as-python
mypython=/g/data/x45/python3.10.4
if [[ -n ${PYTHONPATH} ]] ; then
    export PYTHONPATH=${PYTHONPATH}:${mypython}/lib/python3.10/site-packages/
else
    export PYTHONPATH=${mypython}/lib/python3.10/site-packages/
fi
export PATH=${PATH}:${mypython}/bin

# runs
nrun=$(find ${rundir} -maxdepth 1 -type d -name run\* | wc -l)

for otype in ${otypes} ; do
    echo ""
    echo ${otype}

    file0="${rundir}/run[1-9]/outputs/cru_out_${otype}_${oend}.nc"
    file9="${rundir}/run[1-9][0-9]/outputs/cru_out_${otype}_${oend}.nc"
    file99="${rundir}/run[1-9][0-9][0-9]/outputs/cru_out_${otype}_${oend}.nc"
    file999="${rundir}/run[1-9][0-9][0-9][0-9]/outputs/cru_out_${otype}_${oend}.nc"

    if [[ ${nrun} -gt 999 ]] ; then
        ifiles="${file0} ${file9} ${file99} ${file999}"
    elif [[ ${nrun} -gt 99 ]] ; then
        ifiles="${file0} ${file9} ${file99}"
    elif [[ ${nrun} -gt 9 ]] ; then
        ifiles="${file0} ${file9}"
    else
        ifiles="${file0}"
    fi

    # merge
    ofile="${outdir}/cru_out_${otype}_${oend}-merged1.nc"
    python merge_to_output2d.py -z -v -o ${ofile} ${ifiles}

    # sum patch
    if [[ "${otype}" == "cable" ]] ; then
        pfile="${outdir}/cru_out_${otype}_${oend}-sum_patch1.nc"
        python sum_patchfrac.py -v -z -o ${pfile} ${ofile}
    fi
done
