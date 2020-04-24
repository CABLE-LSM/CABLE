#!/bin/bash
#
# compare_two_cable_runs.sh [-h] dir1 dir2
# Compare restart and output files of two Cable-POP runs performed with script run_cable-pop.sh.
#
# Input
#     dir1        First directory with Cable-POP output, containing subdirectories restart/ and outputs/.
#     dir2        Second directory with Cable-POP output.
#
# Options
#     -h          Prints this help screen.
#
# Example
#     compare_two_cable_runs.sh  dir1 dir2
#
# Copyright (c) 2020 Matthias Cuntz - mc (at) macu (dot) de
#
#
# Example
# -------
#   ./compare_two_cable_runs.sh ~/prog/cable/HarvardForest/run_20200202 ~/prog/cable/HarvardForest/run_20200216
#
#
# History
# -------
# Written  Matthias Cuntz, Mar 2020
# Modified Matthias Cuntz, Mar 2020 - Changed check with ncdump to be correct on more than two land points
#          Matthias Cuntz, Apr 2020 - Check Restart and Output files
#
set -e

prog=$0
pprog=$(basename ${prog})
dprog=$(dirname ${prog})
isdir=${PWD}
pid=$$

# -------------------------------------------------------------------------------------------------
# Functions
#
function usage () {
    printf "${pprog} [-h] dir1 dir2\n"
    printf "Compare restart and output files of two Cable-POP runs performed with script run_cable-pop.sh.\n"
    printf "\n"
    printf "Input\n"
    printf "    dir1        First directory with Cable-POP output, containing subdirectories restart/ and outputs/.\n"
    printf "    dir2        Second directory with Cable-POP output.\n"
    printf "\n"
    printf "Options\n"
    printf "    -h          Prints this help screen.\n"
    printf "\n"
    printf "Example\n"
    printf "    ${pprog}  dir1 dir2\n"
    printf "\n"
    printf "Copyright (c) 2020 Matthias Cuntz - mc (at) macu (dot) de\n"
}

# absolute path
function abspath()
{
    idir=${PWD}
    cd ${1}
    odir=${PWD}
    cd ${idir}
    echo "${odir}"
}

# -------------------------------------------------------------------------------------------------
# Analyse command line
#

# Switches
while getopts "h" Option ; do
    case ${Option} in
        h) usage; exit;;
        *) printf "Error ${pprog}: unimplemented option.\n\n" 1>&2;  usage 1>&2; exit 1;;
    esac
done
shift $((${OPTIND} - 1))

# Check args
if [[ $# -lt 2 ]] ; then 
    printf "Error ${pprog}: two Cable-POP output directories must be given.\n\n" 1>&2
    usage 1>&2
    exit 1
elif [[ $# -gt 2 ]] ; then 
    printf "Error ${pprog}: only two Cable-POP output directories can be given.\n\n" 1>&2
    usage 1>&2
    exit 1
fi

# Input directories
dir1=${1}
dir2=${2}
iout="restart outputs"
for idir in ${dir1} ${dir2} ; do
    if [[ ! -d ${idir} ]] ; then
	printf "Error ${pprog}: Input is not a directory: ${idir}.\n\n" 1>&2
	usage 1>&2
	exit 1
    fi
    for iidir in ${iout} ; do
	if [[ ! -d ${idir}/${iidir} ]] ; then
	    printf "Error ${pprog}: There is no subdirectory ${iidir} in input directory ${idir}.\n\n" 1>&2
	    usage 1>&2
	    exit 1
	fi
    done
done
adir1=$(abspath ${dir1})
adir2=$(abspath ${dir2})

# -------------------------------------------------------------------------------------------------
# Compare restart/output files
#

# Explor
if [[ -d ${HOME}/local/bin ]] ; then export PATH=${PATH}:${HOME}/local/bin ; fi

# Other clusters
set +e
xcdo=$(which cdo 2> /dev/null)
if [[ -z ${xcdo} ]] ; then module load cdo ; fi
xnco=$(which ncks 2> /dev/null)
if [[ -z ${xnco} ]] ; then module load nco ; fi
xncdump=$(which ncdump 2> /dev/null)
if [[ -z ${xncdump} ]] ; then module load netcdf ; fi
set -e

# File suffixes to check
suff1="_climate_restart"
suff2="_zero_biomass _spinup_limit_labile _spinup_analytic_limit_labile _spinup _spinup_analytic _1580_1699 _1700_1899 _1900_2017"

echo ""
echo ""
echo "Compare ${adir1} vs. ${adir2}"

if [[ true ]] ; then
    for idir in ${iout} ; do
        echo ""
        echo "cdo diff in ${idir}"
        cd ${adir1}/${idir}/
        for suff in ${suff1} ${suff2} ; do
            echo '  Step ' ${suff}
            for i in *${suff}.nc ; do
		if [[ -f ${i} ]] ; then
                    echo '    ' ${i}
                    set +e
                    cdo -s diffv ${i} ${adir2}/${idir}/${i} 2> /dev/null
                    set -e
		fi
            done
        done
    done
fi

if [[ true ]] ; then
    idir=restart
    echo ""
    echo "ncdiff pop_cru_ini in ${idir}"
    cd ${adir1}/${idir}/
    for suff in ${suff2} ; do
        i=pop_cru_ini${suff}.nc
	if [[ -f ${i} ]] ; then
            echo '    ' ${i}
            ncdiff ${i} ${adir2}/${idir}/${i} tmp.${pid}.nc
            set +e
            iout=$(ncdump tmp.${pid}.nc | sed -e '/^netcdf/,/^variables:/d' -e '/ latitude/,/;$/d' -e '/ longitude/,/;$/d' -e '/[:{}]/d' -e '/^$/d' -e '/) ;$/d' -e 's/ ;/,/' -e 's/ 0,//g' | sed -e '/^[[:blank:]]*$/d' -e '/=[[:blank:]]*$/d' -e '/[tT]ime =/d')
            set -e
            if [[ -n ${iout} ]] ; then
		echo "     Check ${i} in ${PWD}"
		echo "         ncdiff -O ${i} ${adir2}/${idir}/${i} tmp.nc ; ncdump tmp.nc | sed -e '/^netcdf/,/^variables:/d' -e '/ latitude/,/;$/d' -e '/ longitude/,/;$/d' -e '/[:{}]/d' -e '/^$/d' -e '/) ;$/d' -e 's/ ;/,/' -e 's/ 0,//g' | sed -e '/^[[:blank:]]*$/d' -e '/=[[:blank:]]*$/d' -e '/[tT]ime =/d'"
            fi
            rm tmp.${pid}.nc
	fi
    done
fi

cd ${isdir}

exit
