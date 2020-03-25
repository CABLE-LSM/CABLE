#!/bin/bash
#
# compare_two_cable_runs.sh [-h] dir1 dir2
# Compare restart files of two Cable-POP runs performed with script run_cable-pop.sh.
#
# Input
#     dir1        First directory with Cable-POP output, containing subdirectory restart/.
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
    printf "Compare restart files of two Cable-POP runs performed with script run_cable-pop.sh.\n"
    printf "\n"
    printf "Input\n"
    printf "    dir1        First directory with Cable-POP output, containing subdirectory restart/.\n"
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
if [[ ! -d ${dir1} ]] ; then
    printf "Error ${pprog}: First input is not a directory: ${dir1}.\n\n" 1>&2
    usage 1>&2
    exit 1
fi    
if [[ ! -d ${dir2} ]] ; then
    printf "Error ${pprog}: Second input is not a directory: ${dir2}.\n\n" 1>&2
    usage 1>&2
    exit 1
fi
if [[ ! -d ${dir1}/restart ]] ; then
    printf "Error ${pprog}: First input directory has no subdirectory restart/: ${dir1}.\n\n" 1>&2
    usage 1>&2
    exit 1
fi    
if [[ ! -d ${dir2}/restart ]] ; then
    printf "Error ${pprog}: Second input directory has no subdirectory restart/: ${dir1}.\n\n" 1>&2
    usage 1>&2
    exit 1
fi    
adir1=$(abspath ${dir1})
adir2=$(abspath ${dir2})

# -------------------------------------------------------------------------------------------------
# Compare restart files
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

echo "${adir1} vs. ${adir2}"

cd ${adir1}/restart/

if [[ true ]] ; then
    for suff in ${suff1} ${suff2} ; do
	echo 'Step ' ${suff}
	for i in *${suff}.nc ; do
	    echo ' ' ${i}
	    set +e
	    cdo -s diffv ${i} ${adir2}/restart/${i} 2> /dev/null
	    set -e
	done
    done
fi

if [[ true ]] ; then
    for suff in ${suff2} ; do
	i=pop_cru_ini${suff}.nc
	ncdiff ${i} ${adir2}/restart/${i} tmp.${pid}.nc
	set +e
	# iout=$(ncdump tmp.${pid}.nc | sed '/ latitude/,/;$/d' | sed '/ longitude/,/;$/d' | grep -v '0, 0' | sed -e '/[:;{}]/d' -e '/^$/d' | sed -e '/=[[:blank:]]$/d')
	iout=$(ncdump tmp.${pid}.nc | sed -e '/^netcdf/,/^variables:/d' -e '/ latitude/,/;$/d' -e '/ longitude/,/;$/d' -e '/[:{}]/d' -e '/^$/d' -e '/) ;$/d' -e 's/ ;/,/' -e 's/ 0,//g' | sed -e '/^[[:blank:]]*$/d' -e '/=[[:blank:]]*$/d' -e '/[tT]ime =/d')
	set -e
	if [[ -n ${iout} ]] ; then
	    echo " Check ${i} in ${PWD}"
	    # echo "     ncdiff -O ${i} ${adir2}/restart/${i} tmp.nc ; ncdump tmp.nc | sed '/ latitude/,/;$/d' | sed '/ longitude/,/;$/d' | grep -v '0, 0' | sed -e '/[:;{}]/d' -e '/^$/d'"
	    echo "     ncdiff -O ${i} ${adir2}/restart/${i} tmp.nc ; ncdump tmp.nc | sed -e '/^netcdf/,/^variables:/d' -e '/ latitude/,/;$/d' -e '/ longitude/,/;$/d' -e '/[:{}]/d' -e '/^$/d' -e '/) ;$/d' -e 's/ ;/,/' -e 's/ 0,//g' | sed -e '/^[[:blank:]]*$/d' -e '/=[[:blank:]]*$/d' -e '/[tT]ime =/d'"
	fi
	rm tmp.${pid}.nc
    done
fi

cd ${isdir}

exit
