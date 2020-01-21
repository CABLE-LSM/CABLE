#!/bin/bash
#
# Explor / Pearcey
# https://slurm.schedmd.com/sbatch.html
# Name
#SBATCH -J harvard10
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.out
# Explor partitions (sinfo): std (2x16, parallel), sky (2x16, parallel, AVX512), hf (2x4, serial),
#                            P100 (2x16, GPU), GTX (2x16, GPU), ivy (2x8, parallel), k20 (2x8, GPU)
#SBATCH -p hf
# Nodes / tasks
#SBATCH -N 1
#SBATCH -n 5
#SBATCH --ntasks-per-node=5
# Check memory on *nix with /usr/bin/time -v ./prog
# time (day-hh:mm:ss) / memory (optional, units K,M,G,T)
#SBATCH -t 00:59:59
#SBATCH --mem=4G
# notify: Valid type values are NONE,BEGIN,END,FAIL,REQUEUE,ALL,STAGE_OUT,TIME_LIMIT,TIME_LIMIT_90/80/50,ARRAY_TASKS
#SBATCH --mail-type=FAIL,STAGE_OUT,TIME_LIMIT
#SBATCH --mail-user=matthias.cuntz@inrae.fr
#
# # Raijin
# # https://opus.nci.org.au/display/Help/How+to+submit+a+job
# #PBS -P Harvard
# #PBS -q normal
# #PBS -l walltime=00:19:59
# #PBS -l mem=4GB
# #PBS -l ncpus=1
# #PBS -l jobfs=1GB
# #PBS -l software=netCDF:MPI:Intel:GNU
# #PBS -r y
# #PBS -l wd

system=cuntz@explor # cuntz@explor, cuntz@mcinra, knauer@pearcey, jk8585 or vxh599@raijin

# MPI run or single processor run
# nproc should fit with job tasks 
dompi=1   # 0: normal run: ./cable
          # 1: MPI run: mpiexec -n ${nproc} ./cable_mpi
nproc=3   # Number of cores for MPI runs
          # must be same as above: SBATCH -n nproc or PBS -l ncpus=nproc
ignore_mpi_err=1 # 0/1: continue even if mpi run failed

# --------------------------------------------------------------------
#
# Full Cable run on a single grid cell with biomass spinup, POP, land-use change, etc.
#
# This script uses CRU-JRA forcing.
# The single site met, LUH2, forcing and mask can be extracted from the global data sets
# in step 0.
#
# The run sequence is as follows:
#   1. Create a climate restart file using Cable's default vegetation distribution.
#   2. First phase of spinup with static land use, fixed atmospheric CO2 and
#      N deposition from 1700, and 30 years of repeated meteorology.
#      Zero initial biomass stocks, POP and climate.
#   3. Bring biomass stocks into equilibrium, restricting labile P and mineral N pools.
#      Repeat a and b several times.
#      a) Start from restart files.
#      b) Use dump files for the biophysics. Quasi-equilibrium of soil and litter pools
#         using an analytic steady-state solution.
#   4. Same as 3 but without any restriction on labile P and mineral N pools.
#      Repeat a and b several times.
#   5. Second phase of spinup with dynamic land use, atmospheric CO2 and
#      N deposition.
#      a) Dynamic land use from 1580 to 1699, using still fixed atmospheric CO2 and
#         N deposition from 1700, and 30 years of repeated meteorology.
#      b) Run from 1700 to 1899 with dynmic land use, varying atmospheric CO2 and N deposition,
#         but still with 30 years of repeated meteorology.
#   6. Final run, everything dynamic from 1900 to 2017.
#
# Written, Matthias Cuntz, August-December 2019, following the run scripts and namelists provided by Vanessa Haverd
#
# --------------------------------------------------------------------

set -e

trap cleanup 1 2 3 6

pid=$$
isdir="${PWD}"
prog=$0
pprog=$(basename ${prog})
pdir=$(dirname ${prog})
tmp=${TMPDIR:-"/tmp"}
#
system=$(echo ${system} | tr A-Z a-z)
sys=${system#*@}
user=${system%@*}
# Special things on specific computer system such as loading modules
export mpiexecdir=
if [[ "${sys}" == "explor" ]] ; then
    # prog is slurm_script
    pdir=${isdir}
    # # INTELMPI - load mpi module first, otherwise intel module will not pre-pend LD_LIBRARY_PATH
    # module load intelmpi/2018.5.274
    # module load intel/2018.5
    # export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HOME}/local/lib:${HOME}/local/netcdf-fortran-4.4.4-ifort2018.0/lib
    # export mpiexecdir=/soft/env/soft/all/intel/2018.3/compilers_and_libraries_2018.5.274/linux/mpi/intel64/bin
    # INTEL / OpenMPI - load mpi module first, otherwise intel module will not pre-pend LD_LIBRARY_PATH
    module load openmpi/3.0.0/intel18
    module load intel/2018.5
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HOME}/local/lib:${HOME}/local/netcdf-fortran-4.4.4-ifort2018.0/lib
    export mpiexecdir=/opt/soft/hf/openmpi-3.0.0-intel18/bin
    # # GNU / OpenMPI
    # module load gcc/6.3.0
    # module load openmpi/3.0.1/gcc/6.3.0
    # export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HOME}/local/lib:${HOME}/local/netcdf-fortran-4.4.4-gfortran63/lib
    # export mpiexecdir=/opt/soft/hf/openmpi/3.0.1/gcc/6.3.0/bin
elif [[ "${sys}" == "mcinra" ]] ; then
    # exe="${cablehome}/branches/NESP2pt9_BLAZE/offline/cable-mpi-gfortran"
    export mpiexecdir=/usr/local/openmpi-3.1.4-gfortran/bin
    # # exe="${cablehome}/branches/NESP2pt9_BLAZE/offline/cable-mpi-ifort"
    # export mpiexecdir=/usr/local/openmpi-3.1.5-ifort/bin
elif [[ "${sys}" == "pearcey" ]] ; then
    # prog is slurm_script
    pdir=${isdir}
    module del intel-cc intel-fc
    module add intel-cc/16.0.1.150 intel-fc/16.0.1.150
    module unload intel-mpi/5.0.1.035
    module add netcdf/4.3.3.1 openmpi/1.10.2
elif [[ "${sys}" == "raijin" ]] ; then
    module del intel-cc intel-fc
    module add intel-cc/16.0.1.150 intel-fc/16.0.1.150
    module add netcdf/4.3.3.1
fi
if [[ ! -z ${mpiexecdir} ]] ; then export mpiexecdir="${mpiexecdir}/" ; fi

# --------------------------------------------------------------------
# Sequence switches
#
imeteo=1       # 0: Use global meteo, land use and mask
                # 1: Use local meteo, land use and mask (doextractsite=1)
                # 2: Use global meteo and land use, and local mask (doextractsite=2)
# Step 0
doextractsite=0 # 0: Do not extract meteo, land use and mask at specific site
                # 1: Do extract meteo, land use and mask at specific site (imeteo=1)
                # 2: Do extract mask at specific site, using then global meteo and land use (imeteo=2)
  sitename=HarvardForest10
  # latlon=42.536875,-72.172602 # lat,lon  or  latmin,latmax,lonmin,lonmax  # must have . in numbers otherwise indexes taken
  latlon=42.536875,42.536875,-74.,-72.
# Step 1
doclimate=0     # 1/0: Do/Do not create climate restart file
# Step 2
dofromzero=0    # 1/0: Do/Do not first spinup phase from zero biomass stocks
# Step 3
doequi1=0       # 1/0: Do/Do not bring biomass stocks into quasi-equilibrium with restricted P and N pools
 nequi1=3       #      number of times to repeat steps in doequi1
# Step 4
doequi2=0       # 1/0: Do/Do not bring biomass stocks into quasi-equilibrium with unrestricted P and N pools
 nequi2=3       #      number of times to repeat steps in doequi2
# Step 5
doiniluc=1      # 1/0: Do/Do not spinup with dynamic land use (5a)
doinidyn=1      # 1/0: Do/Do not full dynamic spinup from 1700 to 1899 (5b)
# Step 6
dofinal=1       # 1/0: Do/Do not final run from 1900 to 2017

# --------------------------------------------------------------------
# Other switches
#
# Cable
doc13o2=1           # 1/0: Do/Do not calculate 13C
c13o2_simple_disc=0 # 1/0: simple or full 13C leaf discrimination
explicit_gm=0       # 1/0: explicit (finite) or implicit mesophyll conductance

# --------------------------------------------------------------------
# Setup
#
# Relative directories must be relative to the directory of this script,
#   not relative to the directory from which this script is launched (if different)
#   nor relative to the run path.
#
if [[ "${system}" == "cuntz@explor" ]] ; then
    # Run directory: runpath="${sitepath}/run_xxx"
    sitepath="/home/oqx29/zzy20/prog/cable/single_sites/${sitename}"
    cablehome="/home/oqx29/zzy20/prog/cable"
    # Cable executable
    if [[ ${dompi} -eq 1 ]] ; then
	exe="${cablehome}/branches/NESP2pt9_BLAZE/offline/cable-mpi"
    else
	exe="${cablehome}/branches/NESP2pt9_BLAZE/offline/cable"
    fi
    # CABLE-AUX directory (uses offline/gridinfo_CSIRO_1x1.nc and offline/modis_phenology_csiro.txt)
    aux="${cablehome}/CABLE-AUX"
    # Global Mask
    GlobalLandMaskFile="/home/oqx29/zzy20/data/crujra/daily_1deg/glob_ipsl_1x1.nc"
    # Global CRU
    GlobalMetPath="/home/oqx29/zzy20/data/crujra/daily_1deg"
    # Global LUC
    GlobalTransitionFilePath="/OSM/CBR/OA_GLOBALCABLE/work/LUH2/v3/1deg"
elif [[ "${system}" == "cuntz@mcinra" ]] ; then
    # Run directory: runpath="${sitepath}/run_xxx"
    sitepath="/Users/cuntz/prog/vanessa/cable/single_sites/${sitename}"
    cablehome="/Users/cuntz/prog/vanessa/cable"
    # Cable executable
    if [[ ${dompi} -eq 1 ]] ; then
	# exe="${cablehome}/branches/NESP2pt9_BLAZE/offline/cable-mpi-ifort"
	exe="${cablehome}/branches/NESP2pt9_BLAZE/offline/cable-mpi-gfortran"
    else
	exe="${cablehome}/branches/NESP2pt9_BLAZE/offline/cable"
    fi
    # CABLE-AUX directory (uses offline/gridinfo_CSIRO_1x1.nc and offline/modis_phenology_csiro.txt)
    aux="${cablehome}/CABLE-AUX"
    # Global Mask, CRU, LUC
    GlobalLandMaskFile=
    GlobalMetPath=
    GlobalTransitionFilePath=
elif [[ "${system}" == "knauer@pearcey" ]] ; then
    # Run directory: runpath="${sitepath}/run_xxx"
    sitepath="/OSM/CBR/OA_GLOBALCABLE/work/Juergen/CABLE_run/${sitename}"
    cablehome="/OSM/CBR/OA_GLOBALCABLE/work/Juergen/CABLE_code"
    # Cable executable
    if [[ ${dompi} -eq 1 ]] ; then
	exe="${cablehome}/NESP2pt9_BLAZE/offline/cable-mpi"
    else
	exe="${cablehome}/NESP2pt9_BLAZE/offline/cable"
    fi
    # CABLE-AUX directory (uses offline/gridinfo_CSIRO_1x1.nc and offline/modis_phenology_csiro.txt)
    aux="/OSM/CBR/OA_GLOBALCABLE/work/Vanessa/CABLE-AUX"
    # Global Mask
    GlobalLandMaskFile="/OSM/CBR/OA_GLOBALCABLE/work/Vanessa/MASKS/glob_ipsl_1x1.nc"
    # Global CRU
    GlobalMetPath="/OSM/CBR/OA_GLOBALCABLE/work/CRU-JRA55/crujra/daily_1deg"
    # Global LUC
    GlobalTransitionFilePath="/OSM/CBR/OA_GLOBALCABLE/work/LUH2/v3/1deg"
else
    echo "System not known."
    exit 1
fi
# Run directory
runpath="${sitepath}/run_20200117"

# Cable parameters
namelistpath="../namelists"
filename_veg="../params/def_veg_params.txt"
filename_soil="../params/def_soil_params.txt"
casafile_cnpbiome="../params/pftlookup.csv"
# Mask
LandMaskFile="${sitepath}/mask/glob_ipsl_1x1_${sitename}.nc"
# CRU
MetPath="${sitepath}/met/cru_jra_1deg"
ClimateFile="${sitepath}/met/cru_climate_rst.nc"
# LUC
TransitionFilePath="${sitepath}/LUH2/v3/1deg"
# 13C
filename_d13c_atm="../params/graven_et_al_gmd_2017-table_s1-delta_13c-1700-2025.txt"

# --------------------------------------------------------------------
# Start Script
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Helper functions
#
# usage of script
function usage()
{
    printf "${pprog} [-h]\n"
    printf "Runs Cable on a single grid cell with spinup, POP, land-use change, etc.\n"
    printf "Behaviour of the script is controlled by switches at the top of the script (ca. line 101ff).\n"
    printf "\n"
    printf "Options\n"
    printf "    -h    Prints this help screen.\n"
}
#
# cleanup at end or at trap
function cleanup()
{
    \rm -f ${tmp}/*.${pid}*
    exit 1
}
#
# absolute path
function abspath()
{
    idir=${PWD}
    cd ${1}
    odir=${PWD}
    cd ${idir}
    echo "${odir}"
}
#
# filename with absolute path
function absfile()
{
    f=$(basename ${1})
    d=$(dirname ${1})
    d=$(abspath ${d})
    echo "${d}/${f}"
}
#
# rm file if present
function irm()
{
    for i in "$@" ; do
	if [[ -f ${i} ]] ; then rm ${i} ; fi
    done
}
#
# sed command from comma-separated list: var1=str1,var2=str2
#   no , but = allowed in str
function csed()
{
    com=""
    for i in $(echo ${1} | tr ',' ' ') ; do
        v=$(echo ${i} | cut -d '=' -f 1)
        s=$(echo ${i} | cut -d '=' -f 2-)
        # com="${com} -e '/^[[:blank:]]*${v}[[:blank:]]*=/s/${v}[[:blank:]]*=.*/${v} = ${s}/'"
        # com="${com} -e 's/${v}[[:blank:]]*=.*/${v} = ${s}/'"
        com="${com} -e s|${v}[[:blank:]]*=.*|${v}=${s}|"
    done
    printf "%s" "${com}"
}
#
# get list of variables in file with nco (cdo omits the dimension variables)
function ncvarlist()
{
    out=$(ncks --trd -m ${1} | grep -E ': type' | cut -f 1 -d ' ' | sed 's/://' | sort)
    echo ${out}
}
#
# check if first argument string is in rest of arguments
# returns first argument if present otherwise empty string (check with -z)
function isin()
{
    tofind=${1}
    shift
    out=""
    for i in $@ ; do
	if [[ ${i} == ${tofind} ]] ; then out=${i} ; fi
    done
    echo ${out}
}
#
# returns argument to extract lat and lon with ncks
function nckslatlon()
{
    vars=$(ncvarlist ${1})
    if [[ -z $(isin latitude ${vars}) ]] ; then ilat='lat' ; else ilat='latitude' ; fi
    if [[ -z $(isin longitude ${vars}) ]] ; then ilon='lon' ; else ilon='longitude' ; fi
    if [[ -z $(echo ${2} | cut -f 3 -d ',') || -z $(echo ${2} | cut -f 4 -d ',') ]] ; then
	iilat=$(echo ${2} | cut -f 1 -d ',')
	iilon=$(echo ${2} | cut -f 2 -d ',')
	echo "-d ${ilat},${iilat} -d ${ilon},${iilon}"
    else
	iilat1=$(echo ${2} | cut -f 1 -d ',')
	iilat2=$(echo ${2} | cut -f 2 -d ',')
	iilon1=$(echo ${2} | cut -f 3 -d ',')
	iilon2=$(echo ${2} | cut -f 4 -d ',')
	echo "-d ${ilat},${iilat1},${iilat2} -d ${ilon},${iilon1},${iilon2}"
    fi
}
#
# copy files adding first argument to filenames
function copyid()
{
    rid=${1}
    shift 1
    for i in $@ ; do
	cp ${i} ${i%.*}_${rid}.${i##*.}
    done
}
#
# rename files adding first argument to filenames
function renameid()
{
    rid=${1}
    shift 1
    for i in $@ ; do
	mv ${i} ${i%.*}_${rid}.${i##*.}
    done
}

# --------------------------------------------------------------------------------------------------
# Preparation
#
# Get options
while getopts "h" option ; do
    case ${option} in
        h) usage; exit;;
        *) printf "Error ${pprog}: unimplemented option.\n\n" 1>&2;  usage 1>&2; exit 1;;
    esac
done
shift $((${OPTIND} - 1))

#
# get directories
pdir=$(abspath ${pdir})
cd ${pdir}
adir=$(abspath ${aux})
exe=$(absfile ${exe})
mkdir -p ${runpath}
rdir=$(abspath ${runpath})
ndir=$(abspath ${namelistpath})

#
# prepare run directory
cd ${rdir}
mkdir -p logs
mkdir -p outputs
mkdir -p restart
ln -sf ${adir}
ln -sf ${exe}
iexe=$(basename ${exe})
cd ${pdir}

#
# set stacksize to unlimited if permitted, otherwise to 15 bit if possible
set +e
ulimit -s unlimited 2> /dev/null || ulimit -s 32768
set -e

# --------------------------------------------------------------------
# Sequence
#
t1=$(date +%s)
printf "Started at %s\n" "$(date)"

# 0. Extract meteo, land use and mask for one specific site from global files
if [[ ${doextractsite} -ge 1 ]] ; then
    # xcdo=$(which cdo)
    # if [[ -z ${xcdo} ]] ; then module load cdo ; fi
    xnco=$(which ncks)
    if [[ -z ${xnco} ]] ; then module load nco ; fi
fi
if [[ ${doextractsite} -eq 1 ]] ; then
    cd ${pdir}
    # meteorology
    met_list="pre pres dlwrf dswrf spfh tmax tmin ugrd vgrd ndep"
    mkdir -p ${MetPath}
    MetPath=$(abspath ${MetPath})
    for mm in ${met_list} ; do
        echo ${GlobalMetPath}/${mm}
        mkdir -p ${MetPath}/${mm}
        for nc in ${GlobalMetPath}/${mm}/*.nc ; do
	    ff=$(basename ${nc})
            echo "    ${ff}"
            # cdo -s -f nc4 -z zip sellonlatbox,-72.5,-72.0,42.5,43.0 ${nc} ${MetPath}/${mm}/${ff}
	    ncks -O $(nckslatlon ${nc} ${latlon}) ${nc} ${MetPath}/${mm}/${ff}
        done
    done
    echo ${GlobalMetPath}/co2/*.csv
    mkdir -p ${MetPath}/co2
    cp ${GlobalMetPath}/co2/*.csv ${MetPath}/co2/

    # land use
    mkdir -p ${TransitionFilePath}
    TransitionFilePath=$(abspath ${TransitionFilePath})
    echo ${GlobalTransitionFilePath}
    for nc in ${GlobalTransitionFilePath}/*.nc ; do
	ff=$(basename ${nc})
        echo "    ${ff}"
        # cdo -s -f nc4 -z zip sellonlatbox,-72.5,-72.0,42.5,43.0 ${nc} ${TransitionFilePath}/${ff}
	ncks -O $(nckslatlon ${nc} ${latlon}) ${nc} ${TransitionFilePath}/${ff}
    done

    # mask
    LandMaskFilePath=$(dirname ${LandMaskFile})
    mkdir -p ${LandMaskFilePath}
    LandMaskFile=$(absfile ${LandMaskFile})
    echo $(basename ${LandMaskFile})
    # cdo -s -f nc4 -z zip sellonlatbox,-72.5,-72.0,42.5,43.0 ${GlobalLandMaskFile} ${LandMaskFile}
    ncks -O $(nckslatlon ${GlobalLandMaskFile} ${latlon}) ${GlobalLandMaskFile} ${LandMaskFile}
fi
# ToDo: set land points in mask using Python script
if [[ ${doextractsite} -eq 2 ]] ; then
    cd ${pdir}
    # mask
    LandMaskFilePath=$(dirname ${LandMaskFile})
    mkdir -p ${LandMaskFilePath}
    LandMaskFile=$(absfile ${LandMaskFile})
    echo $(basename ${LandMaskFile})
    # cdo -s -f nc4 -z zip sellonlatbox,-72.5,-72.0,42.5,43.0 ${GlobalLandMaskFile} ${LandMaskFile}
    ncks -O $(nckslatlon ${GlobalLandMaskFile} ${latlon}) ${GlobalLandMaskFile} ${LandMaskFile}
fi

# Choose meteo, land use and mask directories and files
if [[ ${imeteo} -eq 0 ]] ; then
    MetPath=$(abspath ${GlobalMetPath})
    TransitionFilePath=$(abspath ${GlobalTransitionFilePath})
    LandMaskFile=$(absfile ${LandMaskFile})
elif [[ ${imeteo} -eq 1 ]] ; then
    MetPath=$(abspath ${MetPath})
    TransitionFilePath=$(abspath ${TransitionFilePath})
    LandMaskFile=$(absfile ${LandMaskFile})
elif [[ ${imeteo} -eq 2 ]] ; then
    MetPath=$(abspath ${GlobalMetPath})
    TransitionFilePath=$(abspath ${GlobalTransitionFilePath})
    LandMaskFile=$(absfile ${LandMaskFile})
else
    printf "Error ${pprog}: imeteo option unknown: ${imeteo}.\n\n"
    exit 1
fi
# absolute pathes of other parameter files
ClimateFile=$(absfile ${ClimateFile})
filename_veg=$(absfile ${filename_veg})
filename_soil=$(absfile ${filename_soil})
casafile_cnpbiome=$(absfile ${casafile_cnpbiome})
filename_d13c_atm=$(absfile ${filename_d13c_atm})


# 1. Create climate restart file
if [[ ${doclimate} -eq 1 ]] ; then
    echo "1. Create climate restart file"
    rid="climate_restart"
    # CRU
    irm ${rdir}/cru.nml
    com=$(csed "BasePath=\"${MetPath}\",MetPath=\"${MetPath}\",LandMaskFile=\"${LandMaskFile}\"")
    com=${com}$(csed "Run=\"S0_TRENDY\"")
    sed ${com} ${ndir}/cru.nml > ${rdir}/cru.nml
    # LUC
    irm ${rdir}/LUC.nml
    com=$(csed "TransitionFilePath=\"${TransitionFilePath}\",ClimateFile=\"${ClimateFile}\"")
    com=${com}$(csed "YearStart=1700")
    com=${com}$(csed "YearEnd=2017")
    sed ${com} ${ndir}/LUC.nml > ${rdir}/LUC.nml
    # CABLE
    irm ${rdir}/cable.nml
    com=$(csed "filename%veg=\"${filename_veg}\",filename%soil=\"${filename_soil}\"")
    com=${com}$(csed "casafile%cnpbiome=\"${casafile_cnpbiome}\"")
    com=${com}$(csed "filename%restart_in=\"\"")
    com=${com}$(csed "cable_user%CLIMATE_fromZero=.true.")
    com=${com}$(csed "cable_user%YearStart=1860")
    com=${com}$(csed "cable_user%YearEnd=1889")
    com=${com}$(csed "icycle=2")
    com=${com}$(csed "spincasa=.false.")
    com=${com}$(csed "cable_user%CASA_fromZero=.true.")
    com=${com}$(csed "cable_user%CASA_DUMP_READ=.false.")
    com=${com}$(csed "cable_user%CASA_DUMP_WRITE=.true.")
    com=${com}$(csed "cable_user%CASA_OUT_FREQ=\"monthly\"")
    com=${com}$(csed "cable_user%CASA_SPIN_STARTYEAR=1860")
    com=${com}$(csed "cable_user%CASA_SPIN_ENDYEAR=1869")
    com=${com}$(csed "cable_user%limit_labile=.true.")
    com=${com}$(csed "casafile%out=\"outputs/cru_out_casa.nc\"")
    com=${com}$(csed "casafile%cnpipool=\"\"")
    com=${com}$(csed "cable_user%POP_fromZero=.true.")
    com=${com}$(csed "cable_user%POP_out=\"ini\"")
    com=${com}$(csed "cable_user%POP_restart_in=\"\"")
    com=${com}$(csed "cable_user%POPLUC=.false.")
    com=${com}$(csed "cable_user%POPLUC_RunType=\"static\"")
    com=${com}$(csed "cable_user%LUC_outfile=\"outputs/cru_out_LUC.nc\"")
    com=${com}$(csed "cable_user%LUC_restart_in=\"restart/cru_LUC_rst.nc\"")
    com=${com}$(csed "cable_user%LUC_restart_out=\"restart/cru_LUC_rst.nc\"")
    if [[ ${explicit_gm} -eq 1 ]] ; then
        com=${com}$(csed "cable_user%explicit_gm=.true.")
    else
        com=${com}$(csed "cable_user%explicit_gm=.false.")
    fi
    if [[ ${doc13o2} -eq 1 ]] ; then
        com=${com}$(csed "cable_user%c13o2=.true.")
        com=${com}$(csed "cable_user%c13o2_delta_atm_file=\"${filename_d13c_atm}\"")
        com=${com}$(csed "cable_user%c13o2_outfile=\"outputs/cru_out_casa_c13o2.nc\"")
        com=${com}$(csed "cable_user%c13o2_restart_in_flux=\"\"")
        com=${com}$(csed "cable_user%c13o2_restart_out_flux=\"restart/cru_c13o2_flux_rst.nc\"")
        com=${com}$(csed "cable_user%c13o2_restart_in_pools=\"\"")
        com=${com}$(csed "cable_user%c13o2_restart_out_pools=\"restart/cru_c13o2_pools_rst.nc\"")
        if [[ ${c13o2_simple_disc} -eq 1 ]] ; then
            com=${com}$(csed "cable_user%c13o2_simple_disc=.true.")
        else
            com=${com}$(csed "cable_user%c13o2_simple_disc=.false.")
        fi
    fi
    sed ${com} ${ndir}/cable.nml > ${rdir}/cable.nml
    # run model
    cd ${rdir}
    irm logs/log_cable.txt logs/log_out_cable.txt
    if [[ ${dompi} -eq 1 ]] ; then
	if [[ ${ignore_mpi_err} -eq 1 ]] ; then set +e ; fi
	${mpiexecdir}mpiexec -n ${nproc} ./${iexe} > logs/log_out_cable.txt
	if [[ ${ignore_mpi_err} -eq 1 ]] ; then set -e ; fi
    else
	./${iexe} > logs/log_out_cable.txt
    fi
    # save output
    copyid ${rid} cable.nml cru.nml LUC.nml
    mv *_${rid}.nml restart/
    cd logs
    renameid ${rid} log_cable.txt log_out_cable.txt
    cd ../restart
    copyid ${rid} pop_cru_ini.nc cru_climate_rst.nc cru_casa_rst.nc cru_cable_rst.nc
    cp cru_climate_rst.nc ${ClimateFile}
    if [[ ${doc13o2} -eq 1 ]] ; then copyid ${rid} cru_c13o2_flux_rst.nc cru_c13o2_pools_rst.nc ; fi
    cd ../outputs
    renameid ${rid} cru_out_cable.nc cru_out_casa.nc
    if [[ ${doc13o2} -eq 1 ]] ; then copyid ${rid} cru_out_casa_c13o2.nc ; fi
    cd ..
    cd ${pdir}
fi


# 2. First spinup phase from zero biomass
if [[ ${dofromzero} -eq 1 ]] ; then
    echo "2. First spinup from zero biomass"
    rid="zero_biomass"
    # CRU
    irm ${rdir}/cru.nml
    com=$(csed "BasePath=\"${MetPath}\",MetPath=\"${MetPath}\",LandMaskFile=\"${LandMaskFile}\"")
    com=${com}$(csed "Run=\"S0_TRENDY\"")
    sed ${com} ${ndir}/cru.nml > ${rdir}/cru.nml
    # LUC
    irm ${rdir}/LUC.nml
    com=$(csed "TransitionFilePath=\"${TransitionFilePath}\",ClimateFile=\"${ClimateFile}\"")
    com=${com}$(csed "YearStart=1700")
    com=${com}$(csed "YearEnd=2017")
    sed ${com} ${ndir}/LUC.nml > ${rdir}/LUC.nml
    # CABLE
    irm ${rdir}/cable.nml
    com=$(csed "filename%veg=\"${filename_veg}\",filename%soil=\"${filename_soil}\"")
    com=${com}$(csed "casafile%cnpbiome=\"${casafile_cnpbiome}\"")
    com=${com}$(csed "filename%restart_in=\"\"")
    com=${com}$(csed "cable_user%CLIMATE_fromZero=.true.")
    com=${com}$(csed "cable_user%YearStart=1860")
    com=${com}$(csed "cable_user%YearEnd=1889")
    com=${com}$(csed "icycle=2")
    com=${com}$(csed "spincasa=.false.")
    com=${com}$(csed "cable_user%CASA_fromZero=.true.")
    com=${com}$(csed "cable_user%CASA_DUMP_READ=.false.")
    com=${com}$(csed "cable_user%CASA_DUMP_WRITE=.true.")
    com=${com}$(csed "cable_user%CASA_OUT_FREQ=\"monthly\"")
    com=${com}$(csed "cable_user%CASA_SPIN_STARTYEAR=1860")
    com=${com}$(csed "cable_user%CASA_SPIN_ENDYEAR=1869")
    com=${com}$(csed "cable_user%limit_labile=.true.")
    com=${com}$(csed "casafile%out=\"outputs/cru_out_casa.nc\"")
    com=${com}$(csed "casafile%cnpipool=\"\"")
    com=${com}$(csed "cable_user%POP_fromZero=.true.")
    com=${com}$(csed "cable_user%POP_out=\"ini\"")
    com=${com}$(csed "cable_user%POP_restart_in=\"\"")
    com=${com}$(csed "cable_user%POPLUC=.true.")
    com=${com}$(csed "cable_user%POPLUC_RunType=\"static\"")
    com=${com}$(csed "cable_user%LUC_outfile=\"outputs/cru_out_LUC.nc\"")
    com=${com}$(csed "cable_user%LUC_restart_in=\"restart/cru_LUC_rst.nc\"")
    com=${com}$(csed "cable_user%LUC_restart_out=\"restart/cru_LUC_rst.nc\"")
    if [[ ${explicit_gm} -eq 1 ]] ; then
        com=${com}$(csed "cable_user%explicit_gm=.true.")
    else
        com=${com}$(csed "cable_user%explicit_gm=.false.")
    fi
    if [[ ${doc13o2} -eq 1 ]] ; then
        com=${com}$(csed "cable_user%c13o2=.true.")
        com=${com}$(csed "cable_user%c13o2_delta_atm_file=\"${filename_d13c_atm}\"")
        com=${com}$(csed "cable_user%c13o2_outfile=\"outputs/cru_out_casa_c13o2.nc\"")
        com=${com}$(csed "cable_user%c13o2_restart_in_flux=\"\"")
        com=${com}$(csed "cable_user%c13o2_restart_out_flux=\"restart/cru_c13o2_flux_rst.nc\"")
        com=${com}$(csed "cable_user%c13o2_restart_in_pools=\"\"")
        com=${com}$(csed "cable_user%c13o2_restart_out_pools=\"restart/cru_c13o2_pools_rst.nc\"")
        if [[ ${c13o2_simple_disc} -eq 1 ]] ; then
            com=${com}$(csed "cable_user%c13o2_simple_disc=.true.")
        else
            com=${com}$(csed "cable_user%c13o2_simple_disc=.false.")
        fi
    fi
    sed ${com} ${ndir}/cable.nml > ${rdir}/cable.nml
    # run model
    cd ${rdir}
    irm logs/log_cable.txt logs/log_out_cable.txt
    if [[ ${dompi} -eq 1 ]] ; then
	if [[ ${ignore_mpi_err} -eq 1 ]] ; then set +e ; fi
	${mpiexecdir}mpiexec -n ${nproc} ./${iexe} > logs/log_out_cable.txt
	if [[ ${ignore_mpi_err} -eq 1 ]] ; then set -e ; fi
    else
	./${iexe} > logs/log_out_cable.txt
    fi
    # save output
    copyid ${rid} cable.nml cru.nml LUC.nml
    mv *_${rid}.nml restart/
    cd logs
    renameid ${rid} log_cable.txt log_out_cable.txt
    cd ../restart
    copyid ${rid} pop_cru_ini.nc cru_climate_rst.nc cru_casa_rst.nc cru_cable_rst.nc
    if [[ ${doc13o2} -eq 1 ]] ; then copyid ${rid} cru_c13o2_flux_rst.nc cru_c13o2_pools_rst.nc ; fi
    cd ../outputs
    renameid ${rid} cru_out_cable.nc cru_out_casa.nc
    if [[ ${doc13o2} -eq 1 ]] ; then copyid ${rid} cru_out_casa_c13o2.nc ; fi
    cd ..
    cd ${pdir}
fi


# 3. Biomass into quasi-equilibrium with restricted N and P pools
if [[ ${doequi1} -eq 1 ]] ; then
    echo "3. Bring biomass into quasi-equilibrium with restricted N and P pools"
    for ((iequi1=1; iequi1<=${nequi1}; iequi1++)) ; do
	rid="spinup_limit_labile"
	# rid="spinup_limit_labile${iequi}"
        if [[ 1 -eq 1 ]] ; then
            # 3a. 30 year run starting from restart files
            echo "   3a. 30 year spinup from accumulated biomass; iequi1=${iequi1}/${nequi1}"
            # CRU
            irm ${rdir}/cru.nml
            com=$(csed "BasePath=\"${MetPath}\",MetPath=\"${MetPath}\",LandMaskFile=\"${LandMaskFile}\"")
            com=${com}$(csed "Run=\"S0_TRENDY\"")
            sed ${com} ${ndir}/cru.nml > ${rdir}/cru.nml
            # LUC
            irm ${rdir}/LUC.nml
            com=$(csed "TransitionFilePath=\"${TransitionFilePath}\",ClimateFile=\"${ClimateFile}\"")
            com=${com}$(csed "YearStart=1700")
            com=${com}$(csed "YearEnd=2017")
            sed ${com} ${ndir}/LUC.nml > ${rdir}/LUC.nml
            # CABLE
            irm ${rdir}/cable.nml
            com=$(csed "filename%veg=\"${filename_veg}\",filename%soil=\"${filename_soil}\"")
            com=${com}$(csed "casafile%cnpbiome=\"${casafile_cnpbiome}\"")
            com=${com}$(csed "filename%restart_in=\"restart/cru_cable_rst.nc\"")
            com=${com}$(csed "cable_user%CLIMATE_fromZero=.false.")
            com=${com}$(csed "cable_user%YearStart=1840")
            com=${com}$(csed "cable_user%YearEnd=1859")
            com=${com}$(csed "icycle=2")
            com=${com}$(csed "spincasa=.false.")
            com=${com}$(csed "cable_user%CASA_fromZero=.false.")
            com=${com}$(csed "cable_user%CASA_DUMP_READ=.false.")
            com=${com}$(csed "cable_user%CASA_DUMP_WRITE=.true.")
            com=${com}$(csed "cable_user%CASA_OUT_FREQ=\"monthly\"")
            com=${com}$(csed "cable_user%CASA_SPIN_STARTYEAR=1860")
            com=${com}$(csed "cable_user%CASA_SPIN_ENDYEAR=1869")
            com=${com}$(csed "cable_user%limit_labile=.true.")
            com=${com}$(csed "casafile%out=\"outputs/cru_out_casa.nc\"")
            com=${com}$(csed "casafile%cnpipool=\"restart/cru_casa_rst.nc\"")
            com=${com}$(csed "cable_user%POP_fromZero=.false.")
            com=${com}$(csed "cable_user%POP_out=\"ini\"")
            com=${com}$(csed "cable_user%POP_restart_in=\"restart/pop_cru_ini.nc\"")
            com=${com}$(csed "cable_user%POPLUC=.true.")
            com=${com}$(csed "cable_user%POPLUC_RunType=\"static\"")
            com=${com}$(csed "cable_user%LUC_outfile=\"outputs/cru_out_LUC.nc\"")
            com=${com}$(csed "cable_user%LUC_restart_in=\"restart/cru_LUC_rst.nc\"")
            com=${com}$(csed "cable_user%LUC_restart_out=\"restart/cru_LUC_rst.nc\"")
	    if [[ ${explicit_gm} -eq 1 ]] ; then
                com=${com}$(csed "cable_user%explicit_gm=.true.")
            else
                com=${com}$(csed "cable_user%explicit_gm=.false.")
            fi
            if [[ ${doc13o2} -eq 1 ]] ; then
                com=${com}$(csed "cable_user%c13o2=.true.")
		com=${com}$(csed "cable_user%c13o2_delta_atm_file=\"${filename_d13c_atm}\"")
                com=${com}$(csed "cable_user%c13o2_outfile=\"outputs/cru_out_casa_c13o2.nc\"")
                com=${com}$(csed "cable_user%c13o2_restart_in_flux=\"restart/cru_c13o2_flux_rst.nc\"")
                com=${com}$(csed "cable_user%c13o2_restart_out_flux=\"restart/cru_c13o2_flux_rst.nc\"")
                com=${com}$(csed "cable_user%c13o2_restart_in_pools=\"restart/cru_c13o2_pools_rst.nc\"")
                com=${com}$(csed "cable_user%c13o2_restart_out_pools=\"restart/cru_c13o2_pools_rst.nc\"")
                if [[ ${c13o2_simple_disc} -eq 1 ]] ; then
                    com=${com}$(csed "cable_user%c13o2_simple_disc=.true.")
                else
                    com=${com}$(csed "cable_user%c13o2_simple_disc=.false.")
                fi
            fi
            sed ${com} ${ndir}/cable.nml > ${rdir}/cable.nml
            # run model
            cd ${rdir}
	    irm logs/log_cable.txt logs/log_out_cable.txt
	    if [[ ${dompi} -eq 1 ]] ; then
		if [[ ${ignore_mpi_err} -eq 1 ]] ; then set +e ; fi
		${mpiexecdir}mpiexec -n ${nproc} ./${iexe} > logs/log_out_cable.txt
		if [[ ${ignore_mpi_err} -eq 1 ]] ; then set -e ; fi
	    else
		./${iexe} > logs/log_out_cable.txt
	    fi
            # save output
	    copyid ${rid} cable.nml cru.nml LUC.nml
	    mv *_${rid}.nml restart/
    	    cd logs
    	    renameid ${rid} log_cable.txt log_out_cable.txt
    	    cd ../restart
    	    copyid ${rid} pop_cru_ini.nc cru_climate_rst.nc cru_casa_rst.nc cru_cable_rst.nc
    	    if [[ ${doc13o2} -eq 1 ]] ; then copyid ${rid} cru_c13o2_flux_rst.nc cru_c13o2_pools_rst.nc ; fi
    	    cd ../outputs
    	    renameid ${rid} cru_out_cable.nc cru_out_casa.nc
    	    if [[ ${doc13o2} -eq 1 ]] ; then copyid ${rid} cru_out_casa_c13o2.nc ; fi
    	    cd ..
    	    cd ${pdir}
        fi
        #
        if [[ 1 -eq 1 ]] ; then
            # 3b. analytic quasi-equilibrium of biomass pools
            echo "   3b. Analytic solution of biomass pools"
	    rid="spinup_analytic_limit_labile"
	    # rid="spin_casa_limit_labile${iequi}"
            # CRU
            irm ${rdir}/cru.nml
            com=$(csed "BasePath=\"${MetPath}\",MetPath=\"${MetPath}\",LandMaskFile=\"${LandMaskFile}\"")
            com=${com}$(csed "Run=\"S0_TRENDY\"")
            sed ${com} ${ndir}/cru.nml > ${rdir}/cru.nml
            # LUC
            irm ${rdir}/LUC.nml
            com=$(csed "TransitionFilePath=\"${TransitionFilePath}\",ClimateFile=\"${ClimateFile}\"")
            com=${com}$(csed "YearStart=1700")
            com=${com}$(csed "YearEnd=2017")
            sed ${com} ${ndir}/LUC.nml > ${rdir}/LUC.nml
            # CABLE
            irm ${rdir}/cable.nml
            com=$(csed "filename%veg=\"${filename_veg}\",filename%soil=\"${filename_soil}\"")
            com=${com}$(csed "casafile%cnpbiome=\"${casafile_cnpbiome}\"")
            com=${com}$(csed "filename%restart_in=\"restart/cru_cable_rst.nc\"")
            com=${com}$(csed "cable_user%CLIMATE_fromZero=.false.")
            com=${com}$(csed "cable_user%YearStart=1840")
            com=${com}$(csed "cable_user%YearEnd=1859")
            com=${com}$(csed "icycle=12")
            com=${com}$(csed "spincasa=.true.")
            com=${com}$(csed "cable_user%CASA_fromZero=.false.")
            com=${com}$(csed "cable_user%CASA_DUMP_READ=.true.")
            com=${com}$(csed "cable_user%CASA_DUMP_WRITE=.false.")
            com=${com}$(csed "cable_user%CASA_OUT_FREQ=\"monthly\"")
            com=${com}$(csed "cable_user%CASA_SPIN_STARTYEAR=1840")
            com=${com}$(csed "cable_user%CASA_SPIN_ENDYEAR=1859")
            com=${com}$(csed "cable_user%limit_labile=.true.")
            com=${com}$(csed "casafile%out=\"outputs/cru_out_casa.nc\"")
            com=${com}$(csed "casafile%cnpipool=\"restart/cru_casa_rst.nc\"")
            com=${com}$(csed "cable_user%POP_fromZero=.false.")
            com=${com}$(csed "cable_user%POP_out=\"ini\"")
            com=${com}$(csed "cable_user%POP_restart_in=\"restart/pop_cru_ini.nc\"")
            com=${com}$(csed "cable_user%POPLUC=.true.")
            com=${com}$(csed "cable_user%POPLUC_RunType=\"static\"")
            com=${com}$(csed "cable_user%LUC_outfile=\"outputs/cru_out_LUC.nc\"")
            com=${com}$(csed "cable_user%LUC_restart_in=\"restart/cru_LUC_rst.nc\"")
            com=${com}$(csed "cable_user%LUC_restart_out=\"restart/cru_LUC_rst.nc\"")
	    if [[ ${explicit_gm} -eq 1 ]] ; then
                com=${com}$(csed "cable_user%explicit_gm=.true.")
            else
                com=${com}$(csed "cable_user%explicit_gm=.false.")
            fi
            if [[ ${doc13o2} -eq 1 ]] ; then
                com=${com}$(csed "cable_user%c13o2=.true.")
		com=${com}$(csed "cable_user%c13o2_delta_atm_file=\"${filename_d13c_atm}\"")
                com=${com}$(csed "cable_user%c13o2_outfile=\"outputs/cru_out_casa_c13o2.nc\"")
                com=${com}$(csed "cable_user%c13o2_restart_in_flux=\"restart/cru_c13o2_flux_rst.nc\"")
                com=${com}$(csed "cable_user%c13o2_restart_out_flux=\"restart/cru_c13o2_flux_rst.nc\"")
                com=${com}$(csed "cable_user%c13o2_restart_in_pools=\"restart/cru_c13o2_pools_rst.nc\"")
                com=${com}$(csed "cable_user%c13o2_restart_out_pools=\"restart/cru_c13o2_pools_rst.nc\"")
                if [[ ${c13o2_simple_disc} -eq 1 ]] ; then
                    com=${com}$(csed "cable_user%c13o2_simple_disc=.true.")
                else
                    com=${com}$(csed "cable_user%c13o2_simple_disc=.false.")
                fi
            fi
            sed ${com} ${ndir}/cable.nml > ${rdir}/cable.nml
            # run model
            cd ${rdir}
	    irm logs/log_cable.txt logs/log_out_cable.txt
	    if [[ ${dompi} -eq 1 ]] ; then
		if [[ ${ignore_mpi_err} -eq 1 ]] ; then set +e ; fi
		${mpiexecdir}mpiexec -n ${nproc} ./${iexe} > logs/log_out_cable.txt
		if [[ ${ignore_mpi_err} -eq 1 ]] ; then set -e ; fi
	    else
		./${iexe} > logs/log_out_cable.txt
	    fi
            # save output
	    copyid ${rid} cable.nml cru.nml LUC.nml
	    mv *_${rid}.nml restart/
    	    cd logs
    	    renameid ${rid} log_cable.txt log_out_cable.txt
    	    cd ../restart
    	    copyid ${rid} pop_cru_ini.nc cru_casa_rst.nc
    	    if [[ ${doc13o2} -eq 1 ]] ; then copyid ${rid} cru_c13o2_flux_rst.nc cru_c13o2_pools_rst.nc ; fi
	    if [[ ${dompi} -eq 0 ]] ; then # no output only restart if MPI
    		cd ../outputs
    		renameid ${rid} cru_out_casa.nc
    		if [[ ${doc13o2} -eq 1 ]] ; then copyid ${rid} cru_out_casa_c13o2.nc ; fi
    		cd ..
	    fi
    	    cd ${pdir}
        fi
    done
fi


# 4. Biomass into quasi-equilibrium without restricted N and P pools
if [[ ${doequi2} -eq 1 ]] ; then
    echo "4. Bring biomass into quasi-equilibrium without restricted N and P pools"
    for ((iequi2=1; iequi2<=${nequi2}; iequi2++)) ; do
	rid="spinup"
	# rid="spinup${iequi}"
        if [[ 1 -eq 1 ]] ; then
            # 4a. 30 year run starting from restart files
            echo "   4a. 30 year spinup from accumulated biomass; iequi2=${iequi2}/${nequi2}"
            # CRU
            irm ${rdir}/cru.nml
            com=$(csed "BasePath=\"${MetPath}\",MetPath=\"${MetPath}\",LandMaskFile=\"${LandMaskFile}\"")
            com=${com}$(csed "Run=\"S0_TRENDY\"")
            sed ${com} ${ndir}/cru.nml > ${rdir}/cru.nml
            # LUC
            irm ${rdir}/LUC.nml
            com=$(csed "TransitionFilePath=\"${TransitionFilePath}\",ClimateFile=\"${ClimateFile}\"")
            com=${com}$(csed "YearStart=1700")
            com=${com}$(csed "YearEnd=2017")
            sed ${com} ${ndir}/LUC.nml > ${rdir}/LUC.nml
            # CABLE
            irm ${rdir}/cable.nml
            com=$(csed "filename%veg=\"${filename_veg}\",filename%soil=\"${filename_soil}\"")
            com=${com}$(csed "casafile%cnpbiome=\"${casafile_cnpbiome}\"")
            com=${com}$(csed "filename%restart_in=\"restart/cru_cable_rst.nc\"")
            com=${com}$(csed "cable_user%CLIMATE_fromZero=.false.")
            com=${com}$(csed "cable_user%YearStart=1840")
            com=${com}$(csed "cable_user%YearEnd=1859")
            com=${com}$(csed "icycle=2")
            com=${com}$(csed "spincasa=.false.")
            com=${com}$(csed "cable_user%CASA_fromZero=.false.")
            com=${com}$(csed "cable_user%CASA_DUMP_READ=.false.")
            com=${com}$(csed "cable_user%CASA_DUMP_WRITE=.true.")
            com=${com}$(csed "cable_user%CASA_OUT_FREQ=\"monthly\"")
            com=${com}$(csed "cable_user%CASA_SPIN_STARTYEAR=1860")
            com=${com}$(csed "cable_user%CASA_SPIN_ENDYEAR=1869")
            com=${com}$(csed "cable_user%limit_labile=.false.")
            com=${com}$(csed "casafile%out=\"outputs/cru_out_casa.nc\"")
            com=${com}$(csed "casafile%cnpipool=\"restart/cru_casa_rst.nc\"")
            com=${com}$(csed "cable_user%POP_fromZero=.false.")
            com=${com}$(csed "cable_user%POP_out=\"ini\"")
            com=${com}$(csed "cable_user%POP_restart_in=\"restart/pop_cru_ini.nc\"")
            com=${com}$(csed "cable_user%POPLUC=.true.")
            com=${com}$(csed "cable_user%POPLUC_RunType=\"static\"")
            com=${com}$(csed "cable_user%LUC_outfile=\"outputs/cru_out_LUC.nc\"")
            com=${com}$(csed "cable_user%LUC_restart_in=\"restart/cru_LUC_rst.nc\"")
            com=${com}$(csed "cable_user%LUC_restart_out=\"restart/cru_LUC_rst.nc\"")
	    if [[ ${explicit_gm} -eq 1 ]] ; then
                com=${com}$(csed "cable_user%explicit_gm=.true.")
            else
                com=${com}$(csed "cable_user%explicit_gm=.false.")
            fi
            if [[ ${doc13o2} -eq 1 ]] ; then
                com=${com}$(csed "cable_user%c13o2=.true.")
		com=${com}$(csed "cable_user%c13o2_delta_atm_file=\"${filename_d13c_atm}\"")
                com=${com}$(csed "cable_user%c13o2_outfile=\"outputs/cru_out_casa_c13o2.nc\"")
                com=${com}$(csed "cable_user%c13o2_restart_in_flux=\"restart/cru_c13o2_flux_rst.nc\"")
                com=${com}$(csed "cable_user%c13o2_restart_out_flux=\"restart/cru_c13o2_flux_rst.nc\"")
                com=${com}$(csed "cable_user%c13o2_restart_in_pools=\"restart/cru_c13o2_pools_rst.nc\"")
                com=${com}$(csed "cable_user%c13o2_restart_out_pools=\"restart/cru_c13o2_pools_rst.nc\"")
                if [[ ${c13o2_simple_disc} -eq 1 ]] ; then
                    com=${com}$(csed "cable_user%c13o2_simple_disc=.true.")
                else
                    com=${com}$(csed "cable_user%c13o2_simple_disc=.false.")
                fi
            fi
            sed ${com} ${ndir}/cable.nml > ${rdir}/cable.nml
            # run model
            cd ${rdir}
	    irm logs/log_cable.txt logs/log_out_cable.txt
	    if [[ ${dompi} -eq 1 ]] ; then
		if [[ ${ignore_mpi_err} -eq 1 ]] ; then set +e ; fi
		${mpiexecdir}mpiexec -n ${nproc} ./${iexe} > logs/log_out_cable.txt
		if [[ ${ignore_mpi_err} -eq 1 ]] ; then set -e ; fi
	    else
		./${iexe} > logs/log_out_cable.txt
	    fi
            # save output
	    copyid ${rid} cable.nml cru.nml LUC.nml
	    mv *_${rid}.nml restart/
    	    cd logs
    	    renameid ${rid} log_cable.txt log_out_cable.txt
    	    cd ../restart
    	    copyid ${rid} pop_cru_ini.nc cru_climate_rst.nc cru_casa_rst.nc cru_cable_rst.nc
    	    if [[ ${doc13o2} -eq 1 ]] ; then copyid ${rid} cru_c13o2_flux_rst.nc cru_c13o2_pools_rst.nc ; fi
    	    cd ../outputs
    	    renameid ${rid} cru_out_cable.nc cru_out_casa.nc
    	    if [[ ${doc13o2} -eq 1 ]] ; then copyid ${rid} cru_out_casa_c13o2.nc ; fi
    	    cd ..
    	    cd ${pdir}
        fi
        #
        if [[ 1 -eq 1 ]] ; then
            # 4b. analytic quasi-equilibrium of biomass pools
            echo "   4b. Analytic solution of biomass pools"
	    rid="spinup_analytic"
	    # rid="spin_casa${iequi}"
            # CRU
            irm ${rdir}/cru.nml
            com=$(csed "BasePath=\"${MetPath}\",MetPath=\"${MetPath}\",LandMaskFile=\"${LandMaskFile}\"")
            com=${com}$(csed "Run=\"S0_TRENDY\"")
            sed ${com} ${ndir}/cru.nml > ${rdir}/cru.nml
            # LUC
            irm ${rdir}/LUC.nml
            com=$(csed "TransitionFilePath=\"${TransitionFilePath}\",ClimateFile=\"${ClimateFile}\"")
            com=${com}$(csed "YearStart=1700")
            com=${com}$(csed "YearEnd=2017")
            sed ${com} ${ndir}/LUC.nml > ${rdir}/LUC.nml
            # CABLE
            irm ${rdir}/cable.nml
            com=$(csed "filename%veg=\"${filename_veg}\",filename%soil=\"${filename_soil}\"")
            com=${com}$(csed "casafile%cnpbiome=\"${casafile_cnpbiome}\"")
            com=${com}$(csed "filename%restart_in=\"restart/cru_cable_rst.nc\"")
            com=${com}$(csed "cable_user%CLIMATE_fromZero=.false.")
            com=${com}$(csed "cable_user%YearStart=1840")
            com=${com}$(csed "cable_user%YearEnd=1859")
            com=${com}$(csed "icycle=12")
            com=${com}$(csed "spincasa=.true.")
            com=${com}$(csed "cable_user%CASA_fromZero=.false.")
            com=${com}$(csed "cable_user%CASA_DUMP_READ=.true.")
            com=${com}$(csed "cable_user%CASA_DUMP_WRITE=.false.")
            com=${com}$(csed "cable_user%CASA_OUT_FREQ=\"monthly\"")
            com=${com}$(csed "cable_user%CASA_SPIN_STARTYEAR=1840")
            com=${com}$(csed "cable_user%CASA_SPIN_ENDYEAR=1859")
            com=${com}$(csed "cable_user%limit_labile=.false.")
            com=${com}$(csed "casafile%out=\"outputs/cru_out_casa.nc\"")
            com=${com}$(csed "casafile%cnpipool=\"restart/cru_casa_rst.nc\"")
            com=${com}$(csed "cable_user%POP_fromZero=.false.")
            com=${com}$(csed "cable_user%POP_out=\"ini\"")
            com=${com}$(csed "cable_user%POP_restart_in=\"restart/pop_cru_ini.nc\"")
            com=${com}$(csed "cable_user%POPLUC=.true.")
            com=${com}$(csed "cable_user%POPLUC_RunType=\"static\"")
            com=${com}$(csed "cable_user%LUC_outfile=\"outputs/cru_out_LUC.nc\"")
            com=${com}$(csed "cable_user%LUC_restart_in=\"restart/cru_LUC_rst.nc\"")
            com=${com}$(csed "cable_user%LUC_restart_out=\"restart/cru_LUC_rst.nc\"")
	    if [[ ${explicit_gm} -eq 1 ]] ; then
                com=${com}$(csed "cable_user%explicit_gm=.true.")
            else
                com=${com}$(csed "cable_user%explicit_gm=.false.")
            fi
            if [[ ${doc13o2} -eq 1 ]] ; then
                com=${com}$(csed "cable_user%c13o2=.true.")
		com=${com}$(csed "cable_user%c13o2_delta_atm_file=\"${filename_d13c_atm}\"")
                com=${com}$(csed "cable_user%c13o2_outfile=\"outputs/cru_out_casa_c13o2.nc\"")
                com=${com}$(csed "cable_user%c13o2_restart_in_flux=\"restart/cru_c13o2_flux_rst.nc\"")
                com=${com}$(csed "cable_user%c13o2_restart_out_flux=\"restart/cru_c13o2_flux_rst.nc\"")
                com=${com}$(csed "cable_user%c13o2_restart_in_pools=\"restart/cru_c13o2_pools_rst.nc\"")
                com=${com}$(csed "cable_user%c13o2_restart_out_pools=\"restart/cru_c13o2_pools_rst.nc\"")
                if [[ ${c13o2_simple_disc} -eq 1 ]] ; then
                    com=${com}$(csed "cable_user%c13o2_simple_disc=.true.")
                else
                    com=${com}$(csed "cable_user%c13o2_simple_disc=.false.")
                fi
            fi
            sed ${com} ${ndir}/cable.nml > ${rdir}/cable.nml
            # run model
            cd ${rdir}
	    irm logs/log_cable.txt logs/log_out_cable.txt
	    if [[ ${dompi} -eq 1 ]] ; then
		if [[ ${ignore_mpi_err} -eq 1 ]] ; then set +e ; fi
		${mpiexecdir}mpiexec -n ${nproc} ./${iexe} > logs/log_out_cable.txt
		if [[ ${ignore_mpi_err} -eq 1 ]] ; then set -e ; fi
	    else
		./${iexe} > logs/log_out_cable.txt
	    fi
            # save output
	    copyid ${rid} cable.nml cru.nml LUC.nml
	    mv *_${rid}.nml restart/
    	    cd logs
    	    renameid ${rid} log_cable.txt log_out_cable.txt
    	    cd ../restart
    	    copyid ${rid} pop_cru_ini.nc cru_casa_rst.nc
    	    if [[ ${doc13o2} -eq 1 ]] ; then copyid ${rid} cru_c13o2_flux_rst.nc cru_c13o2_pools_rst.nc ; fi
	    if [[ ${dompi} -eq 0 ]] ; then # no output only restart if MPI
    		cd ../outputs
    		renameid ${rid} cru_out_casa.nc
    		if [[ ${doc13o2} -eq 1 ]] ; then copyid ${rid} cru_out_casa_c13o2.nc ; fi
    		cd ..
	    fi
    	    cd ${pdir}
        fi
    done
fi

# 5a. First dynamic land use
if [[ ${doiniluc} -eq 1 ]] ; then
    echo "5a. First dynamic land use"
    rid="1580_1699"
    # CRU
    irm ${rdir}/cru.nml
    com=$(csed "BasePath=\"${MetPath}\",MetPath=\"${MetPath}\",LandMaskFile=\"${LandMaskFile}\"")
    com=${com}$(csed "Run=\"S0_TRENDY\"")
    sed ${com} ${ndir}/cru.nml > ${rdir}/cru.nml
    # LUC
    irm ${rdir}/LUC.nml
    com=$(csed "TransitionFilePath=\"${TransitionFilePath}\",ClimateFile=\"${ClimateFile}\"")
    com=${com}$(csed "YearStart=1580")
    com=${com}$(csed "YearEnd=1699")
    sed ${com} ${ndir}/LUC.nml > ${rdir}/LUC.nml
    # CABLE
    irm ${rdir}/cable.nml
    com=$(csed "filename%veg=\"${filename_veg}\",filename%soil=\"${filename_soil}\"")
    com=${com}$(csed "casafile%cnpbiome=\"${casafile_cnpbiome}\"")
    com=${com}$(csed "filename%restart_in=\"restart/cru_cable_rst.nc\"")
    com=${com}$(csed "cable_user%CLIMATE_fromZero=.false.")
    com=${com}$(csed "cable_user%YearStart=1580")
    com=${com}$(csed "cable_user%YearEnd=1699")
    com=${com}$(csed "icycle=12")
    com=${com}$(csed "spincasa=.false.")
    com=${com}$(csed "cable_user%CASA_fromZero=.false.")
    com=${com}$(csed "cable_user%CASA_DUMP_READ=.true.")
    com=${com}$(csed "cable_user%CASA_DUMP_WRITE=.false.")
    com=${com}$(csed "cable_user%CASA_OUT_FREQ=\"annually\"")
    com=${com}$(csed "cable_user%CASA_SPIN_STARTYEAR=1840")
    com=${com}$(csed "cable_user%CASA_SPIN_ENDYEAR=1859")
    com=${com}$(csed "cable_user%limit_labile=.false.")
    com=${com}$(csed "casafile%out=\"outputs/cru_out_casa.nc\"")
    com=${com}$(csed "casafile%cnpipool=\"restart/cru_casa_rst.nc\"")
    com=${com}$(csed "cable_user%POP_fromZero=.false.")
    com=${com}$(csed "cable_user%POP_out=\"ini\"")
    com=${com}$(csed "cable_user%POP_restart_in=\"restart/pop_cru_ini.nc\"")
    com=${com}$(csed "cable_user%POPLUC=.true.")
    com=${com}$(csed "cable_user%POPLUC_RunType=\"init\"")
    com=${com}$(csed "cable_user%LUC_outfile=\"outputs/cru_out_LUC.nc\"")
    com=${com}$(csed "cable_user%LUC_restart_in=\"\"")
    com=${com}$(csed "cable_user%LUC_restart_out=\"restart/cru_LUC_rst.nc\"")
    if [[ ${explicit_gm} -eq 1 ]] ; then
        com=${com}$(csed "cable_user%explicit_gm=.true.")
    else
        com=${com}$(csed "cable_user%explicit_gm=.false.")
    fi
    if [[ ${doc13o2} -eq 1 ]] ; then
        com=${com}$(csed "cable_user%c13o2=.true.")
        com=${com}$(csed "cable_user%c13o2_delta_atm_file=\"${filename_d13c_atm}\"")
        com=${com}$(csed "cable_user%c13o2_outfile=\"outputs/cru_out_casa_c13o2.nc\"")
        com=${com}$(csed "cable_user%c13o2_restart_in_flux=\"restart/cru_c13o2_flux_rst.nc\"")
        com=${com}$(csed "cable_user%c13o2_restart_out_flux=\"restart/cru_c13o2_flux_rst.nc\"")
        com=${com}$(csed "cable_user%c13o2_restart_in_pools=\"restart/cru_c13o2_pools_rst.nc\"")
        com=${com}$(csed "cable_user%c13o2_restart_out_pools=\"restart/cru_c13o2_pools_rst.nc\"")
        if [[ ${c13o2_simple_disc} -eq 1 ]] ; then
            com=${com}$(csed "cable_user%c13o2_simple_disc=.true.")
        else
            com=${com}$(csed "cable_user%c13o2_simple_disc=.false.")
        fi
    fi
    sed ${com} ${ndir}/cable.nml > ${rdir}/cable.nml
    # run model
    cd ${rdir}
    irm logs/log_cable.txt logs/log_out_cable.txt
    if [[ ${dompi} -eq 1 ]] ; then
	if [[ ${ignore_mpi_err} -eq 1 ]] ; then set +e ; fi
	${mpiexecdir}mpiexec -n ${nproc} ./${iexe} > logs/log_out_cable.txt
	if [[ ${ignore_mpi_err} -eq 1 ]] ; then set -e ; fi
    else
	./${iexe} > logs/log_out_cable.txt
    fi
    # save output
    copyid ${rid} cable.nml cru.nml LUC.nml
    mv *_${rid}.nml restart/
    cd logs
    renameid ${rid} log_cable.txt log_out_cable.txt
    cd ../restart
    copyid ${rid} pop_cru_ini.nc cru_climate_rst.nc cru_casa_rst.nc cru_cable_rst.nc cru_LUC_rst.nc
    if [[ ${doc13o2} -eq 1 ]] ; then copyid ${rid} cru_c13o2_flux_rst.nc cru_c13o2_pools_rst.nc cru_c13o2_luc_rst.nc ; fi
    cd ../outputs
    # MC - Question2VH: cru_out_cable.nc is not produced
    renameid ${rid} cru_out_casa.nc cru_out_LUC.nc
    if [[ ${doc13o2} -eq 1 ]] ; then copyid ${rid} cru_out_casa_c13o2.nc ; fi
    cd ..
    cd ${pdir}
fi


# 5b. Second full dynamic spinup
if [[ ${doinidyn} -eq 1 ]] ; then
    echo "5b. Full dynamic spinup"
    rid="1700_1899"
    # CRU
    irm ${rdir}/cru.nml
    com=$(csed "BasePath=\"${MetPath}\",MetPath=\"${MetPath}\",LandMaskFile=\"${LandMaskFile}\"")
    com=${com}$(csed "Run=\"S1_TRENDY\"")
    sed ${com} ${ndir}/cru.nml > ${rdir}/cru.nml
    # LUC
    irm ${rdir}/LUC.nml
    com=$(csed "TransitionFilePath=\"${TransitionFilePath}\",ClimateFile=\"${ClimateFile}\"")
    com=${com}$(csed "YearStart=1700")
    com=${com}$(csed "YearEnd=1899")
    sed ${com} ${ndir}/LUC.nml > ${rdir}/LUC.nml
    # CABLE
    irm ${rdir}/cable.nml
    com=$(csed "filename%veg=\"${filename_veg}\",filename%soil=\"${filename_soil}\"")
    com=${com}$(csed "casafile%cnpbiome=\"${casafile_cnpbiome}\"")
    com=${com}$(csed "filename%restart_in=\"restart/cru_cable_rst.nc\"")
    com=${com}$(csed "cable_user%CLIMATE_fromZero=.false.")
    com=${com}$(csed "cable_user%YearStart=1700")
    com=${com}$(csed "cable_user%YearEnd=1899")
    com=${com}$(csed "icycle=2")
    com=${com}$(csed "spincasa=.false.")
    com=${com}$(csed "cable_user%CASA_fromZero=.false.")
    com=${com}$(csed "cable_user%CASA_DUMP_READ=.false.")
    com=${com}$(csed "cable_user%CASA_DUMP_WRITE=.true.")
    com=${com}$(csed "cable_user%CASA_OUT_FREQ=\"monthly\"")
    com=${com}$(csed "cable_user%CASA_SPIN_STARTYEAR=1850")
    com=${com}$(csed "cable_user%CASA_SPIN_ENDYEAR=1859")
    com=${com}$(csed "cable_user%limit_labile=.false.")
    com=${com}$(csed "casafile%out=\"outputs/cru_out_casa.nc\"")
    com=${com}$(csed "casafile%cnpipool=\"restart/cru_casa_rst.nc\"")
    com=${com}$(csed "cable_user%POP_fromZero=.false.")
    com=${com}$(csed "cable_user%POP_out=\"ini\"")
    com=${com}$(csed "cable_user%POP_restart_in=\"restart/pop_cru_ini.nc\"")
    com=${com}$(csed "cable_user%POPLUC=.true.")
    com=${com}$(csed "cable_user%POPLUC_RunType=\"restart\"")
    com=${com}$(csed "cable_user%LUC_outfile=\"outputs/cru_out_LUC.nc\"")
    com=${com}$(csed "cable_user%LUC_restart_in=\"restart/cru_LUC_rst.nc\"")
    com=${com}$(csed "cable_user%LUC_restart_out=\"restart/cru_LUC_rst.nc\"")
    if [[ ${explicit_gm} -eq 1 ]] ; then
        com=${com}$(csed "cable_user%explicit_gm=.true.")
    else
        com=${com}$(csed "cable_user%explicit_gm=.false.")
    fi
    if [[ ${doc13o2} -eq 1 ]] ; then
        com=${com}$(csed "cable_user%c13o2=.true.")
        com=${com}$(csed "cable_user%c13o2_delta_atm_file=\"${filename_d13c_atm}\"")
        com=${com}$(csed "cable_user%c13o2_outfile=\"outputs/cru_out_casa_c13o2.nc\"")
        com=${com}$(csed "cable_user%c13o2_restart_in_flux=\"restart/cru_c13o2_flux_rst.nc\"")
        com=${com}$(csed "cable_user%c13o2_restart_out_flux=\"restart/cru_c13o2_flux_rst.nc\"")
        com=${com}$(csed "cable_user%c13o2_restart_in_pools=\"restart/cru_c13o2_pools_rst.nc\"")
        com=${com}$(csed "cable_user%c13o2_restart_out_pools=\"restart/cru_c13o2_pools_rst.nc\"")
        com=${com}$(csed "cable_user%c13o2_restart_in_luc=\"restart/cru_c13o2_luc_rst.nc\"")
        com=${com}$(csed "cable_user%c13o2_restart_out_luc=\"restart/cru_c13o2_luc_rst.nc\"")
        if [[ ${c13o2_simple_disc} -eq 1 ]] ; then
            com=${com}$(csed "cable_user%c13o2_simple_disc=.true.")
        else
            com=${com}$(csed "cable_user%c13o2_simple_disc=.false.")
        fi
    fi
    sed ${com} ${ndir}/cable.nml > ${rdir}/cable.nml
    # run model
    cd ${rdir}
    irm logs/log_cable.txt logs/log_out_cable.txt
    if [[ ${dompi} -eq 1 ]] ; then
	if [[ ${ignore_mpi_err} -eq 1 ]] ; then set +e ; fi
	${mpiexecdir}mpiexec -n ${nproc} ./${iexe} > logs/log_out_cable.txt
	if [[ ${ignore_mpi_err} -eq 1 ]] ; then set -e ; fi
    else
	./${iexe} > logs/log_out_cable.txt
    fi
    # save output
    copyid ${rid} cable.nml cru.nml LUC.nml
    mv *_${rid}.nml restart/
    cd logs
    renameid ${rid} log_cable.txt log_out_cable.txt
    cd ../restart
    copyid ${rid} pop_cru_ini.nc cru_climate_rst.nc cru_casa_rst.nc cru_cable_rst.nc cru_LUC_rst.nc
    if [[ ${doc13o2} -eq 1 ]] ; then copyid ${rid} cru_c13o2_flux_rst.nc cru_c13o2_pools_rst.nc cru_c13o2_luc_rst.nc ; fi
    cd ../outputs
    renameid ${rid} cru_out_cable.nc cru_out_casa.nc cru_out_LUC.nc
    if [[ ${doc13o2} -eq 1 ]] ; then copyid ${rid} cru_out_casa_c13o2.nc ; fi
    cd ..
    cd ${pdir}
fi


# 6. Final centennial run
if [[ ${dofinal} -eq 1 ]] ; then
    echo "6. Final centennial run"
    rid="1900_2017"
    # CRU
    irm ${rdir}/cru.nml
    com=$(csed "BasePath=\"${MetPath}\",MetPath=\"${MetPath}\",LandMaskFile=\"${LandMaskFile}\"")
    com=${com}$(csed "Run=\"S2_TRENDY\"")
    sed ${com} ${ndir}/cru.nml > ${rdir}/cru.nml
    # LUC
    irm ${rdir}/LUC.nml
    com=$(csed "TransitionFilePath=\"${TransitionFilePath}\",ClimateFile=\"${ClimateFile}\"")
    com=${com}$(csed "YearStart=1900")
    com=${com}$(csed "YearEnd=2017")
    sed ${com} ${ndir}/LUC.nml > ${rdir}/LUC.nml
    # CABLE
    irm ${rdir}/cable.nml
    com=$(csed "filename%veg=\"${filename_veg}\",filename%soil=\"${filename_soil}\"")
    com=${com}$(csed "casafile%cnpbiome=\"${casafile_cnpbiome}\"")
    com=${com}$(csed "filename%restart_in=\"restart/cru_cable_rst.nc\"")
    com=${com}$(csed "cable_user%CLIMATE_fromZero=.false.")
    com=${com}$(csed "cable_user%YearStart=1901")
    com=${com}$(csed "cable_user%YearEnd=2017")
    com=${com}$(csed "icycle=2")
    com=${com}$(csed "spincasa=.false.")
    com=${com}$(csed "cable_user%CASA_fromZero=.false.")
    com=${com}$(csed "cable_user%CASA_DUMP_READ=.false.")
    com=${com}$(csed "cable_user%CASA_DUMP_WRITE=.false.")
    com=${com}$(csed "cable_user%CASA_OUT_FREQ=\"monthly\"")
    com=${com}$(csed "cable_user%CASA_SPIN_STARTYEAR=1850")
    com=${com}$(csed "cable_user%CASA_SPIN_ENDYEAR=1859")
    com=${com}$(csed "cable_user%limit_labile=.false.")
    com=${com}$(csed "casafile%out=\"outputs/cru_out_casa.nc\"")
    com=${com}$(csed "casafile%cnpipool=\"restart/cru_casa_rst.nc\"")
    com=${com}$(csed "cable_user%POP_fromZero=.false.")
    com=${com}$(csed "cable_user%POP_out=\"ini\"")
    com=${com}$(csed "cable_user%POP_restart_in=\"restart/pop_cru_ini.nc\"")
    com=${com}$(csed "cable_user%POPLUC=.true.")
    com=${com}$(csed "cable_user%POPLUC_RunType=\"restart\"")
    com=${com}$(csed "cable_user%LUC_outfile=\"outputs/cru_out_LUC.nc\"")
    com=${com}$(csed "cable_user%LUC_restart_in=\"restart/cru_LUC_rst.nc\"")
    com=${com}$(csed "cable_user%LUC_restart_out=\"restart/cru_LUC_rst.nc\"")
    if [[ ${explicit_gm} -eq 1 ]] ; then
        com=${com}$(csed "cable_user%explicit_gm=.true.")
    else
        com=${com}$(csed "cable_user%explicit_gm=.false.")
    fi
    if [[ ${doc13o2} -eq 1 ]] ; then
        com=${com}$(csed "cable_user%c13o2=.true.")
        com=${com}$(csed "cable_user%c13o2_delta_atm_file=\"${filename_d13c_atm}\"")
        com=${com}$(csed "cable_user%c13o2_outfile=\"outputs/cru_out_casa_c13o2.nc\"")
        com=${com}$(csed "cable_user%c13o2_restart_in_flux=\"restart/cru_c13o2_flux_rst.nc\"")
        com=${com}$(csed "cable_user%c13o2_restart_out_flux=\"restart/cru_c13o2_flux_rst.nc\"")
        com=${com}$(csed "cable_user%c13o2_restart_in_pools=\"restart/cru_c13o2_pools_rst.nc\"")
        com=${com}$(csed "cable_user%c13o2_restart_out_pools=\"restart/cru_c13o2_pools_rst.nc\"")
        com=${com}$(csed "cable_user%c13o2_restart_in_luc=\"restart/cru_c13o2_luc_rst.nc\"")
        com=${com}$(csed "cable_user%c13o2_restart_out_luc=\"restart/cru_c13o2_luc_rst.nc\"")
        if [[ ${c13o2_simple_disc} -eq 1 ]] ; then
            com=${com}$(csed "cable_user%c13o2_simple_disc=.true.")
        else
            com=${com}$(csed "cable_user%c13o2_simple_disc=.false.")
        fi
    fi
    sed ${com} ${ndir}/cable.nml > ${rdir}/cable.nml
    # run model
    cd ${rdir}
    irm logs/log_cable.txt logs/log_out_cable.txt
    if [[ ${dompi} -eq 1 ]] ; then
	if [[ ${ignore_mpi_err} -eq 1 ]] ; then set +e ; fi
	${mpiexecdir}mpiexec -n ${nproc} ./${iexe} > logs/log_out_cable.txt
	if [[ ${ignore_mpi_err} -eq 1 ]] ; then set -e ; fi
    else
	./${iexe} > logs/log_out_cable.txt
    fi
    # save output
    copyid ${rid} cable.nml cru.nml LUC.nml
    mv *_${rid}.nml restart/
    cd logs
    renameid ${rid} log_cable.txt log_out_cable.txt
    cd ../restart
    copyid ${rid} pop_cru_ini.nc cru_climate_rst.nc cru_casa_rst.nc cru_cable_rst.nc cru_LUC_rst.nc
    if [[ ${doc13o2} -eq 1 ]] ; then copyid ${rid} cru_c13o2_flux_rst.nc cru_c13o2_pools_rst.nc cru_c13o2_luc_rst.nc ; fi
    cd ../outputs
    renameid ${rid} cru_out_cable.nc cru_out_casa.nc cru_out_LUC.nc
    if [[ ${doc13o2} -eq 1 ]] ; then copyid ${rid} cru_out_casa_c13o2.nc ; fi
    cd ..
    cd ${pdir}
fi


# --------------------------------------------------------------------
# Finish
#
cd ${isdir}

t2=$(date +%s)
dt=$((t2-t1))
if [[ ${dt} -lt 60 ]] ; then
    printf "Finished at %s   in %i seconds.\n" "$(date)" ${dt}
else
    dm=$(echo "(${t2}-${t1})/60." | bc -l)
    printf "Finished at %s   in %.2f minutes.\n" "$(date)" ${dm}
fi

exit 0

# grep -Ei '(filename%restart_in|cable_user%CLIMATE_fromZero|cable_user%YearStart|cable_user%YearEnd|icycle|spincasa|cable_user%CASA_fromZero|cable_user%CASA_DUMP_READ|cable_user%CASA_DUMP_WRITE|cable_user%CASA_OUT_FREQ|cable_user%CASA_SPIN_STARTYEAR|cable_user%CASA_SPIN_ENDYEAR|cable_user%limit_labile|casafile%out|casafile%cnpipool|cable_user%POP_fromZero|cable_user%POP_out|cable_user%POP_restart_in|cable_user%POPLUC|cable_user%POPLUC_RunType|cable_user%LUC_outfile|cable_user%LUC_restart_in|cable_user%LUC_restart_out)' ../driver_files/cable.nml.1900_2017 | pbcopy
