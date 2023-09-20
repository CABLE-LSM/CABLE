#!/usr/bin/env bash

# Gadi
# https://opus.nci.org.au/display/Help/How+to+submit+a+job
#PBS -N TRENDY_S0
#PBS -P vl59
#PBS -q normal
#PBS -l walltime=48:00:00
#PBS -l mem=192GB
#PBS -l ncpus=96
#PBS -l storage=gdata/vl59+gdata/x45
#PBS -l software=netCDF:MPI:Intel:GNU
#PBS -r y
#PBS -l wd
#PBS -j oe
#PBS -S /bin/bash

# cuntz@explor, cuntz@mc16, cuntz@mcinra, moc801@gadi cuntz@gadi
# jk8585@gadi knauer@gadi
system=jk8585@gadi

# MPI run or single processor run
# nproc should fit with job tasks
dompi=1      # 0: normal run: ./cable
             # 1: MPI run: mpiexec -n ${nproc} ./cable_mpi
nproc=96     # Number of cores for MPI runs
             # must be same as above: PBS -l ncpus=nproc
run_type=0   # 0: Full run (global)
             # 1: Test Suite run on 1000pts
             # 2: Full test run on smaller spatial domain

# --------------------------------------------------------------------
#
# Full Cable run with biomass spinup, POP, land-use change, etc.
#
# This script uses CRU-JRA forcing.
#
# Global meteo and land-use change data can be used with a mask giving land points.
# In step 0, the land mask can be extracted for a single point, an area
# or a number of random points can be chosen.
# Alternatively, single site met, LUH2, forcing and mask can be extracted from the global data sets.
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
#   5. Second phase of spinup with dynamic land use, atmospheric CO2 and N deposition.
#      a) Dynamic land use from 1580 to 1699, using still fixed atmospheric CO2 and
#         N deposition from 1700, and 30 years of repeated meteorology.
#      b) Run from 1700 to 1899 with dynmic land use, varying atmospheric CO2 and N deposition,
#         but still with 30 years of repeated meteorology.
#   6. Final historical run, everything dynamic from 1900 to 2020.
#
# Written,  Matthias Cuntz, Aug 2019, following the run scripts and namelists provided by Vanessa Haverd
# Modified, Jurgen Knauer, 2020      - gm_explicit, coordination, acclimation
#                                    - bios, plume, future runs
#           Matthias Cuntz, Mar 2021 - functions into run_cable-pop_lib.sh
#
# --------------------------------------------------------------------


#---------------------------------------------------------------------
# Spatial extent
# --------------------------------------------------------------------
# global runs:   imeteo=0, doextractsite=0
# regional runs: imeteo=1, doextractsite=1, latlon or randompoints must be specified

imeteo=0        # 0: Use global meteo, land use and mask
                # 1: Use local mask, but global meteo and land use (doextractsite=1)
                # 2: Use local meteo, land use and mask (doextractsite=2)
doextractsite=0 # 0: Do not extract local meteo, land use nor mask
                # 1: Do extract only mask at specific site/region (imeteo=1)
                # 2: Do extract meteo, land use and mask at specific site/region (imeteo=2)
                #    Does not work with randompoints /= 0 but with latlon
randompoints=0  # <0: use -1*randompoints from file ${LandMaskFilePath}/${experiment}_points.csv if existing
                # 0:  use latlon
                # >0: generate and use randompoints random grid points from GlobalLandMaskFile
# lat,lon  or  latmin,latmax,lonmin,lonmax   # must have . in numbers otherwise indexes taken
purge_restart=1 # 1/0: Do/Do not delete all restart files (completely new run, e.g. if settings changed)
purge_spinup=0  # 1/0: Do/Do not delete all spinup outputs (only used if run_type -eq 1)


# --------------------------------------------------------------------
# TRENDY experiment
# --------------------------------------------------------------------
## for TRENDY, the steps run depend on the experiment
# Steps:
# 1) doclimate   1/0  Do/Do not create climate restart file
# 2) dofromzero  1/0  Do/Do not first spinup phase from zero biomass stocks
# 3) doequi1     1/0  Do/Do not bring biomass stocks into quasi-equilibrium with unrestricted P and N pools
#    nequi1           number of times to repeat steps in doequi1
# 4) doequi2     1/0  Do/Do not bring biomass stocks into quasi-equilibrium with restricted P and N pools
#    nequi2           number of times to repeat steps in doequi2
# 5) doiniluc    1/0  Do/Do not spinup with dynamic land use (initialise land use)
# 6) doinidyn    1/0  Do/Do not full dynamic spinup (transient run) from 1701 to 1900
# 7) dofinal     1/0  Do/Do not final run from 1901 to 2019


# v10: TRENDY runs must be conducted in 3 steps:
# 1) run spinup
# 2) run transient run (1700-1900) with restart files copied from spinup directory
#    note: restart files are always taken from the 'spinup' folder.
#          If restart files are taken from a different folder, those must be copied to the 'spinup' folder first.
# 3) run historic run (1901-2020)

# v11: TRENDY conducted in 5 steps, see 'run_step' options below.

experiment="S0"        # TRENDY experiment (spinup, S0, S1, S2, S3, S4, S5, S6)     
run_step="1901-2021"   # if runtype -eq 0 (full run), the following steps are implemented:
                       # spinup1
                       # spinup2
                       # spinup3
                       # 1700-1800
                       # 1801-1900
                       # 1901-2021
name_ext=""            # extension to experiment name (only used for $sitepath,
                                   # i.e. for name of output folder. Should begin with a "_"
                                   # or leave empty "". Only used if run_type .gt. 0

if [[ $run_type -eq 1 ]] ; then # do all the steps for a test suite run on 1000 pts
    
    experiment_name=${experiment}_1000pts${name_ext}
    purge_restart=1
    
    doclimate=1     
    dofromzero=1
    doequi1=1      
    nequi1=4        
    doequi2=1       
    nequi2=5
    if [[ "${experiment}" == "S3" ]] ; then
	doiniluc=1
    else
	doiniluc=0
    fi
    doinidyn=1      
    dofinal=1
    
elif [[ $run_type -eq 2 ]] ; then # full test run on a smaller spatial domain
    
    experiment_name=${experiment}_test${name_ext}
    imeteo=1
    doextractsite=1
    #latlon=-34.5,-33.5,149.0,156.5
    #latlon=-44.0,-10.0,110.0,155.0  # Australia
    #latlon=-44.0,-10.0,145.0,155.0  # East Australia
    latlon=44.0,50.0,0.0,6.0   # Middle Europe
    
    doclimate=1     # 1/0: Do/Do not create climate restart file
    dofromzero=1
    doequi1=1       # 1/0: Do/Do not bring biomass stocks into quasi-equilibrium with restricted P and N pools
    nequi1=2        #      number of times to repeat steps in doequi1
    doequi2=1       # 1/0: Do/Do not bring biomass stocks into quasi-equilibrium with unrestricted P and N pools
    nequi2=4        #      number of times to repeat steps in doequi2
    if [[ "${experiment}" == "S3" ]] ; then
	doiniluc=1      # 1/0: Do/Do not spinup with dynamic land use (initialise land use)
    else
	doiniluc=0
    fi
    doinidyn=1      # 1/0: Do/Do not full dynamic spinup (transient run) from 1701 to 1899
    dofinal=1       # 1/0: Do/Do not final run from 1901 to 2021
    
elif [[ $run_type -eq 0 ]] ; then # full run
    
    experiment_name=${experiment}
    
    if [[ "${run_step}" == "spinup1" ]] ; then
       purge_restart=1 # 1/0: Do/Do not delete all restart files

       doclimate=1     
       dofromzero=1
       doequi1=1       
       nequi1=2        
       doequi2=0       
       nequi2=0       
       doiniluc=0   
       doinidyn=0      
       dofinal=0
    elif [[ "${run_step}" == "spinup2" ]] ; then
       purge_restart=0 # 1/0: Do/Do not delete all restart files

       doclimate=0     
       dofromzero=0
       doequi1=1       
       nequi1=4  # but starting at iequi1=3        
       doequi2=0       
       nequi2=0       
       doiniluc=0   
       doinidyn=0      
       dofinal=0
    elif [[ "${run_step}" == "spinup3" ]] ; then
       purge_restart=0 # 1/0: Do/Do not delete all restart files

       doclimate=0     
       dofromzero=0
       doequi1=0       
       nequi1=0        
       doequi2=1       
       nequi2=5       
       doiniluc=0   
       doinidyn=0      
       dofinal=0
    elif [[ "${run_step}" == "1700-1800" ]] ; then
       purge_restart=0 # 1/0: Do/Do not delete all restart files

       doclimate=0     
       dofromzero=0
       doequi1=0       
       nequi1=0        
       doequi2=0       
       nequi2=0
       if [[ "${experiment}" == "S3" ]] ; then
	  doiniluc=1
       else
	  doiniluc=0
       fi   
       doinidyn=1      
       dofinal=0
    elif [[ "${run_step}" == "1801-1900" ]] ; then
       purge_restart=0 # 1/0: Do/Do not delete all restart files

       doclimate=0     
       dofromzero=0
       doequi1=0       
       nequi1=0        
       doequi2=0       
       nequi2=0
       doiniluc=0
       doinidyn=1      
       dofinal=0
    elif [[ "${run_step}" == "1901-2021" ]] ; then
       purge_restart=0 # 1/0: Do/Do not delete all restart files

       doclimate=0     
       dofromzero=0
       doequi1=0       
       nequi1=0        
       doequi2=0       
       nequi2=0
       doiniluc=0
       doinidyn=0      
       dofinal=1
    fi # run_step
fi # $run_type


# MetType
mettype="cru"       # "cru", "plume", "bios"

# Cable
read_fdiff=1        # 1/0: do/do not read in diffuse radiation fraction
call_blaze=0        # 1/0: do/do not call BLAZE
explicit_gm=0       # 1/0: explicit (finite) or implicit mesophyll conductance
use_LUTgm=1         # 1/0: Do/Do not use lookup table for parameter conversion accounting for gm (only used if explicit_gm=1)
Rubisco_params="Bernacchi_2002"   # "Bernacchi_2002" or "Walker_2013"
coordinate_photosyn=1 # 1/0: Do/Do not coordinate photosynthesis
coord=F               # T/F: version of photosyn. optimisation (optimised(F) or forced (T))
acclimate_photosyn=1  # 1/0: Do/Do not acclimate photosynthesis
call_pop=1          # 1/0: Do/Do not use POP population dynamics model, coupled to CASA
doc13o2=0           # 1/0: Do/Do not calculate 13C
c13o2_simple_disc=0 # 1/0: simple or full 13C leaf discrimination

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

system=$(echo ${system} | tr A-Z a-z)
sys=${system#*@}
user=${system%@*}

#
# Special things on specific computer system such as loading modules
#

export mpiexecdir=
if [[ "${sys}" == "explor" ]] ; then
    # prog is slurm_script
    pdir=${isdir}
    # # INTELMPI - load mpi module first, otherwise intel module will not prepend LD_LIBRARY_PATH
    # module load intelmpi/2018.5.274
    # module load intel/2018.5
    # unset I_MPI_PMI_LIBRARY
    # export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HOME}/local/lib:${HOME}/local/netcdf-fortran-4.4.4-ifort2018.0/lib
    # export mpiexecdir=/soft/env/soft/all/intel/2018.3/compilers_and_libraries_2018.5.274/linux/mpi/intel64/bin
    # # INTEL / OpenMPI - load mpi module first, otherwise intel module will not pre-pend LD_LIBRARY_PATH
    # module load openmpi/3.0.0/intel18
    # module load intel/2018.5
    # export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HOME}/local/lib:${HOME}/local/netcdf-fortran-4.4.4-ifort2018.0/lib
    # export mpiexecdir=/opt/soft/hf/openmpi-3.0.0-intel18/bin
    # GNU / OpenMPI
    module load gcc/6.3.0
    module load openmpi/3.0.1/gcc/6.3.0
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HOME}/local/lib:${HOME}/local/netcdf-fortran-4.4.4-gfortran63/lib
    export mpiexecdir=/opt/soft/hf/openmpi/3.0.1/gcc/6.3.0/bin
    if [[ ${doextractsite} -ge 1 ]] ; then module load python/intel/2019/3 ; fi
elif [[ "${sys}" == "mc16" ]] ; then
    # export mpiexecdir=/usr/local/openmpi-4.0.4-gfortran/bin
    export mpiexecdir=/usr/local/openmpi-4.0.5-ifort/bin
elif [[ "${sys}" == "mcinra" ]] ; then
    export mpiexecdir=/usr/local/openmpi-3.1.4-gfortran/bin
    # export mpiexecdir=/usr/local/openmpi-3.1.5-ifort/bin
elif [[ "${sys}" == "pearcey" ]] ; then
    # prog is slurm_script
    pdir=${isdir}
    module del intel-cc intel-fc
    module add intel-cc/16.0.1.150 intel-fc/16.0.1.150
    module unload intel-mpi/5.0.1.035
    module add netcdf/4.3.3.1 openmpi/1.8.8
elif [[ "${sys}" == "gadi" ]] ; then
    pdir=${isdir}
    . /etc/bashrc
    module purge
    #module load intel-compiler/2019.5.281
    #module load intel-mpi/2019.5.281
    #module load netcdf/4.6.3
    module load intel-compiler/2021.5.0
    module load intel-mpi/2021.5.1
    module load netcdf/4.8.0
    # module load hdf5/1.10.5
    if [[ ${doextractsite} -ge 1 ]] ; then
        module load python3/3.7.4
	unset PYTHONPATH
	export PYTHONPATH=/g/data/x45/python/lib/python3.7/site-packages:/g/data/x45/intelpython/lib/python3.7/site-packages
        #export PYTHONPATH=${PYTHONPATH}:/g/data/x45/python/lib/python3.7/site-packages  # TODO: check python avail. on x45
    fi
    if [[ ${randompoints} -eq 0 ]] ; then module load nco/4.9.2 ; fi  # needed for cropping outputs
    #export mpiexecdir=/apps/intel-mpi/2019.5.281/intel64/bin
    export mpiexecdir=/apps/intel-mpi/2021.5.1/bin
fi
if [[ ! -z ${mpiexecdir} ]] ; then export mpiexecdir="${mpiexecdir}/" ; fi

#
# Directories of things
#

#
# Relative directories must be relative to the directory of this script,
#   not relative to the directory from which this script is launched (if different)
#   nor relative to the run path.
#
if [[ "${system}" == "cuntz@explor" ]] ; then
    # Run directory: runpath="${sitepath}/run_xxx"
    cablebase="/home/oqx29/zzy20/prog/cable"
    sitepath="${cablebase}/runs/single_sites/${experiment}"
    cablehome="${cablebase}/branches/NESP2pt9_BLAZE"
    # Cable executable
    if [[ ${dompi} -eq 1 ]] ; then
        exe="${cablehome}/offline/cable-mpi"
    else
        exe="${cablehome}/offline/cable"
    fi
    # CABLE-AUX directory (uses offline/gridinfo_CSIRO_1x1.nc and offline/modis_phenology_csiro.txt)
    aux="${cablebase}/CABLE-AUX"
    # Global Mask - for create_landmask.py
    GlobalLandMaskFile="${aux}/offline/gridinfo_CSIRO_1x1.nc"
    SurfaceFile="${aux}/offline/gridinfo_CSIRO_1x1.nc"   # note that SurfaceFile does not need subsetting
    # # Global Mask - for global run
    # GlobalLandMaskFile="/home/oqx29/zzy20/data/crujra/daily_1deg/glob_ipsl_1x1.nc"
    # Global CRU
    GlobalMetPath="/home/oqx29/zzy20/data/crujra/daily_1deg"
    # Global LUC
    GlobalTransitionFilePath="/home/oqx29/zzy20/data/LUH2_v3_1deg/"
elif [[ "${system}" == "cuntz@mc16" || "${system}" == "cuntz@mcinra" ]] ; then
    # Run directory: runpath="${sitepath}/run_xxx"
    cablebase="/Users/cuntz/prog/vanessa/cable"
    sitepath="${cablebase}/runs/single_sites/${experiment}"
    cablehome="${cablebase}/branches/NESP2pt9_BLAZE"
    # Cable executable
    if [[ ${dompi} -eq 1 ]] ; then
        # exe="${cablehome}/branches/NESP2pt9_BLAZE/offline/cable-mpi-gfortran"
        exe="${cablehome}/offline/cable-mpi-ifort"
    else
        exe="${cablehome}/offline/cable"
    fi
    # CABLE-AUX directory (uses offline/gridinfo_CSIRO_1x1.nc and offline/modis_phenology_csiro.txt)
    aux="${cablebase}/CABLE-AUX"
    # Global Mask, CRU, LUC
    GlobalLandMaskFile="${aux}/offline/gridinfo_CSIRO_1x1.nc"
    SurfaceFile="${aux}/offline/gridinfo_CSIRO_1x1.nc"   # note that SurfaceFile does not need subsetting
    GlobalMetPath=
    GlobalTransitionFilePath=
elif [[ "${system}" == "moc801@gadi" || "${system}" == "cuntz@gadi" ]] ; then
    # Run directory: runpath="${sitepath}/run_xxx"
    # sitepath="/home/801/moc801/prog/cable/runs/single_sites/${experiment}"
    sitepath="/scratch/x45/moc801/cable/c13"
    cablehome="/home/801/moc801/prog/cable/branches/NESP2pt9_BLAZE"
    # Cable executable
    if [[ ${dompi} -eq 1 ]] ; then
        exe="${cablehome}/offline/cable-mpi"
    else
        exe="${cablehome}/offline/cable"
    fi
    # CABLE-AUX directory (uses offline/gridinfo_CSIRO_1x1.nc and offline/modis_phenology_csiro.txt)
    aux="/g/data/x45/CABLE-AUX"
    # Global Mask
    # GlobalLandMaskFile="${aux}/offline/gridinfo_CSIRO_1x1.nc"
    GlobalLandMaskFile="/g/data/x45/ipbes/masks/glob_ipsl_1x1.nc"
    # Global CRU
    GlobalMetPath="/g/data/x45/crujra/daily_1deg"
    # Global LUC
    GlobalTransitionFilePath="/g/data/x45/LUH2/GCB_2018/1deg/EXTRACT"
elif [[ "${system}" == "jk8585@gadi" || "${system}" == "knauer@gadi" ]] ; then
    # Run directory: runpath="${sitepath}/run"
    sitepath="/g/data/vl59/TRENDY_v11/${experiment_name}"
    cablecode="/home/599/jk8585/CABLE_code/SHARE/CABLE-POP_TRENDY"
    cablehome="/home/599/jk8585/CABLE_run/TRENDY_v11"
    # Cable executable
    if [[ ${dompi} -eq 1 ]] ; then
	exe=${cablehome}/exes/cable-mpi
        #exe="${cablecode}/offline/cable-mpi"
    else
	exe=${cablehome}/exes/cable
        #exe="${cablecode}/offline/cable"
    fi
    # CABLE-AUX directory (uses offline/gridinfo_CSIRO_1x1.nc and offline/modis_phenology_csiro.txt)
    aux="/g/data/vl59/TRENDY_v11/aux"
    # Global Mask
    SurfaceFile="${aux}/gridinfo_CSIRO_1x1.nc"   # note that SurfaceFile does not need subsetting
    # Global Met
    GlobalMetPath="/g/data/x45/CRUJRA2022/daily_1deg"
    # Global Land Mask
    if [[ ${run_type} -eq 1 ]] ; then
	GlobalLandMaskFile="${aux}/landmasks/landmask_1000pts.nc"
    else
	GlobalLandMaskFile="${aux}/landmasks/glob_ipsl_1x1.nc"
    fi
    # Global LUC
    #GlobalTransitionFilePath="/g/data/x45/LUH2/GCB_2021/1deg/EXTRACT"
    GlobalTransitionFilePath="/g/data/x45/LUH2/GCB_2022/1deg/EXTRACT"
else
    echo "System not known."
    exit 1
fi

# Run directory
runpath="${sitepath}"

# Cable parameters
namelistpath="${cablehome}/namelists"
filename_veg="${cablehome}/params/def_veg_params.txt"
filename_soil="${cablehome}/params/def_soil_params.txt"
casafile_cnpbiome="${cablehome}/params/pftlookup.csv"
# Other scripts
ScriptsPath="${cablehome}/scripts"
# Mask (not used unless imeteo > 0)
LandMaskFile="${sitepath}/mask/${experiment}_landmask.nc"
# CRU (not used unless imeteo > 0)
MetPath="${sitepath}/met/cru_jra_1deg"

# changes for TRENDY >= v11: ClimateFile always created!
# ClimateFile="/g/data/x45/ipbes/cable_climate/ipsl_climate_rst_glob_1deg.nc"
ClimateFile="$(dirname ${sitepath})/climate_restart/cru_climate_rst.nc"

# LUC
#TransitionFilePath="${sitepath}/LUH2/v3/1deg"
# gm lookup tables
gm_lut_bernacchi_2002="${cablehome}/params/gm_LUT_351x3601x7_1pt8245_Bernacchi2002.nc"
gm_lut_walker_2013="${cablehome}/params/gm_LUT_351x3601x7_1pt8245_Walker2013.nc"
# 13C
filename_d13c_atm="${cablehome}/params/graven_et_al_gmd_2017-table_s1-delta_13c-1700-2025.txt"

# --------------------------------------------------------------------
# Start Script
# --------------------------------------------------------------------
echo $MetPath
# --------------------------------------------------------------------
# Helper functions, most functions are in plumber_cable-pop_lib.sh
#

source ${pdir}/run_cable-pop_lib.sh

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

# cleanup at end or at trap
function cleanup()
{
    \rm -f ${tmp}/*.${pid}*
    exit 1
}

# returns argument to extract lat and lon with ncks
function nckslatlon()
{
    vars=$(ncvarlist ${1})
    if [[ -z $(isin latitude ${vars}) ]] ; then ilat="lat" ; else ilat="latitude" ; fi
    if [[ -z $(isin longitude ${vars}) ]] ; then ilon="lon" ; else ilon="longitude" ; fi
    if [[ -z $(echo ${2} | cut -f 3 -d ",") || -z $(echo ${2} | cut -f 4 -d ",") ]] ; then
        iilat=$(echo ${2} | cut -f 1 -d ",")
        iilon=$(echo ${2} | cut -f 2 -d ",")
        echo "-d ${ilat},${iilat} -d ${ilon},${iilon}"
    else
        iilat1=$(echo ${2} | cut -f 1 -d ",")
        iilat2=$(echo ${2} | cut -f 2 -d ",")
        iilon1=$(echo ${2} | cut -f 3 -d ",")
        iilon2=$(echo ${2} | cut -f 4 -d ",")
        echo "-d ${ilat},${iilat1},${iilat2} -d ${ilon},${iilon1},${iilon2}"
    fi
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
mkdir -p ${sitepath}/mask
pdir=$(abspath ${pdir})
cd ${pdir}
adir=$(abspath ${aux})
exe=$(absfile ${exe})
mkdir -p ${runpath}
rdir=$(abspath ${runpath})
ndir=$(abspath ${namelistpath})
sdir=$(abspath ${ScriptsPath})

#
# prepare run directory
cd ${rdir}
mkdir -p logs
mkdir -p outputs
mkdir -p restart
ln -sf ${adir}
# ln -sf ${exe}
cp ${exe} ./
iexe=$(basename ${exe})
cd ${pdir}

#
# set stacksize to unlimited if permitted, otherwise to 15 bit if possible
set +e
ulimit -s unlimited 2> /dev/null || ulimit -s 32768
set -e

# --------------------------------------------------------------------
# Info
#
t1=$(date +%s)
printf "Started at %s\n" "$(date)"

printf "\nSetup\n"
printf "    Serial / Parallel\n"
printf "        dompi=${dompi}\n"
printf "            nproc=${nproc}\n"
printf "\n"
printf "    Sequence\n"
printf "        imeteo=${imeteo}\n"
printf "        doextractsite=${doextractsite}\n"
printf "            experiment=${experiment}\n"
printf "            randompoints=${randompoints}\n"
printf "            latlon=${latlon}\n"
printf "        doclimate=${doclimate}\n"
printf "        dofromzero=${dofromzero}\n"
printf "        doequi1=${doequi1}\n"
printf "            nequi1=${nequi1}\n"
printf "        doequi2=${doequi2}\n"
printf "            nequi2=${nequi2}\n"
printf "        doiniluc=${doiniluc}\n"
printf "        doinidyn=${doinidyn}\n"
printf "        dofinal=${dofinal}\n"
printf "        dofuture=${dofuture}\n"
printf "\n"
printf "    Options\n"
printf "        mettype=${mettype}\n"
printf "        metmodel=${metmodel}\n"
printf "        RCP=${RCP}\n"
printf "        explicit_gm=${explicit_gm}\n"
printf "        use_LUTgm=${use_LUTgm}\n"
printf "        Rubisco_params=${Rubisco_params}\n"
printf "        coordinate_photosyn=${coordinate_photosyn}\n"
printf "        coord=${coord}\n"
printf "        acclimate_photosyn=${acclimate_photosyn}\n"
printf "        call_pop=${call_pop}\n"
printf "        doc13o2=${doc13o2}\n"
printf "        c13o2_simple_disc=${c13o2_simple_disc}\n"
printf "\n"
printf "    Directories\n"
printf "        sitepath=${sitepath}\n"
printf "        cablehome=${cablehome}\n"
printf "        exe=${exe}\n"
printf "        aux=${aux}\n"
printf "        GlobalLandMaskFile=${GlobalLandMaskFile}\n"
printf "        SurfaceFile=${SurfaceFile}\n"
printf "        GlobalMetPath=${GlobalMetPath}\n"
printf "        GlobalTransitionFilePath=${GlobalTransitionFilePath}\n"
printf "        runpath=${runpath}\n"
printf "        namelistpath=${namelistpath}\n"
printf "        filename_veg=${filename_veg}\n"
printf "        filename_soil=${filename_soil}\n"
printf "        casafile_cnpbiome=${casafile_cnpbiome}\n"
printf "        LandMaskFile=${LandMaskFile}\n"
printf "        MetPath=${MetPath}\n"
printf "        ClimateFile=${ClimateFile}\n"
printf "        TransitionFilePath=${TransitionFilePath}\n"
printf "        gm_lut_bernacchi_2002=${gm_lut_bernacchi_2002}\n"
printf "        gm_lut_walker_2013=${gm_lut_walker_2013}\n"
printf "        filename_d13c_atm=${filename_d13c_atm}\n"
printf "\n"

# --------------------------------------------------------------------
# Prep input
#

# 0. Extract meteo, land use and mask for one specific site from global files
if [[ ${doextractsite} -ge 2 ]] ; then
    set +e
    xnco=$(which ncks 2> /dev/null)
    set -e
    if [[ -z ${xnco} ]] ; then module load nco ; fi
fi

if [[ ${doextractsite} -eq 1 ]] ; then
    echo "0. Set local mask"
    cd ${pdir}
    # mask
    LandMaskFilePath=$(dirname ${LandMaskFile})
    mkdir -p ${LandMaskFilePath}
    LandMaskFile=$(absfile ${LandMaskFile})
    # echo $(basename ${LandMaskFile})
    # generate random points if ${randompoints} > 0
    if [[ ${randompoints} -ne 0 ]] ; then
        dogeneraterandom=1
        rpoints=${randompoints}
        if [[ ${randompoints} -lt 0 ]] ; then
            if [[ -f ${LandMaskFilePath}/${experiment}_points.csv ]] ; then dogeneraterandom = 0 ; fi
            rpoints=$(( ${randompoints} * -1 ))
        fi
        # generate random points
        cat > ${tmp}/sedtmp.${pid} << EOF
            basepath      = "${sitepath}"
            gridinfo_file = "${GlobalLandMaskFile}"
            outname       = "${LandMaskFilePath}/${experiment}_points.csv"
EOF
        applysed ${tmp}/sedtmp.${pid} ${sdir}/generate_latlonlist.py ${LandMaskFilePath}/generate_latlonlist.py
        python3 ${LandMaskFilePath}/generate_latlonlist.py ${rpoints}

        # set mask to generated random points
        cat > ${tmp}/sedtmp.${pid} << EOF
            path          = "${LandMaskFilePath}"
            maskfname     = "${LandMaskFile}"
            latlonfile    = "${LandMaskFilePath}/${experiment}_points.csv"
            gridinfo_file = "${GlobalLandMaskFile}"
EOF
	      if [[ "${mettype}" == "bios" ]] ; then
            echo "res = 0.25" >> ${tmp}/sedtmp.${pid}
	      fi
        applysed ${tmp}/sedtmp.${pid} ${sdir}/create_landmask.py ${LandMaskFilePath}/create_landmask.py
        sed -i -e "s|from lnutils.*|sys.path.insert(1,'${sdir}'); from lnutils import latlon2ixjy|" ${LandMaskFilePath}/create_landmask.py
        python3 ${LandMaskFilePath}/create_landmask.py
    else  # no random points
        # # cdo -s -f nc4 -z zip sellonlatbox,-72.5,-72.0,42.5,43.0 ${GlobalLandMaskFile} ${LandMaskFile}
        # ncks -O $(nckslatlon ${GlobalLandMaskFile} ${latlon}) ${GlobalLandMaskFile} ${LandMaskFile}
        cat > ${tmp}/sedtmp.${pid} << EOF
            path          = "${LandMaskFilePath}"
            maskfname     = "${LandMaskFile}"
            gridinfo_file = "${GlobalLandMaskFile}"
EOF
	      if [[ "${mettype}" == "bios" ]] ; then
            echo "res = 0.25" >> ${tmp}/sedtmp.${pid}
	      fi
        applysed ${tmp}/sedtmp.${pid} ${sdir}/create_landmask.py ${LandMaskFilePath}/create_landmask.py
        sed -i -e "s|from lnutils.*|sys.path.insert(1,'${sdir}'); from lnutils import latlon2ixjy|" ${LandMaskFilePath}/create_landmask.py
        python3 ${LandMaskFilePath}/create_landmask.py ${latlon}
    fi
elif [[ ${doextractsite} -eq 2 ]] ; then
    echo "0. Extract local meteo and mask"
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


# --------------------------------------------------------------------
# Prepare sequence
#


# Choose meteo, land use and mask directories and files
if [[ ${imeteo} -eq 0 ]] ; then
    MetPath=$(abspath ${GlobalMetPath})
    TransitionFilePath=$(abspath ${GlobalTransitionFilePath})
    LandMaskFile=$(absfile ${GlobalLandMaskFile})
elif [[ ${imeteo} -eq 1 ]] ; then
    MetPath=$(abspath ${GlobalMetPath})
    TransitionFilePath=$(abspath ${GlobalTransitionFilePath})
    LandMaskFile=$(absfile ${LandMaskFile})
elif [[ ${imeteo} -eq 2 ]] ; then
    MetPath=$(abspath ${MetPath})
    TransitionFilePath=$(abspath ${TransitionFilePath})
    LandMaskFile=$(absfile ${LandMaskFile})
else
    printf "Error ${pprog}: imeteo option unknown: ${imeteo}.\n\n"
    exit 1
fi

# absolute paths of other parameter files
ClimateFile=$(absfile ${ClimateFile})
filename_veg=$(absfile ${filename_veg})
filename_soil=$(absfile ${filename_soil})
casafile_cnpbiome=$(absfile ${casafile_cnpbiome})
gm_lut_bernacchi_2002=$(absfile ${gm_lut_bernacchi_2002})
gm_lut_walker_2013=$(absfile ${gm_lut_walker_2013})
filename_d13c_atm=$(absfile ${filename_d13c_atm})
if [[ "${Rubisco_params}" == "Bernacchi_2002" ]] ; then
    filename_gm_lut=${gm_lut_bernacchi_2002}
elif [[ "${Rubisco_params}" == "Walker_2013" ]] ; then
    filename_gm_lut=${gm_lut_walker_2013}
else
    filename_gm_lut=""
fi


# delete all restart files if required
if [[ ${purge_restart} -eq 1 ]] ; then
    rm -f ${rdir}/restart/*
    rm -f ${ClimateFile}
fi


# ------------------------------------------------------------------------------------
# Write global namelists
# ------------------------------------------------------------------------------------
# Write standard namelists with options that are common to all steps of the sequence.
# They can, however, be overwritten in later steps.


# global meteo namelist file
if [[ ${read_fdiff} -eq 1 ]] ; then
    fdiff_bool=.true.
else
    fdiff_bool=.false.
fi
cat > ${tmp}/sedtmp.${pid} << EOF
    BasePath     = "${MetPath}"
    MetPath      = "${MetPath}"
    ReadDiffFrac = ${fdiff_bool}
    LandMaskFile = "${LandMaskFile}"
    Run          = "S0_TRENDY"
EOF
#if [[ ${read_fdiff} -eq 1 ]] ; then
#    sed -i -e "/ReadDiffFrac/s/=.*/= .true./" ${tmp}/sedtmp.${pid}
#fi
applysed ${tmp}/sedtmp.${pid} ${ndir}/cru.nml ${rdir}/cru_${experiment}.nml


# global landuse change namelist
#if [[ "${experiment}" == "S3" || "${experiment}" == "test_S3" ]] ; then
#    YearStart=1580
#else
#    YearStart=1700
#fi  
cat > ${tmp}/sedtmp.${pid} << EOF
    TransitionFilePath = "${TransitionFilePath}"
    ClimateFile        = "${ClimateFile}"
    YearStart          = 1700
    YearEnd            = 2021
EOF
applysed ${tmp}/sedtmp.${pid} ${ndir}/luc.nml ${rdir}/luc_${experiment}.nml



# global Cable namelist
cat > ${tmp}/sedtmp.${pid} << EOF
    filename%met                       = ""
    filename%veg                       = "${filename_veg}"
    filename%soil                      = "${filename_soil}"
    filename%type                      = "${SurfaceFile}"
    filename%out                       = "outputs/${mettype}_out_cable.nc"
    filename%restart_in                = "restart/${mettype}_cable_rst.nc"
    filename%restart_out               = "restart/${mettype}_cable_rst.nc"
    casafile%cnpbiome                  = "${casafile_cnpbiome}"
    casafile%out                       = "outputs/${mettype}_out_casa.nc"
    casafile%cnpipool                  = "restart/${mettype}_casa"
    casafile%cnpepool                  = "restart/${mettype}_casa"
    cable_user%CASA_OUT_FREQ           = "monthly"
    cable_user%POP_restart_in          = "restart/pop_${mettype}_ini.nc"
    cable_user%POP_restart_out         = "restart/pop_${mettype}_ini.nc"
    cable_user%LUC_restart_in          = "restart/${mettype}_LUC_rst.nc"
    cable_user%LUC_restart_out         = "restart/${mettype}_LUC_rst.nc"
    cable_user%LUC_outfile             = "outputs/${mettype}_out_LUC.nc"
    cable_user%climate_restart_in      = "restart/${mettype}_climate_rst.nc"
    cable_user%climate_restart_out     = "restart/${mettype}_climate_rst.nc"
    cable_user%RunIden                 = "${mettype}"
    cable_user%MetType                 = "${mettype}"
    output%averaging                   = "monthly"
    output%grid                        = "mask"
    leaps                              = .false.
    cable_user%SOIL_STRUC              = "default"
    cable_user%Rubisco_parameters      = "${Rubisco_params}"
    cable_user%CALL_POP                = .false.
    cable_user%coordinate_photosyn     = .false.
    cable_user%acclimate_photosyn      = .false.
    cable_user%explicit_gm             = .false.
    cable_user%gm_LUT_file             = "${filename_gm_lut}"
    cable_user%CALL_BLAZE              = .false.
    cable_user%c13o2                   = .false.
    cable_user%c13o2_simple_disc       = .false.
    cable_user%c13o2_delta_atm_file    = "${filename_d13c_atm}"
    cable_user%c13o2_outfile           = "outputs/${mettype}_out_casa_c13o2.nc"
    cable_user%c13o2_restart_in_flux   = "restart/${mettype}_c13o2_flux_rst.nc"
    cable_user%c13o2_restart_out_flux  = "restart/${mettype}_c13o2_flux_rst.nc"
    cable_user%c13o2_restart_in_pools  = "restart/${mettype}_c13o2_pools_rst.nc"
    cable_user%c13o2_restart_out_pools = "restart/${mettype}_c13o2_pools_rst.nc"
    cable_user%c13o2_restart_in_luc    = "restart/${mettype}_c13o2_luc_rst.nc"
    cable_user%c13o2_restart_out_luc   = "restart/${mettype}_c13o2_luc_rst.nc"
EOF
if [[ ${call_pop} -eq 1 ]] ; then
    sed -i -e "/cable_user%CALL_POP/s/=.*/= .true./" ${tmp}/sedtmp.${pid}
fi
if [[ ${coordinate_photosyn} -eq 1 ]] ; then
    sed -i -e "/cable_user%coordinate_photosyn/s/=.*/= .true./" ${tmp}/sedtmp.${pid}
fi
if [[ ${acclimate_photosyn} -eq 1 ]] ; then
    sed -i -e "/cable_user%acclimate_photosyn/s/=.*/= .true./" ${tmp}/sedtmp.${pid}
fi
if [[ ${explicit_gm} -eq 1 ]] ; then
    sed -i -e "/cable_user%explicit_gm/s/=.*/= .true./" ${tmp}/sedtmp.${pid}
fi
if [[ ${doc13o2} -eq 1 ]] ; then
    sed -i -e "/cable_user%c13o2/s/=.*/= .true./" ${tmp}/sedtmp.${pid}
    if [[ ${c13o2_simple_disc} -eq 1 ]] ; then
        sed -i -e "/cable_user%c13o2_simple_disc/s/=.*/= .true./" ${tmp}/sedtmp.${pid}
    fi
fi
if [[ ${call_blaze} -eq 1 ]] ; then
    sed -i -e "cable_user%CALL_BLAZE/s/=.*/= .true./" ${tmp}/sedtmp.${pid}
fi
applysed ${tmp}/sedtmp.${pid} ${ndir}/cable.nml ${rdir}/cable_${experiment}.nml


## BLAZE (no changes at the moment)
if [[ ${call_blaze} -eq 1 ]] ; then
    cp ${ndir}/blaze.nml ${rdir}/blaze.nml
fi


# --------------------------------------------------------------------
# Sequence
#


# --------------------------------------------------------------------
# 1. Create climate restart file
if [[ ${doclimate} -eq 1 ]] ; then
    echo "1. Create climate restart file"
    rid="climate_restart"
    
    # Met forcing
    cp ${rdir}/cru_${experiment}.nml ${rdir}/cru.nml
	
    # LUC
    cp ${rdir}/luc_${experiment}.nml ${rdir}/luc.nml
    # Cable
    #   do not calculate 13C because there is no 13C in the climate restart file
    #MCTEST
    # cable_user%YearEnd = 1889
    # cable_user%CASA_SPIN_ENDYEAR = 1869
    #MCTEST
    cat > ${tmp}/sedtmp.${pid} << EOF
        filename%restart_in            = ""
        cable_user%CLIMATE_fromZero    = .true.
        cable_user%YearStart           = 1860
        cable_user%YearEnd             = 1889
        icycle                         = 2
        spincasa                       = .false.
        cable_user%CASA_fromZero       = .true.
        cable_user%CASA_DUMP_READ      = .false.
        cable_user%CASA_DUMP_WRITE     = .true.
        cable_user%CASA_SPIN_STARTYEAR = 1860
        cable_user%CASA_SPIN_ENDYEAR   = 1869
        cable_user%limit_labile        = .true.
        casafile%cnpipool              = ""
        cable_user%POP_fromZero        = .true.
        cable_user%POP_out             = "ini"
        cable_user%POP_restart_in      = ""
        cable_user%POPLUC              = .false.
        cable_user%POPLUC_RunType      = "static"
        cable_user%c13o2               = .false.
EOF
    applysed ${tmp}/sedtmp.${pid} ${rdir}/cable_${experiment}.nml ${rdir}/cable.nml
    # run model
    cd ${rdir}
    irm logs/log_cable.txt logs/log_out_cable.txt
    if [[ ${dompi} -eq 1 ]] ; then
        ${mpiexecdir}mpiexec -n ${nproc} ./${iexe} > logs/log_out_cable.txt
    else
        ./${iexe} > logs/log_out_cable.txt
    fi
    # save output
    renameid ${rid} ${mettype}.nml luc.nml cable.nml
    imv *_${rid}.nml restart/
    cd logs
    renameid ${rid} log_cable.txt log_out_cable.txt
    cd ../restart
    copyid ${rid} ${mettype}_climate_rst.nc
    cp ${mettype}_climate_rst.nc ${ClimateFile}  # new for TRENDY >= v11   
    cd ../outputs
    renameid ${rid} ${mettype}_out_cable.nc ${mettype}_out_casa.nc
    cd ..
    cd ${pdir}
fi


# --------------------------------------------------------------------
# 2. First spinup phase from zero biomass
if [[ ${dofromzero} -eq 1 ]] ; then
    echo "2. First spinup from zero biomass"
    rid="zero_biomass"

    # Met forcing
    cp ${rdir}/cru_${experiment}.nml ${rdir}/cru.nml

    # LUC
    cp ${rdir}/luc_${experiment}.nml ${rdir}/luc.nml
    # Cable
    #MCTEST
    # cable_user%YearEnd = 1889
    # cable_user%CASA_SPIN_ENDYEAR = 1869
    # remove
    #   cable_user%CASA_OUT_FREQ
    #   output%averaging
    #MCTEST
    cat > ${tmp}/sedtmp.${pid} << EOF
        filename%restart_in               = ""
        cable_user%CLIMATE_fromZero       = .true.
        cable_user%YearStart              = 1860
        cable_user%YearEnd                = 1889
        icycle                            = 2
        spincasa                          = .false.
        cable_user%CASA_OUT_FREQ          = "monthly"
        cable_user%CASA_fromZero          = .true.
        cable_user%CASA_DUMP_READ         = .false.
        cable_user%CASA_DUMP_WRITE        = .true.
        cable_user%CASA_SPIN_STARTYEAR    = 1850
        cable_user%CASA_SPIN_ENDYEAR      = 1859
        cable_user%limit_labile           = .true.
        casafile%cnpipool                 = ""
        cable_user%POP_fromZero           = .true.
        cable_user%POP_out                = "ini"
        cable_user%POP_restart_in         = ""
        cable_user%POPLUC                 = .true.
        cable_user%POPLUC_RunType         = "static"
        cable_user%c13o2_restart_in_flux  = ""
        cable_user%c13o2_restart_in_pools = ""
EOF
    applysed ${tmp}/sedtmp.${pid} ${rdir}/cable_${experiment}.nml ${rdir}/cable.nml

    # run model
    cd ${rdir}
    irm logs/log_cable.txt logs/log_out_cable.txt
    if [[ ${dompi} -eq 1 ]] ; then
        ${mpiexecdir}mpiexec -n ${nproc} ./${iexe} > logs/log_out_cable.txt
    else
        ./${iexe} > logs/log_out_cable.txt
    fi
    # save output
    renameid ${rid} ${mettype}.nml luc.nml cable.nml
    imv *_${rid}.nml restart/
    cd logs
    renameid ${rid} log_cable.txt log_out_cable.txt
    cd ../restart
    copyid ${rid} ${mettype}_climate_rst.nc ${mettype}_cable_rst.nc pop_${mettype}_ini.nc
    copyid ${rid} ${mettype}_casa_biome_rst.nc ${mettype}_casa_met_rst.nc ${mettype}_casa_pool_rst.nc ${mettype}_casa_phen_rst.nc
    copyid ${rid} ${mettype}_casa_flux_rst.nc ${mettype}_casa_bal_rst.nc
    if [[ ${doc13o2} -eq 1 ]] ; then copyid ${rid} ${mettype}_c13o2_flux_rst.nc ${mettype}_c13o2_pools_rst.nc ; fi
    cd ../outputs
    renameid ${rid} ${mettype}_out_cable.nc ${mettype}_out_casa.nc
    if [[ ${doc13o2} -eq 1 ]] ; then renameid ${rid} ${mettype}_out_casa_c13o2.nc ; fi
    cd ..
    cd ${pdir}
fi


# --------------------------------------------------------------------
# 3. Biomass into quasi-equilibrium with restricted N and P pools
if [[ ${doequi1} -eq 1 ]] ; then
    echo "3. Bring biomass into quasi-equilibrium with unrestricted N and P pools"
    if [[ "${run_step}" == "spinup2" ]] ; then
	start_iequi1=3
    else
	start_iequi1=1
    fi
    
    for ((iequi1=${start_iequi1}; iequi1<=${nequi1}; iequi1++)) ; do
        # 3a. 30 year run starting from restart files
        echo "   3a. 30 year spinup from accumulated biomass; iequi1=${iequi1}/${nequi1}"
        rid="spinup_limit_labile_${iequi1}"
        # rid="spinup_limit_labile${iequi}"

	# Met forcing
        cp ${rdir}/cru_${experiment}.nml ${rdir}/cru.nml

        # LUC
        cp ${rdir}/luc_${experiment}.nml ${rdir}/luc.nml
        # Cable
        #MCTEST
        # cable_user%YearEnd = 1859
        # cable_user%CASA_SPIN_ENDYEAR = 1869
        #MCTEST
        cat > ${tmp}/sedtmp.${pid} << EOF
            cable_user%CLIMATE_fromZero    = .false.
            cable_user%YearStart           = 1840
            cable_user%YearEnd             = 1859
            icycle                         = 2
            spincasa                       = .false.
            cable_user%CASA_fromZero       = .false.
            cable_user%CASA_DUMP_READ      = .false.
            cable_user%CASA_DUMP_WRITE     = .true.
            cable_user%CASA_SPIN_STARTYEAR = 1850
            cable_user%CASA_SPIN_ENDYEAR   = 1859
            cable_user%limit_labile        = .true.
            cable_user%POP_fromZero        = .false.
            cable_user%POP_out             = "ini"
            cable_user%POPLUC              = .true.
            cable_user%POPLUC_RunType      = "static"
EOF
        applysed ${tmp}/sedtmp.${pid} ${rdir}/cable_${experiment}.nml ${rdir}/cable.nml
        # run model
        cd ${rdir}
        irm logs/log_cable.txt logs/log_out_cable.txt
        if [[ ${dompi} -eq 1 ]] ; then
            ${mpiexecdir}mpiexec -n ${nproc} ./${iexe} > logs/log_out_cable.txt
        else
            ./${iexe} > logs/log_out_cable.txt
        fi
        # save output
        renameid ${rid} ${mettype}.nml luc.nml cable.nml
        mv *_${rid}.nml restart/
        cd logs
        renameid ${rid} log_cable.txt log_out_cable.txt
        cd ../restart
        copyid ${rid} ${mettype}_climate_rst.nc ${mettype}_cable_rst.nc pop_${mettype}_ini.nc
	copyid ${rid} ${mettype}_casa_biome_rst.nc ${mettype}_casa_met_rst.nc ${mettype}_casa_pool_rst.nc ${mettype}_casa_phen_rst.nc
        copyid ${rid} ${mettype}_casa_flux_rst.nc ${mettype}_casa_bal_rst.nc
        if [[ ${doc13o2} -eq 1 ]] ; then copyid ${rid} ${mettype}_c13o2_flux_rst.nc ${mettype}_c13o2_pools_rst.nc ; fi
        cd ../outputs
        renameid ${rid} ${mettype}_out_cable.nc ${mettype}_out_casa.nc
	if [[ ${doc13o2} -eq 1 ]] ; then renameid ${rid} ${mettype}_out_casa_c13o2.nc ; fi
        cd ..
        cd ${pdir}
	
        #
        # 3b. analytic quasi-equilibrium of biomass pools
	echo "   3b. Analytic solution of biomass pools"
        rid="spinup_analytic_limit_labile"
        # rid="spin_casa_limit_labile${iequi}"

        # Met forcing
        cp ${rdir}/cru_${experiment}.nml ${rdir}/cru.nml

        # LUC
        cp ${rdir}/luc_${experiment}.nml ${rdir}/luc.nml

        # Cable
        cat > ${tmp}/sedtmp.${pid} << EOF
            cable_user%CLIMATE_fromZero    = .false.
            cable_user%YearStart           = 1840
            cable_user%YearEnd             = 1859
            icycle                         = 12
            spincasa                       = .true.
            cable_user%CASA_fromZero       = .false.
            cable_user%CASA_DUMP_READ      = .true.
            cable_user%CASA_DUMP_WRITE     = .false.
            cable_user%CASA_SPIN_STARTYEAR = 1840
            cable_user%CASA_SPIN_ENDYEAR   = 1859
            cable_user%limit_labile        = .true.
            cable_user%POP_fromZero        = .false.
            cable_user%POP_out             = "ini"
            cable_user%POPLUC              = .true.
            cable_user%POPLUC_RunType      = "static"
EOF
        applysed ${tmp}/sedtmp.${pid} ${rdir}/cable_${experiment}.nml ${rdir}/cable.nml
        # run model
        cd ${rdir}
        irm logs/log_cable.txt logs/log_out_cable.txt
        if [[ ${dompi} -eq 1 ]] ; then
             ${mpiexecdir}mpiexec -n ${nproc} ./${iexe} > logs/log_out_cable.txt
        else
            ./${iexe} > logs/log_out_cable.txt
        fi
        # save output
        renameid ${rid} ${mettype}.nml luc.nml cable.nml
        mv *_${rid}.nml restart/
        cd logs
        renameid ${rid} log_cable.txt log_out_cable.txt
        cd ../restart
        copyid ${rid} pop_${mettype}_ini.nc
	copyid ${rid} ${mettype}_casa_biome_rst.nc ${mettype}_casa_met_rst.nc ${mettype}_casa_pool_rst.nc ${mettype}_casa_phen_rst.nc
        copyid ${rid} ${mettype}_casa_flux_rst.nc ${mettype}_casa_bal_rst.nc
        if [[ ${doc13o2} -eq 1 ]] ; then copyid ${rid} ${mettype}_c13o2_flux_rst.nc ${mettype}_c13o2_pools_rst.nc ; fi
        if [[ ${dompi} -eq 0 ]] ; then # no output only restart if MPI
            cd ../outputs
            renameid ${rid} ${mettype}_out_casa.nc
            if [[ ${doc13o2} -eq 1 ]] ; then renameid ${rid} ${mettype}_out_casa_c13o2.nc ; fi
            cd ..
        fi
        cd ${pdir}
    done
fi



# --------------------------------------------------------------------
# 4. Biomass into quasi-equilibrium without restricted N and P pools
if [[ ${doequi2} -eq 1 ]] ; then
    echo "4. Bring biomass into quasi-equilibrium with restricted N and P pools"
    for ((iequi2=1; iequi2<=${nequi2}; iequi2++)) ; do
        # 4a. 30 year run starting from restart files
        echo "   4a. 30 year spinup from accumulated biomass; iequi2=${iequi2}/${nequi2}"
        #rid="spinup_nutrient_limited"
        rid="spinup_nutrient_limited_${iequi2}"

	# Met forcing
        cp ${rdir}/cru_${experiment}.nml ${rdir}/cru.nml

	# LUC
	cp ${rdir}/luc_${experiment}.nml ${rdir}/luc.nml

	# Cable
        cat > ${tmp}/sedtmp.${pid} << EOF
            cable_user%CLIMATE_fromZero    = .false.
            cable_user%YearStart           = 1840
            cable_user%YearEnd             = 1859
            icycle                         = 2
            spincasa                       = .false.
            cable_user%CASA_fromZero       = .false.
            cable_user%CASA_DUMP_READ      = .false.
            cable_user%CASA_DUMP_WRITE     = .true.
            cable_user%CASA_SPIN_STARTYEAR = 1850
            cable_user%CASA_SPIN_ENDYEAR   = 1859
            cable_user%limit_labile        = .false.
            cable_user%POP_fromZero        = .false.
            cable_user%POP_out             = "ini"
            cable_user%POPLUC              = .true.
            cable_user%POPLUC_RunType      = "static"
EOF
        applysed ${tmp}/sedtmp.${pid} ${rdir}/cable_${experiment}.nml ${rdir}/cable.nml
        # run model
        cd ${rdir}
        irm logs/log_cable.txt logs/log_out_cable.txt
        if [[ ${dompi} -eq 1 ]] ; then
            ${mpiexecdir}mpiexec -n ${nproc} ./${iexe} > logs/log_out_cable.txt
        else
            ./${iexe} > logs/log_out_cable.txt
        fi
        # save output
        renameid ${rid} ${mettype}.nml luc.nml cable.nml
        mv *_${rid}.nml restart/
        cd logs
        renameid ${rid} log_cable.txt log_out_cable.txt
        cd ../restart
        copyid ${rid} ${mettype}_climate_rst.nc ${mettype}_cable_rst.nc pop_${mettype}_ini.nc
	copyid ${rid} ${mettype}_casa_biome_rst.nc ${mettype}_casa_met_rst.nc ${mettype}_casa_pool_rst.nc ${mettype}_casa_phen_rst.nc
        copyid ${rid} ${mettype}_casa_flux_rst.nc ${mettype}_casa_bal_rst.nc
        if [[ ${doc13o2} -eq 1 ]] ; then copyid ${rid} ${mettype}_c13o2_flux_rst.nc ${mettype}_c13o2_pools_rst.nc ; fi
        cd ../outputs
        renameid ${rid} ${mettype}_out_cable.nc ${mettype}_out_casa.nc
	if [[ ${doc13o2} -eq 1 ]] ; then renameid ${rid} ${mettype}_out_casa_c13o2.nc ; fi
        cd ..
        cd ${pdir}
	
        #
        # 4b. analytic quasi-equilibrium of biomass pools
        echo "   4b. Analytic solution of biomass pools"
        #rid="spinup_analytic"
        rid="spinup_analytic_${iequi2}"
        # Met forcing
        cp ${rdir}/cru_${experiment}.nml ${rdir}/cru.nml

        # LUC
        cp ${rdir}/luc_${experiment}.nml ${rdir}/luc.nml
        #applysed ${tmp}/sedtmp.${pid} ${rdir}/LUC_${experiment}.nml ${rdir}/LUC.nml


        # Cable
        #MCTEST
        # cable_user%YearEnd = 1859
        # cable_user%CASA_SPIN_ENDYEAR = 1859
        #MCTEST
        cat > ${tmp}/sedtmp.${pid} << EOF
            cable_user%CLIMATE_fromZero    = .false.
            cable_user%YearStart           = 1840
            cable_user%YearEnd             = 1859
            icycle                         = 12
            spincasa                       = .true.
            cable_user%CASA_fromZero       = .false.
            cable_user%CASA_DUMP_READ      = .true.
            cable_user%CASA_DUMP_WRITE     = .false.
            cable_user%CASA_SPIN_STARTYEAR = 1840
            cable_user%CASA_SPIN_ENDYEAR   = 1859
            cable_user%limit_labile        = .false.
            cable_user%POP_fromZero        = .false.
            cable_user%POP_out             = "ini"
            cable_user%POPLUC              = .true.
            cable_user%POPLUC_RunType      = "static"
EOF
        applysed ${tmp}/sedtmp.${pid} ${rdir}/cable_${experiment}.nml ${rdir}/cable.nml
        # run model
        cd ${rdir}
        irm logs/log_cable.txt logs/log_out_cable.txt
        if [[ ${dompi} -eq 1 ]] ; then
            ${mpiexecdir}mpiexec -n ${nproc} ./${iexe} > logs/log_out_cable.txt
        else
            ./${iexe} > logs/log_out_cable.txt
        fi
        # save output
        renameid ${rid} ${mettype}.nml luc.nml cable.nml
        mv *_${rid}.nml restart/
        cd logs
        renameid ${rid} log_cable.txt log_out_cable.txt
        cd ../restart
        copyid ${rid} ${mettype}_cable_rst.nc pop_${mettype}_ini.nc
	copyid ${rid} ${mettype}_casa_biome_rst.nc ${mettype}_casa_met_rst.nc ${mettype}_casa_pool_rst.nc ${mettype}_casa_phen_rst.nc
        copyid ${rid} ${mettype}_casa_flux_rst.nc ${mettype}_casa_bal_rst.nc
        if [[ ${doc13o2} -eq 1 ]] ; then copyid ${rid} ${mettype}_c13o2_flux_rst.nc ${mettype}_c13o2_pools_rst.nc ; fi
        if [[ ${dompi} -eq 0 ]] ; then # no output only restart if MPI
            cd ../outputs
            renameid ${rid} ${mettype}_out_casa.nc
	    if [[ ${doc13o2} -eq 1 ]] ; then renameid ${rid} ${mettype}_out_casa_c13o2.nc ; fi
            cd ..
        fi
        cd ${pdir}
    done
fi



# -------------------------------------------------------------------------------------------
### End of spinup here. Copy spinup restart files to the actual folders (global runs only)
# -------------------------------------------------------------------------------------------
if [[ ${run_type} -eq 0 && "${run_step}" == "1700-1800" ]] ; then
    
    echo "copying spinup restart files to experiment folder"

    if [[ ! -f $(dirname ${rdir})/spinup/restart/cru_cable_rst.nc ]] ; then
	echo "No existing restart files! Run run_step 'spinup1' first!"
	exit 1
    fi 
    
    cp $(dirname ${rdir})/spinup/restart/cru_climate_rst.nc     ${rdir}/restart/cru_climate_rst.nc 
    cp $(dirname ${rdir})/spinup/restart/cru_cable_rst.nc       ${rdir}/restart/cru_cable_rst.nc
    cp $(dirname ${rdir})/spinup/restart/pop_cru_ini.nc         ${rdir}/restart/pop_cru_ini.nc
    cp $(dirname ${rdir})/spinup/restart/cru_casa_biome_rst.nc  ${rdir}/restart/cru_casa_biome_rst.nc
    cp $(dirname ${rdir})/spinup/restart/cru_casa_met_rst.nc    ${rdir}/restart/cru_casa_met_rst.nc
    cp $(dirname ${rdir})/spinup/restart/cru_casa_pool_rst.nc   ${rdir}/restart/cru_casa_pool_rst.nc
    cp $(dirname ${rdir})/spinup/restart/cru_casa_phen_rst.nc   ${rdir}/restart/cru_casa_phen_rst.nc
    cp $(dirname ${rdir})/spinup/restart/cru_casa_flux_rst.nc   ${rdir}/restart/cru_casa_flux_rst.nc
    cp $(dirname ${rdir})/spinup/restart/cru_casa_bal_rst.nc    ${rdir}/restart/cru_casa_bal_rst.nc 

    ## for runs S3 upwards also copy dump files:
    if [[ "${experiment}" == "S3" ]] ; then
        cp $(dirname ${rdir})/spinup/*dump.nc ${rdir}/.
    fi
fi



# --------------------------------------------------------------------
# 5. First dynamic land use (Initialise land use)
if [[ ${doiniluc} -eq 1 ]] ; then
    echo "5. First dynamic land use (initialise land use)"

    # Met forcing
    YearStart=1580  # should be the same as in the global luc.nml file!!
    YearEnd=1699
    cp ${rdir}/cru_${experiment}.nml ${rdir}/cru.nml
    rid="init_land_use"
    
    # LUC
    cat > ${tmp}/sedtmp.${pid} << EOF
         YearStart = ${YearStart}
         YearEnd   = ${YearEnd}
EOF
    applysed ${tmp}/sedtmp.${pid} ${rdir}/luc_${experiment}.nml ${rdir}/luc.nml

    # Cable
    cat > ${tmp}/sedtmp.${pid} << EOF
        cable_user%CLIMATE_fromZero     = .false.
        cable_user%YearStart            = ${YearStart}
        cable_user%YearEnd              = ${YearEnd}
        icycle                          = 12
        spincasa                        = .false.
        cable_user%CASA_OUT_FREQ        = "annually"
        cable_user%CASA_fromZero        = .false.
        cable_user%CASA_DUMP_READ       = .true.
        cable_user%CASA_DUMP_WRITE      = .false.
        cable_user%CASA_SPIN_STARTYEAR  = 1840
        cable_user%CASA_SPIN_ENDYEAR    = 1859
        cable_user%limit_labile         = .false.
        cable_user%POP_fromZero         = .false.
        cable_user%POP_out              = "ini"
        cable_user%POPLUC               = .true.
        cable_user%POPLUC_RunType       = "init"
        cable_user%LUC_restart_in       = ""
        cable_user%c13o2_restart_in_luc = ""
EOF
    applysed ${tmp}/sedtmp.${pid} ${rdir}/cable_${experiment}.nml ${rdir}/cable.nml
    # run model
    cd ${rdir}
    irm logs/log_cable.txt logs/log_out_cable.txt
    if [[ ${dompi} -eq 1 ]] ; then
        ${mpiexecdir}mpiexec -n ${nproc} ./${iexe} > logs/log_out_cable.txt
    else
        ./${iexe} > logs/log_out_cable.txt
    fi
    # save output
    renameid ${rid} ${mettype}.nml luc.nml cable.nml
    mv *_${rid}.nml restart/
    cd logs
    renameid ${rid} log_cable.txt log_out_cable.txt
    cd ../restart
    #copyid ${rid} ${mettype}_cable_rst.nc ${mettype}_casa_rst.nc ${mettype}_LUC_rst.nc pop_${mettype}_ini.nc
    copyid ${rid} ${mettype}_LUC_rst.nc pop_${mettype}_ini.nc
    copyid ${rid} ${mettype}_casa_biome_rst.nc ${mettype}_casa_met_rst.nc ${mettype}_casa_pool_rst.nc ${mettype}_casa_phen_rst.nc
    copyid ${rid} ${mettype}_casa_flux_rst.nc ${mettype}_casa_bal_rst.nc
    if [[ ${doc13o2} -eq 1 ]] ; then copyid ${rid} ${mettype}_c13o2_pools_rst.nc ${mettype}_c13o2_luc_rst.nc ; fi
    cd ../outputs
    #renameid ${rid} ${mettype}_out_LUC.nc ${mettype}_out_casa.nc ${mettype}_out_cable.nc
    renameid ${rid} ${mettype}_out_LUC.nc ${mettype}_out_casa.nc
    if [[ ${doc13o2} -eq 1 ]] ; then renameid ${rid} ${mettype}_out_casa_c13o2.nc ; fi
    cd ..
    cd ${pdir}
fi



# --------------------------------------------------------------------
# 6. Transient run
if [[ ${doinidyn} -eq 1 ]] ; then
    echo "6. Transient run (full dynamic spinup)"

    # Met forcing
    if [[ ${run_type} -eq 0 ]] ; then   # global run
       if [[ "${run_step}" == "1700-1800" ]] ; then
          YearStart=1700
          YearEnd=1800
       elif [[ "${run_step}" == "1801-1900" ]] ; then
          YearStart=1801
          YearEnd=1900
          if [[ ! -f ${rdir}/restart/cru_cable_rst.nc ]] ; then
	     echo "No exisiting restart files! Run run_step '1700-1800' first!"
	     exit 1
	  fi
       else
          echo "invalid option for 'run_step'!! Exiting..."
          exit 1
       fi
    else
	YearStart=1700
	YearEnd=1900
    fi
	  
	  

       
    if [[ "${experiment}" == "S0" ]] ; then
        cat > ${tmp}/sedtmp.${pid} << EOF
            Run = "S0_TRENDY"
EOF
    elif [[ "${experiment}" == "S1" || "${experiment}" == "S2" || "${experiment}" == "S3" ]] ; then
        cat > ${tmp}/sedtmp.${pid} << EOF
            Run = "S1_TRENDY"
EOF
    fi	
    applysed ${tmp}/sedtmp.${pid} ${rdir}/cru_${experiment}.nml ${rdir}/cru.nml
    rid=${YearStart}_${YearEnd}

    
    # LUC
    if [[ "${experiment}" == "S3" ]] ; then
       cat > ${tmp}/sedtmp.${pid} << EOF
           YearStart = ${YearStart}
           YearEnd   = ${YearEnd}
EOF
    else
       cat > ${tmp}/sedtmp.${pid} << EOF
           YearStart = 1700
           YearEnd   = ${YearEnd}
EOF
    fi
    applysed ${tmp}/sedtmp.${pid} ${rdir}/luc_${experiment}.nml ${rdir}/luc.nml

    
    # Cable
    if [[ "${experiment}" == "S3" ]] ; then
	POPLUC_RunType="restart"
    else
        POPLUC_RunType="static"
    fi
    
    cat > ${tmp}/sedtmp.${pid} << EOF
        cable_user%CLIMATE_fromZero    = .false.
        cable_user%YearStart           = ${YearStart}
        cable_user%YearEnd             = ${YearEnd}
        icycle                         = 2
        spincasa                       = .false.
        cable_user%CASA_fromZero       = .false.
        cable_user%CASA_DUMP_READ      = .false.
        cable_user%CASA_DUMP_WRITE     = .false.
        cable_user%CASA_SPIN_STARTYEAR = 1850
        cable_user%CASA_SPIN_ENDYEAR   = 1859
        cable_user%limit_labile        = .false.
        cable_user%POP_fromZero        = .false.
        cable_user%POP_out             = "ini"
        cable_user%POPLUC              = .true.
        cable_user%POPLUC_RunType      = "${POPLUC_RunType}"
EOF
    applysed ${tmp}/sedtmp.${pid} ${rdir}/cable_${experiment}.nml ${rdir}/cable.nml
    # run model
    cd ${rdir}
    irm logs/log_cable.txt logs/log_out_cable.txt
    if [[ ${dompi} -eq 1 ]] ; then
        ${mpiexecdir}mpiexec -n ${nproc} ./${iexe} > logs/log_out_cable.txt
    else
        ./${iexe} > logs/log_out_cable.txt
    fi
    # save output
    renameid ${rid} ${mettype}.nml luc.nml cable.nml
    mv *_${rid}.nml restart/
    cd logs
    renameid ${rid} log_cable.txt log_out_cable.txt
    cd ../restart
    copyid ${rid} ${mettype}_climate_rst.nc ${mettype}_cable_rst.nc
    copyid ${rid} ${mettype}_LUC_rst.nc pop_${mettype}_ini.nc
    copyid ${rid} ${mettype}_casa_biome_rst.nc ${mettype}_casa_met_rst.nc ${mettype}_casa_pool_rst.nc ${mettype}_casa_phen_rst.nc
    copyid ${rid} ${mettype}_casa_flux_rst.nc ${mettype}_casa_bal_rst.nc
    if [[ ${doc13o2} -eq 1 ]] ; then copyid ${rid} ${mettype}_c13o2_flux_rst.nc ${mettype}_c13o2_pools_rst.nc ${mettype}_c13o2_luc_rst.nc ; fi
    cd ../outputs
    renameid ${rid} ${mettype}_out_cable.nc ${mettype}_out_casa.nc ${mettype}_out_LUC.nc
    if [[ ${doc13o2} -eq 1 ]] ; then renameid ${rid} ${mettype}_out_casa_c13o2.nc ; fi
    cd ..
    cd ${pdir}
fi


# --------------------------------------------------------------------
# 7. Final centennial run
if [[ ${dofinal} -eq 1 ]] ; then
    echo "7. Final centennial run"

    # Met forcing
    YearStart=1901
    YearEnd=2021
    if [[ "${experiment}" == "S0" ]] ; then
        cat > ${tmp}/sedtmp.${pid} << EOF
            Run = "S0_TRENDY"
EOF
    elif [[ "${experiment}" == "S1" ]] ; then
        cat > ${tmp}/sedtmp.${pid} << EOF
            Run = "S1_TRENDY"
EOF
    elif [[ "${experiment}" == "S2" || "${experiment}" == "S3" ]] ; then
        cat > ${tmp}/sedtmp.${pid} << EOF
            Run = "S2_TRENDY"
EOF
    fi
    applysed ${tmp}/sedtmp.${pid} ${rdir}/cru_${experiment}.nml ${rdir}/cru.nml
    rid=${YearStart}_${YearEnd}

    
    # LUC
    if [[ "${experiment}" == "S3" ]] ; then
       cat > ${tmp}/sedtmp.${pid} << EOF
           YearStart = ${YearStart}
           YearEnd   = ${YearEnd}
EOF
    else
       cat > ${tmp}/sedtmp.${pid} << EOF
           YearStart = 1700
           YearEnd   = ${YearEnd}
EOF
    fi
    applysed ${tmp}/sedtmp.${pid} ${rdir}/luc_${experiment}.nml ${rdir}/luc.nml

    
    # Cable
    if [[ "${experiment}" == "S3" ]] ; then
	POPLUC_RunType="restart"
    else
        POPLUC_RunType="static"
    fi
    
    cat > ${tmp}/sedtmp.${pid} << EOF
        cable_user%CLIMATE_fromZero    = .false.
        cable_user%YearStart           = ${YearStart}
        cable_user%YearEnd             = ${YearEnd}
        icycle                         = 2
        spincasa                       = .false.
        cable_user%CASA_fromZero       = .false.
        cable_user%CASA_DUMP_READ      = .false.
        cable_user%CASA_DUMP_WRITE     = .false.
        cable_user%CASA_SPIN_STARTYEAR = 1850
        cable_user%CASA_SPIN_ENDYEAR   = 1859
        cable_user%limit_labile        = .false.
        cable_user%POP_fromZero        = .false.
        cable_user%POP_out             = "ini"
        cable_user%POPLUC              = .true.
        cable_user%POPLUC_RunType      = "${POPLUC_RunType}"
EOF
    applysed ${tmp}/sedtmp.${pid} ${rdir}/cable_${experiment}.nml ${rdir}/cable.nml
    # run model
    cd ${rdir}
    irm logs/log_cable.txt logs/log_out_cable.txt
    if [[ ${dompi} -eq 1 ]] ; then
        ${mpiexecdir}mpiexec -n ${nproc} ./${iexe} > logs/log_out_cable.txt
    else
        ./${iexe} > logs/log_out_cable.txt
    fi
    # save output
    renameid ${rid} ${mettype}.nml luc.nml cable.nml
    mv *_${rid}.nml restart/
    cd logs
    renameid ${rid} log_cable.txt log_out_cable.txt
    cd ../restart
    copyid ${rid} ${mettype}_climate_rst.nc ${mettype}_cable_rst.nc
    copyid ${rid} ${mettype}_LUC_rst.nc pop_${mettype}_ini.nc
    copyid ${rid} ${mettype}_casa_biome_rst.nc ${mettype}_casa_met_rst.nc ${mettype}_casa_pool_rst.nc ${mettype}_casa_phen_rst.nc
    copyid ${rid} ${mettype}_casa_flux_rst.nc ${mettype}_casa_bal_rst.nc
    if [[ ${doc13o2} -eq 1 ]] ; then copyid ${rid} ${mettype}_c13o2_flux_rst.nc ${mettype}_c13o2_pools_rst.nc ${mettype}_c13o2_luc_rst.nc ; fi
    cd ../outputs
    renameid ${rid} ${mettype}_out_cable.nc ${mettype}_out_casa.nc ${mettype}_out_LUC.nc
    if [[ ${doc13o2} -eq 1 ]] ; then renameid ${rid} ${mettype}_out_casa_c13o2.nc ; fi
    cd ..
    cd ${pdir}
fi



## Remove spinup files where they aren't needed (e.g. for 1000pt runs)
if [[ ${run_type} -eq 1 && ${purge_spinup} -eq 1 ]] ; then

    cd ${rdir}/outputs
    rm *_spinup_*
    rm *_zero_biomass*
    rm *_1700_1900.nc
    
fi



# --------------------------------------------------------------------
# Finish
# --------------------------------------------------------------------

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
