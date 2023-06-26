#!/usr/bin/env bash

# Gadi
# https://opus.nci.org.au/display/Help/How+to+submit+a+job
#PBS -N S3_test8_serial_global
#PBS -P x45
#PBS -q normal
#PBS -p 400
#PBS -l walltime=12:00:00
#PBS -l mem=4GB
#PBS -l ncpus=1
#PBS -l storage=gdata/x45+scratch/pt17+gdata/vl59
#PBS -l software=netCDF:MPI:Intel:GNU:scorep
#PBS -r y
#PBS -l wd
#PBS -j oe
#PBS -S /bin/bash

#export SCOREP_ENABLE_PROFILING=1
#export SCOREP_MEMORY_RECORDING=1
#export SCOREP_ENABLE_TRACING=1
#export SCOREP_TOTAL_MEMORY=4GB
#export SCOREP_EXPERIMENT_DIRECTORY=/scratch/pt17/jk8585/CABLE_tests/scorep_profiler_memory
#export SCOREP_FILTERING_FILE=/scratch/pt17/jk8585/profile.filter

# --------------------------------------------------------------------
#
# Full Cable run with biomass spinup, POP, land-use change, etc.
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
#           Jurgen Knauer, Mar 2023  - modified to run CABLE in serial with smaller land masks
#
# --------------------------------------------------------------------


# --------------------------------------------------------------------
# Load Modules
# --------------------------------------------------------------------
pdir=${isdir}
. /etc/bashrc
module purge
module load intel-compiler/2021.5.0
module load intel-mpi/2021.5.1
module load netcdf/4.8.0
export mpiexecdir=/apps/intel-mpi/2021.5.1/bin
if [[ ! -z ${mpiexecdir} ]] ; then export mpiexecdir="${mpiexecdir}/" ; fi



## ------------------------------------------------------------------
## Basic settings (parsed through from wrapper script)
## ------------------------------------------------------------------
# TRENDY experiment (S0, S1, S2, S3, S4, S5, S6):     
experiment='S3'
# Name of the experiment (= name of output folder)     
experiment_name='S3_test8_serial_global'
# Code directory
cablecode='/home/599/jk8585/CABLE_code/SHARE/CABLE-POP_TRENDY'
# Script directory
cablehome='/home/599/jk8585/CABLE_run/TRENDY_v11'
# Cable executable
exe='/home/599/jk8585/CABLE_code/SHARE/CABLE-POP_TRENDY/offline/cable'
# CABLE-AUX directory (uses offline/gridinfo_CSIRO_1x1.nc and offline/modis_phenology_csiro.txt)
aux='/g/data/vl59/TRENDY_v11/aux'
# Global Meteorology
MetPath='/g/data/x45/CRUJRA2022/daily_1deg'
# Global LUC
TransitionFilePath='/g/data/x45/LUH2/GCB_2022/1deg/EXTRACT'
# Global Surface file 
SurfaceFile='/g/data/vl59/TRENDY_v11/aux/gridinfo_CSIRO_1x1.nc'
# Output directory of the run
runpath='/scratch/pt17/jk8585/CABLE_tests/S3_test8_serial_global/run100'
# Land Mask used for this run
LandMaskFile='/scratch/pt17/jk8585/CABLE_tests/S3_test8_serial_global/run100/landmask/landmask100.nc'


## ----------------------------------------------------------------
## Run Sequence
## ----------------------------------------------------------------
doclimate=1     # 1/0: Do/Do not create climate restart file
dofromzero=1    # 1/0  Do/Do not first spinup phase from zero biomass stocks
doequi1=1       # 1/0: Do/Do not bring biomass stocks into quasi-equilibrium with restricted P and N pools
nequi1=2        #      number of times to repeat steps in doequi1
doequi2=1       # 1/0: Do/Do not bring biomass stocks into quasi-equilibrium with unrestricted P and N pools
nequi2=2        #      number of times to repeat steps in doequi2
if [[ "${experiment}" == "S3" ]] ; then
    doiniluc=1      # 1/0: Do/Do not spinup with dynamic land use (initialise land use)
else
    doiniluc=0
fi
doinidyn=1      # 1/0: Do/Do not full dynamic spinup (transient run) from 1701 to 1900
dofinal=1       # 1/0: Do/Do not final run from 1901 to 2021
    
purge_restart=0  # Delete all restart files?

## ----------------------------------------------------------------
## CABLE Settings
## ----------------------------------------------------------------
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
# Parameter files
namelistpath="${cablehome}/namelists"
filename_veg="${cablehome}/params/def_veg_params.txt"
filename_soil="${cablehome}/params/def_soil_params.txt"
casafile_cnpbiome="${cablehome}/params/pftlookup.csv"
# Climate restart file 
# changes for TRENDY >= v11: ClimateFile always created!
# ClimateFile="/g/data/x45/ipbes/cable_climate/ipsl_climate_rst_glob_1deg.nc"
#ClimateFile="$(dirname ${runpath})/climate_restart/cru_climate_rst.nc"
ClimateFile="${runpath}/cru_climate_rst.nc"
# gm lookup tables
gm_lut_bernacchi_2002="${cablehome}/params/gm_LUT_351x3601x7_1pt8245_Bernacchi2002.nc"
gm_lut_walker_2013="${cablehome}/params/gm_LUT_351x3601x7_1pt8245_Walker2013.nc"
# 13C
filename_d13c_atm="${cablehome}/params/graven_et_al_gmd_2017-table_s1-delta_13c-1700-2025.txt"


# --------------------------------------------------------------------
# Setup and Functions
# --------------------------------------------------------------------

set -e

trap cleanup 1 2 3 6

pid=$$
isdir="${PWD}"
prog=$0
pprog=$(basename ${prog})
pdir=$(dirname ${prog})
tmp=${TMPDIR:-"/tmp"}

# Helper functions, most functions are in plumber_cable-pop_lib.sh
source ${isdir}/run_cable-pop_lib.sh

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


# ---------------------------------------------------------------------------------
# Preparation
# ---------------------------------------------------------------------------------
# Get options
while getopts "h" option ; do
    case ${option} in
        h) usage; exit;;
        *) printf "Error ${pprog}: unimplemented option.\n\n" 1>&2;  usage 1>&2; exit 1;;
    esac
done
shift $((${OPTIND} - 1))

# get directories
pdir=$(abspath ${pdir})
cd ${pdir}
adir=$(abspath ${aux})
exe='/home/599/jk8585/CABLE_code/SHARE/CABLE-POP_TRENDY/offline/cable'
mkdir -p ${runpath}
rdir=$(abspath ${runpath})
ndir=$(abspath ${namelistpath})

# prepare run directory
cd ${rdir}
mkdir -p logs
mkdir -p outputs
mkdir -p restart
ln -sf ${adir}
cp ${exe} ./
iexe=$(basename ${exe})
cd ${pdir}

# set stacksize to unlimited if permitted, otherwise to 15 bit if possible
set +e
#ulimit -s unlimited 2> /dev/null || ulimit -s 32768
set -e

# --------------------------------------------------------------------
# Print Info
# --------------------------------------------------------------------
t1=$(date +%s)
printf "Started at %s\n" "$(date)"

printf "\nSetup\n"
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
printf "        LandMaskFile=${LandMaskFile}\n"
printf "        SurfaceFile=${SurfaceFile}\n"
printf "        runpath=${runpath}\n"
printf "        namelistpath=${namelistpath}\n"
printf "        filename_veg=${filename_veg}\n"
printf "        filename_soil=${filename_soil}\n"
printf "        casafile_cnpbiome=${casafile_cnpbiome}\n"
printf "        MetPath=${MetPath}\n"
printf "        ClimateFile=${ClimateFile}\n"
printf "        TransitionFilePath=${TransitionFilePath}\n"
printf "        gm_lut_bernacchi_2002=${gm_lut_bernacchi_2002}\n"
printf "        gm_lut_walker_2013=${gm_lut_walker_2013}\n"
printf "        filename_d13c_atm=${filename_d13c_atm}\n"
printf "\n"



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
applysed ${tmp}/sedtmp.${pid} ${ndir}/cru.nml ${rdir}/cru_${experiment}.nml


# global landuse change namelist
cat > ${tmp}/sedtmp.${pid} << EOF
    TransitionFilePath = "${TransitionFilePath}"
    ClimateFile        = "${ClimateFile}"
    YearStart          = 1700
    YearEnd            = 2021
EOF
applysed ${tmp}/sedtmp.${pid} ${ndir}/luc.nml ${rdir}/luc_${experiment}.nml


# global Cable namelist
if [[ "${Rubisco_params}" == "Bernacchi_2002" ]] ; then
    filename_gm_lut=${gm_lut_bernacchi_2002}
elif [[ "${Rubisco_params}" == "Walker_2013" ]] ; then
    filename_gm_lut=${gm_lut_walker_2013}
else
    filename_gm_lut=""
fi

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
    output%grid                        = "land"
    output%vars5D                      = .FALSE.
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


# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# Start Runs according to sequence specified above
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------

# delete all restart files if required
if [[ ${purge_restart} -eq 1 ]] ; then
    rm -f ${rdir}/restart/*
    #rm -f ${ClimateFile}
fi

# --------------------------------------------------------------------
# 1. Create climate restart file
# --------------------------------------------------------------------
if [[ ${doclimate} -eq 1 ]] ; then
    echo "1. Create climate restart file"
    rid="climate_restart"
    
    # Met forcing
    cp ${rdir}/cru_${experiment}.nml ${rdir}/cru.nml
	
    # LUC
    cp ${rdir}/luc_${experiment}.nml ${rdir}/luc.nml
    
    # Cable
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
    #valgrind --tool=massif --xtree-memory=full --pages-as-heap=yes ./${iexe} > logs/log_out_cable.txt
    ./${iexe} > logs/log_out_cable.txt
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
# --------------------------------------------------------------------
if [[ ${dofromzero} -eq 1 ]] ; then
    echo "2. First spinup from zero biomass"
    rid="zero_biomass"

    # Met forcing
    cp ${rdir}/cru_${experiment}.nml ${rdir}/cru.nml

    # LUC
    cp ${rdir}/luc_${experiment}.nml ${rdir}/luc.nml
    # Cable
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
    #valgrind --tool=massif --xtree-memory=full --pages-as-heap=yes ./${iexe} > logs/log_out_cable.txt
    ./${iexe} > logs/log_out_cable.txt
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
# --------------------------------------------------------------------
if [[ ${doequi1} -eq 1 ]] ; then
    echo "3. Bring biomass into quasi-equilibrium with unrestricted N and P pools"
    for ((iequi1=1; iequi1<=${nequi1}; iequi1++)) ; do
        # 3a. 30 year run starting from restart files
        echo "   3a. 30 year spinup from accumulated biomass; iequi1=${iequi1}/${nequi1}"
        rid="spinup_limit_labile_${iequi1}"
        # rid="spinup_limit_labile${iequi}"

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
        ./${iexe} > logs/log_out_cable.txt
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
        ./${iexe} > logs/log_out_cable.txt
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
        #cd ../outputs
        #renameid ${rid} ${mettype}_out_casa.nc
        #if [[ ${doc13o2} -eq 1 ]] ; then renameid ${rid} ${mettype}_out_casa_c13o2.nc ; fi
        #cd ..
        cd ${pdir}
    done
fi



# ------------------------------------------------------------------------
# 4. Biomass into quasi-equilibrium without restricted N and P pools
# ------------------------------------------------------------------------
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
        ./${iexe} > logs/log_out_cable.txt
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
        ./${iexe} > logs/log_out_cable.txt
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
        #if [[ ${dompi} -eq 0 ]] ; then # no output only restart if MPI
        #    cd ../outputs
        #    renameid ${rid} ${mettype}_out_casa.nc
	    #if [[ ${doc13o2} -eq 1 ]] ; then renameid ${rid} ${mettype}_out_casa_c13o2.nc ; fi
        #    cd ..
        #fi
        cd ${pdir}
    done
fi



# --------------------------------------------------------------------
# 5. First dynamic land use (Initialise land use)
# --------------------------------------------------------------------
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
    ./${iexe} > logs/log_out_cable.txt
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
# --------------------------------------------------------------------
if [[ ${doinidyn} -eq 1 ]] ; then
    echo "6. Transient run (full dynamic spinup)"

    # Met forcing
	YearStart=1700
	YearEnd=1900
	#YearEnd=1709 
	         
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
    ./${iexe} > logs/log_out_cable.txt
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
# --------------------------------------------------------------------
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
    ./${iexe} > logs/log_out_cable.txt
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
