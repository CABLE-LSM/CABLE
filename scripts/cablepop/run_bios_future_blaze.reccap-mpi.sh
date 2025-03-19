#!/usr/bin/env bash

# Gadi
# https://opus.nci.org.au/display/Help/How+to+submit+a+job
 #PBS -N BIOS_BLAZE_tests_mpi_reccap
 #PBS -P x45
 # express / normal / copyq (2x24, cascadelake)
 #PBS -q normal
 # Typical for global or Aust continent at 0.25, 192 GB memory and 48 cpus,
 # maybe 12 hours walltime
# Typical for small runs, requires fewer cpus than pixels
 #PBS -l walltime=24:00:00
 #PBS -l mem=96GB
 #PBS -l ncpus=48
 # #PBS -l jobfs=1GB
 #PBS -l storage=gdata/x45
 #PBS -l software=netCDF:MPI:Intel:GNU
 #PBS -r y
 #PBS -l wd
 #PBS -j oe
 #PBS -S /bin/bash
 #PBS -M ian.harman@csiro.au
 #PBS -m ae

# script varied from original run_cable-pop.sh script for BIOS3 future runs - most of systems/other users stuff stripped out
# inh599@gadi harman@gadi 2023 

system=inh599@gadi

# MPI run or single processor run
# nproc should fit with job tasks
dompi=1   # 0: normal run: ./cable
          # 1: MPI run: mpiexec -n ${nproc} ./cable_mpi
nproc=48  # Number of cores for MPI runs
          # must be same as above: SBATCH -n nproc or PBS -l ncpus=nproc
          # 4 for acctest9 and 48 for other cases 

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
#      b) Run from 1700 to 1950 with dynamic land use, varying atmospheric CO2 and N deposition,
#         but still with 30 years of repeated meteorology.
#   6. Final historical run
#          a) using AGCD forcing from 1951 to GCMstart, everything dynamic
#          b) CCAM-derived historical data from GCMstart to GCMswitch, everything dynamic
#   7. Future run, everything dynamic using CCAM-derived forcing
#
# Written,  Matthias Cuntz, Aug 2019, following the run scripts,namelists provided by V Haverd
# Modified, Jurgen Knauer, 2020      - gm_explicit, coordination, acclimation
#                                    - bios, plume, future runs
#           Matthias Cuntz, Mar 2021 - functions into run_cable-pop_lib.sh
#           Ian Harman, October 2023 - streamlined for BIOS binary only
#                                    - and CCAM-derived future runs.
#
# --------------------------------------------------------------------
# Sequence switches
#
# imeteo no longer active - removed
# doextract no longer active - removed
# randompoints no longer active - removed

experiment=BLAZE25_mpi_newswitch_reccap  #experiment name

# Step 0
purge_restart=0 # 1/0: Do/Do not delete all restart files (completely new run, e.g. if settings changed)
# Step 1
doclimate=1     # 1/0: Do/Do not create climate restart file
# Step 2
dofromzero=1    # 1/0: Do/Do not first spinup phase from zero biomass stocks
# Step 3
doequi1=1       # 1/0: Do/Do not bring biomass stocks into quasi-equilibrium with restricted P and N pools
nequi1=1       #      number of times to repeat steps in doequi1
# Step 4
doequi2=1       # 1/0: Do/Do not bring biomass stocks into quasi-equilibrium with unrestricted P and N pools
nequi2=1        #      number of times to repeat steps in doequi2
# Step 5a
doiniluc=0      # 1/0: Do/Do not spinup with dynamic land use (5a)
# Step 5b
doinidyn=1      # 1/0: Do/Do not full dynamic spinup from 1700 to 1899 (5b)
# Step 6
dofinal1=1     # 1/0: Do/Do not final run from 1900 to GCMstart
dofinal2=0      # 1/0: Do/Do not final run from GCMstart to GCMswitch
# Step 7
dofuture=0      # 1/0: Do/Do not future runs GCMswtich to GCMend

# --------------------------------------------------------------------
# Other switches
restarttype='None' # 'None' 'AGCD_1951' 'AGCD_1978' 'GCM_2015'
landmasktype='land' # 'mask' 'land' # determines the GlobalLandMask.  
                    # Must match output%grid: mask = gridded, land = points, i.e act9 and reccap

# MetType
mettype='bios'        # 'bios' only
domain='reccap1000pts' #  bios domain - 'acttest9','reccap1000pts','aust_0.25_pts','australia'
GCM=''  # 'NCC-NorESM2-MM' 'ECMWF-ERA5' 'CNRM-ESM2-2'
RCP='historical'           # 'historical', 'ssp126', 'ssp370', 'evaluation' - note no future run if historical
GCMstart=2024         #  start year for CCAM derived meteorology - must be between 1951-2014 (GCM), 1979-2021 (ERA5), 2024 for AGCD runs
GCMsw=            #  start year of 'future' in CCAM derived meteorology (2015 for CCAM met, 2022 for ERA5)
GCMend=           #  end year for CCAM derived meteorology - must be less than 2100

# Cable science switches
explicit_gm=0       # 1/0: explicit (finite) or implicit mesophyll conductance
use_LUTgm=1         # 1/0: Do/Do not use lookup table for parameter conversion accounting for gm (only used if explicit_gm=1)
Rubisco_params="Bernacchi_2002"   # "Bernacchi_2002" or "Walker_2013"
coordinate_photosyn=1 # 1/0: Do/Do not coordinate photosynthesis
coord=T               # T/F: version of photosyn. optimisation (optimised(F) or forced (T))
acclimate_photosyn=1  # 1/0: Do/Do not acclimate photosynthesis
call_pop=1          # 1/0: Do/Do not use POP population dynamics model, coupled to CASA
call_blaze=1        # 1/0 Do/Do not use BLAZE fire model.
doc13o2=0           # 1/0: Do/Do not calculate 13C
c13o2_simple_disc=0 # 1/0: simple or full 13C leaf discrimination

#INH need to add a top level switch for S0, S1, S2, S3 type experiments? 

#switch dependent stuff
if [[ "${domain}" == "aust_0.25_pts" ]] ; then
    degrees=0.25
else
    degrees=0.05
fi

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
# Special things on specific computer system such as loading modules ----------------
#

export mpiexecdir=

#if [[ "${sys}" == "pearcey" ]] ; then
#    # prog is slurm_script
#    pdir=${isdir}
#    module del intel-cc intel-fc
#    module add intel-cc/16.0.1.150 intel-fc/16.0.1.150
#    module unload intel-mpi/5.0.1.035
#    module add netcdf/4.3.3.1 openmpi/1.8.8
if [[ "${sys}" == "gadi" ]] ; then
# INH I couldn't get this if/fi to work - ?due to not running through PBS
# ACB: got it working running through PBS.
    pdir=${isdir}
    #. /etc/bashrc
    module purge
    # module load intel-compiler/2019.5.281
    # module load intel-mpi/2019.5.281
    # module load netcdf/4.6.3
    # module load intel-compiler/2021.5.0
    # module load intel-mpi/2021.5.1
    # module load netcdf/4.8.0
    # # module load hdf5/1.10.5
    module load intel-compiler-llvm/2023.0.0
    module load intel-mpi/2021.8.0
    module load netcdf/4.9.2
    #if [[ ${randompoints} -eq 0 ]] ; then module load nco/4.9.2 ; fi  # needed for cropping outputs
    export mpiexecdir=/apps/intel-mpi/2019.5.281/intel64/bin
fi

if [[ ! -z ${mpiexecdir} ]] ; then export mpiexecdir="${mpiexecdir}/" ; fi

#
# Directories of things ----------------------------------------------------------
#
#   Relative directories must be relative to the directory of this script,
#   not relative to the directory from which this script is launched (if different)
#   nor relative to the run path.  EDIT: JUST HARD_WIRE THE FULL PATHS
#
if [[ "${system}" == "inh599@gadi" || "${system}" == "harman@gadi" ]] ; then
    # Run directory: runpath="${sitepath}/run"
    #sitepath="/g/data/x45/BIOS3_output/${experiment}" # Results
    sitepath="/scratch/x45/inh599/BIOStests/${experiment}" # Results
    workpath="/home/599/inh599/BLAZE25/runs/"              # run directory
    #workpath="/home/563/ab7412/CABLE_run/BIOS/CCAM"       # other scripts and params directory
    cablehome="/home/599/inh599/BLAZE25/CABLE/"            # source code and exe
elif [[ "${system}" == "ab7412@gadi" ]] ; then
     # Run directory: runpath="${sitepath}/run"
    #sitepath="/g/data/x45/BIOS3_output/${experiment}" # Results
    sitepath="/scratch/x45/ab7412/CABLE_BIOS/${experiment}" # Results
    workpath="/home/563/ab7412/CABLE_run/BIOS/CCAM" # run directory
    cablehome="/home/563/ab7412/cable_code/Github/CABLE"
else
    echo "System not known."
    exit 1
fi

# Cable executable
if [[ ${dompi} -eq 1 ]] ; then
        exe="${cablehome}/offline/cable-mpi"
  else
        exe="${cablehome}/offline/cable"
fi

# CABLE-AUX directory (uses offline/gridinfo_CSIRO_1x1.nc and offline/modis_phenology_csiro.txt)
  aux="/g/data/x45/CABLE-AUX"
  BlazeDataPath="/g/data/x45/Data_BLAZE"

# Global Mask
  SurfaceFile="${aux}/offline/gridinfo_CSIRO_CRU05x05_4tiles.nc"   # note that SurfaceFile does not need subsetting

# Met is assumed to be BIOS
  if [[ "${mettype}" == "bios" ]] ; then
   	GlobalMetPath="/g/data/x45/BIOS3_forcing/${domain}/met/"                # last slash is needed | ACB: metpath for spinup
    ParamPath="/g/data/x45/BIOS3_forcing/${domain}/params/"                 # only in bios.nml
    GlobalTransitionFilePath="/g/data/x45/LUH2/v3h/${degrees}deg_aust/EXTRACT" #LUC information

    #set landmask dependent up on land or mask switch
    if [[ "${landmasktype}" == "mask" ]] ; then
        degrees_msk=${degrees/0./} # remove 0. from degrees
        GlobalLandMaskFile="/g/data/x45/BIOS3_forcing/${domain}/australia_op_maskv2ctr${degrees_msk}"      # no file extension | ACB: Landmask file. 
    elif [[ "${landmasktype}" == "land" ]] ; then
        GlobalLandMaskFile="/g/data/x45/BIOS3_forcing/${domain}/${domain}"      # no file extension | ACB: Landmask file. 
else
        echo "Landmasktype not known. Check Landmasktype switch is defined correctly."
        exit 1
    fi
   
  fi


# Run directory
runpath="${sitepath}/run"

# Cable parameters
if [[ "${mettype}" == 'bios' ]] ; then
    namelistpath="${workpath}/namelists_bios"
    filename_veg="${workpath}/params_bios/def_veg_params.txt"
    filename_soil="${workpath}/params_bios/def_soil_params.txt"
    casafile_cnpbiome="${workpath}/params_bios/trendy_v11_pftlookup.csv"
fi

# Other scripts
ScriptsPath="${cablehome}/scripts"

# Mask (AB: I think this can be removed as it's overwritten further down)
LandMaskFile="${sitepath}/mask/${experiment}_landmask.nc"

# Met
if [[ "${mettype}" == 'bios' ]] ; then
    MetPath="${sitepath}/met/bios_${degrees}deg"
    ClimateFile="${sitepath}/mask/bios_climate_rst.nc"
    if [[ "${sys}" == "gadi" ]]; then
        if [[ (${doclimate} -eq 0) && (! -f ${ClimateFile}) ]] ; then
            ClimateFile="/g/data/x45/BIOS3_output/bio_climate_acttest9/bios_climate_rst.nc"
        fi
    fi
fi

# LUC (AB - I think this can be removed as it is overwritten further down)
TransitionFilePath="${sitepath}/LUH2/v3/1deg"

# gm lookup tables
gm_lut_bernacchi_2002="${cablehome}/params/gm_LUT_351x3601x7_1pt8245_Bernacchi2002.nc"
gm_lut_walker_2013="${cablehome}/params/gm_LUT_351x3601x7_1pt8245_Walker2013.nc"

# 13C
filename_d13c_atm="${cablehome}/params/graven_et_al_gmd_2017-table_s1-delta_13c-1700-2025.txt"

# --------------------------------------------------------------------
# Start Script
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Helper functions, most functions are in run_cable-pop_lib.sh
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
printf "            experiment=${experiment}\n"
#printf "            randompoints=${randompoints}\n"
printf "            latlon=${latlon}\n"
printf "        doclimate=${doclimate}\n"
printf "        dofromzero=${dofromzero}\n"
printf "        doequi1=${doequi1}\n"
printf "            nequi1=${nequi1}\n"
printf "        doequi2=${doequi2}\n"
printf "            nequi2=${nequi2}\n"
printf "        doiniluc=${doiniluc}\n"
printf "        doinidyn=${doinidyn}\n"
printf "        dofinal1=${dofinal1}\n"
printf "        dofinal2=${dofinal2}\n"
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
# Prep input - all done in preprocessing
#
# --------------------------------------------------------------------
# Prepare sequence
#
#BIOS configuration paths - historical meteorology and landmask
MetPath=$(abspath ${GlobalMetPath})
TransitionFilePath=$(abspath ${GlobalTransitionFilePath})
LandMaskFile=$(absfile ${GlobalLandMaskFile})

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
fi

# make a copy of the run script and save into the log directory

rs_name=$(basename "$0")
new_name="run_script_${rs_name%.sh}.sh"
copy_path="${rdir}/logs/${new_name}"
cp "$0" "$copy_path"
chmod 755 "$copy_path"
echo "run script copied to:" $copy_path

#make a copy of the parameter files and save into the log directory
copy_path="${rdir}/logs"
cp ${filename_veg} "$copy_path"
cp ${filename_soil} "$copy_path"
cp ${casafile_cnpbiome} "$copy_path"
echo "parameter file " ${filename_veg} " copied to: " $copy_path
echo "parameter file " ${filename_soil} " copied to: " $copy_path
echo "parameter file " ${casafile_cnpbiome} " copied to: " $copy_path

# Write standard namelists with options that are common to all steps of the sequence.
# They can, however, be overwritten in later steps.

# global meteo namelist file - start with the historical information
if [[ "${mettype}" == "bios" ]] ; then
    cat > ${tmp}/sedtmp.${pid} << EOF
        met_path         = "${MetPath}/"
        param_path       = "${ParamPath}"
        landmaskflt_file = "${GlobalLandMaskFile}.flt"
        landmaskhdr_file = "${GlobalLandMaskFile}.hdr"
        rain_file        = "1900010120231231_rain_recal_b2405.bin"
        swdown_file      = "1900010120231231_rad_b2405.bin"
        tairmax_file     = "1900010120231231_tmax_noclim_b2405.bin"
        tairmin_file     = "1900010120231231_tmin_noclim_b2405.bin"
        wind_file        = "1900010120231231_windspeed_ms_b2405.bin"
        vp0900_file      = "1900010120231231_vph09_b2405.bin"
        vp1500_file      = "1900010120231231_vph15_b2405.bin"
        co2_file         = "1700_2023_trendy_global_co2_ann.bin"
EOF
    applysed ${tmp}/sedtmp.${pid} ${ndir}/bios.nml ${rdir}/bios_${experiment}.nml
fi

# global landuse change namelist
cat > ${tmp}/sedtmp.${pid} << EOF
    TransitionFilePath = "${TransitionFilePath}"
    ClimateFile        = "${ClimateFile}"
    YearStart          = 1700
    YearEnd            = 2017
EOF
applysed ${tmp}/sedtmp.${pid} ${ndir}/LUC.nml ${rdir}/LUC_${experiment}.nml

# Blaze namelist !CLN CHECK - NB blazeTStep is overwritten in code
cat > ${tmp}/sedtmp.${pid} << EOF
    blazeTStep       = "annually"  ! Call frequency ("daily", "monthly", "annually")
    BurnedAreaSource = "SIMFIRE"   ! Burnt Area ("PRESCRIBED", "SIMFIRE", "GFED4")
    BurnedAreaFile   = "${BlazeDataPath}/BA_Aust_2001-2019.nc"  ! used for Prescribed fires !CLN not available for now!
    SIMFIRE_REGION   = "ANZ"       ! ("ANZ", "EUROPE", "GLOBAL")
    HydePath         = "${BlazeDataPath}/HYDE3.1"  ! Path to Hyde3.1 population density data
    BurnedAreaClimatologyFile = "${BlazeDataPath}/simfire_monthly_ba.nc"  ! BA climatology file (needed when blazeTStep!="annually")
EOF
applysed ${tmp}/sedtmp.${pid} ${ndir}/blaze.nml ${rdir}/blaze_${experiment}.nml
cp ${rdir}/blaze_${experiment}.nml ${rdir}/blaze.nml
#if BLAZE active then set rblaze=1 for use in all steps (else=0)
rblaze=0
if [[ ${call_blaze} -eq 1 ]] ; then
    rblaze=1
fi

# global Cable namelist
cat > ${tmp}/sedtmp.${pid} << EOF
    filename%met                       = "${mettype}"
    filename%veg                       = "${filename_veg}"
    filename%soil                      = "${filename_soil}"
    filename%type                      = "${SurfaceFile}"
    filename%out                       = "outputs/${mettype}_out_cable.nc"
    filename%restart_in                = "restart/${mettype}_cable_rst.nc"
    filename%restart_out               = "restart/${mettype}_cable_rst.nc"
    casafile%cnpbiome                  = "${casafile_cnpbiome}"
    casafile%out                       = "outputs/${mettype}_out_casa.nc"
    casafile%cnpipool                  = "restart/${mettype}_casa_rst.nc"
    casafile%cnpepool                  = "restart/${mettype}_casa_rst.nc"
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
    leaps                              = .false.
    cable_user%SOIL_STRUC              = "sli"
    cable_user%Rubisco_parameters      = "${Rubisco_params}"
    cable_user%CALL_POP                = .false.
    cable_user%coordinate_photosyn     = .false.
    cable_user%acclimate_photosyn      = .false.
    cable_user%explicit_gm             = .false.
    cable_user%gm_LUT_file             = "${filename_gm_lut}"
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
    cable_user%CALL_BLAZE              = ${rblaze}
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
applysed ${tmp}/sedtmp.${pid} ${ndir}/cable.nml ${rdir}/cable_${experiment}.nml


# --------------------------------------------------------------------
# Sequence
#


# --------------------------------------------------------------------
# 1. Create climate restart file
if [[ ${doclimate} -eq 1 ]] ; then
    echo "1. Create climate restart file"
    rid="climate_restart"
    # Met forcing
    if [[ "${mettype}" == "bios" ]] ; then
        cat > ${tmp}/sedtmp.${pid} << EOF
	          Run = "spinup"
EOF
        applysed ${tmp}/sedtmp.${pid} ${rdir}/bios_${experiment}.nml ${rdir}/bios.nml
    fi
    # LUC
    cp ${rdir}/LUC_${experiment}.nml ${rdir}/LUC.nml
    # Cable
    #   do not calculate 13C because there is no 13C in the climate restart file
    cat > ${tmp}/sedtmp.${pid} << EOF
        filename%restart_in            = ""
        cable_user%CLIMATE_fromZero    = .true.
        cable_user%YearStart           = 1951
        cable_user%YearEnd             = 1980
        icycle                         = 2
        spincasa                       = .false.
        cable_user%CASA_fromZero       = .true.
        cable_user%CASA_DUMP_READ      = .false.
        cable_user%CASA_DUMP_WRITE     = .true.
        cable_user%CASA_SPIN_STARTYEAR = 1951
        cable_user%CASA_SPIN_ENDYEAR   = 1960
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
    renameid ${rid} ${mettype}.nml LUC.nml cable.nml
    imv *_${rid}.nml restart/
    cd logs
    renameid ${rid} log_cable.txt log_out_cable.txt
    cd ../restart
    copyid ${rid} ${mettype}_climate_rst.nc
    cp ${mettype}_climate_rst.nc ${ClimateFile}
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
    if [[ "${mettype}" == "bios" ]] ; then
        cat > ${tmp}/sedtmp.${pid} << EOF
	          Run = "spinup"
EOF
        applysed ${tmp}/sedtmp.${pid} ${rdir}/bios_${experiment}.nml ${rdir}/bios.nml
    fi
    # LUC
    cp ${rdir}/LUC_${experiment}.nml ${rdir}/LUC.nml
    # Cable namelist
    cat > ${tmp}/sedtmp.${pid} << EOF
        filename%restart_in               = ""
        cable_user%CLIMATE_fromZero       = .true.
        cable_user%YearStart              = 1951
        cable_user%YearEnd                = 1980
        icycle                            = 2
        spincasa                          = .false.
        cable_user%CASA_OUT_FREQ          = "monthly"
        cable_user%CASA_fromZero          = .true.
        cable_user%CASA_DUMP_READ         = .false.
        cable_user%CASA_DUMP_WRITE        = .true.
        output%averaging                  = "monthly"
        cable_user%CASA_SPIN_STARTYEAR    = 1951
        cable_user%CASA_SPIN_ENDYEAR      = 1960
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
    renameid ${rid} ${mettype}.nml LUC.nml cable.nml
    imv *_${rid}.nml restart/
    cd logs
    renameid ${rid} log_cable.txt log_out_cable.txt
    cd ../restart
    copyid ${rid} ${mettype}_climate_rst.nc ${mettype}_casa_rst.nc ${mettype}_cable_rst.nc pop_${mettype}_ini.nc
    copyid ${rid} ${mettype}_c13o2_flux_rst.nc ${mettype}_c13o2_pools_rst.nc
    cd ../outputs
    renameid ${rid} ${mettype}_out_cable.nc ${mettype}_out_casa.nc ${mettype}_out_casa_c13o2.nc
    cd ..
    cd ${pdir}
fi


# --------------------------------------------------------------------
# 3. Biomass into quasi-equilibrium with restricted N and P pools
if [[ ${doequi1} -eq 1 ]] ; then
    echo "3. Bring biomass into quasi-equilibrium with restricted N and P pools"
    for ((iequi1=1; iequi1<=${nequi1}; iequi1++)) ; do
        # 3a. 30 year run starting from restart files
        echo "   3a. 30 year spinup from accumulated biomass; iequi1=${iequi1}/${nequi1}"
        #rid="spinup_limit_labile"
        rid="spinup_limit_labile${iequi1}"
        # Met forcing
        if [[ "${mettype}" == "bios" ]] ; then
            cat > ${tmp}/sedtmp.${pid} << EOF
    	          Run = "spinup"
EOF
            applysed ${tmp}/sedtmp.${pid} ${rdir}/bios_${experiment}.nml ${rdir}/bios.nml
        fi
        # LUC
        cp ${rdir}/LUC_${experiment}.nml ${rdir}/LUC.nml
        # Cable namelist
        cat > ${tmp}/sedtmp.${pid} << EOF
            cable_user%CLIMATE_fromZero    = .false.
            cable_user%YearStart           = 1951
            cable_user%YearEnd             = 1980
            icycle                         = 2
            spincasa                       = .false.
            cable_user%CASA_fromZero       = .false.
            cable_user%CASA_DUMP_READ      = .false.
            cable_user%CASA_DUMP_WRITE     = .true.
            cable_user%CASA_SPIN_STARTYEAR = 1951
            cable_user%CASA_SPIN_ENDYEAR   = 1960
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
        renameid ${rid} ${mettype}.nml LUC.nml cable.nml
        mv *_${rid}.nml restart/
        cd logs
        renameid ${rid} log_cable.txt log_out_cable.txt
        cd ../restart
        copyid ${rid} ${mettype}_climate_rst.nc ${mettype}_casa_rst.nc ${mettype}_cable_rst.nc pop_${mettype}_ini.nc
        copyid ${rid} ${mettype}_c13o2_flux_rst.nc ${mettype}_c13o2_pools_rst.nc
        cd ../outputs
        renameid ${rid} ${mettype}_out_cable.nc ${mettype}_out_casa.nc ${mettype}_out_casa_c13o2.nc
        cd ..
        cd ${pdir}

	#
        # 3b. analytic quasi-equilibrium of biomass pools
        echo "   3b. Analytic solution of biomass pools"
        rid="spin_casa_limit_labile${iequi1}"
        # Met forcing
        if [[ "${mettype}" == "bios" ]] ; then
            cat > ${tmp}/sedtmp.${pid} << EOF
    	          Run = "spinup"
EOF
            applysed ${tmp}/sedtmp.${pid} ${rdir}/bios_${experiment}.nml ${rdir}/bios.nml
        fi
        # LUC
        cp ${rdir}/LUC_${experiment}.nml ${rdir}/LUC.nml
        # Cable
        cat > ${tmp}/sedtmp.${pid} << EOF
            cable_user%CLIMATE_fromZero    = .false.
            cable_user%YearStart           = 1951
            cable_user%YearEnd             = 1980
            icycle                         = 12
            spincasa                       = .true.
            cable_user%CASA_fromZero       = .false.
            cable_user%CASA_DUMP_READ      = .true.
            cable_user%CASA_DUMP_WRITE     = .false.
            cable_user%CASA_SPIN_STARTYEAR = 1951
            cable_user%CASA_SPIN_ENDYEAR   = 1960
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
        renameid ${rid} ${mettype}.nml LUC.nml cable.nml
        mv *_${rid}.nml restart/
        cd logs
        renameid ${rid} log_cable.txt log_out_cable.txt
        cd ../restart
        copyid ${rid} ${mettype}_casa_rst.nc pop_${mettype}_ini.nc
        copyid ${rid} ${mettype}_c13o2_flux_rst.nc ${mettype}_c13o2_pools_rst.nc
        if [[ ${dompi} -eq 0 ]] ; then # no output only restart if MPI
            cd ../outputs
            renameid ${rid} ${mettype}_out_casa.nc ${mettype}_out_casa_c13o2.nc
            cd ..
        fi
        cd ${pdir}
    done
fi


# --------------------------------------------------------------------
# 4. Biomass into quasi-equilibrium without restricted N and P pools
# if BLAZE active then for doequi2 onwards switch to fully active (=3)
if [[ ${call_blaze} -eq 1 ]] ; then
    rblaze=3
fi
if [[ ${doequi2} -eq 1 ]] ; then
    echo "4. Bring biomass into quasi-equilibrium without restricted N and P pools"
    for ((iequi2=1; iequi2<=${nequi2}; iequi2++)) ; do
        # 4a. 30 year run starting from restart files
        echo "   4a. 30 year spinup from accumulated biomass; iequi2=${iequi2}/${nequi2}"
        
        rid="spinup_nutrient_limited${iequi2}"
        # Met forcing
        if [[ "${mettype}" == "bios" ]] ; then
            cat > ${tmp}/sedtmp.${pid} << EOF
    	          Run = "spinup"
EOF
            applysed ${tmp}/sedtmp.${pid} ${rdir}/bios_${experiment}.nml ${rdir}/bios.nml
        fi
        # LUC
        cp ${rdir}/LUC_${experiment}.nml ${rdir}/LUC.nml
        # Cable
        cat > ${tmp}/sedtmp.${pid} << EOF
            cable_user%CLIMATE_fromZero    = .false.
            cable_user%YearStart           = 1951
            cable_user%YearEnd             = 1980
            icycle                         = 2
            spincasa                       = .false.
            cable_user%CASA_fromZero       = .false.
            cable_user%CASA_DUMP_READ      = .false.
            cable_user%CASA_DUMP_WRITE     = .true.
            cable_user%CASA_SPIN_STARTYEAR = 1951
            cable_user%CASA_SPIN_ENDYEAR   = 1960
            cable_user%limit_labile        = .false.
            cable_user%POP_fromZero        = .false.
            cable_user%POP_out             = "ini"
            cable_user%POPLUC              = .true.
            cable_user%POPLUC_RunType      = "static"
            cable_user%CALL_BLAZE          = ${rblaze}
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
        renameid ${rid} ${mettype}.nml LUC.nml cable.nml
        mv *_${rid}.nml restart/
        cd logs
        renameid ${rid} log_cable.txt log_out_cable.txt
        cd ../restart
        copyid ${rid} ${mettype}_climate_rst.nc ${mettype}_casa_rst.nc ${mettype}_cable_rst.nc pop_${mettype}_ini.nc
        copyid ${rid} ${mettype}_c13o2_flux_rst.nc ${mettype}_c13o2_pools_rst.nc
        cd ../outputs
        renameid ${rid} ${mettype}_out_cable.nc ${mettype}_out_casa.nc ${mettype}_out_casa_c13o2.nc
        cd ..
        cd ${pdir}
        #
        # 4b. analytic quasi-equilibrium of biomass pools
        echo "   4b. Analytic solution of biomass pools"
        rid="spin_casa_nutrient_limited${iequi2}"
	
        # Met forcing
        if [[ "${mettype}" == "bios" ]] ; then
            cat > ${tmp}/sedtmp.${pid} << EOF
    	          Run = "spinup"
EOF
            applysed ${tmp}/sedtmp.${pid} ${rdir}/bios_${experiment}.nml ${rdir}/bios.nml
        fi
        # LUC
        cp ${rdir}/LUC_${experiment}.nml ${rdir}/LUC.nml
        # Cable
        cat > ${tmp}/sedtmp.${pid} << EOF
            cable_user%CLIMATE_fromZero    = .false.
            cable_user%YearStart           = 1951
            cable_user%YearEnd             = 1980
            icycle                         = 12
            spincasa                       = .true.
            cable_user%CASA_fromZero       = .false.
            cable_user%CASA_DUMP_READ      = .true.
            cable_user%CASA_DUMP_WRITE     = .false.
            cable_user%CASA_SPIN_STARTYEAR = 1951
            cable_user%CASA_SPIN_ENDYEAR   = 1960
            cable_user%limit_labile        = .false.
            cable_user%POP_fromZero        = .false.
            cable_user%POP_out             = "ini"
            cable_user%POPLUC              = .true.
            cable_user%POPLUC_RunType      = "static"
            cable_user%CALL_BLAZE          = ${rblaze}
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
        renameid ${rid} ${mettype}.nml LUC.nml cable.nml
        mv *_${rid}.nml restart/
        cd logs
        renameid ${rid} log_cable.txt log_out_cable.txt
        cd ../restart
        copyid ${rid} ${mettype}_casa_rst.nc pop_${mettype}_ini.nc
        copyid ${rid} ${mettype}_c13o2_flux_rst.nc ${mettype}_c13o2_pools_rst.nc
        if [[ ${dompi} -eq 0 ]] ; then # no output only restart if MPI
            cd ../outputs
            #renameid ${rid} ${mettype}_out_casa.nc ${mettype}_out_casa_c13o2.nc
            cd ..
        fi
        cd ${pdir}
    done
fi

# --------------------------------------------------------------------
# 5a. First dynamic land use
if [[ ${doiniluc} -eq 1 ]] ; then
    echo "5a. First dynamic land use"
    # Met forcing
    if [[ "${mettype}" == "bios" ]] ; then
        YearStart=1580
        YearEnd=1699
        cat > ${tmp}/sedtmp.${pid} << EOF
	          Run = "spinup"
EOF
        applysed ${tmp}/sedtmp.${pid} ${rdir}/bios_${experiment}.nml ${rdir}/bios.nml
    fi
    rid=${YearStart}_${YearEnd}
    # LUC
    cat > ${tmp}/sedtmp.${pid} << EOF
         YearStart = ${YearStart}
         YearEnd   = ${YearEnd}
EOF
    applysed ${tmp}/sedtmp.${pid} ${rdir}/LUC_${experiment}.nml ${rdir}/LUC.nml

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
        cable_user%CASA_SPIN_STARTYEAR  = 1951
        cable_user%CASA_SPIN_ENDYEAR    = 1960
        cable_user%limit_labile         = .false.
        cable_user%POP_fromZero         = .false.
        cable_user%POP_out              = "ini"
        cable_user%POPLUC               = .true.
        cable_user%POPLUC_RunType       = "init"
        cable_user%LUC_restart_in       = ""
        cable_user%c13o2_restart_in_luc = ""
        cable_user%CALL_BLAZE          = ${rblaze}
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
    renameid ${rid} ${mettype}.nml LUC.nml cable.nml
    mv *_${rid}.nml restart/
    cd logs
    renameid ${rid} log_cable.txt log_out_cable.txt
    cd ../restart
    copyid ${rid} ${mettype}_casa_rst.nc ${mettype}_LUC_rst.nc pop_${mettype}_ini.nc
    copyid ${rid} ${mettype}_c13o2_pools_rst.nc ${mettype}_c13o2_luc_rst.nc
     cd ../outputs
     renameid ${rid} ${mettype}_out_LUC.nc
     cd ..
    cd ${pdir}
fi


# --------------------------------------------------------------------
# 5b. Second full dynamic spinup
if [[ ${doinidyn} -eq 1 ]] ; then
    echo "5b. Full dynamic spinup"
    # Met forcing
    if [[ "${mettype}" == "bios" ]] ; then
        YearStart=1700
        YearEnd=1950
        cat > ${tmp}/sedtmp.${pid} << EOF
	          Run = "premet"
EOF
        applysed ${tmp}/sedtmp.${pid} ${rdir}/bios_${experiment}.nml ${rdir}/bios.nml
    fi
    rid=${YearStart}_${YearEnd}
    # LUC
    cat > ${tmp}/sedtmp.${pid} << EOF
         YearStart = 1900
         YearEnd   = 1900
EOF
    applysed ${tmp}/sedtmp.${pid} ${rdir}/LUC_${experiment}.nml ${rdir}/LUC.nml
    # Cable
    cat > ${tmp}/sedtmp.${pid} << EOF
        cable_user%CLIMATE_fromZero    = .false.
        cable_user%YearStart           = ${YearStart}
        cable_user%YearEnd             = ${YearEnd}
        icycle                         = 2
        spincasa                       = .false.
        cable_user%CASA_fromZero       = .false.
        cable_user%CASA_DUMP_READ      = .false.
        cable_user%CASA_DUMP_WRITE     = .false.
        cable_user%CASA_SPIN_STARTYEAR = 1951
        cable_user%CASA_SPIN_ENDYEAR   = 1960
        cable_user%limit_labile        = .false.
        cable_user%POP_fromZero        = .false.
        cable_user%POP_out             = "ini"
        cable_user%POPLUC              = .true.
        cable_user%POPLUC_RunType      = "static"
        cable_user%CALL_BLAZE          = ${rblaze}
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
    renameid ${rid} ${mettype}.nml LUC.nml cable.nml
    mv *_${rid}.nml restart/
    cd logs
    renameid ${rid} log_cable.txt log_out_cable.txt
    cd ../restart
    copyid ${rid} ${mettype}_climate_rst.nc ${mettype}_casa_rst.nc ${mettype}_cable_rst.nc
    copyid ${rid} ${mettype}_LUC_rst.nc pop_${mettype}_ini.nc
    copyid ${rid} ${mettype}_c13o2_flux_rst.nc ${mettype}_c13o2_pools_rst.nc ${mettype}_c13o2_luc_rst.nc
    cd ../outputs
    renameid ${rid} ${mettype}_out_cable.nc ${mettype}_out_casa.nc ${mettype}_out_LUC.nc
    renameid ${rid} ${mettype}_out_casa_c13o2.nc
    cd ..
    cd ${pdir}
fi


# --------------------------------------------------------------------
# 6a. Final run - 1951 to GCMstart
if [[ ${dofinal1} -eq 1 ]] ; then

      # check restart switch and copy restarts if required
    if [[ "${restarttype}" == "AGCD_1951" ]] ; then
    restart_path="/g/data/x45/BIOS3_forcing/CCAM/restart_files/${restarttype}/${domain}/"
    cd ${restart_path}
    cp *.nc ${rdir}/restart
    fi


    echo "6. Final run part 1"
    # Met forcing
    if [[ "${mettype}" == "bios" ]] ; then
        YearStart=1951
	value=$(expr $GCMstart - 1)
        YearEnd=${value}
        cat > ${tmp}/sedtmp.${pid} << EOF
	          Run = "standard"
EOF
        applysed ${tmp}/sedtmp.${pid} ${rdir}/bios_${experiment}.nml ${rdir}/bios.nml
    fi
    rid=${YearStart}_${YearEnd}
    # LUC
    cat > ${tmp}/sedtmp.${pid} << EOF
         YearStart = 1900
         YearEnd   = 1900
EOF
    applysed ${tmp}/sedtmp.${pid} ${rdir}/LUC_${experiment}.nml ${rdir}/LUC.nml
    # Cable
    cat > ${tmp}/sedtmp.${pid} << EOF
        cable_user%CLIMATE_fromZero    = .false.
        cable_user%YearStart           = ${YearStart}
        cable_user%YearEnd             = ${YearEnd}
        icycle                         = 2
        spincasa                       = .false.
        cable_user%CASA_fromZero       = .false.
        cable_user%CASA_DUMP_READ      = .false.
        cable_user%CASA_DUMP_WRITE     = .false.
        cable_user%CASA_SPIN_STARTYEAR = 1951
        cable_user%CASA_SPIN_ENDYEAR   = 1960
        cable_user%limit_labile        = .false.
        cable_user%POP_fromZero        = .false.
        cable_user%POP_out             = "ini"
        cable_user%POPLUC              = .true.
        cable_user%POPLUC_RunType      = "static"
        cable_user%CALL_BLAZE          = ${rblaze}
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
    rid=${YearStart}_${YearEnd}
    # save output
    renameid ${rid} ${mettype}.nml LUC.nml cable.nml
    mv *_${rid}.nml restart/
    cd logs
    renameid ${rid} log_cable.txt log_out_cable.txt
    cd ../restart
    copyid ${rid} ${mettype}_climate_rst.nc ${mettype}_casa_rst.nc ${mettype}_cable_rst.nc
    copyid ${rid} ${mettype}_LUC_rst.nc pop_${mettype}_ini.nc
    copyid ${rid} ${mettype}_c13o2_flux_rst.nc ${mettype}_c13o2_pools_rst.nc ${mettype}_c13o2_luc_rst.nc
    cd ../outputs
    renameid ${rid} ${mettype}_out_cable.nc ${mettype}_out_casa.nc ${mettype}_out_LUC.nc
    renameid ${rid} ${mettype}_out_casa_c13o2.nc
    cd ..
    cd ${pdir}
fi

#-----
#update the path to the CCAM derived meteorology
if [[ "${domain}" == "australia" ]] ; then
    domain="aust_0.05_pts"
fi
if [[ "${mettype}" == "bios" ]] ; then
    MetPath="/g/data/x45/BIOS3_forcing/CCAM/${GCM}/historical/${domain}/met/"      # last slash is needed # need to update path if evaluation
fi
#-----

# 6b. Final run - GCMstart to GCMswitch


if [[ ${dofinal2} -eq 1 ]] ; then
    echo "6. Final run - part 2"
    # check restart switch and copy restarts if required
    if [[ "${restarttype}" == "AGCD_1951" ]] ; then
    restart_path="/g/data/x45/BIOS3_forcing/CCAM/restart_files/${restarttype}/${domain}/"
    cd ${restart_path}
    cp *.nc ${rdir}/restart
    fi


    # Met forcing - update to use GCM information
    if [[ "${mettype}" == "bios" ]] ; then
        YearStart=${GCMstart}
	    value=$(expr $GCMsw - 1)
        YearEnd=${value}

	# update files
	cat > ${tmp}/sedtmp.${pid} << EOF
            met_path         = "${MetPath}"
            rain_file        = "pr_Adjust_${GCM}_historical_${YearStart}_${YearEnd}.bin" 
            swdown_file      = "solar_Adjust_${GCM}_historical_${YearStart}_${YearEnd}.bin"
            tairmax_file     = "tasmax_Adjust_${GCM}_historical_${YearStart}_${YearEnd}.bin"
            tairmin_file     = "tasmin_Adjust_${GCM}_historical_${YearStart}_${YearEnd}.bin"
            wind_file        = "sfcWind_Adjust_${GCM}_historical_${YearStart}_${YearEnd}.bin"
            vp0900_file      = "vph09_Adjust_${GCM}_historical_${YearStart}_${YearEnd}.bin"
            vp1500_file      = "vph15_Adjust_${GCM}_historical_${YearStart}_${YearEnd}.bin"
            co2_file         = "0000_2014_CO2_time_series_ccam_historical.bin"
	        Run              = "standard"
EOF
        applysed ${tmp}/sedtmp.${pid} ${rdir}/bios_${experiment}.nml ${rdir}/bios.nml
    fi
    rid=${YearStart}_${YearEnd}
    # LUC
    cat > ${tmp}/sedtmp.${pid} << EOF
         YearStart = 1900
         YearEnd   = 1900
EOF
    applysed ${tmp}/sedtmp.${pid} ${rdir}/LUC_${experiment}.nml ${rdir}/LUC.nml
    # Cable
    cat > ${tmp}/sedtmp.${pid} << EOF
        cable_user%CLIMATE_fromZero    = .false.
        cable_user%YearStart           = ${YearStart}
        cable_user%YearEnd             = ${YearEnd}
        icycle                         = 2
        spincasa                       = .false.
        cable_user%CASA_fromZero       = .false.
        cable_user%CASA_DUMP_READ      = .false.
        cable_user%CASA_DUMP_WRITE     = .false.
        cable_user%CASA_SPIN_STARTYEAR = 1951
        cable_user%CASA_SPIN_ENDYEAR   = 1960
        cable_user%limit_labile        = .false.
        cable_user%POP_fromZero        = .false.
        cable_user%POP_out             = "ini"
        cable_user%POPLUC              = .true.
        cable_user%POPLUC_RunType      = "static"
        cable_user%CALL_BLAZE          = ${rblaze}
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
    renameid ${rid} ${mettype}.nml LUC.nml cable.nml
    mv *_${rid}.nml restart/
    cd logs
    renameid ${rid} log_cable.txt log_out_cable.txt
    cd ../restart
    copyid ${rid} ${mettype}_climate_rst.nc ${mettype}_casa_rst.nc ${mettype}_cable_rst.nc
    copyid ${rid} ${mettype}_LUC_rst.nc pop_${mettype}_ini.nc
    copyid ${rid} ${mettype}_c13o2_flux_rst.nc ${mettype}_c13o2_pools_rst.nc ${mettype}_c13o2_luc_rst.nc
    cd ../outputs
    renameid ${rid} ${mettype}_out_cable.nc ${mettype}_out_casa.nc ${mettype}_out_LUC.nc
    renameid ${rid} ${mettype}_out_casa_c13o2.nc
    cd ..
    cd ${pdir}
fi



# --------------------------------------------------------------------
# 7. Future run
# update the path to the CCAM derived meteorology
if [[ "${domain}" == "australia" ]] ; then
    domain="aust_0.05_pts"
fi
if [[ "${mettype}" == "bios" ]] ; then
    MetPath="/g/data/x45/BIOS3_forcing/CCAM/${GCM}/${RCP}/${domain}/met/"               # last slash is needed
fi

if [[ (${dofuture} -eq 1) && ("${RCP}" != "historical") ]] ; then
    echo "7. Future run"
    rcpnd=$(echo ${RCP} | sed "s|\.||")
    rcpd=$(echo ${RCP} | sed "s|\.|p|")
 
    if [[ "${restarttype}" == "GCM_2014" ]] ; then
    restart_path="/g/data/x45/BIOS3_forcing/CCAM/restart_files/${GCM}/${domain}/"
    cd ${restart_path}
    cp *.nc ${rdir}/restart
    fi

    # Met forcing
    YearStart=${GCMsw}
    YearEnd=${GCMend}
    cat > ${tmp}/sedtmp.${pid} << EOF
        met_path         = "${MetPath}"
        rain_file        = "pr_Adjust_${GCM}_${RCP}_${GCMsw}_${GCMend}.bin"
        swdown_file      = "solar_Adjust_${GCM}_${RCP}_${GCMsw}_${GCMend}.bin"
        tairmax_file     = "tasmax_Adjust_${GCM}_${RCP}_${GCMsw}_${GCMend}.bin"
        tairmin_file     = "tasmin_Adjust_${GCM}_${RCP}_${GCMsw}_${GCMend}.bin"
        wind_file        = "sfcWind_Adjust_${GCM}_${RCP}_${GCMsw}_${GCMend}.bin"
        vp0900_file      = "vph09_Adjust_${GCM}_${RCP}_${GCMsw}_${GCMend}.bin"
        vp1500_file      = "vph15_Adjust_${GCM}_${RCP}_${GCMsw}_${GCMend}.bin"
	    Run              = "standard"
        RCP              = "${RCP}"
        CO2              = "varying"
        co2_file         = "2015_2500_CO2_time_series_ccam_${RCP}.bin"
        NDEP             = "varying"
        NDEPfile         = "${NdepPath}/RCP${rcpnd}/ndep_total_2000-2109_1.0x1.0_FD.nc"
EOF
    applysed ${tmp}/sedtmp.${pid} ${rdir}/bios_${experiment}.nml ${rdir}/bios.nml
    rid=${YearStart}_${YearEnd}_${rcpd}
    # LUC
    cat > ${tmp}/sedtmp.${pid} << EOF
         YearStart = 1900
         YearEnd   = 1900
EOF
    applysed ${tmp}/sedtmp.${pid} ${rdir}/LUC_${experiment}.nml ${rdir}/LUC.nml
    # Cable
    cat > ${tmp}/sedtmp.${pid} << EOF
        cable_user%CLIMATE_fromZero    = .false.
        cable_user%YearStart           = ${YearStart}
        cable_user%YearEnd             = ${YearEnd}
        icycle                         = 2
        spincasa                       = .false.
        cable_user%CASA_fromZero       = .false.
        cable_user%CASA_DUMP_READ      = .false.
        cable_user%CASA_DUMP_WRITE     = .false.
        cable_user%CASA_SPIN_STARTYEAR = 1951
        cable_user%CASA_SPIN_ENDYEAR   = 1960
        cable_user%limit_labile        = .false.
        cable_user%POP_fromZero        = .false.
        cable_user%POP_out             = "ini"
        cable_user%POPLUC              = .true.
        cable_user%POPLUC_RunType      = "static"
        cable_user%CALL_BLAZE          = ${rblaze}
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
    renameid ${rid} ${mettype}.nml LUC.nml cable.nml
    mv *_${rid}.nml restart/
    cd logs
    renameid ${rid} log_cable.txt log_out_cable.txt
    cd ../restart
    copyid ${rid} ${mettype}_climate_rst.nc ${mettype}_casa_rst.nc ${mettype}_cable_rst.nc pop_${mettype}_ini.nc
    copyid ${rid} ${mettype}_c13o2_flux_rst.nc ${mettype}_c13o2_pools_rst.nc
    cd ../outputs
    renameid ${rid} ${mettype}_out_cable.nc ${mettype}_out_casa.nc
    renameid ${rid} ${mettype}_out_casa_c13o2.nc
    cd ..
    cd ${pdir}
fi

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
