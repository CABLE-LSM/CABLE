#!/usr/bin/env bash

# Explor / Pearcey (launch with:  sbatch --ignore-pbs)
# https://slurm.schedmd.com/sbatch.html
# Name - 8 letters and digits
#SBATCH -J x0003-s0
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.out
# Explor partitions (sinfo): std (2x16, parallel), sky (2x16, parallel, AVX512), hf (2x4, serial),
#                            P100 (2x16, GPU), GTX (2x16, GPU), ivy (2x8, parallel), k20 (2x8, GPU)
#SBATCH -p std
# -N Nodes / -n tasks (mpiexec, srun, ...) / -c cpus_per_task (OpenMP, make-j, ...)
#SBATCH -N 1
#SBATCH -n 1
# #SBATCH --ntasks-per-node=32
# Check memory on *nix with /usr/bin/time -v ./prog
# time (day-hh:mm:ss) / memory (optional, units K,M,G,T)
#SBATCH -t 03:59:59
#SBATCH --mem=100G
# notify: Valid type values are NONE,BEGIN,END,FAIL,REQUEUE,ALL,STAGE_OUT,TIME_LIMIT,TIME_LIMIT_90/80/50,ARRAY_TASKS
#SBATCH --mail-type=FAIL,STAGE_OUT,TIME_LIMIT
#SBATCH --mail-user=matthias.cuntz@inrae.fr

# Gadi
# https://opus.nci.org.au/display/Help/4.+PBS+Jobs
#PBS -N x0001-s1
#PBS -P x45
# express / normal / copyq (2x24, cascadelake)
# expressbw / normalbw (2x14, broadwell) / normalsl (2x16, skylake)- ex-Raijin nodes
#PBS -q normal
# Global run or Australia at 0.25: 192 GB memory and 48 cpus, ca. 12 hours walltime
#PBS -l walltime=04:30:00
#PBS -l mem=100GB
#PBS -l ncpus=1
# #PBS -l jobfs=1GB
#PBS -l storage=gdata/x45
#PBS -l software=netCDF:MPI:Intel:GNU
#PBS -r y
#PBS -l wd
#PBS -j oe
#PBS -S /bin/bash
#PBS -M matthias.cuntz@inrae.fr
#PBS -m ae

# cuntz@explor, cuntz@mc16, cuntz@mcinra, moc801@gadi cuntz@gadi
# kna016@pearcey knauer@pearcey, jk8585@gadi knauer@gadi
# bri220@pearcey pcb599@gadi, briggs@gadi
# vil029@percey yc3714@gadi villalobos@gadi
# nieradzik@aurora
system=cuntz@gadi

# MPI run or single processor run
# nproc should fit with job tasks
dompi=0   # 0: normal run: ./cable
          # 1: MPI run: mpiexec -n ${nproc} ./cable_mpi
nproc=1   # Number of cores for MPI runs
          # must be same as above: SBATCH -n nproc or PBS -l ncpus=nproc

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
#   6. Final historical run, everything dynamic from 1900 to 2017.
#   7. Future run, everything dynamic (not all met types).
#
# Written,  Matthias Cuntz, Aug 2019, following the run scripts and namelists provided by Vanessa Haverd
# Modified, Jurgen Knauer, 2020      - gm_explicit, coordination, acclimation
#                                    - bios, plume, future runs
#           Matthias Cuntz, Mar 2021 - functions into run_cable-pop_lib.sh
#
# --------------------------------------------------------------------

#ASKJK - changes in comparison to gm_acclim_coord
# 1. Removed switch Test_parallel=1 -> used [[ dompi -eq 1 ]]
# 2. ScriptsPath="${workpath}/scripts"
#    but
#    namelistpath="$(dirname ${workpath})/namelists"
#    with
#      workpath="/home/599/jk8585/CABLE_run/gm_acclim_coord/global_runs"
#      cablehome="/home/599/jk8585/CABLE_code"
#    -> changed to similar of ScriptsPath="$(dirname ${workpath})/scripts"
# 3. Need plume.nml, bios.nml, gm_LUT_*.nc
# 4. output%grid = "mask" (in gm_acclim_coord) or "land" (default before)
# 5. What should be for Run in bios.nml after 1. Climate restart?
# 6. Why is YearEnd different for plume (1849) compared to cru (1699) in 5a. First dynamic land use?
# 7. Do we need chunking in 5a, 5b, 6, and 7?
# 8. Do we need cropping output to latlon region at the end: is this not in step 0 with ${doextractsite} -eq 1?
#ASKJK - changes in comparison to gm_acclim_coord

# --------------------------------------------------------------------
# Sequence switches
#

imeteo=0        # 0: Use global meteo, land use and mask
                # 1: Use local mask, but global meteo and land use (doextractsite=1)
                # 2: Use local meteo, land use and mask (doextractsite=2)
# Step 0
purge_restart=0 # 1/0: Do/Do not delete all restart files (completely new run, e.g. if settings changed)
doextractsite=0 # 0: Do not extract local meteo, land use nor mask
                # 1: Do extract only mask at specific site/region (imeteo=1)
                # 2: Do extract meteo, land use and mask at specific site/region (imeteo=2)
                #    Does not work with randompoints /= 0 but with latlon
    experiment=HarvardForest
    randompoints=0   # <0: use -1*randompoints from file ${LandMaskFilePath}/${experiment}_points.csv if existing
                     # 0:  use latlon
                     # >0: generate and use randompoints random grid points from GlobalLandMaskFile
    # lat,lon  or  latmin,latmax,lonmin,lonmax   # must have . in numbers otherwise indexes taken
    latlon=42.536875,-72.172602
    # latlon=-34.5,-33.5,149.5,156.5
    # latlon=42.5,43.5,109.5,110.5
    # latlon=-44.0,-10.0,110.0,155.0  # Australia
# Step 1
doclimate=1     # 1/0: Do/Do not create climate restart file
# Step 2
dofromzero=1    # 1/0: Do/Do not first spinup phase from zero biomass stocks
# Step 3
doequi1=1       # 1/0: Do/Do not bring biomass stocks into quasi-equilibrium with restricted P and N pools
nequi1=4        #      number of times to repeat steps in doequi1
# Step 4
doequi2=1       # 1/0: Do/Do not bring biomass stocks into quasi-equilibrium with unrestricted P and N pools
nequi2=2        #      number of times to repeat steps in doequi2
# Step 5a
doiniluc=1      # 1/0: Do/Do not spinup with dynamic land use (5a)
# Step 5b
doinidyn=1      # 1/0: Do/Do not full dynamic spinup from 1700 to 1899 (5b)
# Step 6
dofinal=1       # 1/0: Do/Do not final run from 1900 to 2017
# Step 7
dofuture=0      # 1/0: Do/Do not future runs (plume only)

# --------------------------------------------------------------------
# Other switches
#

# MetType
mettype="cru"       # "cru", "plume", "bios"
degrees=0.05        #  bios only: resolution of met and LUC (0.05 or 0.25)
metmodel="hadgem2"  # "hadgem2", "ipsl" (only used if mettype="plume")
RCP="hist"          # "hist", "2.6", "4.5", "6.0", "8.5" (no future runs if RCP="hist")

# Cable
# Accounting for mesophyll conductance
# 1: explicit (finite)
# 0: implicit
explicit_gm=1
# Parameter conversion accounting for gm (only used if explicit_gm=1)
# 1: Use lookup table for parameter conversion accounting for gm
# 0: Do not use lookup table
use_LUTgm=1
# Parameters for RubisCO temperature function
# "Bernacchi_2002" or "Walker_2013"
Rubisco_params="Bernacchi_2002"
# Coordination of photosynthesis
# 1: Coordinate photosynthetic rates
# 0: Do not coordinate
coordinate_photosyn=1
# Coordination type (only used if coordinate_photosyn=1)
# T: force coordinated RubisCO- and energy-limited rates
# F: optimise photosynthesis
coord=F
# Acclimation of photosynthesis
# 1: Acclimate photosynthesis
# 0: Do not acclimate
acclimate_photosyn=1
# POP
# 1: Use POP, population dynamics model coupled to CASA
# 0: Do not use POP
call_pop=1
# Soil and snow model
# "default" or "sli"
soil_struc="default"
# 13C
# 1: Calculate 13C
# 0: no 13C
doc13o2=0
# 13C leaf discrimination
# 1: use simple 13C leaf discrimination (a+(b-a)*Ci/Ca)
# 0: calculate full 13C leaf discrimination
c13o2_simple_disc=0

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
    # export mpiexecdir=/usr/local/openmpi-4.1.1-gfortran/bin
    export mpiexecdir=/usr/local/openmpi-4.1.1-ifort/bin
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
        export PYTHONPATH=${PYTHONPATH}:/g/data/x45/python/lib/python3.7/site-packages
    fi
    if [[ ${randompoints} -eq 0 ]] ; then module load nco/4.9.2 ; fi  # needed for cropping outputs
    export mpiexecdir=/apps/intel-mpi/2019.5.281/intel64/bin
elif [[ "${sys}" == "aurora" ]] ; then
    pdir=${isdir}
    module purge
    module load intel/2020a
    module load netCDF-Fortran/4.5.2
    if [[ ${doextractsite} -ge 1 ]] ; then
        module load GCCcore/10.3.0 Python/3.9.5
        export PYTHONPATH=${PYTHONPATH}  #CLN not copied!#:/g/data/x45/python/lib/python3.7/site-packages
	      echo 'Extraction not possible at the moment. Set Py-Path!'
	      exit -666
    fi
    export mpiexecdir=""
fi
if [[ ! -z ${mpiexecdir} ]] ; then export mpiexecdir="${mpiexecdir}/" ; fi

#
# Directories of things
#

# Naming of experiments for gm_acclim_coord experiments of Jurgen
if [[ ("${system}" == "jk8585@gadi") || ("${system}" == "knauer@gadi") ]] ; then
    experiment=test_${mettype}
    if [[ "${mettype}" == "plume" ]] ; then
        # (TEST_)(exp|imp)_(no)coord(T|F)(_PSacclim)_(B2002|W2013)
        if [[ ${explicit_gm} -eq 1 ]] ; then
            if [[ ${coordinate_photosyn} -eq 1 ]] ; then
                if [[ ${acclimate_photosyn} -eq 1 ]] ; then
                    experiment=exp_coord${coord}_PSacclim
                else
	                  experiment=exp_coord${coord}
                fi
            else
	              if [[ ${acclimate_photosyn} -eq 1 ]] ; then
	                  experiment=exp_nocoord_PSacclim
	              else
	                  experiment=exp_nocoord
	              fi
	          fi
        else
	          if [[ ${coordinate_photosyn} -eq 1 ]] ; then
	              if [[ ${acclimate_photosyn} -eq 1 ]] ; then
		                experiment=imp_coord${coord}_PSacclim
	              else
		                experiment=imp_coord${coord}
	              fi
	          else
	              if [[ ${acclimate_photosyn} -eq 1 ]] ; then
		                experiment=imp_nocoord_PSacclim
	              else
		                experiment=imp_nocoord
	              fi
            fi
        fi

        if [[ "${Rubisco_params}" == "Bernacchi_2002" ]] ; then
            experiment=${experiment}_B2002
        elif [[ "${Rubisco_params}" == "Walker_2013" ]] ; then
            experiment=${experiment}_W2013
        else
	          echo "unknown option for 'Rubisco_params', exiting..."
	          exit 1
        fi

        if [[ ${dompi} -ne 0 ]] ; then
	          experiment=TEST_${experiment}
        fi
    fi
fi

#
# Relative directories must be relative to the directory of this script,
#   not relative to the directory from which this script is launched (if different)
#   nor relative to the run path.
#
if [[ "${system}" == "cuntz@explor" ]] ; then
    # Run directory: runpath="${sitepath}/run_xxx"
    cablebase="/home/oqx29/zzy20/prog/cable"
    sitepath="${cablebase}/runs/single_sites/${experiment}"
    cablehome="${cablebase}/branches/CABLE-POP_TRENDY"
    workpath=${cablehome}
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
    MetVersion="CRUJRA_2018"
    # Global LUC
    GlobalTransitionFilePath="/home/oqx29/zzy20/data/LUH2_v3_1deg/"
elif [[ "${system}" == "cuntz@mc16" || "${system}" == "cuntz@mcinra" ]] ; then
    # Run directory: runpath="${sitepath}/run_xxx"
    cablebase="/Users/cuntz/prog/vanessa/cable"
    sitepath="${cablebase}/runs/single_sites/${experiment}"
    cablehome="${cablebase}/branches/CABLE-POP_TRENDY"
    workpath=${cablehome}
    # Cable executable
    if [[ ${dompi} -eq 1 ]] ; then
        # exe="${cablehome}/offline/cable-mpi-gfortran"
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
    MetVersion="CRUJRA_2018"
    GlobalTransitionFilePath=
elif [[ "${system}" == "moc801@gadi" || "${system}" == "cuntz@gadi" ]] ; then
    # Run directory: runpath="${sitepath}/run_xxx"
    # sitepath="/home/801/moc801/prog/cable/runs/single_sites/${experiment}"
    sitepath="/scratch/x45/moc801/cable/c13"
    cablehome="/home/801/moc801/prog/cable/branches/CABLE-POP_TRENDY"
    workpath=${cablehome}
    # Cable executable
    if [[ ${dompi} -eq 1 ]] ; then
        exe="${cablehome}/offline/cable-mpi"
    else
        exe="${cablehome}/offline/cable"
    fi
    # CABLE-AUX directory (uses offline/gridinfo_CSIRO_1x1.nc and offline/modis_phenology_csiro.txt)
    aux="/g/data/x45/CABLE-AUX"
    # BLAZE directory
    BlazeDataPath="/g/data/x45/Data_BLAZE"
    # Global MET, MASK, LUH
    if [[ "${mettype}" == "cru" ]] ; then
        # Global MASK
        SurfaceFile="${aux}/offline/gridinfo_CSIRO_1x1.nc"   # note that SurfaceFile does not need subsetting
        # Global MET
        # GlobalLandMaskFile="${aux}/offline/gridinfo_CSIRO_1x1.nc"
        GlobalLandMaskFile="/g/data/x45/ipbes/masks/glob_ipsl_1x1.nc"
        MetVersion="CRUJRA_2021"
        GlobalMetPath="/g/data/x45/CRUJRA2021/daily_1deg"
        # Global LUC
        GlobalTransitionFilePath="/g/data/x45/LUH2/GCB_2019/1deg/EXTRACT"
    elif [[ "${mettype}" == "plume" ]] ; then
        # Global MASK
        SurfaceFile="${aux}/offline/gridinfo_CSIRO_1x1.nc"   # note that SurfaceFile does not need subsetting
        # Global MET
	      # GlobalLandMaskFile="/g/data/x45/ipbes/masks/glob_ipsl_1x1.nc"
	      GlobalLandMaskFile="/g/data/x45/ipbes/masks/gridinfo_CSIRO_1x1.nc"
        MetVersion=
	      GlobalMetPath="/g/data/x45/ipbes/${metmodel}/1deg"
        # only in plume.nml
        CO2Path="/g/data/x45/ipbes/co2"
	      NdepPath="/g/data/x45/ipbes/ndep"
        # Global LUC
        GlobalTransitionFilePath="/g/data/x45/LUH2/GCB_2019/1deg/EXTRACT"
    elif [[ "${mettype}" == "bios" ]] ; then
        # Global MASK
        SurfaceFile="${aux}/offline/gridinfo_CSIRO_CRU05x05_4tiles.nc"   # note that SurfaceFile does not need subsetting
        # Global MET
        GlobalLandMaskFile="/g/data/x45/BIOS3_forcing/acttest9/acttest9"  # no file extension
        MetVersion=
        GlobalMetPath="/g/data/x45/BIOS3_forcing/acttest9/met/"  # last slash is needed
	      ParamPath="/g/data/x45/BIOS3_forcing/acttest9/params/"  # only in bios.nml
        # Global LUC
        GlobalTransitionFilePath="/g/data/x45/LUH2/v3h/${degrees}deg_aust/EXTRACT"
    fi
elif [[ "${system}" == "kna016@pearcey" || "${system}" == "knauer@pearcey" ]] ; then
    # Run directory: runpath="${sitepath}/run_xxx"
    sitepath="/OSM/CBR/OA_GLOBALCABLE/work/Juergen/CABLE_run/parallel_runs/${experiment}"
    cablehome="/OSM/CBR/OA_GLOBALCABLE/work/Juergen/CABLE_code/CABLE-POP_TRENDY"
    workpath=${cablehome}
    # Cable executable
    if [[ ${dompi} -eq 1 ]] ; then
        exe="${cablehome}/offline/cable-mpi"
    else
        exe="${cablehome}/offline/cable"
    fi
    # CABLE-AUX directory (uses offline/gridinfo_CSIRO_1x1.nc and offline/modis_phenology_csiro.txt)
    aux="/OSM/CBR/OA_GLOBALCABLE/work/Vanessa/CABLE-AUX"
    # Global Mask
    GlobalLandMaskFile="${aux}/offline/gridinfo_CSIRO_1x1.nc"
    SurfaceFile="${aux}/offline/gridinfo_CSIRO_1x1.nc"   # note that SurfaceFile does not need subsetting
    # Global CRU
    GlobalMetPath="/OSM/CBR/OA_GLOBALCABLE/work/CRU-JRA55/crujra/daily_1deg"
    MetVersion="CRUJRA_2021"
    # Global LUC
    GlobalTransitionFilePath="/OSM/CBR/OA_GLOBALCABLE/work/LUH2/v3/1deg"
elif [[ "${system}" == "jk8585@gadi" || "${system}" == "knauer@gadi" ]] ; then
    # Run directory: runpath="${sitepath}/run"
    sitepath="/g/data/x45/jk8585/global_runs/gm_acclim_coord/${experiment}"
    cablecode="/home/599/jk8585/CABLE_code/gm_acclim_coord"
    cablehome="/home/599/jk8585/CABLE_run/gm_acclim_coord"
    workpath=${cablehome}
    # Cable executable
    if [[ ${dompi} -eq 1 ]] ; then
        exe="${cablecode}/offline/cable-mpi"
    else
        exe="${cablecode}/offline/cable"
    fi
    # CABLE-AUX directory (uses offline/gridinfo_CSIRO_1x1.nc and offline/modis_phenology_csiro.txt)
    aux="/g/data/x45/CABLE-AUX"
    # Global Mask
    SurfaceFile="${aux}/offline/gridinfo_CSIRO_1x1.nc"   # note that SurfaceFile does not need subsetting
    # Global Met
    if [[ "${mettype}" == "cru" ]] ; then
	      GlobalLandMaskFile="/g/data/x45/ipbes/masks/glob_ipsl_1x1.nc"
        GlobalMetPath="/g/data/x45/CRUJRA2020/daily_1deg"
        MetVersion="CRUJRA_2021"
    elif [[ "${mettype}" == "plume" ]] ; then
	      # GlobalLandMaskFile="/g/data/x45/ipbes/masks/glob_ipsl_1x1.nc"
	      GlobalLandMaskFile="/g/data/x45/ipbes/masks/gridinfo_CSIRO_1x1.nc"
	      GlobalMetPath="/g/data/x45/ipbes/${metmodel}/1deg"
        MetVersion=
        # only in plume.nml
        CO2Path="/g/data/x45/ipbes/co2"
	      NdepPath="/g/data/x45/ipbes/ndep"
    elif [[ "${mettype}" == "bios" ]] ; then
	      GlobalLandMaskFile="/g/data/x45/BIOS3_forcing/acttest9/acttest9" # no file extension
        GlobalMetPath="/g/data/x45/BIOS3_forcing/acttest9/met/"  # last slash is needed
        MetVersion=
        # only in bios.nml
	      ParamPath="/g/data/x45/BIOS3_forcing/acttest9/params/"
    fi
    # Global LUC
    GlobalTransitionFilePath="/g/data/x45/LUH2/GCB_2019/1deg/EXTRACT"
elif [[ "${system}" == "yc3714@gadi" || "${system}" == "villalobos@gadi" ]] ; then
    # Run directory: runpath="${sitepath}/run"
    sitepath="/g/data/x45/BIOS3_output/${experiment}" # Results 
    cablehome="/home/563/yc3714/CSIRO/CABLE_code/9011/NESP2pt9_BLAZE" # model home
    workpath="/home/563/yc3714/CSIRO/CABLE_run/BIOS3_blaze" # run directory
    # Cable executable
    if [[ ${dompi} -eq 1 ]] ; then
        exe="${cablehome}/offline/cable-mpi"
    else
        exe="${cablehome}/offline/cable"
    fi
    # CABLE-AUX directory (uses offline/gridinfo_CSIRO_1x1.nc and offline/modis_phenology_csiro.txt)
    aux="/g/data/x45/CABLE-AUX"
    # Global Mask
    SurfaceFile="${aux}/offline/gridinfo_CSIRO_CRU05x05_4tiles.nc"   # note that SurfaceFile does not need subsetting
    # Global Met
    if [[ "${mettype}" == "cru" ]] ; then
	      GlobalLandMaskFile="/g/data/x45/ipbes/masks/glob_ipsl_1x1.nc"
        MetVersion="CRUJRA_2021"
        GlobalMetPath="/g/data/x45/CRUJRA2021/daily_1deg"
    elif [[ "${mettype}" == "plume" ]] ; then
	      GlobalLandMaskFile="/g/data/x45/ipbes/masks/gridinfo_CSIRO_1x1.nc"
        MetVersion=
	      GlobalMetPath="/g/data/x45/ipbes/${metmodel}/1deg"
        # only in plume.nml
        CO2Path="/g/data/x45/ipbes/co2"
	      NdepPath="/g/data/x45/ipbes/ndep"
    elif [[ "${mettype}" == "bios" ]] ; then
    	  GlobalLandMaskFile="/g/data/x45/BIOS3_forcing/acttest9/acttest9" # no file extension
        MetVersion=
        GlobalMetPath="/g/data/x45/BIOS3_forcing/acttest9/met/"          # last slash is needed
	      ParamPath="/g/data/x45/BIOS3_forcing/acttest9/params/"           # only in bios.nml
        GlobalTransitionFilePath="/g/data/x45/LUH2/v3h/${degrees}deg_aust/EXTRACT"
    fi
    # Global LUC
    # GlobalTransitionFilePath="/g/data/x45/LUH2/GCB_2019/1deg/EXTRACT"
elif [[ "${system}" == "bri220@pearcey" || "${system}" == "briggs@pearcey" ]] ; then
    # Run directory: runpath="${sitepath}/run_xxx"
    sitepath="/OSM/CBR/OA_GLOBALCABLE/work/Peter/CABLE_run/parallel_runs/${experiment}"
    cablehome="/OSM/CBR/OA_GLOBALCABLE/work/Peter/CABLE_code"
    workpath=${cablehome}
    # Cable executable
    if [[ ${dompi} -eq 1 ]] ; then
        exe="${cablehome}/NESP2pt9_BLAZE/offline/cable-mpi"
    else
        exe="${cablehome}/NESP2pt9_BLAZE/offline/cable"
    fi
    # CABLE-AUX directory (uses offline/gridinfo_CSIRO_1x1.nc and offline/modis_phenology_csiro.txt)
    aux="/OSM/CBR/OA_GLOBALCABLE/work/Vanessa/CABLE-AUX"
    # Global Mask
    GlobalLandMaskFile="${aux}/offline/gridinfo_CSIRO_1x1.nc"
    # Global CRU
    MetVersion="CRUJRA_2021"
    GlobalMetPath="/OSM/CBR/OA_GLOBALCABLE/work/CRU-JRA55/crujra/daily_1deg"
    # Global LUC
    GlobalTransitionFilePath="/OSM/CBR/OA_GLOBALCABLE/work/LUH2/v3/1deg"
elif [[ "${system}" == "pcb599@gadi" || "${system}" == "briggs@gadi" ]] ; then
    # Run directory: runpath="${sitepath}/run_xxx"
    sitepath="/g/data/x45/BIOS3_output/${experiment}"
    cablehome="/home/599/pcb599/CABLE_code/NESP2pt9_BLAZE"
    workpath="/home/599/pcb599/CABLE_run/BIOS3/run"
    # Cable executable
    if [[ ${dompi} -eq 1 ]] ; then
	      exe="${cablehome}/offline/cable-mpi"
    else
	      exe="${cablehome}/offline/cable"
    fi
    # CABLE-AUX directory (uses offline/gridinfo_CSIRO_1x1.nc and offline/modis_phenology_csiro.txt)
    aux="/g/data/x45/CABLE-AUX"
    # Global Mask
    GlobalLandMaskFile="${aux}/offline/gridinfo_CSIRO_1x1.nc"
    SurfaceFile="${aux}/offline/gridinfo_CSIRO_1x1.nc"   # note that SurfaceFile does not need subsetting
    # Global LUC
    GlobalTransitionFilePath="/g/data/x45/LUH2/GCB_2019/1deg/EXTRACT"
    # Global Met
    if [[ "${mettype}" == 'plume' ]] ; then
	      GlobalMetPath="/g/data/x45/ipbes/${metmodel}/1deg"
        CO2Path="/g/data/x45/ipbes/co2"
	      NdepPath="/g/data/x45/ipbes/ndep"
        MetVersion=
	      GlobalLandMaskFile="/g/data/x45/ipbes/masks/glob_ipsl_1x1.nc"
	      GlobalLandMaskFile="/g/data/x45/ipbes/masks/gridinfo_CSIRO_1x1.nc"
    elif [[ "${mettype}" == 'bios' ]] ; then
	      GlobalLandMaskFile="/g/data/x45/BIOS3_forcing/acttest9/acttest9" # no file extension
        GlobalMetPath="/g/data/x45/BIOS3_forcing/acttest9/met/"  # last slash is needed
        MetVersion=
	      ParamPath="/g/data/x45/BIOS3_forcing/acttest9/params/"
        GlobalTransitionFilePath="/g/data/x45/LUH2/v3h/${degrees}deg_aust/EXTRACT"
    elif [[ "${mettype}" == 'cru' ]] ; then
	      GlobalLandMaskFile="/g/data/x45/ipbes/masks/glob_ipsl_1x1.nc"
        MetVersion="CRUJRA_2021"
        GlobalMetPath="/g/data/x45/CRUJRA201/daily_1deg"
    fi
elif [[ "${system}" == "nieradzik@aurora" ]] ; then
    # Run directory: runpath="${sitepath}/run_xxx"
    sitepath="/home/x_larni/STOREDIR/RUNDIR/CABLE/${experiment}"
    cablehome="/home/x_larni/SRC/CABLE/NESP2pt9_BLAZE"
    workpath=${cablehome}
    # Cable executable
    if [[ ${dompi} -eq 1 ]] ; then
	      exe="${cablehome}/offline/cable-mpi"
    else
	      exe="${cablehome}/offline/cable"
    fi
    # CABLE-AUX directory (uses offline/gridinfo_CSIRO_1x1.nc and offline/modis_phenology_csiro.txt)
    aux="/home/x_larni/STOREDIR/DATA/CABLE-AUX"
    # Global Mask
    GlobalLandMaskFile="${aux}/offline/gridinfo_CSIRO_1x1.nc"
    SurfaceFile="${aux}/offline/gridinfo_CSIRO_1x1.nc"   # note that SurfaceFile does not need subsetting
    # Global LUC
    GlobalTransitionFilePath="/g/data/x45/LUH2/GCB_2019/1deg/EXTRACT"
    # Global Met
    if [[ "${mettype}" == 'plume' ]] ; then
        MetVersion=
	      GlobalMetPath="/g/data/x45/ipbes/${metmodel}/1deg"
        CO2Path="/g/data/x45/ipbes/co2"
	      NdepPath="/g/data/x45/ipbes/ndep"
	      GlobalLandMaskFile="/g/data/x45/ipbes/masks/glob_ipsl_1x1.nc"
	      GlobalLandMaskFile="/g/data/x45/ipbes/masks/gridinfo_CSIRO_1x1.nc"
    elif [[ "${mettype}" == 'bios' ]] ; then
	      GlobalLandMaskFile="/home/x_larni/STOREDIR/DATA/CABLE_INPUT/acttest9/acttest9" # no file extension
        MetVersion=
        GlobalMetPath="/home/x_larni/STOREDIR/DATA/CABLE_INPUT/acttest9/met/"  # last slash is needed
	      ParamPath="/home/x_larni/STOREDIR/DATA/CABLE_INPUT/acttest9/params/"
        GlobalTransitionFilePath="/home/x_larni/STOREDIR/DATA/CABLE_INPUT/LUH2/v3h/${degrees}deg_aust/EXTRACT"
    elif [[ "${mettype}" == 'cru' ]] ; then
	      GlobalLandMaskFile="/g/data/x45/ipbes/masks/glob_ipsl_1x1.nc"
        MetVersion="CRUJRA_2021"
        GlobalMetPath="/g/data/x45/CRUJRA2021/daily_1deg"
    fi
else
    echo "System not known."
    exit 1
fi
# Run directory
# runpath="${sitepath}/run_20220224"
runpath="${sitepath}/run"

# Cable parameters
if [[ "${mettype}" == 'bios' ]] ; then
    namelistpath="${workpath}/namelists_bios"
    filename_veg="${workpath}/params_bios/def_veg_params.txt"
    filename_soil="${workpath}/params_bios/def_soil_params.txt"
    casafile_cnpbiome="${workpath}/params_bios/pftlookup.csv"
else
    namelistpath="${workpath}/namelists"
    filename_veg="${workpath}/params/def_veg_params.txt"
    filename_soil="${workpath}/params/def_soil_params.txt"
    casafile_cnpbiome="${workpath}/params/pftlookup.csv"
fi
# Other scripts
ScriptsPath="${cablehome}/scripts"
# Mask (not used unless imeteo > 0)
LandMaskFile="${sitepath}/mask/${experiment}_landmask.nc"
# Met (not used unless imeteo > 0)
if [[ "${mettype}" == 'bios' ]] ; then
    MetPath="${sitepath}/met/bios_${degrees}deg"
    ClimateFile="${sitepath}/mask/bios_climate_rst.nc"
    if [[ "${sys}" == "gadi" ]]; then
        if [[ (${doclimate} -eq 0) && (! -f ${ClimateFile}) ]] ; then
            ClimateFile="/g/data/x45/BIOS3_output/bio_climate_acttest9/bios_climate_rst.nc"
        fi
    fi
else
    MetPath="${sitepath}/met/cru_jra_1deg"
    ClimateFile="${sitepath}/mask/cru_climate_rst.nc"
    if [[ "${sys}" == "gadi" ]]; then
        if [[ (${doclimate} -eq 0) && (! -f ${ClimateFile}) ]] ; then
            ClimateFile="/g/data/x45/ipbes/cable_climate/ipsl_climate_rst_glob_1deg.nc"
        fi
    fi
fi
# LUC
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
printf "        soil_struc=${soil_struc}\n"
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


# absolute pathes of other parameter files
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
    # rm -f ${ClimateFile}
fi


# Write standard namelists with options that are common to all steps of the sequence.
# They can, however, be overwritten in later steps.

# global meteo namelist file
if [[ "${mettype}" == "cru" ]] ; then
    cat > ${tmp}/sedtmp.${pid} << EOF
        BasePath     = "${MetPath}"
        MetPath      = "${MetPath}"
        MetVersion   = "${MetVersion}"
        LandMaskFile = "${LandMaskFile}"
        Run          = "S0_TRENDY"
EOF
    applysed ${tmp}/sedtmp.${pid} ${ndir}/cru.nml ${rdir}/cru_${experiment}.nml
elif [[ "${mettype}" == "plume" ]] ; then
    #TODO: NDEPfile, CO2file
    cat > ${tmp}/sedtmp.${pid} << EOF
        BasePath     = "${MetPath}"
        Forcing      = "hadgem2-es"
    	  LandMaskFile = "${LandMaskFile}"
    	  Run          = "spinup"
    	  RCP          = "hist"
    	  CO2          = "static1850"
    	  CO2file      = "${CO2Path}/co2_1850_2005_hist.dat"
    	  NDEP         = "static1850"
        NDEPfile     = "${NdepPath}/NOy_plus_NHx_dry_plus_wet_deposition_hist_1850_2015_annual_1deg.nc"
EOF
    if [[ "${metmodel}" == "ipsl" ]] ; then
        sed -i -e "/Forcing/s/=.*/= 'ipsl-cm5a-lr'/" ${tmp}/sedtmp.${pid}
    fi
    applysed ${tmp}/sedtmp.${pid} ${ndir}/plume.nml ${rdir}/plume_${experiment}.nml
elif [[ "${mettype}" == "bios" ]] ; then
    cat > ${tmp}/sedtmp.${pid} << EOF
        met_path         = "${MetPath}/"
        param_path       = "${ParamPath}"
        landmaskflt_file = "${GlobalLandMaskFile}.flt"
        landmaskhdr_file = "${GlobalLandMaskFile}.hdr"
        rain_file        = "1900010120201231_rain_b2003.bin"
        swdown_file      = "1900010120201231_rad_b2003.bin"
        tairmax_file     = "1900010120201231_tmax_noclim_b2003.bin"
        tairmin_file     = "1900010120201231_tmin_noclim_b2003.bin"
        wind_file        = "1900010120201231_windspeed_ms_b2003.bin"
        vp0900_file      = "1900010120201231_vph09_b2003.bin"
        vp1500_file      = "1900010120201231_vph15_b2003.bin"
        co2_file         = "1750_2020_globalCO2_time_series.bin"
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

# Blaze namelist !CLN CHECK
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

# casafile%cnpipool                  = "restart/${mettype}_casa_rst.nc"
# casafile%cnpepool                  = "restart/${mettype}_casa_rst.nc"
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
    leaps                              = .false.
    cable_user%SOIL_STRUC              = "${soil_struc}"
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
    if [[ "${mettype}" == "cru" ]] ; then
        cp ${rdir}/cru_${experiment}.nml ${rdir}/cru.nml
    elif [[ "${mettype}" == "plume" ]] ; then
        cp ${rdir}/plume_${experiment}.nml ${rdir}/plume.nml
    elif [[ "${mettype}" == "bios" ]] ; then
        cat > ${tmp}/sedtmp.${pid} << EOF
	          Run = "spinup"
EOF
        applysed ${tmp}/sedtmp.${pid} ${rdir}/bios_${experiment}.nml ${rdir}/bios.nml
    fi
    # LUC
    cp ${rdir}/LUC_${experiment}.nml ${rdir}/LUC.nml
    # Cable
    #   do not calculate 13C because there is no 13C in the climate restart file
    #MCTEST
    # cable_user%YearEnd             = 1861
    # cable_user%CASA_SPIN_ENDYEAR   = 1861
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
    if [[ "${mettype}" == "cru" ]] ; then
        cp ${rdir}/cru_${experiment}.nml ${rdir}/cru.nml
    elif [[ "${mettype}" == "plume" ]] ; then
        cp ${rdir}/plume_${experiment}.nml ${rdir}/plume.nml
    elif [[ "${mettype}" == "bios" ]] ; then
        cat > ${tmp}/sedtmp.${pid} << EOF
	          Run = "spinup"
EOF
        applysed ${tmp}/sedtmp.${pid} ${rdir}/bios_${experiment}.nml ${rdir}/bios.nml
    fi
    # LUC
    cp ${rdir}/LUC_${experiment}.nml ${rdir}/LUC.nml
    # Cable
    #MCTEST
    # cable_user%YearEnd                = 1861
    # cable_user%CASA_SPIN_ENDYEAR      = 1861
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
        output%averaging                  = "all"
        cable_user%CASA_SPIN_STARTYEAR    = 1860
        cable_user%CASA_SPIN_ENDYEAR      = 1869
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
    # copyid ${rid} ${mettype}_climate_rst.nc ${mettype}_casa_rst.nc ${mettype}_cable_rst.nc pop_${mettype}_ini.nc
    copyid ${rid} ${mettype}_climate_rst.nc ${mettype}_casa_met_rst.nc ${mettype}_casa_bal_rst.nc ${mettype}_casa_biome_rst.nc ${mettype}_casa_pool_rst.nc ${mettype}_casa_flux_rst.nc ${mettype}_casa_phen_rst.nc ${mettype}_cable_rst.nc pop_${mettype}_ini.nc
    copyid ${rid} ${mettype}_c13o2_flux_rst.nc ${mettype}_c13o2_pools_rst.nc
    cd ../outputs
    renameid ${rid} ${mettype}_out_cable.nc ${mettype}_out_casa.nc ${mettype}_out_casa_c13o2.nc
    cd ..
    cd ${pdir}
fi


# --------------------------------------------------------------------
# 3. Biomass into quasi-equilibrium with unrestricted N and P pools
if [[ ${doequi1} -eq 1 ]] ; then
    echo "3. Bring biomass into quasi-equilibrium with unrestricted N and P pools"
    for ((iequi1=1; iequi1<=${nequi1}; iequi1++)) ; do
        # 3a. 30 year run starting from restart files
        echo "   3a. 30 year spinup from accumulated biomass; iequi1=${iequi1}/${nequi1}"
        rid="spinup_limit_labile"
        # rid="spinup_limit_labile${iequi}"
        # Met forcing
        if [[ "${mettype}" == "cru" ]] ; then
            cp ${rdir}/cru_${experiment}.nml ${rdir}/cru.nml
        elif [[ "${mettype}" == "plume" ]] ; then
            cp ${rdir}/plume_${experiment}.nml ${rdir}/plume.nml
        elif [[ "${mettype}" == "bios" ]] ; then
            cat > ${tmp}/sedtmp.${pid} << EOF
    	          Run = "spinup"
EOF
            applysed ${tmp}/sedtmp.${pid} ${rdir}/bios_${experiment}.nml ${rdir}/bios.nml
        fi
        # LUC
        cp ${rdir}/LUC_${experiment}.nml ${rdir}/LUC.nml
        # Cable
        #MCTEST
        # cable_user%YearEnd             = 1841
        # cable_user%CASA_SPIN_ENDYEAR   = 1861
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
            cable_user%CASA_SPIN_STARTYEAR = 1860
            cable_user%CASA_SPIN_ENDYEAR   = 1869
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
        # copyid ${rid} ${mettype}_climate_rst.nc ${mettype}_casa_rst.nc ${mettype}_cable_rst.nc pop_${mettype}_ini.nc
        copyid ${rid} ${mettype}_climate_rst.nc ${mettype}_casa_met_rst.nc ${mettype}_casa_bal_rst.nc ${mettype}_casa_biome_rst.nc ${mettype}_casa_pool_rst.nc ${mettype}_casa_flux_rst.nc ${mettype}_casa_phen_rst.nc ${mettype}_cable_rst.nc pop_${mettype}_ini.nc
        copyid ${rid} ${mettype}_c13o2_flux_rst.nc ${mettype}_c13o2_pools_rst.nc
        cd ../outputs
        renameid ${rid} ${mettype}_out_cable.nc ${mettype}_out_casa.nc ${mettype}_out_casa_c13o2.nc
        cd ..
        cd ${pdir}
        #
        # 3b. analytic quasi-equilibrium of biomass pools
        echo "   3b. Analytic solution of biomass pools"
        rid="spinup_analytic_limit_labile"
        # rid="spin_casa_limit_labile${iequi}"
        # Met forcing
        if [[ "${mettype}" == "cru" ]] ; then
            cp ${rdir}/cru_${experiment}.nml ${rdir}/cru.nml
        elif [[ "${mettype}" == "plume" ]] ; then
            cp ${rdir}/plume_${experiment}.nml ${rdir}/plume.nml
        elif [[ "${mettype}" == "bios" ]] ; then
            cat > ${tmp}/sedtmp.${pid} << EOF
    	          Run = "spinup"
EOF
            applysed ${tmp}/sedtmp.${pid} ${rdir}/bios_${experiment}.nml ${rdir}/bios.nml
        fi
        # LUC
        cp ${rdir}/LUC_${experiment}.nml ${rdir}/LUC.nml
        # Cable
        #MCTEST
        # cable_user%YearEnd             = 1841
        # cable_user%CASA_SPIN_ENDYEAR   = 1841
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
        # copyid ${rid} ${mettype}_casa_rst.nc pop_${mettype}_ini.nc
        copyid ${rid} ${mettype}_casa_met_rst.nc ${mettype}_casa_bal_rst.nc ${mettype}_casa_biome_rst.nc ${mettype}_casa_pool_rst.nc ${mettype}_casa_flux_rst.nc ${mettype}_casa_phen_rst.nc pop_${mettype}_ini.nc
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
# 4. Biomass into quasi-equilibrium with restricted N and P pools
if [[ ${doequi2} -eq 1 ]] ; then
    echo "4. Bring biomass into quasi-equilibrium with restricted N and P pools"
    for ((iequi2=1; iequi2<=${nequi2}; iequi2++)) ; do
        # 4a. 30 year run starting from restart files
        echo "   4a. 30 year spinup from accumulated biomass; iequi2=${iequi2}/${nequi2}"
        rid="spinup"
        # rid="spinup${iequi}"
        # Met forcing
        if [[ "${mettype}" == "cru" ]] ; then
            cp ${rdir}/cru_${experiment}.nml ${rdir}/cru.nml
        elif [[ "${mettype}" == "plume" ]] ; then
            cp ${rdir}/plume_${experiment}.nml ${rdir}/plume.nml
        elif [[ "${mettype}" == "bios" ]] ; then
            cat > ${tmp}/sedtmp.${pid} << EOF
    	          Run = "spinup"
EOF
            applysed ${tmp}/sedtmp.${pid} ${rdir}/bios_${experiment}.nml ${rdir}/bios.nml
        fi
        # LUC
        cp ${rdir}/LUC_${experiment}.nml ${rdir}/LUC.nml
        # Cable
        #MCTEST
        # cable_user%YearEnd             = 1841
        # cable_user%CASA_SPIN_ENDYEAR   = 1861
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
            cable_user%CASA_SPIN_STARTYEAR = 1860
            cable_user%CASA_SPIN_ENDYEAR   = 1869
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
        renameid ${rid} ${mettype}.nml LUC.nml cable.nml
        mv *_${rid}.nml restart/
        cd logs
        renameid ${rid} log_cable.txt log_out_cable.txt
        cd ../restart
        # copyid ${rid} ${mettype}_climate_rst.nc ${mettype}_casa_rst.nc ${mettype}_cable_rst.nc pop_${mettype}_ini.nc
        copyid ${rid} ${mettype}_climate_rst.nc ${mettype}_casa_met_rst.nc ${mettype}_casa_bal_rst.nc ${mettype}_casa_biome_rst.nc ${mettype}_casa_pool_rst.nc ${mettype}_casa_flux_rst.nc ${mettype}_casa_phen_rst.nc ${mettype}_cable_rst.nc pop_${mettype}_ini.nc
        copyid ${rid} ${mettype}_c13o2_flux_rst.nc ${mettype}_c13o2_pools_rst.nc
        cd ../outputs
        renameid ${rid} ${mettype}_out_cable.nc ${mettype}_out_casa.nc ${mettype}_out_casa_c13o2.nc
        cd ..
        cd ${pdir}
        #
        # 4b. analytic quasi-equilibrium of biomass pools
        echo "   4b. Analytic solution of biomass pools"
        rid="spinup_analytic"
        # rid="spin_casa${iequi}"
        # Met forcing
        if [[ "${mettype}" == "cru" ]] ; then
            cp ${rdir}/cru_${experiment}.nml ${rdir}/cru.nml
        elif [[ "${mettype}" == "plume" ]] ; then
            cp ${rdir}/plume_${experiment}.nml ${rdir}/plume.nml
        elif [[ "${mettype}" == "bios" ]] ; then
            cat > ${tmp}/sedtmp.${pid} << EOF
    	          Run = "spinup"
EOF
            applysed ${tmp}/sedtmp.${pid} ${rdir}/bios_${experiment}.nml ${rdir}/bios.nml
        fi
        # LUC
        cp ${rdir}/LUC_${experiment}.nml ${rdir}/LUC.nml
        # Cable
        #MCTEST
        # cable_user%YearEnd             = 1841
        # cable_user%CASA_SPIN_ENDYEAR   = 1841
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
        renameid ${rid} ${mettype}.nml LUC.nml cable.nml
        mv *_${rid}.nml restart/
        cd logs
        renameid ${rid} log_cable.txt log_out_cable.txt
        cd ../restart
        # copyid ${rid} ${mettype}_casa_rst.nc pop_${mettype}_ini.nc
        copyid ${rid} ${mettype}_casa_met_rst.nc ${mettype}_casa_bal_rst.nc ${mettype}_casa_biome_rst.nc ${mettype}_casa_pool_rst.nc ${mettype}_casa_flux_rst.nc ${mettype}_casa_phen_rst.nc pop_${mettype}_ini.nc
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
# 5a. First dynamic land use (initialise land use)
if [[ ${doiniluc} -eq 1 ]] ; then
    echo "5a. First dynamic land use (initialise land use)"
    # Met forcing
    if [[ "${mettype}" == "cru" ]] ; then
        YearStart=1580
        YearEnd=1699
        cp ${rdir}/cru_${experiment}.nml ${rdir}/cru.nml
    elif [[ "${mettype}" == "plume" ]] ; then
	      YearStart=1580
        YearEnd=1849
        cp ${rdir}/plume_${experiment}.nml ${rdir}/plume.nml
    elif [[ "${mettype}" == "bios" ]] ; then
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
    #MCTEST
    # cable_user%CASA_SPIN_ENDYEAR    = 1859
    # cable_user%YearEnd              = ${YearEnd}
    # cable_user%CASA_SPIN_ENDYEAR    = 1841
    # cable_user%YearEnd              = $(( ${YearStart} + 1 ))
    #MCTEST
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
    renameid ${rid} ${mettype}.nml LUC.nml cable.nml
    mv *_${rid}.nml restart/
    cd logs
    renameid ${rid} log_cable.txt log_out_cable.txt
    cd ../restart
    # copyid ${rid} ${mettype}_casa_rst.nc ${mettype}_LUC_rst.nc pop_${mettype}_ini.nc
    copyid ${rid} ${mettype}_casa_met_rst.nc ${mettype}_casa_bal_rst.nc ${mettype}_casa_biome_rst.nc ${mettype}_casa_pool_rst.nc ${mettype}_casa_flux_rst.nc ${mettype}_casa_phen_rst.nc ${mettype}_LUC_rst.nc pop_${mettype}_ini.nc
    copyid ${rid} ${mettype}_c13o2_pools_rst.nc ${mettype}_c13o2_luc_rst.nc
    # cd ../outputs
    # renameid ${rid} ${mettype}_out_LUC.nc
    # cd ..
    cd ${pdir}
fi


# --------------------------------------------------------------------
# 5b. Second full dynamic spinup
if [[ ${doinidyn} -eq 1 ]] ; then
    echo "5b. Full dynamic spinup"
    # Met forcing
    if [[ "${mettype}" == "cru" ]] ; then
        YearStart=1700
        YearEnd=1899
        cat > ${tmp}/sedtmp.${pid} << EOF
            Run = "S1_TRENDY"
EOF
        applysed ${tmp}/sedtmp.${pid} ${rdir}/cru_${experiment}.nml ${rdir}/cru.nml
    elif [[ "${mettype}" == "plume" ]] ; then
	      YearStart=1850
        YearEnd=1900
        cat > ${tmp}/sedtmp.${pid} << EOF
    	      Run  = "${YearStart}_${YearEnd}"
    	      CO2  = "varying"
    	      NDEP = "varying"
EOF
        applysed ${tmp}/sedtmp.${pid} ${rdir}/plume_${experiment}.nml ${rdir}/plume.nml
    elif [[ "${mettype}" == "bios" ]] ; then
        YearStart=1700
        YearEnd=1899
        cat > ${tmp}/sedtmp.${pid} << EOF
	          Run = "premet"
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
    #MCTEST
    # cable_user%CASA_SPIN_ENDYEAR   = 1851
    # cable_user%YearEnd             = $(( ${YearStart} + 1 ))
    # cable_user%CASA_SPIN_ENDYEAR = 1859
    # cable_user%YearEnd             = ${YearEnd}
    #MCTEST
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
        cable_user%POPLUC_RunType      = "restart"
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
    # copyid ${rid} ${mettype}_climate_rst.nc ${mettype}_casa_rst.nc ${mettype}_cable_rst.nc
    copyid ${rid} ${mettype}_climate_rst.nc ${mettype}_casa_met_rst.nc ${mettype}_casa_bal_rst.nc ${mettype}_casa_biome_rst.nc ${mettype}_casa_pool_rst.nc ${mettype}_casa_flux_rst.nc ${mettype}_casa_phen_rst.nc ${mettype}_cable_rst.nc
    copyid ${rid} ${mettype}_LUC_rst.nc pop_${mettype}_ini.nc
    copyid ${rid} ${mettype}_c13o2_flux_rst.nc ${mettype}_c13o2_pools_rst.nc ${mettype}_c13o2_luc_rst.nc
    cd ../outputs
    renameid ${rid} ${mettype}_out_cable.nc ${mettype}_out_casa.nc ${mettype}_out_LUC.nc
    renameid ${rid} ${mettype}_out_casa_c13o2.nc
    cd ..
    cd ${pdir}
fi


# --------------------------------------------------------------------
# 6. Final centennial run
if [[ ${dofinal} -eq 1 ]] ; then
    echo "6. Final centennial run"
    # Met forcing
    if [[ "${mettype}" == "cru" ]] ; then
        YearStart=1900
        YearEnd=2017
        cat > ${tmp}/sedtmp.${pid} << EOF
            Run = "S2_TRENDY"
EOF
        applysed ${tmp}/sedtmp.${pid} ${rdir}/cru_${experiment}.nml ${rdir}/cru.nml
    elif [[ "${mettype}" == "plume" ]] ; then
	      YearStart=1901
        YearEnd=2005
        cat > ${tmp}/sedtmp.${pid} << EOF
    	      Run  = "${YearStart}_${YearEnd}"
    	      CO2  = "varying"
    	      NDEP = "varying"
EOF
        applysed ${tmp}/sedtmp.${pid} ${rdir}/plume_${experiment}.nml ${rdir}/plume.nml
    elif [[ "${mettype}" == "bios" ]] ; then
        YearStart=1900
        YearEnd=2019
        cat > ${tmp}/sedtmp.${pid} << EOF
	          Run = "standard"
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
    #MCTEST
    # cable_user%CASA_SPIN_ENDYEAR   = 1851
    # cable_user%YearStart           = $(( ${YearStart} + 1 ))
    # cable_user%YearEnd             = $(( ${YearStart} + 2 ))
    # cable_user%CASA_SPIN_ENDYEAR = 1859
    # cable_user%YearStart           = ${YearStart}
    # cable_user%YearEnd             = ${YearEnd}
    #MCTEST
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
        cable_user%POPLUC_RunType      = "restart"
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
    # copyid ${rid} ${mettype}_climate_rst.nc ${mettype}_casa_rst.nc ${mettype}_cable_rst.nc
    copyid ${rid} ${mettype}_climate_rst.nc ${mettype}_casa_met_rst.nc ${mettype}_casa_bal_rst.nc ${mettype}_casa_biome_rst.nc ${mettype}_casa_pool_rst.nc ${mettype}_casa_flux_rst.nc ${mettype}_casa_phen_rst.nc ${mettype}_cable_rst.nc
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
if [[ (${dofuture} -eq 1) && ("${mettype}" == "plume") && ("${RCP}" != "hist") ]] ; then
    echo "7. Future run"
    rcpnd=$(echo ${RCP} | sed "s|\.||")
    rcpd=$(echo ${RCP} | sed "s|\.|p|")
    # Met forcing
	  YearStart=2006
    YearEnd=2099
    cat > ${tmp}/sedtmp.${pid} << EOF
        Run          = "${YearStart}_${YearEnd}"
        RCP          = "${RCP}"
        CO2          = "varying"
        CO2file      = "${CO2Path}/co2_2006_2100_rcp${rcpnd}.dat"
        NDEP         = "varying"
        NDEPfile     = "${NdepPath}/RCP${rcpnd}/ndep_total_2000-2109_1.0x1.0_FD.nc"
EOF
    applysed ${tmp}/sedtmp.${pid} ${rdir}/plume_${experiment}.nml ${rdir}/plume.nml
    rid=${YearStart}_${YearEnd}_${rcpd}
    # LUC
    cat > ${tmp}/sedtmp.${pid} << EOF
         YearStart = ${YearStart}
         YearEnd   = ${YearEnd}
EOF
    applysed ${tmp}/sedtmp.${pid} ${rdir}/LUC_${experiment}.nml ${rdir}/LUC.nml
    # Cable
    #MCTEST
    # cable_user%CASA_SPIN_ENDYEAR = 1859
    # cable_user%YearEnd             = ${YearEnd}
    #MCTEST
    cat > ${tmp}/sedtmp.${pid} << EOF
        cable_user%CLIMATE_fromZero    = .false.
        cable_user%YearStart           = ${YearStart}
        cable_user%YearEnd             = $(( ${YearStart} + 1 ))
        icycle                         = 2
        spincasa                       = .false.
        cable_user%CASA_fromZero       = .false.
        cable_user%CASA_DUMP_READ      = .false.
        cable_user%CASA_DUMP_WRITE     = .false.
        cable_user%CASA_SPIN_STARTYEAR = 1850
        cable_user%CASA_SPIN_ENDYEAR   = 1851
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
    renameid ${rid} ${mettype}.nml LUC.nml cable.nml
    mv *_${rid}.nml restart/
    cd logs
    renameid ${rid} log_cable.txt log_out_cable.txt
    cd ../restart
    # copyid ${rid} ${mettype}_climate_rst.nc ${mettype}_casa_rst.nc ${mettype}_cable_rst.nc pop_${mettype}_ini.nc
    copyid ${rid} ${mettype}_climate_rst.nc ${mettype}_casa_met_rst.nc ${mettype}_casa_bal_rst.nc ${mettype}_casa_biome_rst.nc ${mettype}_casa_pool_rst.nc ${mettype}_casa_flux_rst.nc ${mettype}_casa_phen_rst.nc ${mettype}_cable_rst.nc pop_${mettype}_ini.nc
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
