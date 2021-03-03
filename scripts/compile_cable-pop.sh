#!/usr/bin/env bash
# --------------------------------------------------------------------
#
# Compile Cable-POP using korn-shell scripts build(_mpi).ksh.
# This script can be submitted to compile with -xHost on the host machine.
#
# --------------------------------------------------------------------
# Explor / Pearcey (launch with: sbatch --ignore-pbs)
# https://slurm.schedmd.com/sbatch.html
# Name - 8 letters and digits
#SBATCH -J ccable
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.out
# Explor partitions (sinfo): std (2x16, parallel), sky (2x16, parallel, AVX512), hf (2x4, serial),
#                            P100 (2x16, GPU), GTX (2x16, GPU), ivy (2x8, parallel), k20 (2x8, GPU)
#SBATCH -p std
# -N Nodes / -n tasks (mpiexec, srun, ...) / -c cpus_per_task (OpenMP, make-j, ...)
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 4
# time (day-hh:mm:ss) / memory (optional, units K,M,G,T)
#SBATCH -t 00:19:59
#SBATCH --mem=4G
# notify: Valid type values are NONE,BEGIN,END,FAIL,REQUEUE,ALL,STAGE_OUT,TIME_LIMIT,TIME_LIMIT_90/80/50,ARRAY_TASKS
#SBATCH --mail-type=FAIL,STAGE_OUT,TIME_LIMIT
#SBATCH --mail-user=matthias.cuntz@inrae.fr
#
# --------------------------------------------------------------------
# Gadi
# https://opus.nci.org.au/display/Help/How+to+submit+a+job
#PBS -N ccable
#PBS -P x45
# express / normal / copyq (2x24, cascadelake)
# expressbw / normalbw (2x14, broadwell) / normalsl (2x16, skylake)- ex-Raijin nodes
#PBS -q normal
#PBS -l walltime=00:19:59
#PBS -l mem=4GB
#PBS -l ncpus=4
# #PBS -l jobfs=1GB
#PBS -l software=netCDF:MPI:Intel:GNU
#PBS -r y
#PBS -l wd
#PBS -j oe
#PBS -S /bin/bash
#PBS -M matthias.cuntz@inrae.fr

# Compile MPI version
dompi=1   # 0: serial code; 1: parallel code
doclean=1 # 0: make; 1: clean existing compile directory first
mflag=4   # Numbers of jobs of make to run simultaneously. Set to "" to not limit number of jobs.
          # Could be 2*${SLURM_CPUS_PER_TASK}+1 in case of Slurm

# --------------------------------------------------------------------

set -e

trap cleanup 1 2 3 6

# cleanup at trap
function cleanup()
{
    rm -f build*.${pid}
    exit 1
}

pid=$$
isdir="${PWD}"
prog=$0
pprog=$(basename ${prog})
pdir=$(dirname ${prog})
tmp=${TMPDIR:-"/tmp"}
t1=$(date +%s)

# --------------------------------------------------------------------

# Check environment
# echo "uname -n: "$(uname -n)
# echo "uname -n 4: "$(uname -n | cut -c 1-4 | tr - _)
# echo "uname -a: "$(uname -a)
# echo "hostname --domain: "$(hostname --domain)

# jobid=""
# jobname=""
# jobdir=${isdir}
# jobuser=""
# if [[ -n ${PBS_ENVIRONMENT} ]] ; then
#     jobid=${PBS_JOBID}
#     jobname=${PBS_JOBNAME}
#     jobdir=${PBS_O_WORKDIR}
#     jobuser=${PBS_USER}
# elif [[ -n ${SLURM_JOB_ID} ]] ; then
#     jobid=${SLURM_JOB_ID}
#     jobname=${SLURM_JOB_NAME}
#     jobdir=${SLURM_SUBMIT_DIR}
#     jobuser=${SLURM_JOB_USER}
#     # mflag=$(( 2 * ${SLURM_CPUS_PER_TASK} + 1 ))
# fi

# Get offline path
if [[ ${dompi} -eq 1 ]] ; then
    build_script="build_mpi.ksh"
else
    build_script="build.ksh"
fi
bfile=$(find . -name ${build_script} | sed -e '/UM\/build.ksh/d')
if [[ -z ${bfile} ]] ; then
    cd ..
    bfile=$(find . -name ${build_script} | sed -e '/UM\/build.ksh/d')
fi  
bdir=$(dirname ${bfile})
bfile=$(basename ${bfile})

# Change number of tasks in build script and compile
cd ${bdir}
sed \
    -e '/OPTFLAG=/s/OPTFLAG=.*/OPTFLAG="-xHost"/' \
    -e '/export MFLAGS/s/MFLAGS=.*/MFLAGS="-j '${mflag}'"/' \
    ${bfile} > ${bfile}.${pid}
if [[ ${doclean} -eq 1 ]] ; then
    bash ${bfile}.${pid} cclean
else
    bash ${bfile}.${pid}
fi
rm ${bfile}.${pid}

# --------------------------------------------------------------------

t2=$(date +%s)
dt=$(( t2 - t1 ))
printf "\n"
if [[ ${dt} -lt 60 ]] ; then
    printf "Compiled in %i seconds.\n" ${dt}
else
    dm=$(echo "(${t2}-${t1})/60." | bc -l)
    printf "Compiled in %.2f minutes.\n" ${dm}
fi

exit 0
