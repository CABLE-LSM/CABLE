#!/bin/ksh

export dosvn=0 # 1/0: do/do not check svn

# so that script can be called by 'bash build.ksh' if no ksh installed
if [ -z ${PS3} ] ; then # https://stackoverflow.com/questions/3327013/how-to-determine-the-current-shell-im-working-on
    eval 'function print(){ printf "$@\n"; }'
fi

known_hosts()
{
    if [ -z ${PS3} ] ; then
        kh=(kh gadi pear mcin vm_o)
    else
        set -A kh gadi pear mcin vm_o
    fi
}

known_domains()
{
    if [ -z ${PS3} ] ; then
        kd=(kd nci.org.au pear local explor)
    else
        set -A kd nci.org.au pear local explor
    fi
}


## gadi.nci.org.au
host_gadi()
{
   if [ -z ${PS3} ] ; then
      . /etc/bashrc
   else
      . /etc/kshrc
   fi
   module purge
   module load intel-compiler/2019.5.281
   module load intel-mpi/2019.5.281
   module load netcdf/4.6.3

   export FC=mpif90
   export NCMOD=${NETCDF_ROOT}"/include/Intel"
   if [[ ${1} == "debug" ]]; then
       # debug
       # export CFLAGS='-O0 -fpp -traceback -g -fp-model precise -ftz -fpe0'
       export CFLAGS="-fpp -O0 -debug extended -traceback -g -check all,noarg_temp_created -warn all -fp-stack-check -nofixed -assume byterecl -fp-model precise -diag-disable=10382 -fpe0" # -fpe-all=0 -no-ftz -ftrapuv"
       export LDFLAGS="-O0"
       OPTFLAG=""
   else
       # release
       # export CFLAGS='-O2 -fpp -fp-model precise'
       export CFLAGS="-fpp -O3 -nofixed -assume byterecl -fp-model precise -ip -diag-disable=10382"
       export LDFLAGS="-O3"
       OPTFLAG="-xCASCADELAKE"
       # OPTFLAG="-xCORE-AVX2 -axSKYLAKE-AVX512,CASCADELAKE" # given in user training: does not work
       # OPTFLAG="-xCASCADELAKE" # or -xCORE-AVX512;                           queues: express / normal
       # OPTFLAG="-xBROADWELL"   # or -xCORE-AVX512;                           queues: expressbw / normalbw
       # OPTFLAG="-xSKYLAKE"     # or -xSKYLAKE-AVX512 depends on performance; queues: normalsl
   fi
   export CFLAGS="${CFLAGS} ${OPTFLAG}"
   export CFLAGS="${CFLAGS} -D__MPI__"
   export CFLAGS="${CFLAGS} -D__CRU2017__"
   export CFLAGS="${CFLAGS} -D__NETCDF3__"
   export LDFLAGS="-L"${NCDIR}" "${LDFLAGS}
   export LD="-lnetcdf -lnetcdff"
   export MFLAGS="-j 8"
   build_build
   cd ../
   build_status
}


## pearcey.hpsc.csiro.au
host_pear()
{
    #    export LD_PRELOAD=/apps/netcdf/4.3.3/lib/libnetcdf.so
    #    export LD_PRELOAD=/apps/openmpi/1.8.4/lib/libopen-rte.so.7:/apps/openmpi/1.8.4/lib/libopen-pal.so.6
    if [ -z ${PS3} ] ; then
        . /apps/modules/Modules/default/init/ksh
    fi

    #   module add netcdf/4.3.3.1 openmpi/1.7.5
    #   module add netcdf/4.3.3.1 openmpi/1.8.8

    module del intel-cc intel-fc
    module add intel-cc/16.0.1.150 intel-fc/16.0.1.150
    module add netcdf/4.3.3.1 openmpi/1.8.8

    export NCDIR=$NETCDF_ROOT'/lib/'
    export NCMOD=$NETCDF_ROOT'/include/'
    export FC='mpifort' #'mpif90'
    export CFLAGS='-O0 -fp-model precise -fpp'
    #   export CFLAGS='-O0 -C'
    # best settings for debugging
    #   export CFLAGS='-O0 -C -g -debug all -traceback -check all,noarg_temp_created, -C  '
    #   export CFLAGS='-O0 '

    #export CFLAGS='-O0 -fp-model precise -g -debug -traceback -fpp '
    #export CFLAGS="${CFLAGS} -D__CRU2018__"
    export CFLAGS="${CFLAGS} -D__CRU2017__"
    export CFLAGS="${CFLAGS} -D__NETCDF3__"

    #   export CFLAGS='-O0 -fp-model precise -g -debug all -traceback -fpe0 '
       #export CFLAGS='  -g -debug -traceback -fp-stack-check -O0 -debug -fpe0 -no-ftz -ftrapuv'

    # best debug flags
    #   export LDFLAGS='-g -L'$NCDIR  #'-L'$NCDIR' -O2'
    export CFLAGS="${CFLAGS} -D__MPI__"
    export LDFLAGS='-O0 -L'$NCDIR''
    export MFLAGS='-j 8'
    export LD='-lnetcdf -lnetcdff'
    build_build
    cd ../
    build_status
}


## MatthiasCuntz@INRAE
host_mcin()
{
    idebug=0
    iintel=1
    ignu=0
    inag=0
    np=$#
    for ((i=0; i<${np}; i++)) ; do
        if [[ "${1}" == "debug" ]] ; then
            idebug=1
            shift 1
        elif [[ "${1}" == "ifort" || "${1}" == "intel" ]] ; then
            iintel=1
            ignu=0
            inag=0
            shift 1
        elif [[ "${1}" == "gfortran" || "${1}" == "gnu" ]] ; then
            iintel=0
            ignu=1
            inag=0
            shift 1
        elif [[ "${1}" == "nagfor" || "${1}" == "nag" ]] ; then
            iintel=0
            ignu=0
            inag=1
            shift 1
        else
            echo "Error: command line option not known: " ${1}
            exit 1
        fi
    done
    if [[ ${iintel} -eq 1 ]] ;  then
        # INTEL
        /opt/intel/compilers_and_libraries/mac/bin/compilervars.sh intel64
        export FC=/usr/local/openmpi-3.1.5-ifort/bin/mpifort
        # release
        export CFLAGS="-fpp -O3 -nofixed -assume byterecl -fp-model precise -ip -diag-disable=10382"
        export LDFLAGS="-O3"
	OPTFLAG="-xHost"
        if [[ ${idebug} -eq 1 ]] ; then
            # debug
            export CFLAGS="-fpp -O0 -debug extended -traceback -g -check all,noarg_temp_created -warn all -fp-stack-check -nofixed -assume byterecl -fp-model precise -diag-disable=10382 -fpe0" # -fpe-all=0 -no-ftz -ftrapuv -init=arrays,snan
            export LDFLAGS="-O0"
	    OPTFLAG=
        fi
        export CFLAGS="${CFLAGS} -D__INTEL__ -D__INTEL_COMPILER__"
        export LD=""
        export NCROOT="/usr/local/netcdf-fortran-4.4.5-ifort"
        export cdir=".mpitmp-ifort"
        export PROG=cable-mpi-ifort
    elif [[ ${ignu} -eq 1 ]] ;  then
        # GFORTRAN
        export FC=/usr/local/openmpi-3.1.4-gfortran/bin/mpifort
        # release
        export CFLAGS="-cpp -O3 -Wno-aggressive-loop-optimizations -ffree-form -ffixed-line-length-132"
        export LDFLAGS="-O3"
	OPTFLAG="-march=native"
        if [[ ${idebug} -eq 1 ]] ; then
            # debug
            export CFLAGS="-cpp -O -g -pedantic-errors -Wall -W -Wno-maybe-uninitialized -ffree-form -ffixed-line-length-132 -fbacktrace -ffpe-trap=zero,overflow -finit-real=nan" #  -ffpe-trap=zero,overflow,underflow
            export LDFLAGS="-O"
	    OPTFLAG=
        fi
        # export CFLAGS="${CFLAGS} -march=native"
        export CFLAGS="${CFLAGS} -D__GFORTRAN__ -D__gFortran__"
        export LD=""
        export NCROOT="/usr/local/netcdf-fortran-4.4.5-gfortran"
        export cdir=".mpitmp-gfortran"
        export PROG=cable-mpi-gfortran
    elif [[ ${inag} -eq 1 ]] ;  then
        # NAG
        export FC=/usr/local/openmpi-3.1.5-nagfor/bin/mpifort
        # release
        export CFLAGS="-O4"
        export LDFLAGS="-O4"
	OPTFLAG=
        if [[ ${idebug} -eq 1 ]] ; then
            # debug
            # export CFLAGS="-C -C=dangling -g -nan -O0 -strict95 -gline"
            # set runtime environment variables: export NAGFORTRAN_RUNTIME_OPTIONS=show_dangling
            export CFLAGS="-C=alias -C=array -C=bits -C=dangling -C=do -C=intovf -C=present -C=pointer -C=recursion -g -nan -O0 -strict95 -gline"
            export LDFLAGS="-O0"
	    OPTFLAG=
        fi
        export CFLAGS="${CFLAGS} -fpp -colour -unsharedf95 -kind=byte -ideclient -ieee=full -free -not_openmp"
        export CFLAGS="${CFLAGS} -mismatch"
        export CFLAGS="${CFLAGS} -D__NAG__"
        export LD="-ideclient -unsharedrts"
        export NCROOT="/usr/local/netcdf-fortran-4.4.5-nagfor"
        export cdir=".mpitmp-nagfor"
        export PROG=cable-mpi-nagfor
    fi

    # All compilers
    export CFLAGS="${CFLAGS} ${OPTFLAG}"
    export CFLAGS="${CFLAGS} -D__MPI__"
    export CFLAGS="${CFLAGS} -D__CRU2017__"
    export CFLAGS="${CFLAGS} -D__NETCDF3__"
    # export CFLAGS="${CFLAGS} -D__C13DEBUG__"

    export NCCROOT="/usr/local"
    export NCCLIB=${NCCROOT}"/lib"
    export NCLIB=${NCROOT}"/lib"
    export NCMOD=${NCROOT}"/include"
    export LDFLAGS="-L${NCLIB} -lnetcdff -L${NCCLIB} -lnetcdf -lhdf5_hl -lhdf5 -lsz -lz "${LDFLAGS}
    export dosvn=0
    export MFLAGS="-j 8"
    build_build ${cdir} ${PROG}
    cd ../
    build_status ${cdir} ${PROG}
}


## MatthiasCuntz@Explor
host_vm_o()
{
    idebug=0
    iintel=1
    ignu=0
    np=$#
    for ((i=0; i<${np}; i++)) ; do
        if [[ "${1}" == "debug" ]] ; then
            idebug=1
            shift 1
        elif [[ "${1}" == "ifort" || "${1}" == "intel" ]] ; then
            iintel=1
            ignu=0
            shift 1
        elif [[ "${1}" == "gfortran" || "${1}" == "gnu" ]] ; then
            iintel=0
            ignu=1
            shift 1
        else
            echo "Error: command line option not known: " ${1}
            exit 1
        fi
    done
    if [[ ${iintel} -eq 1 ]] ;  then
        # INTEL - load mpi module first, otherwise intel module will not pre-pend LD_LIBRARY_PATH
        # module load intelmpi/2018.5.274
        # module load intel/2018.5
        # export FC=mpiifort
        module load openmpi/3.0.0/intel18
        module load intel/2018.5
        export FC=mpifort
        # release
        export CFLAGS="-fpp -O3 -nofixed -assume byterecl -fp-model precise -ip -diag-disable=10382"
        export LDFLAGS="-O3"
	OPTFLAG="-xBROADWELL"
        if [[ ${idebug} -eq 1 ]] ; then
            # debug
            export CFLAGS="-fpp -O0 -debug extended -traceback -g -check all,noarg_temp_created -warn all -fp-stack-check -nofixed -assume byterecl -fp-model precise -diag-disable=10382 -fpe0" # -fpe-all=0 -no-ftz -ftrapuv -init=arrays,snan
            export LDFLAGS="-O0"
	    OPTFLAG=
        fi
        # OPTFLAG="${CFLAGS} -march=broadwell"     # std / hf
        # OPTFLAG="${CFLAGS} -march=core-avx2"     # std / hf
        # OPTFLAG="${CFLAGS} -mtune=broadwell"     # std / hf
        # OPTFLAG="${CFLAGS} -march=skylake-avx512 # sky
        # OPTFLAG="${CFLAGS} -march=ivybridge"     # ivy / k20
        # OPTFLAG="${CFLAGS} -march=avx"           # ivy / k20
        # OPTFLAG="${CFLAGS} -mtune=ivybridge"     # ivy / k20
        export CFLAGS="${CFLAGS} -D__INTEL__ -D__INTEL_COMPILER__"
        export LD=""
        export NCROOT="/home/oqx29/zzy20/local/netcdf-fortran-4.4.4-ifort2018.0"
    else
        # GFORTRAN # 6.3.0 because of netcdf-fortran
        module load gcc/6.3.0
        module load openmpi/3.0.1/gcc/6.3.0
        export FC=mpifort
        # release
        export CFLAGS="-cpp -O3 -Wno-aggressive-loop-optimizations -ffree-form -ffixed-line-length-132"
        export LDFLAGS="-O3"
	OPTFLAG="-march=broadwell"
        if [[ ${idebug} -eq 1 ]] ; then
            # debug
            export CFLAGS="-cpp -O -g -pedantic-errors -Wall -W -Wno-maybe-uninitialized -ffree-form -ffixed-line-length-132 -fbacktrace -ffpe-trap=zero,overflow -finit-real=nan" #  -ffpe-trap=zero,overflow,underflow
            export LDFLAGS="-O"
	    OPTFLAG=
        fi
        # OPTFLAG="${CFLAGS} -march=broadwell"     # std / hf
        # OPTFLAG="${CFLAGS} -mavx2"               # std / hf
        # OPTFLAG="${CFLAGS} -march=skylake-avx512 # sky
        # OPTFLAG="${CFLAGS} -march=ivybridge"     # ivy / k20
        # OPTFLAG="${CFLAGS} -mavx"                # ivy / k20
        export CFLAGS="${CFLAGS} -D__GFORTRAN__ -D__gFortran__"
        export LD=""
        export NCROOT="/home/oqx29/zzy20/local/netcdf-fortran-4.4.4-gfortran63"
    fi

    # All compilers
    export CFLAGS="${CFLAGS} ${OPTFLAG}"
    export CFLAGS="${CFLAGS} -D__MPI__"
    # export CFLAGS="${CFLAGS} -D__C13DEBUG__"
    export CFLAGS="${CFLAGS} -D__CRU2017__"
    export CFLAGS="${CFLAGS} -D__NETCDF3__"

    export NCCROOT="/home/oqx29/zzy20/local"
    export NCCLIB=${NCCROOT}"/lib"
    export NCLIB=${NCROOT}"/lib"
    export NCMOD=${NCROOT}"/include"
    export LDFLAGS="-L${NCCLIB} -L${NCLIB} -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lsz -lz "${LDFLAGS}
    export dosvn=0
    export MFLAGS="-j 8"
    build_build
    cd ../
    build_status
}


clean_ask()
{
    print '\n\tPress Enter too continue building, Control-C to abort now.\n'
    read dummy
}


clean_build()
{
    rm -fr .mpitmp*
}


do_i_no_u()
{
   if [ -z ${PS3} ] ; then
       kmax=${#kh[*]}
       k=0
   else
       integer kmax=${#kh[*]}
       integer k=0
   fi
   typeset -f subr

   while [[ $k -lt $kmax ]] ; do
      if [[ $HOST_MACH = ${kh[$k]} || ${domain} = ${kd[$k]} ]] ; then
         echo 'Host recognized as' $HOST_MACH
         subr=host_${kh[$k]}
         $subr $*
      fi
      (( k = k + 1 ))
   done
}


build_status()
{
   if [[ $# -gt 0 ]] ; then export cdir="${1}" ; else export cdir='.mpitmp' ; fi
   if [[ $# -gt 1 ]] ; then export PROG="${2}" ; else export PROG='cable-mpi' ; fi

   if [[ -f ${cdir}/${PROG} ]]; then
        mv ${cdir}/${PROG} .
        print '\nBUILD OK\n'
   else
      print '\nOooops. Something went wrong\n'
      print '\nKnown build issues:\n'
      print '\nSome systems require additional library. \n'
      print '\nEdit Makefile_offline; add -lnetcdff to LD = ...\n'
   fi

   exit
}


build_build()
{
   if [[ ${dosvn} -eq 1 ]] ; then
       # write file for consumption by Fortran code
       # get SVN revision number
       CABLE_REV=`svn info | grep Revis |cut -c 11-18`
       if [[ $CABLE_REV = "" ]]; then
           echo "this is not an svn checkout"
           CABLE_REV=0
           echo "setting CABLE revision number to " $CABLE_REV
       fi
       print $CABLE_REV > ~/.cable_rev
       # get SVN status
       CABLE_STAT=`svn status`
       print $CABLE_STAT >> ~/.cable_rev
   fi

   if [[ $# -gt 0 ]] ; then export cdir="${1}" ; else export cdir='.mpitmp' ; fi
   if [[ $# -gt 1 ]] ; then export PROG="${2}" ; else export PROG='cable-mpi' ; fi

   if [[ ! -d ${cdir} ]]; then
      mkdir ${cdir}
   fi

   if [[ -f ${PROG} ]]; then
      print '\ncable-mpi executable exists. copying to a dated backup file\n'
      mv ${PROG} ${PROG}.`date +%d.%m.%y`
   fi

   # directories contain source code
   PHYS="../core/biogeophys"
   UTIL="../core/utils"
   DRV="."
   CASA="../core/biogeochem"
   BLAZE="../core/blaze"

   /bin/cp -p $PHYS/*90  ./${cdir}
   /bin/cp -p $UTIL/*90  ./${cdir}
   /bin/cp -p $DRV/*90   ./${cdir}
   /bin/cp -p $CASA/*90  ./${cdir}
   /bin/cp -p $BLAZE/*90 ./${cdir}

   /bin/cp -p Makefile_mpi ./${cdir}

   cd ${cdir}/
   make -f Makefile_mpi ${MFLAGS} PROG=${PROG}
}


#############################################
## build_mpi.ksh - MAIN SCRIPT STARTS HERE ##
#############################################

if [[ $1 = 'cclean' ]]; then
    print '\nCleaning up\n'
    clean_build
    shift 1
fi
if [[ $1 = 'clean' ]]; then
    print '\nCleaning up\n'
    clean_ask
    clean_build
    shift 1
fi

known_hosts
known_domains
HOST_MACH=`uname -n | cut -c 1-4 | tr - _`
if [ `uname -s` = 'Darwin' ] ; then
    domain=`hostname`
else
    domain=`hostname --domain`
fi
domain=${domain#*.}
do_i_no_u $*

echo "Error: host not recognized: "${HOST_MACH}" nor domain: "${domain}
exit 1
