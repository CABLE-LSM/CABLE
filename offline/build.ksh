#!/bin/ksh

export dosvn=0 # 1/0: do/do not check svn

# so that script can be called by 'bash build.ksh' if no ksh installed
if [ -z ${PS3} ] ; then # https://stackoverflow.com/questions/3327013/how-to-determine-the-current-shell-im-working-on
    eval 'function print(){ printf "$@\n"; }'
fi

known_hosts()
{
    if [ -z ${PS3} ] ; then
        kh=(kh gadi pear mcin mc16 vm_o)
    else
        set -A kh gadi pear mcin mc16 vm_o
    fi
}


known_domains()
{
    if [ -z ${PS3} ] ; then
        kd=(kd nci.org.au pear local local explor)
    else
        set -A kd nci.org.au pear local local explor
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
    module load netcdf/4.6.3

    export FC=ifort
    export NCDIR=${NETCDF_ROOT}"/lib/Intel"
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
        # OPTFLAG="${CFLAGS} -xCORE-AVX2 -axSKYLAKE-AVX512,CASCADELAKE" # given in user training: does not work
        # OPTFLAG="${CFLAGS} -xCASCADELAKE" # or -xCORE-AVX512;                           queues: express / normal
        # OPTFLAG="${CFLAGS} -xBROADWELL"   # or -xCORE-AVX512;                           queues: expressbw / normalbw
        # OPTFLAG="${CFLAGS} -xSKYLAKE"     # or -xSKYLAKE-AVX512 depends on performance; queues: normalsl
    fi
    export CFLAGS="${CFLAGS} ${OPTFLAG}"
    # export CFLAGS="${CFLAGS} -D__CRU2017__"
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
    # if [ ! -z ${PS3} ] ; then
    . /apps/modules/Modules/default/init/ksh
    # fi

    module del intel-cc intel-fc
    module add intel-cc/16.0.1.150 intel-fc/16.0.1.150
    module add netcdf/4.3.3.1

    export NCDIR=$NETCDF_ROOT'/lib/'
    export NCMOD=$NETCDF_ROOT'/include/'
    export FC='ifort'
    #export  CFLAGS='-O0 -fp-model precise -fpe0 -fpp -g -debug -traceback -fp-stack-check -no-ftz -ftrapuv -check all,noarg_temp_created -C '
    #export CFLAGS='-O0 -fpe=0 -fpe-all=0 -fpp -g -debug -traceback -fp-stack-check -no-ftz -ftrapuv -check bounds
    export CFLAGS='-O2 -fp-model precise -fpp'
    # export CFLAGS="${CFLAGS} -D__CRU2017__"
    export CFLAGS="${CFLAGS} -D__NETCDF3__"
    export LDFLAGS='-g -L'$NCDIR' -O2'
    export LD='-lnetcdf -lnetcdff'
    export dosvn=0
    export MFLAGS='-j 8'
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
        export FC=ifort
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
    elif [[ ${ignu} -eq 1 ]] ;  then
        # GFORTRAN
        export FC=gfortran
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
    elif [[ ${inag} -eq 1 ]] ;  then
        # NAG
        export FC=nagfor
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
        # export CFLAGS="${CFLAGS} -march=native"
        export CFLAGS="${CFLAGS} -D__NAG__"
        export LD="-ideclient -unsharedrts"
        export NCROOT="/usr/local/netcdf-fortran-4.4.5-nagfor"
    fi

    # All compilers
    export CFLAGS="${CFLAGS} ${OPTFLAG}"
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
    build_build
    cd ../
    build_status
}


host_mc16()
{
    idebug=0
    iintel=0
    ignu=1
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
        export FC=ifort
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
        export NCROOT="/usr/local/netcdf-fortran-4.5.3-ifort"
    elif [[ ${ignu} -eq 1 ]] ;  then
        # GFORTRAN
        export FC=gfortran
        # release
        export CFLAGS="-cpp -O3 -Wno-aggressive-loop-optimizations -ffree-form -ffixed-line-length-132 -frecursive"
        export LDFLAGS="-O3"
        OPTFLAG="-march=native"
        if [[ ${idebug} -eq 1 ]] ; then
            # debug
            export CFLAGS="-cpp -O -g -pedantic-errors -Wall -W -Wno-maybe-uninitialized -ffree-form -ffixed-line-length-132 -frecursive -fbacktrace -ffpe-trap=zero,overflow -finit-real=nan" #  -ffpe-trap=zero,overflow,underflow
            export LDFLAGS="-O"
            OPTFLAG=
        fi
        # export CFLAGS="${CFLAGS} -march=native"
        export CFLAGS="${CFLAGS} -D__GFORTRAN__ -D__gFortran__"
        export LD=""
        export NCROOT="/usr/local/netcdf-fortran-4.5.3-gfortran"
    elif [[ ${inag} -eq 1 ]] ;  then
        # NAG
        export FC=nagfor
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
        # export CFLAGS="${CFLAGS} -march=native"
        export CFLAGS="${CFLAGS} -D__NAG__"
        export LD="-ideclient -unsharedrts"
        export NCROOT="/usr/local/netcdf-fortran-4.5.3-nagfor"
    fi

    # All compilers
    export CFLAGS="${CFLAGS} ${OPTFLAG}"
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
    build_build
    cd ../
    build_status
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
        # INTEL
        module load intel/2018.5
        export FC=ifort
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
        export FC=gfortran
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
        export CFLAGS="${CFLAGS} -D__GFORTRAN__ -D__gFortran__" # -pg # gprof --line ./cable gmon.out
        export LD=""
        export NCROOT="/home/oqx29/zzy20/local/netcdf-fortran-4.4.4-gfortran63"
    fi

    # All compilers
    export CFLAGS="${CFLAGS} ${OPTFLAG}"
    # export CFLAGS="${CFLAGS} -D__C13DEBUG__"
    export CFLAGS="${CFLAGS} -D__CRU2017__"
    export CFLAGS="${CFLAGS} -D__NETCDF3__"

    export NCCROOT="/home/oqx29/zzy20/local"
    export NCCLIB=${NCCROOT}"/lib"
    export NCLIB=${NCROOT}"/lib"
    export NCMOD=${NCROOT}"/include"
    export LDFLAGS="-L${NCCLIB} -L${NCLIB} -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lsz -lz "${LDFLAGS} # " -pg"
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
    rm -fr .tmp
}


# returns the index of the first appearance of
# the first argument string in the rest of the arguments.
# returns -1 if not present.
# Can return 0, i.e. script exists if set -e
iisin()
{
    tofind=${1}
    shift
    ii=0
    lauf=-1
    for i in $@ ; do
        (( lauf = lauf + 1 ))
        if [[ ${i} == ${tofind} ]] ; then
            ii=${lauf}
            break
        fi
    done
    echo ${ii}
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

    ii=`iisin ${HOST_MACH} ${kh[*]}`
    if [[ ${ii} -ge 0 ]] ; then
        echo 'Host recognized as' ${HOST_MACH}
        subr=host_${kh[$ii]}
    else
        ii=`iisin ${domain} ${kd[*]}`
        if [[ ${ii} -ge 0 ]] ; then
            echo 'Domain recognized as' ${domain}
            subr=host_${kd[$ii]}
        else
            echo "Neither host nor domain recognized: host ${HOST_MACH} and domain ${domain}."
            exit
        fi
    fi
    ${subr} $*
}


build_status()
{
    if [[ -f .tmp/cable ]]; then
        mv .tmp/cable .
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
        echo $CABLE_REV > ~/.cable_rev
        # get SVN status
        CABLE_STAT=`svn status`
        echo $CABLE_STAT >> ~/.cable_rev
    fi

    if [[ ! -d .tmp ]]; then
        mkdir .tmp
    fi

    if [[ -f cable ]]; then
        print '\ncable executable exists. copying to dated backup file\n'
        mv cable cable.`date +%d.%m.%y`
    fi

    # directories contain source code
    PHYS="../core/biogeophys"
    UTIL="../core/utils"
    DRV="."
    CASA="../core/biogeochem"
    BLAZE="../core/blaze"

    /bin/cp -p $PHYS/*90  ./.tmp
    /bin/cp -p $UTIL/*90  ./.tmp
    /bin/cp -p $DRV/*90   ./.tmp
    /bin/cp -p $CASA/*90  ./.tmp
    /bin/cp -p $BLAZE/*90 ./.tmp

    /bin/cp -p Makefile_offline ./.tmp

    cd .tmp/
    make -f Makefile_offline ${MFLAGS}
}


###########################################
## build.ksh - MAIN SCRIPT STARTS HERE   ##
###########################################

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
