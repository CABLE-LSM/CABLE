#!/bin/ksh

export dosvn=0 # 1/0: do/do not check svn

# so that script can be called by 'bash build.ksh' if no ksh installed
if [ -z ${PS3} ] ; then # https://stackoverflow.com/questions/3327013/how-to-determine-the-current-shell-im-working-on
    eval 'function print(){ printf "$@\n"; }'
fi

known_hosts()
{
    if [ -z ${PS3} ] ; then
        kh=(kh gadi pear mcin mc16 mcmi zlle zlhp vm_o logi auro bioc)
    else
        set -A kh gadi pear mcin mc16 mcmi zlle zlhp vm_o logi auro bioc
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
    # module load intel-compiler/2019.5.281
    # module load netcdf/4.6.3
    module load intel-compiler-llvm/2023.0.0
    module load netcdf/4.9.2

    export FC=ifort
    export NCDIR=${NETCDF_ROOT}"/lib/Intel"
    export NCMOD=${NETCDF_ROOT}"/include/Intel"
    # release
    # -ip only in CFLAGS, or -ipo in CFLAGS and LDFLAGS
    export CFLAGS="-fpp -O3 -nofixed -assume byterecl -fp-model precise -ip -diag-disable=10382,15009"
    export LDFLAGS="-O3"
    OPTFLAG="-march=broadwell -axSKYLAKE-AVX512,CASCADELAKE,SAPPHIRERAPIDS"
    export CFLAGS="${CFLAGS} ${OPTFLAG}"
    export CFLAGS="${CFLAGS} -D__INTEL__ -D__INTEL_COMPILER__"
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
        export CFLAGS="${CFLAGS} -fpp -colour -unsharedf95 -ideclient -ieee=full -free -not_openmp"
        export CFLAGS="${CFLAGS} -mismatch -mdir ."
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
    ignu=1
    inag=0
    np=$#
    for ((i=0; i<${np}; i++)) ; do
        if [[ "${1}" == "debug" ]] ; then
            idebug=1
            shift 1
        elif [[ "${1}" == "gfortran" || "${1}" == "gnu" ]] ; then
            ignu=1
            inag=0
            shift 1
        elif [[ "${1}" == "nagfor" || "${1}" == "nag" ]] ; then
            ignu=0
            inag=1
            shift 1
        else
            echo "Error: command line option not known: " ${1}
            exit 1
        fi
    done
    if [[ ${ignu} -eq 1 ]] ;  then
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
        export NCROOT="/usr/local/netcdf-fortran-4.6.1-gfortran"
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
        export CFLAGS="${CFLAGS} -fpp -colour -unsharedf95 -ideclient -ieee=full -free -not_openmp"
        export CFLAGS="${CFLAGS} -mismatch -mdir ."
        # export CFLAGS="${CFLAGS} -march=native"
        export CFLAGS="${CFLAGS} -D__NAG__"
        export LD="-ideclient -unsharedrts"
        export LDFLAGS="${LDFLAGS} -Wl,-ld_classic"  # until gcc update of Homebrew
        export NCROOT="/usr/local/netcdf-fortran-4.6.1-nagfor"
    fi

    # All compilers
    export CFLAGS="${CFLAGS} ${OPTFLAG}"
    # export CFLAGS="${CFLAGS} -D__CRU2017__"
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


host_mcmi()
{
    idebug=0
    ignu=1
    inag=0
    np=$#
    for ((i=0; i<${np}; i++)) ; do
        if [[ "${1}" == "debug" ]] ; then
            idebug=1
            shift 1
        elif [[ "${1}" == "gfortran" || "${1}" == "gnu" ]] ; then
            ignu=1
            inag=0
            shift 1
        elif [[ "${1}" == "nagfor" || "${1}" == "nag" ]] ; then
            ignu=0
            inag=1
            shift 1
        else
            echo "Error: command line option not known: " ${1}
            exit 1
        fi
    done
    if [[ ${ignu} -eq 1 ]] ;  then
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
        export NCROOT="/usr/local/netcdf-fortran-4.6.1-gfortran"
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
        # export CFLAGS="${CFLAGS} -fpp -colour -unsharedf95 -kind=byte -ideclient -ieee=full -free -not_openmp"
        export CFLAGS="${CFLAGS} -fpp -colour -unsharedf95 -ideclient -ieee=full -free -not_openmp"
        export CFLAGS="${CFLAGS} -mismatch -mdir ."
        # export CFLAGS="${CFLAGS} -march=native"
        export CFLAGS="${CFLAGS} -D__NAG__"
        export LD="-ideclient -unsharedrts"
        export NCROOT="/usr/local/netcdf-fortran-4.6.1-nagfor"
    fi

    # All compilers
    export CFLAGS="${CFLAGS} ${OPTFLAG}"
    # export CFLAGS="${CFLAGS} -D__CRU2017__"
    export CFLAGS="${CFLAGS} -D__NETCDF3__"
    # export CFLAGS="${CFLAGS} -D__C13DEBUG__"

    export NCCROOT="/opt/homebrew"
    export NCCLIB=${NCCROOT}"/lib"
    export NCLIB=${NCROOT}"/lib"
    export NCMOD=${NCROOT}"/include"
    export LDFLAGS="-L${NCLIB} -lnetcdff -L${NCCLIB} -lnetcdf  -lhdf5 -lhdf5_hl -lsz -lz -ldl -lzstd -lbz2 -lcurl -lxml2 "${LDFLAGS}
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
    module purge
    if [[ ${iintel} -eq 1 ]] ;  then
        # INTEL
        # 2018
        module load intel/2018.5
        # # 2023
        # module load intel/2019.4-full
        # module load hdf5/1.10.1/intelmpi/intel17
        # module load netcdf-c/4.7.2
        # module load netcdf-fortran/4.5.2
        export FC=ifort
        # release
        export CFLAGS="-fpp -O3 -nofixed -assume byterecl -fp-model precise -ip -diag-disable=10382"
        export LDFLAGS="-O3"
        # OPTFLAG="-xBROADWELL"
        OPTFLAG="-xHost"
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
        # 2018
        export NCCROOT="/home/oqx29/shared/local.save"
        export NCROOT="${NCCROOT}/netcdf-fortran-4.4.4-ifort2019.4"
        # # 2023
        # export NCCROOT=`nc-config --prefix`
        # export NCROOT=`nf-config --prefix`
    else
        # 2018
        # GFORTRAN # 6.3.0 because of netcdf-fortran
        module load gcc/6.3.0
        # # 2023
        # module purge
        # module use /opt/modulefiles/shared/mcs_mod
        # # module load softwares/anaconda3/2022.05
        # module load softwares/anaconda3/2023.03
        # source ${HOME_ANACONDA}/anaconda.rc
        # conda activate ${HOME}/.conda/envs/cpystd
        export FC=gfortran
        # release
        export CFLAGS="-cpp -O3 -Wno-aggressive-loop-optimizations -ffree-form -ffixed-line-length-132"
        export LDFLAGS="-O3"
        # OPTFLAG="-march=broadwell"
        OPTFLAG="-march=native"
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
        # 2018
        # export NCCROOT="/home/oqx29/shared/local.gnu"
        # export NCROOT=${NCCROOT}
        export NCCROOT="/home/oqx29/shared/local.save"
        export NCROOT="/home/oqx29/shared/local.save/netcdf-fortran-4.4.4-gfortran63/"
        # export LDFLAGS="-lhdf5_hl -lhdf5 -lsz -lz -lssl -lcrypto ${LDFLAGS}"
        # # 2023
        # export NCCROOT=`nc-config --prefix`
        # export LDFLAGS="-L${NCLIB}/lib -lcurl -L/lib ${LDFLAGS}"  # for libc
    fi

    # All compilers
    export CFLAGS="${CFLAGS} ${OPTFLAG}"
    # export CFLAGS="${CFLAGS} -D__C13DEBUG__"
    # export CFLAGS="${CFLAGS} -D__CRU2017__"
    export CFLAGS="${CFLAGS} -D__NETCDF3__"

    export NCCLIB=${NCCROOT}"/lib"
    export NCLIB=${NCROOT}"/lib"
    export NCMOD=${NCROOT}"/include"
    # 2018
    export LDFLAGS="-L${NCCLIB} -L${NCLIB} -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lsz -lz -lssl -lcrypto "${LDFLAGS} # " -pg"
    # # 2023
    # export LDFLAGS="-L${NCLIB} -lnetcdff -L${NCCLIB} -lnetcdf ${LDFLAGS}" # " -pg"
    export dosvn=0
    export MFLAGS="-j 8"
    build_build
    cd ../
    build_status
}


## MatthiasCuntz@curta
host_logi()
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
        # # INTEL - curta
        # module load intel-oneapi-compilers/2023.1.0/gcc@11.2.0-57kkzxu
        # module load netcdf-c/4.9.2/oneapi@2023.0.0-2rbrj54
        # module load netcdf-fortran/4.6.0/oneapi@2023.0.0-n24bykw
	# module load gcc/11.2.0
        # INTEL - curta2
        module load intel-oneapi-compilers/2023.1.0/gcc@11.2.0-rex53zv
        module load netcdf-c/4.9.2/oneapi@2023.0.0-eyrwrda
        module load netcdf-fortran/4.6.0/oneapi@2023.0.0-2jxietz
	module load gcc/11.2.0
	# # Does not work on curta with the following error:
	# #     ifort: /lib64/libstdc++.so.6: version `GLIBCXX_3.4.21' not found (required by ifort)
	# # Does not work on curta2 because compiler/intel/2020.4.304 does not exist
	# module load compiler/intel/2020.4.304
	# module load netcdf-intel
	# module load netcdf-fortran-intel
        export FC=ifort
        # release
        export CFLAGS="-fpp -O3 -nofixed -assume byterecl -fp-model precise -ip -diag-disable=10382"
        export LDFLAGS="-O3"
        # OPTFLAG="-xBROADWELL"
        OPTFLAG="-xHost"
        if [[ ${idebug} -eq 1 ]] ; then
            # debug
            export CFLAGS="-fpp -O0 -debug extended -traceback -g -check all,noarg_temp_created -warn all -fp-stack-check -nofixed -assume byterecl -fp-model precise -diag-disable=10382 -fpe0" # -fpe-all=0 -no-ftz -ftrapuv -init=arrays,snan
            export LDFLAGS="-O0"
            OPTFLAG=
        fi
        export CFLAGS="${CFLAGS} -D__INTEL__ -D__INTEL_COMPILER__"
        export LD=""
	NCCFLAGS=`pkg-config --cflags netcdf`
	NCFLAGS=`pkg-config --cflags netcdf-fortran`
	export CFLAGS="${CFLAGS} ${NCFLAGS} ${NCCFLAGS}"
    else
        # 2018
        # GFORTRAN # 6.3.0 because of netcdf-fortran
        module load gcc/6.3.0
        # # 2023
        # module purge
        # module use /opt/modulefiles/shared/mcs_mod
        # # module load softwares/anaconda3/2022.05
        # module load softwares/anaconda3/2023.03
        # source ${HOME_ANACONDA}/anaconda.rc
        # conda activate ${HOME}/.conda/envs/cpystd
        export FC=gfortran
        # release
        export CFLAGS="-cpp -O3 -Wno-aggressive-loop-optimizations -ffree-form -ffixed-line-length-132"
        export LDFLAGS="-O3"
        # OPTFLAG="-march=broadwell"
        OPTFLAG="-march=native"
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
        # 2018
        # export NCCROOT="/home/oqx29/shared/local.gnu"
        # export NCROOT=${NCCROOT}
        export NCCROOT="/home/oqx29/shared/local.save"
        export NCROOT="/home/oqx29/shared/local.save/netcdf-fortran-4.4.4-gfortran63/"
        # export LDFLAGS="-lhdf5_hl -lhdf5 -lsz -lz -lssl -lcrypto ${LDFLAGS}"
        # # 2023
        # export NCCROOT=`nc-config --prefix`
        # export LDFLAGS="-L${NCLIB}/lib -lcurl -L/lib ${LDFLAGS}"  # for libc
    fi

    # All compilers
    export CFLAGS="${CFLAGS} ${OPTFLAG}"
    # export CFLAGS="${CFLAGS} -D__C13DEBUG__"
    # export CFLAGS="${CFLAGS} -D__CRU2017__"
    export CFLAGS="${CFLAGS} -D__NETCDF3__"

    # NCCLIBS=`nc-config --static --libs`
    # NCLIBS=`nf-config --flibs`
    # export LDFLAGS="${NCLIBS} ${NCCLIBS} -static-intel ${LDFLAGS}"
    NCCLIBS=`pkg-config --libs netcdf`
    NCLIBS=`pkg-config --libs netcdf-fortran`
    export NCMOD=`nf-config --includedir`
    export LDFLAGS="${NCLIBS} ${NCCLIBS} ${LDFLAGS}"
    export dosvn=0
    export MFLAGS="-j 8"
    build_build
    cd ../
    build_status
}

## Lars Nieradzik @ aurora.lunarc.lu.se
host_auro()
{
    if [ -z ${PS3} ] ; then
        . /etc/bashrc
    else
        . /etc/kshrc
    fi
    module purge
    module load intel/2020a
    module load netCDF-Fortran/4.5.2

    export FC=ifort
    export NETCDF_RT="/sw/easybuild/software/netCDF-Fortran/4.5.2-iimpi-2020a"
    export NCDIR=${NETCDF_RT}"/lib"
    export NCMOD=${NETCDF_RT}"/include"
    if [[ ${1} == "debug" ]]; then
        # debug
        # export CFLAGS='-O0 -fpp -traceback -g -fp-model precise -ftz -fpe0'
        #CLNexport CFLAGS="-fpp -O0 -debug extended -traceback -g -check all,noarg_temp_created -warn all -fp-stack-check -nofixed -assume byterecl -fp-model precise -diag-disable=10382 -fpe0" # -fpe-all=0 -no-ftz -ftrapuv"
        export CFLAGS="-fpp -O0 -debug extended -traceback -g -check all,noarg_temp_created -warn all -fp-stack-check -nofixed -assume byterecl -fp-model precise -diag-disable=10382 -fpe0" # -fpe-all=0 -no-ftz -ftrapuv"export LDFLAGS="-O0"
        OPTFLAG=""
    else
        # release
        export CFLAGS='-O2 -fpp -fp-model precise'
        export CFLAGS="-fpp -O3 -nofixed -assume byterecl -fp-model precise -ip -diag-disable=10382"
        export LDFLAGS="-O3"
        OPTFLAG=""
    fi
    export CFLAGS="${CFLAGS} ${OPTFLAG}"
    export CFLAGS="${CFLAGS} -D__INTEL__ -D__INTEL_COMPILER__"
    export LDFLAGS="-L"${NCDIR}" "${LDFLAGS}
    export LD="-lnetcdf -lnetcdff"
    export MFLAGS="-j 8"
    build_build
    cd ../
    build_status
}


host_zlle()
{
    idebug=0
    ignu=1
    np=$#
    for ((i=0; i<${np}; i++)) ; do
        if [[ "${1}" == "debug" ]] ; then
            idebug=1
            shift 1
        else
            echo "Error: command line option not known: " ${1}
            exit 1
        fi
    done
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
    export NCROOT="/usr"

    # All compilers
    export CFLAGS="${CFLAGS} ${OPTFLAG}"
    # export CFLAGS="${CFLAGS} -D__CRU2017__"
    export CFLAGS="${CFLAGS} -D__NETCDF3__"
    # export CFLAGS="${CFLAGS} -D__C13DEBUG__"

    export NCCROOT="/usr"
    export NCCLIB=${NCCROOT}"/lib/x86_64-linux-gnu"
    export NCLIB=${NCROOT}"/lib/x86_64-linux-gnu"
    export NCMOD=${NCROOT}"/include"
    export LDFLAGS="-Wl,-Bsymbolic-functions -flto=auto -ffat-lto-objects -Wl,-z,relro -Wl,-z,now "${LDFLAGS}
    export LD="-L${NCLIB} -lnetcdff -L${NCCLIB} -lnetcdf -lm"
    export dosvn=0
    export MFLAGS="-j 8"
    build_build
    cd ../
    build_status
}

host_zlhp()
{
    idebug=0
    ignu=1
    np=$#
    for ((i=0; i<${np}; i++)) ; do
        if [[ "${1}" == "debug" ]] ; then
            idebug=1
            shift 1
        else
            echo "Error: command line option not known: " ${1}
            exit 1
        fi
    done
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
    export NCROOT="/usr"

    # All compilers
    export CFLAGS="${CFLAGS} ${OPTFLAG}"
    # export CFLAGS="${CFLAGS} -D__CRU2017__"
    export CFLAGS="${CFLAGS} -D__NETCDF3__"
    # export CFLAGS="${CFLAGS} -D__C13DEBUG__"

    export NCCROOT="/usr"
    export NCCLIB=${NCCROOT}"/lib/x86_64-linux-gnu"
    export NCLIB=${NCROOT}"/lib/x86_64-linux-gnu"
    export NCMOD=${NCROOT}"/include"
    export LDFLAGS="-Wl,-Bsymbolic-functions -flto=auto -ffat-lto-objects -Wl,-z,relro -Wl,-z,now "${LDFLAGS}
    export LD="-L${NCLIB} -lnetcdff -L${NCCLIB} -lnetcdf -lm"
    export dosvn=0
    export MFLAGS="-j 8"
    build_build
    cd ../
    build_status
}

## biocomp
host_bioc()
{
    idebug=0
    np=$#
    for ((i=0; i<${np}; i++)) ; do
        if [[ "${1}" == "debug" ]] ; then
            idebug=1
            shift 1
        else
            echo "Error: command line option not known: " ${1}
            exit 1
        fi
    done
    #
    eval "$(${HOME}/miniconda3/bin/conda 'shell.bash' 'hook' 2> /dev/null)"
    conda activate pystd
    #
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
    export CFLAGS="${CFLAGS} -D__GFORTRAN__ -D__gFortran__" # -pg # gprof --line ./cable gmon.out
    export CFLAGS="${CFLAGS} ${OPTFLAG}"
    # export CFLAGS="${CFLAGS} -D__C13DEBUG__"
    export CFLAGS="${CFLAGS} -D__NETCDF3__"
    #
    export LD=""
    export NCCROOT="${HOME}/miniconda3/envs/pystd"
    export NCROOT="${HOME}/miniconda3/envs/pystd"
    # export LDFLAGS="-lmfhdf -ldf -lhdf5_hl -lhdf5 -lcrypto -lcurl -lpthread -lsz -lz -ldl -lm -lzip -lblosc -lzstd -lbz2 -lxml2 ${LDFLAGS}"
    export NCCLIB=${NCCROOT}"/lib"
    export NCLIB=${NCROOT}"/lib"
    export NCMOD=${NCROOT}"/include"
    export LDFLAGS="-L${NCCLIB} -L${NCLIB} -lnetcdff -lnetcdf "${LDFLAGS} # " -pg"
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

    # # directories contain source code (relative to offline/)
    # PHYS="../core/biogeophys"
    # UTIL="../core/utils"
    # DRV="."
    # CASA="../core/biogeochem"
    # BLAZE="../core/blaze"

    # directories contain source code (relative to offline/old/)
    PHYS="../../core/biogeophys"
    UTIL="../../core/utils"
    DRV=".."
    CASA="../../core/biogeochem"
    BLAZE="../../core/blaze"

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
# host2: for compilation on compute nodes of Explor
host2=`echo ${HOST_MACH} | cut -c1,2`
if [ ${host2} = 'cn' ] ; then HOST_MACH="vm_o" ; fi
if [ `uname -s` = 'Darwin' ] ; then
    domain=`hostname`
else
    domain=`hostname --domain`
fi
domain=${domain#*.}
do_i_no_u $*

echo "Error: host not recognized: "${HOST_MACH}" nor domain: "${domain}
exit 1