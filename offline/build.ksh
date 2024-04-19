#!/bin/ksh

export dosvn=0 # 1/0: do/do not check svn

# so that script can be called by 'bash build.ksh' if no ksh installed
if [ -z ${PS3} ] ; then # https://stackoverflow.com/questions/3327013/how-to-determine-the-current-shell-im-working-on
    eval 'function print(){ printf "$@\n"; }'
fi

known_hosts()
{
    if [ -z ${PS3} ] ; then
        kh=(kh gadi pear mcin mc16 mcmi vm_o auro)
    else
        set -A kh gadi pear mcin mc16 mcmi vm_o auro
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
    export MFLAGS=""
    export CFLAGS='-O0 -traceback -g -fp-model precise -ftz -fpe0 -fpp'
    #export CFLAGS=" -O2 -fpp -fp-model precise"
    #debugPOP#export CFLAGS="-fpp -fpe0 -O0 -fpe-all=0 -g -debug -traceback -fp-stack-check -no-ftz -ftrapuv -check bounds "
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
    UTIL3="../util"
    DRV="."
    CASA="../core/biogeochem"
    BLAZE="../core/blaze"
    SCI="../science/*"
   PAR="../params"

    /bin/cp -p $PHYS/*90  ./.tmp
    /bin/cp -p $UTIL/*90  ./.tmp
    /bin/cp -p $UTIL3/*90  ./.tmp
    /bin/cp -p $DRV/*90   ./.tmp
    #/bin/cp -p $CASA/*90  ./.tmp
    /bin/cp -p $BLAZE/*90 ./.tmp
    /bin/cp -p $SCI/*90   ./.tmp
/bin/cp -p $PAR/*90 ./.tmp

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
