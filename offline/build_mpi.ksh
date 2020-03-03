#!/bin/ksh

export dosvn=0 # 1/0: do/do not check svn

# so that script can be called by 'bash build.ksh' if no ksh installed
if [ -z ${PS3} ] ; then # https://stackoverflow.com/questions/3327013/how-to-determine-the-current-shell-im-working-on
    eval 'function print(){ printf "$@\n"; }'
fi

known_hosts()
{
    if [ -z ${PS3} ] ; then
	kh=(kh cher burn shin raij pear mcin vm_o gadi)
    else
	set -A kh cher burn shin raij pear mcin vm_o gadi
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
   # module load hdf5/1.10.5
   
   export FC=mpif90
   # export NCCLIB="${NETCDF_ROOT}/lib"
   # export NCLIB="${NETCDF_ROOT}/lib/Intel"
   export NCMOD="${NETCDF_ROOT}/include/Intel"
   # export HDF5LIB="${HDF5_ROOT}/lib"
   if [[ ${1} == 'debug' ]]; then
       #trunk export CFLAGS='-O0 -traceback -g -fp-model precise -ftz -fpe0'
       #build.ksh export CFLAGS='-O0 -fpp -traceback -g -fp-model precise -ftz -fpe0'
       export CFLAGS='-g -debug -traceback -fpp -check all,noarg_temp_created -fp-stack-check -O0 -debug -fpe=0 -fpe-all=0 -no-ftz -ftrapuv'
   else
       export CFLAGS='-O2 -fpp -fp-model precise'
       # export CFLAGS="${CFLAGS} -xCORE-AVX2 -axSKYLAKE-AVX512,CASCADELAKE" # given in user training: does not work
       export CFLAGS="${CFLAGS} -xCASCADELAKE" # or -xCORE-AVX512 # queues: express / normal
       # export CFLAGS="${CFLAGS} -xBROADWELL"  # or -xCORE-AVX512; queues: expressbw / normalbw
       # export CFLAGS="${CFLAGS} -xSKYLAKE"    # or -xSKYLAKE-AVX512 depends on performance; queues: normalsl
   fi
   export CFLAGS="${CFLAGS} -D__MPI__"
   export CFLAGS="${CFLAGS} -D__CRU2017__"
   export LDFLAGS="-O2"
   # export LD='-L${NCLIB} -lnetcdff -L${NCCLIB} -lnetcdf -L${HDF5LIB} -lhdf5_hl -lhdf5 -ldl -lz -lcurl -lm'
   export LD='-lnetcdf -lnetcdff'
   export MFLAGS='-j 8'
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
     export CFLAGS='-O2 -fp-model precise -fpp'
    #   export CFLAGS='-O0 -C'
    # best settings for debugging
    #   export CFLAGS='-O0 -C -g -debug all -traceback -check all,noarg_temp_created, -C  '
    #   export CFLAGS='-O0 '

    #export CFLAGS='-O0 -fp-model precise -g -debug -traceback -fpp '
    export CFLAGS="${CFLAGS} -D__CRU2018__"
    #   export CFLAGS='-O0 -fp-model precise -g -debug all -traceback -fpe0 '
    #   export CFLAGS='  -g -debug -traceback -fp-stack-check -O0 -debug -fpe0 -no-ftz -ftrapuv'

    # best debug flags
    #   export LDFLAGS='-g -L'$NCDIR  #'-L'$NCDIR' -O2'
    export CFLAGS="${CFLAGS} -D__MPI__"
    export LDFLAGS='-O2 -L'$NCDIR''
    export MFLAGS='-j 8'
    export LD='-lnetcdf -lnetcdff'
    build_build
    cd ../
    build_status
}


# MatthiasCuntz@INRAE
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
	fi
    done
    if [[ ${iintel} -eq 1 ]] ;  then
	# INTEL
	/opt/intel/compilers_and_libraries/mac/bin/compilervars.sh intel64
	export FC=/usr/local/openmpi-3.1.5-ifort/bin/mpifort
	# release
	# export CFLAGS="-O3 -fpp -nofixed -assume byterecl -fp-model precise -m64 -ip -xHost -diag-disable=10382"
	# Gadi flags
	export CFLAGS='-O2 -fpp -fp-model precise'
	if [[ ${idebug} -eq 1 ]] ; then
	    # debug -fpe0
	    # export CFLAGS="-check all,noarg_temp_created -warn all -g -debug -traceback -fp-stack-check -O0 -debug -fpp -nofixed -assume byterecl -fp-model precise -m64 -ip -xHost -diag-disable=10382"
            export CFLAGS="-check all,noarg_temp_created -warn all -g -debug -traceback -fp-stack-check -O0 -fpp -nofixed -assume byterecl -fp-model precise -diag-disable=10382 -fpe0 -fpe-all=0 -no-ftz -ftrapuv -init=arrays,snan"
	fi
	# export CFLAGS="${CFLAGS} -mtune=corei7"
	# export CFLAGS="${CFLAGS} -march=native"
	export CFLAGS="${CFLAGS} -D__INTEL__ -D__INTEL_COMPILER__"
	export LD=''
	export NCROOT='/usr/local/netcdf-fortran-4.4.5-ifort'
	export cdir='.mpitmp-ifort'
	export PROG=cable-mpi-ifort
    elif [[ ${ignu} -eq 1 ]] ;  then
        # GFORTRAN
	export FC=/usr/local/openmpi-3.1.4-gfortran/bin/mpifort
	# release
	export CFLAGS="-O3 -Wno-aggressive-loop-optimizations -cpp -ffree-form -ffixed-line-length-132"
	if [[ ${idebug} -eq 1 ]] ; then
	    # debug
	    export CFLAGS="-pedantic-errors -Wall -W -O0 -g -Wno-maybe-uninitialized -cpp -ffree-form -ffixed-line-length-132"
	fi
	# export CFLAGS="${CFLAGS} -march=native"
	export CFLAGS="${CFLAGS} -D__GFORTRAN__ -D__gFortran__"
	export LD=''
	export NCROOT='/usr/local/netcdf-fortran-4.4.5-gfortran'
	export cdir='.mpitmp-gfortran'
	export PROG=cable-mpi-gfortran
    elif [[ ${inag} -eq 1 ]] ;  then
        # NAG
	export FC=/usr/local/openmpi-3.1.5-nagfor/bin/mpifort
	# release
        export CFLAGS="-O4"
	if [[ ${idebug} -eq 1 ]] ; then
	    # debug
            #export CFLAGS="-C -C=dangling -g -nan -O0 -strict95 -gline"
	    # set runtime environment variables: export NAGFORTRAN_RUNTIME_OPTIONS=show_dangling
            export CFLAGS="-C=alias -C=array -C=bits -C=dangling -C=do -C=intovf -C=present -C=pointer -C=recursion -g -nan -O0 -strict95 -gline"
	fi
	export CFLAGS="${CFLAGS} -fpp -colour -unsharedf95 -kind=byte -ideclient -ieee=full -free -not_openmp"
	export CFLAGS="${CFLAGS} -mismatch"
	# export CFLAGS="${CFLAGS} -march=native"
	export CFLAGS="${CFLAGS} -D__NAG__"
        export LD='-ideclient -unsharedrts'
        export NCROOT='/usr/local/netcdf-fortran-4.4.5-nagfor'
	export cdir='.mpitmp-nagfor'
	export PROG=cable-mpi-nagfor
    fi
    export CFLAGS="${CFLAGS} -D__MPI__"
    # export CFLAGS="${CFLAGS} -D__C13DEBUG__"
    export CFLAGS="${CFLAGS} -D__CRU2017__"

    # # NAG - Does not work for pop_io.f90
    # export FC=nagfor
    # # release
    # export CFLAGS="-O4 -fpp -colour -unsharedf95 -kind=byte -ideclient -ieee=full -free"
    # if [[ ${1} = 'debug' ]] ; then
    #     # debug
    #     export CFLAGS="-C -C=dangling -g -nan -O0 -strict95 -gline -fpp -colour -unsharedf95 -kind=byte -ideclient -ieee=full -free -D__NAG__"
    # fi
    # export LD='-ideclient -unsharedrts'
    # export NCROOT='/usr/local/netcdf-fortran-4.4.5-nagfor'

    # All compilers
    export NCCROOT='/usr/local'
    export NCCLIB=${NCCROOT}'/lib'
    export NCLIB=${NCROOT}'/lib'
    export NCMOD=${NCROOT}'/include'
    export LDFLAGS="-L${NCLIB} -lnetcdff -L${NCCLIB} -lnetcdf -lhdf5_hl -lhdf5 -lsz -lz"
    # Gadi flags
    export LDFLAGS="-O0 ${LDFLAGS}"
    export dosvn=0
    export MFLAGS='-j 8'
    build_build ${cdir} ${PROG}
    cd ../
    build_status ${cdir} ${PROG}
}


# MatthiasCuntz@Explor
host_vm_o()
{
    idebug=0
    iintel=0
    np=$#
    for ((i=0; i<${np}; i++)) ; do
        if [[ "${1}" == "debug" ]] ; then
            idebug=1
            shift 1
        elif [[ "${1}" == "ifort" || "${1}" == "intel" ]] ; then
            iintel=1
            shift 1
        elif [[ "${1}" == "gfortran" || "${1}" == "gnu" ]] ; then
            iintel=0
            shift 1
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
        export CFLAGS="-O3 -fpp -nofixed -assume byterecl -fp-model precise -m64 -ip -xHost -diag-disable=10382"
        if [[ ${idebug} -eq 1 ]] ; then
            # debug
            export CFLAGS="-check all,noarg_temp_created -warn all -g -debug -traceback -fp-stack-check -O0 -debug -fpp -nofixed -assume byterecl -fp-model precise -m64 -ip -xHost -diag-disable=10382"
        fi
	# export CFLAGS="${CFLAGS} -march=broadwell"     # std / hf
	# export CFLAGS="${CFLAGS} -march=core-avx2"     # std / hf
	# export CFLAGS="${CFLAGS} -mtune=broadwell"     # std / hf
	# export CFLAGS="${CFLAGS} -march=skylake-avx512 # sky
	# export CFLAGS="${CFLAGS} -march=ivybridge"     # ivy / k20
	# export CFLAGS="${CFLAGS} -march=avx"           # ivy / k20
	# export CFLAGS="${CFLAGS} -mtune=ivybridge"     # ivy / k20
	export CFLAGS="${CFLAGS} -D__INTEL__ -D__INTEL_COMPILER__"
        export LD=''
        export NCROOT='/home/oqx29/zzy20/local/netcdf-fortran-4.4.4-ifort2018.0'
    else
        # GFORTRAN # 6.3.0 because of netcdf-fortran
        module load gcc/6.3.0
	module load openmpi/3.0.1/gcc/6.3.0
        export FC=mpifort
        # release
        export CFLAGS="-O3 -Wno-aggressive-loop-optimizations -cpp -ffree-form -ffixed-line-length-132"
        if [[ ${idebug} -eq 1 ]] ; then
            # debug
            export CFLAGS="-pedantic-errors -Wall -W -O -g -Wno-maybe-uninitialized -cpp -ffree-form -ffixed-line-length-132"
        fi
	# export CFLAGS="${CFLAGS} -march=broadwell"     # std / hf
	# export CFLAGS="${CFLAGS} -mavx2"               # std / hf
	# export CFLAGS="${CFLAGS} -march=skylake-avx512 # sky
	# export CFLAGS="${CFLAGS} -march=ivybridge"     # ivy / k20
	# export CFLAGS="${CFLAGS} -mavx"                # ivy / k20
	export CFLAGS="${CFLAGS} -D__GFORTRAN__ -D__gFortran__"
        export LD=''
        export NCROOT='/home/oqx29/zzy20/local/netcdf-fortran-4.4.4-gfortran63'
    fi
    export CFLAGS="${CFLAGS} -D__MPI__"
    # export CFLAGS="${CFLAGS} -D__C13DEBUG__"
    export CFLAGS="${CFLAGS} -D__CRU2017__"

    # # NAG - Does not work for pop_io.f90
    # export FC=nagfor
    # # release
    # export CFLAGS="-O4 -fpp -colour -unsharedf95 -kind=byte -ideclient -ieee=full -free"
    # if [[ ${1} = 'debug' ]] ; then
    #     # debug
    #     export CFLAGS="-C -C=dangling -g -nan -O0 -strict95 -gline -fpp -colour -unsharedf95 -kind=byte -ideclient -ieee=full -free -D__NAG__"
    # fi
    # export LD='-ideclient -unsharedrts'
    # export NCROOT='/usr/local/netcdf-fortran-4.4.5-nagfor'

    # All compilers

    # All compilers
    export NCCROOT='/home/oqx29/zzy20/local'
    export NCCLIB=${NCCROOT}'/lib'
    export NCLIB=${NCROOT}'/lib'
    export NCMOD=${NCROOT}'/include'
    export LDFLAGS="-L${NCCLIB} -L${NCLIB} -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lsz -lz"
    export dosvn=0
    export MFLAGS='-j 8'
    build_build
    cd ../
    build_status
}


## unknown machine, user entering options stdout 
host_read()
{
   print "\n\tWhat is the ROOT path of your NetCDF library" \
         "and .mod file. "
   print "\tRemember these have to be created by the same " \
         "Fortran compiler you" 
   print "\twant to use to build CABLE. e.g./usr/local/intel"
   read NCDF_ROOT
   
   print "\n\tWhat is the path, relative to the above ROOT, of " \
         "your NetCDF library." 
   print "\n\tPress enter for default [lib]."
   read NCDF_DIR
   if [[ $NCDF_DIR == '' ]]; then
      export NCDIR=$NCDF_ROOT/'lib'
   else   
      export NCDIR=$NCDF_ROOT/$NCDF_DIR
   fi

   print "\n\tWhat is the path, relative to the above ROOT, of " \
         "your NetCDF .mod file."
   print "\n\tPress enter for default [include]."
   read NCDF_MOD
   if [[ $NCDF_MOD == '' ]]; then
      export NCMOD=$NCDF_ROOT/'include'
   else   
      export NCMOD=$NCDF_ROOT/$NCDF_MOD
   fi

   print "\n\tWhat is the Fortran compiler you wish to use."
   print "\te.g. ifort, gfortran, nagfor"
   
   print "\n\tPress enter for default [ifort]."
   read FCRESPONSE 
   if [[ $FCRESPONSE == '' ]]; then
      export FC='ifort'
   else   
      export FC=$FCRESPONSE
   fi

   print "\n\tWhat are the approriate compiler options"
   print "\te.g.(ifort) -O2 -fp-model precise "
   print "\n\tPress enter for default [-O2 -fp-model precise]."
   read CFLAGRESPONSE 
   if [[ $CFLAGRESPONSE == '' ]]; then
      export CFLAGS='-O2 -fp-model precise'
   else   
      export CFLAGS=$CFLAGRESPONSE
   fi

   iflags='-L'$NCDIR' -O2'
   export LDFLAGS=$iflags

   print "\n\tWhat are the approriate libraries to link"
   print "\te.g.(most systems) -lnetcdf "
   print "\n\tPress enter for default [-lnetcdf]."
   read LDRESPONSE 
   if [[ $LDRESPONSE == '' ]]; then
      export LD='-lnetcdf'
   else   
      export LD=$LDRESPONSE
   fi
}


host_write()
{
   print '#!/bin/ksh' > junk
   print '' >> junk
   print 'known_hosts()' >> junk
   print '{' >> junk
   print '   set -A kh' ${kh[*]} $HOST_MACH >> junk
   print '}' >> junk
   print '' >> junk
   print '' >> junk
   print '## '$HOST_COMM >> junk
   print 'host_'$HOST_MACH'()' >> junk
   print '{' >> junk
   print '   export NCDIR='"'"$NCDIR"'" >> junk
   print '   export NCMOD='"'"$NCMOD"'" >> junk
   print '   export FC='$FC >> junk
   print '   export CFLAGS='"'"$CFLAGS"'" >> junk
   print '   export LD='"'"$LD"'" >> junk
   print '   export LDFLAGS='"'"$LDFLAGS"'" >> junk
   print '   build_build' >> junk
   print '   cd ../' >> junk
   print '   build_status' >> junk
   print '}' >> junk
   print '' >> junk
   print '' >> junk
}


clean_build()
{
      print '\ncleaning up\n'
      print '\n\tPress Enter too continue buiding, Control-C to abort now.\n'
      read dummy 
      rm -fr .mpitmp*
}


set_up_CABLE_AUX()
{
      print "\n\tYou do not have a ~/CABLE-AUX/ directory. This directory"
      print "\tcontains configuration and data essential to using CABLE."
      print "\tNCI account holders can have this set up for you now (anywhere)."
      print "\tOthers will have to use the tarball available for download at ..."
      print "\n\tDo you want to run set up this directory now? y/[n]"
      print "\n\t B Y P A S S E D by LN"
      #read setup_CABLE_AUX
      setup_CABLE_AUX='n'
      if [[ $setup_CABLE_AUX = 'y' ]]; then
         print "\n\tPlease enter your NCI user ID"
         read NCI_USERID 
         mkdir ~/CABLE-AUX 
         
         fscp1="scp -r "
         fscp2="@vayu.nci.org.au:/projects/access/CABLE-AUX/"
         fscp3="offline "
         fscp4=$HOME"/CABLE-AUX/"
         fscp5=$fscp1$NCI_USERID$fscp2
         fscp=$fscp5$fscp3$fscp4$fscp3
         $fscp
          
         RC=$?
         if [[ $RC > 0 ]];then 
            print "ERROR: scp of ~/CABLE-AUX/offline failed" 
            exit $RC 
         fi
         
         fscp3="core "
         fscp=$fscp5$fscp3$fscp4$fscp3
         $fscp
         
         RC=$?
         if [[ $RC > 0 ]];then 
            print "ERROR: scp of ~/CABLE-AUX/core failed" 
            exit $RC 
         fi
      fi        
}



not_recognized()
{  
   print "\n\n\tThis is not a recognized host for which we " \
         "know the location of the" 
   print "\tnetcdf distribution and correct compiler switches."

   print "\n\tPlease enter these details as prompted, and the " \
         "script will be " 
   print "\tupdated accordingly. " 
   print "\n\tIf this is a common machine for CABLE users, " \
         "please email"
   print "\n\t\t cable_help@nf.nci.org.au "  
   print "\n\talong with your new build_mpi.ksh so that we can " \
         "update the script "
   print "\tfor all users. "
   print "\n\tTo enter compile options for this build press " \
         "enter, otherwise " 
   print "\tControl-C to abort script."           
   
   host_read

   print "\n\tPlease supply a comment include the new build " \
         "script." 
   print "\n\tGenerally the host URL e.g. raijin.nci.org.au "
   read HOST_COMM
   
   build_build
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
   
   while [[ $k -lt $kmax ]]; do
      if [[ $HOST_MACH = ${kh[$k]} ]];then
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

      
i_do_now()
{
      cd ../
      host_write
      tail -n +7 build_mpi.ksh > build_mpi.ksh.tmp
      cat junk build_mpi.ksh.tmp > build_mpi.ksh.new
      mv build_mpi.ksh.new build_mpi.ksh
      chmod u+x build_mpi.ksh 
      rm -f build_mpi.ksh.tmp build_mpi.ksh.new junk 
      build_status
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

if [[ $1 = 'clean' ]]; then
    clean_build
    shift 1
fi

known_hosts
HOST_MACH=`uname -n | cut -c 1-4 | tr - _`
do_i_no_u $*

# only ksh because host_write writes ksh script
not_recognized
i_do_now
