#!/bin/ksh

export dosvn=1 # 1/0: do/do not check svn

known_hosts()
{
   set -A kh vayu cher pear shin jigg nXXX raij ces2 ccrc mael valh
}

## 
host_valh()
{
   # GFORTRAN
   export FC=gfortran
   # debug
  
   export CFLAGS=" -Wall -W -O -g -Wno-maybe-uninitialized -cpp -ffree-form -ffixed-line-length-132 -Wno-tabs -Wconversion" # -DVanessas_common"
   # # release
   # export CFLAGS="-O3 -Wno-aggressive-loop-optimizations -cpp -ffree-form -ffixed-line-length-132"
   export LD=''
   export NCROOT='/usr/local/netcdf-fortran-4.4.2-gfortran/'

   # # NAG
   # export FC=nagfor
   # # debug
   # export CFLAGS="-C -C=dangling -g -nan -O0 -strict95 -gline -fpp -colour -unsharedf95 -kind=byte -ideclient -ieee=full -free -DNAG"
   # # # release
   # # export CFLAGS="-O4 -fpp -colour -unsharedf95 -kind=byte -ideclient -ieee=full -free"
   # export LD='-ideclient -unsharedrts'
   # export NCROOT='/usr/local/netcdf-fortran-4.4.2-nagfor'

   # All compilers
   export NCCROOT='/usr/local'
   export NCCLIB=${NCROOT}'/lib'
   export NCLIB=${NCROOT}'/lib'
   export NCMOD=${NCROOT}'/include'
   export LDFLAGS="-L${NCCLIB} -L${NCLIB} -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lsz -lz"
   export dosvn=0
   build_build
   cd ../
   build_status
}
host_mael()
{
   export NCDIR='/share/apps/netcdf/intel/4.1.3/lib'
   export NCMOD='/share/apps/netcdf/intel/4.1.3/include'
   export FC=ifort
   export CFLAGS='-O2 -fp-model precise  '
   #export CFLAGS='-O3 -fp-model precise  -ipo --parallel '   
   export LDFLAGS='-L/share/apps/intel/Composer/lib/intel64 -L/share/apps/netcdf/intel/4.1.3/lib  -O2'
   if [[ $1 = 'debug' ]]; then
      export CFLAGS='-O0 -traceback -g -fp-model precise  ' 
      export LDFLAGS='-L/share/apps/intel/Composer/lib/intel64 -L/share/apps/netcdf/intel/4.1.3/lib '
   fi
   export LD='-lnetcdf -lnetcdff'
   build_build
   cd ../
   build_status
}



host_ccrc()
{
   export NCDIR='/usr/local/netcdf/intel/4.1.3/lib'
   export NCMOD='/usr/local/netcdf/intel/4.1.3/include'
   export FC=ifort
   export CFLAGS='-O2 -fp-model precise '   #-traceback
   if [[ $1 = 'debug' ]]; then
      export CFLAGS='-O0 -debug -g -ftrapuv -CB -check bounds -diag-enable warn'
# -diag-enable sc2 -diag-enable sc-single-file
   fi
   export LD='-lnetcdf -lnetcdff'
   export LDFLAGS='-L/usr/local/intel/Compiler/11.1/lib/intel64 -L//usr/local/netcdf/intel/4.1.3/lib -O2'
   if [[ $1 = 'debug' ]]; then
      export LDFLAGS='-L/usr/local/intel/Compiler/11.1/lib/intel64 -L//usr/local/netcdf/intel/4.1.3/lib -O0 -debug -g -ftrapuv -diag-enable warn'
# -diag-enable sc2 -diag-enable sc-single-file
   fi
   build_build
   cd ../
   build_status
}


## 
host_ces2()
{
   # GFORTRAN
   export FC=gfortran
   # debug
   export CFLAGS="-pedantic-errors -Wall -W -O -g -Wno-maybe-uninitialized -cpp -ffree-form -ffixed-line-length-132"
   # # release
   # export CFLAGS="-O3 -Wno-aggressive-loop-optimizations -cpp -ffree-form -ffixed-line-length-132"
   export LD=''
   export NCROOT='/usr/local/netcdf-fortran-4.4.1-gfortran'

   # # NAG
   # export FC=nagfor
   # # debug
   # export CFLAGS="-C -C=dangling -g -nan -O0 -strict95 -gline -fpp -colour -unsharedf95 -kind=byte -ideclient -ieee=full -free -DNAG"
   # # # release
   # # export CFLAGS="-O4 -fpp -colour -unsharedf95 -kind=byte -ideclient -ieee=full -free"
   # export LD='-ideclient -unsharedrts'
   # export NCROOT='/usr/local/netcdf-fortran-4.4.1-nagfor'

   # All compilers
   export NCCROOT='/usr/local'
   export NCCLIB=${NCROOT}'/lib'
   export NCLIB=${NCROOT}'/lib'
   export NCMOD=${NCROOT}'/include'
   export LDFLAGS="-L${NCCLIB} -L${NCLIB} -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lsz -lz"
   export dosvn=0
   build_build
   cd ../
   build_status
}




## Interactive Job nXXX@burnet.hpsc.csiro.au  
host_nXXX()
{
   export NCDIR=$NETCDF_ROOT'/lib/'
   export NCMOD=$NETCDF_ROOT'/include/'
   #export FC=$F90
   export FC=ifort
   #vanessa's test options
#   export CFLAGS='  -g -debug -traceback -fp-stack-check -O0 -debug -fpe=0 -fpe-all=0 -no-ftz -ftrapuv'
# export CFLAGS='-warn all,nounused  -check all,noarg_temp_created -g -debug -traceback -fp-stack-check -O0 -debug -fpe1 -no-ftz -ftrapuv'
#   export CFLAGS='-O0 -fp-model precise -debug all -g  '
   export CFLAGS='-O2 -fp-model precise'
   export LDFLAGS='-L'$NCDIR' -O2'
   export LD='-lnetcdf -lnetcdff'
   build_build
   cd ../
   build_status
}

## jiggle
host_jigg()
{
   export NCDIR='/usr/local/lib'
   export NCMOD='/usr/local/include'
   export FC=gfortran
   export CFLAGS='-O2 -x f95-cpp-input'
   export LD='-lnetcdf -lnetcdff'
   export LDFLAGS='-L/usr/local/lib -O2'
   build_build
   cd ../
   build_status
}



## pearcey.hpsc.csiro.au 
host_pear()
{
   . /apps/modules/Modules/default/init/ksh
 
   module del intel-cc intel-fc
   module add intel-cc/16.0.1.150 intel-fc/16.0.1.150
   module add netcdf/4.3.3.1

   export NCDIR=$NETCDF_ROOT'/lib/'
   export NCMOD=$NETCDF_ROOT'/include/'
   export FC='ifort'
   export CFLAGS='-O0 -fp-model precise -g -debug -traceback -fpe=0 -fp-stack-check -no-ftz -ftrapuv   -check all,noarg_temp_created -C '
   #export CFLAGS='-O2 -fp-model precise'
   #export CFLAGS='  -g -debug -traceback -fp-stack-check -O0 -debug -fpe=0 -fpe-all=0 -no-ftz -ftrapuv -check bounds'
   export LDFLAGS='-g -L'$NCDIR' -O0'
   export LD='-lnetcdf -lnetcdff'
   build_build
   cd ../
   build_status
}



## raijin.nci.org.au
host_raij()
{
   module del intel-cc intel-fc
   module add intel-cc/16.0.1.150 intel-fc/16.0.1.150
   module add netcdf/4.3.3.1

   export NCDIR=$NETCDF_ROOT'/lib/Intel'
   export NCMOD=$NETCDF_ROOT'/include/Intel'
   export FC=$F90
   export CFLAGS='-O0 -fp-model precise'
   if [[ $1 = 'debug' ]]; then
      export CFLAGS='-O0 -traceback -g -fp-model precise -ftz -fpe0'
   fi
   export LDFLAGS='-L'$NCDIR' -O0'
   export LD='-lnetcdf -lnetcdff'
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
   print "\te.g. ifort, gfortran"
   
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
   print '   set -A kh' ${kh[*]} $HOST_MACH >> junk ccrc ccrc Mart Mart Mart Mart Mart Mart Mart Mart Mart Mart Mart Mart Mart Mart Mart Mart Mart ccrc ccrc Mart ccrc ccrc ccrc ccrc ccrc ccrc ccrc ccrc
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
      rm -fr .tmp
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
   print "\n\talong with your new build.ksh so that we can " \
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
   integer kmax=${#kh[*]}
   integer k=0
   typeset -f subr

   # for specific nodes on burnet
   ic=`echo $HOST_MACH | cut -c 1`
   in=`echo $HOST_MACH | cut -c 2-4`
   if [[ $ic == 'n' ]]; then
       if [ $in -gt 0 -a $in -lt 1000 ]; then
	   HOST_MACH=nXXX
       fi
   fi
   while [[ $k -lt $kmax ]]; do
      if [[ $HOST_MACH = ${kh[$k]} ]];then
         print 'Host recognized as' $HOST_MACH
         subr=host_${kh[$k]}
         $subr $1
      fi        
      (( k = k + 1 ))
   done 
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


      
i_do_now()
{
      cd ../
      host_write
      tail -n +7 build.ksh > build.ksh.tmp
      cat junk build.ksh.tmp > build.ksh.new
      mv build.ksh.new build.ksh
      chmod u+x build.ksh 
      rm -f build.ksh.tmp build.ksh.new junk 
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
   
   /bin/cp -p $PHYS/*90 ./.tmp
   /bin/cp -p $UTIL/*90 ./.tmp
   /bin/cp -p $DRV/*90 ./.tmp
   /bin/cp -p $CASA/*90 ./.tmp
   
   print "\n\n\tPlease note: CASA-CNP files are included in build only for " 
   print "\ttechnical reasons. Implementation is not officially available with" 
   print "\tthe release of CABLE 2.0\n"
    
   /bin/cp -p Makefile_offline  ./.tmp
   
   cd .tmp/
   make -f Makefile_offline
 
}

###########################################
## build.ksh - MAIN SCRIPT STARTS HERE   ##
###########################################

if [[ $1 = 'clean' ]]; then
   clean_build
fi

   
known_hosts

HOST_MACH=`uname -n | cut -c 1-4`

do_i_no_u $1

not_recognized

i_do_now

