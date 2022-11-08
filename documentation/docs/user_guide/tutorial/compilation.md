# Install and compile
To install cable you need svn.

## To install CABLE

CABLE can be downloaded from the [subverssion repository](https://trac.nci.org.au/svn/cable) hosted at [NCI](https://nci.org.au) once a user has requested admission to the CABLE software group at NCI, as described [here](https://trac.nci.org.au/trac/cable/wiki/CABLE_Registration). We recommend that within the repository you make a copy of the [trunk](https://trac.nci.org.au/svn/cable/trunk) to your own user branch, and then checkout CABLE from there. 

Example user **NCI_ID** creating **MyFirstCABLECheckout** 
    
    svn copy https://trac.nci.org.au/svn/cable/trunk https://trac.nci.org.au/svn/cable/Users/NCI_ID/MyFirstCABLECheckout
    
    svn checkout https://trac.nci.org.au/svn/cable/Users/NCI_ID//MyFirstCABLECheckout




## Building CABLE

In this section we discuss building CABLE offline for single-site and regional/global applications.

In principle, the only requirements for building and running CABLE are a Fortran compiler and a netcdf distribution. However, in reality it is rarely so simple as this makes it sound. 

Whilst we have endeavoured to make CABLE as portable as possible, and indeed we have used CABLE on a variety of platforms, our discussion here is limited to building CABLE on UNIX/Linux platforms only. Specifically, on Gadi@NCI, using Intel fortran. On other HPC systems there are likely modules available for both fortran and netcdf. We reccommend to get advice from your local sytem administrators.

CABLE has the directory structure:

    science/
    util/
    params/
    offline/
    coupled/


All applications use **science/ util/ params/** directories. Offline applications also use the **offline/** directory. Online applications use the **coupled/** directory inatead, however shall not be discussed further here.

CABLE supports both serial and parallel application. 

For single-site investigations, assuming you have a suitable restart file, CABLE usually executes a year of forcing data done in a matter of seconds/minutes (on a single processor). Even in the absence of a restart file, thus requiring spinup of the model, this usually doesn't take more than a few minutes. In this case a serial version of CABLE will suffice. 

For global (or regional) offline runs, CABLE can still be run in serial mode (about 15 minutes/year for GSWP global run at 1x1 degree resolution). However, running on multiple processors speeds up the simulation considerably.

CABLE can be built using the bash script [build3.sh](https://trac.nci.org.au/svn/cable/trunk/offline/build3.sh) in the **offline/** directory. In particular there is a function in the build script, *host_gadi*, which is the appropriate configuration to build on gadi . This may be of use as a template building CABLE on another Linux/HPC system.

Both the serial and parallel executable for CABLE offline are built in the {{{offline/}}} directory by executing 
the build script. To build serial model execute: 


    ./build3.sh

The build script loads modules for Fortran compilers and NETCDF packages, sets compiler flags, and creates a hidden tempoary directory (.tmp) in which it then compiles CABLE.  A [Makefile](https://trac.nci.org.au/svn/cable/trunk/offline/Makefile) compiles CABLE code that is common to both serial and parallel CABLE. This includes all of the **science/**, **util/**, and **params/** code. A script, **serial_cable**, is then called which compiles the serial driver and links everything together to produce an executable, **cable**.

To build parrallel model execute the same build script but witb the argumnent **mpi**.   

    ./build3.sh mpi

A script, **parallel_cable**, is then called which compiles the required MPI driver(s) and links everything together to produce an executable, **cable-mpi**.

The hidden **.tmp/** directory is really only needed for technical reasons and you should never need to go 
into the .''tmp'' directory. However they persist in this directory, if for some reason you want the **.o** object files created by the compiler. Alternatively, it is sometimes useful to verify that the files you want to compile are actually making it into this directory. This is particuarly relevant if you are adding new you files to CABLE. 

One of the features of the build process is that only source files which are
modified are re-built, followed by their dependents. This is possible because the .''tmp'' directory is
overwritten by the build script, preserving timestamps of the source files from their original location.
If you change the files in **./tmp** directly, those changes will not be picked up by the next build
and perhaps worse still, they will be overwritten and lost.


