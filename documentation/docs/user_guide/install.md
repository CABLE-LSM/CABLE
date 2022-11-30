# Installation instructions

To install CABLE you need to have the following already installed on your system:

 - SVN
 - a Fortran compiler
 - the netCDF library

## Getting the CABLE source code

CABLE can be downloaded from the [subversion repository][cable-svn] hosted at [NCI][NCI] once a user has requested admission to the CABLE software group at NCI, as described [here][registration]. To install the latest version of CABLE, you simply need to checkout the `trunk`:
    
    svn checkout https://trac.nci.org.au/svn/cable/trunk


## Building CABLE

In this section, we discuss building CABLE offline for single-site and regional/global applications.

In principle, the only requirements for building and running CABLE are a Fortran compiler and a netCDF distribution. However, in reality it is rarely so simple as this makes it sound. 

Whilst we have endeavoured to make CABLE as portable as possible, and indeed we have used CABLE on a variety of platforms, our discussion here is limited to building CABLE on UNIX/Linux platforms only. Specifically, on Gadi@NCI, using Intel Fortran compiler. On other HPC systems there are likely modules available for both Fortran and netCDF. We reccommend to get advice from your local sytem administrators.

CABLE has the directory structure:

    science/
    util/
    params/
    offline/
    coupled/


All applications use the `science/ util/ params/` directories. Offline applications also use the `offline/` directory. Coupled applications use the `coupled/` directory instead. We will not discuss the coupled applications any further here. Please refer to the documentation for a specific coupled application with CABLE to learn more.

### Launching the build

CABLE supports both serial and parallel applications. 

???+ tip "Serial for single-site"
    On a single processor, for single-site investigations with a suitable restart file, CABLE usually executes a year of forcing data in a matter of seconds/minutes. Even in the absence of a restart file, thus requiring spinup of the model, this usually doesn't take more than a few minutes. In this case a serial version of CABLE will suffice. 

???+ tip "Multiprocessor for global simulations"
    For global (or regional) offline simulations, CABLE can still be run in serial mode (about 15 minutes/year for GSWP global run at 1x1 degree resolution). However, running on multiple processors speeds up the simulation considerably.

CABLE can be built using the BASH script [build3.sh][build3] in the **offline/** directory. 

???+ tip "Build on other HPC"
    The function *host_gadi* in the build script shows the appropriate configuration to build on Gadi.
    This may be of use as a template building CABLE on another Linux/HPC system.

Both the serial and parallel executable for CABLE offline are built in the `offline/` directory by executing 
the build script. 

To build the serial model execute: 

    ./build3.sh

To build the parrallel model execute the same build script but with the argumnent **mpi**.   

    ./build3.sh mpi

???+ warning
    If you need to switch between a serial compilation and a parallel compilation, you need to completely [clean the previous build][clean-build] first.
### Description of the build process

The build script:

1. loads modules for the Fortran compiler and the netCDF library
2. sets the compiler flags
3. creates a hidden tempoary directory (`.tmp/`) in which it then compiles CABLE.  
   
A [Makefile][makefile] compiles the CABLE code that is common to both serial and parallel CABLE. This includes all the files under the `science/`, `util/`, and `params/` directories. 

A serial compilation then calls `serial_cable` to compile the serial driver and link everything together to produce the executable, `cable`.

A parallel compilation then calls `parallel_cable` to compile the required MPI driver(s) and link everything together to produce the executable, `cable-mpi`.

The hidden `.tmp/` directory is really only needed for technical reasons and you should never need to go 
into the `.tmp/` directory. However, if for some reason you want the **.o** object files created by the compiler, they persist in this directory. Alternatively, it is sometimes useful to verify that the files you want to compile are actually making it into this directory. This is particuarly relevant if you are adding new files to CABLE. 

???+ warning "Do not change files under `.tmp/`"
    If you change the files in `.tmp/` directly, the next build will not pick up these changes and will overwrite them.

One of the features of the build process is that only source files which are
modified are re-built, followed by their dependents. This is possible because the `.tmp/` directory is
overwritten by the build script, preserving timestamps of the source files from their original location.

### Cleaning the build

From time to time, it might be useful to clean a previous build completely and restart the build from scratch. This is required 
when switching between serial and parallel builds.

To clean the build, you need to run:

    ./build3.sh clean



[cable-svn]: https://trac.nci.org.au/svn/cable
[NCI]: https://nci.org.au
[registration]: https://trac.nci.org.au/trac/cable/wiki/CABLE_Registration
[build3]: https://trac.nci.org.au/svn/cable/trunk/offline/build3.sh
[makefile]: https://trac.nci.org.au/svn/cable/trunk/offline/Makefile
[clean-build]: install.md/#cleaning-the-build