# Installation instructions

To install CABLE you need to have the following already installed on your system:

- `git`
- `cmake`
- `pkg-config`
- MPI (for parallel applications)
- a Fortran compiler
- the netCDF library

## Flowchart of the installation

``` mermaid
graph TD

    A(Clone CABLE with git):::UserAction -->|Serial?| B(Run `./build.bash`):::UserAction;
    B --> D[load modules and invoke `cmake`];
    D -->|Serial?| E;
    E[executable `cable`]:::Output;
    A -->|Parallel?| C(Run `./build.bash --mpi`):::UserAction;
    C --> D;
    D -->|Parallel?| F;
    F[executable `cable-mpi`]:::Output;
    click A "#getting-the-cable-source-code"
    click B "#launching-the-build"
    click C "#launching-the-build"
    click D "#description-of-the-build-process"
    click E "#description-of-the-build-process"
    click F "#description-of-the-build-process"

    UserAction ---- Automatic ---- Output;

    UserAction(Actions from the user):::UserAction;
    Automatic(Automatic steps of the build script);
    Output(Output of the build script):::Output;

    classDef UserAction fill:#FEFB8E
    classDef Output fill:#cefe8e
```

## Getting the CABLE source code

CABLE can be downloaded from the [GitHub repository][cable-github]. To install the latest version of CABLE, you simply need to clone the repository:

    git clone https://github.com/CABLE-LSM/CABLE.git

!!! tip "Choice of protocol"

    If you have SSH keys installed for GitHub, you can also use the SSH protocol to clone the CABLE repository.

## Building CABLE

In this section, we discuss building CABLE offline for single-site and regional/global applications.

In principle, the only requirements for building and running CABLE are a Fortran compiler and a netCDF distribution. However, in reality it is rarely so simple as this makes it sound.

Whilst we have endeavoured to make CABLE as portable as possible, and indeed we have used CABLE on a variety of platforms, our discussion here is limited to building CABLE on UNIX/Linux platforms only. Specifically, on Gadi@NCI, using Intel Fortran compiler. On other HPC systems there are likely modules available for both Fortran and netCDF. We recommend to get advice from your local system administrators.

See [**Build System**][build-system] for more extensive documentation on the build system.

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

CABLE can be built using the BASH script [build.bash][build.bash] in the project root directory.

???+ tip "Build on other HPC"
    Gadi specific configuration is guarded by a check on the current hostname (see [here][build.bash-hostname-check]). This may be of use as a template for building CABLE on another Linux/HPC system. For advice on issues relating to porting to other HPC systems, please get in touch with ACCESS-NRI either via GitHub or the [ACCESS-Hive Forum][hive-forum-cable].

Executables are built in the `<project_root>/build` directory. Once built successfully, they are then installed in the `<project_root>/bin` directory.

To build the serial model execute:

    ./build.bash

To build the parallel model execute the same build script but with the `--mpi` flag.

    ./build.bash --mpi

### Cleaning the build

From time to time, it might be useful to clean a previous build completely and restart the build from scratch.

To clean the previous build prior to compiling, specify the `--clean` flag to `build.bash`.

    ./build.bash --clean

### Enabling debug mode

The release build is default. To enable debug mode, specify the CMake option `-DCMAKE_BUILD_TYPE=Debug` when invoking `build.bash`.

### Enabling verbose output from Makefile builds

To enable more verbose output from Makefile builds, specify the CMake option `-DCMAKE_VERBOSE_MAKEFILE=ON` when invoking `build.bash`.

???+ tip
    Run `./build.bash --help` for information on supported options.

[cable-github]: https://github.com/CABLE-LSM/cable.git
[NCI]: https://nci.org.au
[registration]: https://trac.nci.org.au/trac/cable/wiki/CABLE_Registration
[build.bash]: https://github.com/CABLE-LSM/CABLE/blob/main/build.bash
[build.bash-hostname-check]: https://github.com/CABLE-LSM/CABLE/blob/main/build.bash#L45-L55
[clean-build]: installation.md/#cleaning-the-build
[build-system]: ../developer_guide/other_resources/build_system.md
[hive-forum-cable]: https://forum.access-hive.org.au/c/land/cable/18