#!/bin/bash

nproc_default=4

script_name=$(basename "${0}")

show_help() {
    cat << EOF
Usage: ./$script_name [OPTIONS]

Build script wrapper around CMake. Supplied arguments that do not match the
options below will be passed to CMake when generating the build system.

Options:
      --clean   Delete build directory before invoking CMake.
      --mpi     Compile MPI executable.
  -h, --help    Show this screen.

Enabling debug mode:

  The release build is default. To enable debug mode, specify the CMake option
  -DCMAKE_BUILD_TYPE=Debug when invoking $script_name.

Enabling verbose output from Makefile builds:

  To enable more verbose output from Makefile builds, specify the CMake option
  -DCMAKE_VERBOSE_MAKEFILE=ON when invoking $script_name.

Parallel compilation:

  By default, the number of parallel jobs used in the compilation is
  $nproc_default. This value can be overwritten by setting the environment
  variable CMAKE_BUILD_PARALLEL_LEVEL.

EOF
}

cmake_args=(-DCMAKE_BUILD_TYPE=Release)

# Argument parsing adapted and stolen from http://mywiki.wooledge.org/BashFAQ/035#Complex_nonstandard_add-on_utilities
while [ $# -gt 0 ]; do
    case $1 in
        --clean)
            rm -r build
            ;;
        --mpi)
            mpi=1
            cmake_args+=(-DCABLE_MPI="ON")
            ;;
        -h|--help)
            show_help
            exit
            ;;
        ?*)
            cmake_args+=("$1")
            ;;
    esac
    shift
done

if hostname -f | grep gadi.nci.org.au > /dev/null; then
    . /etc/bashrc
    module purge
    module add cmake/3.24.2
    module add intel-compiler/2019.5.281
    module add netcdf/4.6.3
    # This is required so that the netcdf-fortran library is discoverable by
    # pkg-config:
    prepend_path PKG_CONFIG_PATH "${NETCDF_BASE}/lib/Intel/pkgconfig"

    if [[ -n $mpi ]]; then
        module add intel-mpi/2019.5.281
    fi

    if module is-loaded openmpi; then
        # This is required so that the openmpi MPI libraries are discoverable
        # via CMake's `find_package` mechanism:
        prepend_path CMAKE_PREFIX_PATH "${OPENMPI_BASE}/include/Intel"
    fi
fi

export CMAKE_BUILD_PARALLEL_LEVEL="${CMAKE_BUILD_PARALLEL_LEVEL:=$nproc_default}"

cmake -S . -B build "${cmake_args[@]}" &&\
cmake --build build &&\
cmake --install build --prefix .
