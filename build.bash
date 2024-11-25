#!/usr/bin/env bash

# Exit immediately on error
set -e

ncpus_default=4

script_name=$(basename "${0}")

show_help() {
    cat << EOF
Usage: ./${script_name} [OPTIONS]

Build script wrapper around CMake. Supplied arguments that do not match the
options below will be passed to CMake when generating the build system.

Options:
  -c, --clean   Delete build directory before invoking CMake.
  -m, --mpi     Compile MPI executable.
  -C, --compiler <compiler>
                Specify the compiler to use.
  -n, --ncpus <ncpus>
                Specify the number of parallel jobs in the compilation. By
                default this value is set to ${ncpus_default}.
  -l, --library
                Build just CABLE science library (libscable_science.a)
  -h, --help    Show this screen.

Enabling debug mode:

  The release build is default. To enable debug mode, specify the CMake option
  -DCMAKE_BUILD_TYPE=Debug when invoking ${script_name}.

Enabling verbose output from Makefile builds:

  To enable more verbose output from Makefile builds, specify the CMake option
  -DCMAKE_VERBOSE_MAKEFILE=ON when invoking ${script_name}.

EOF
}

# DEFAULTS

# Configure
cmake_args=(-DCMAKE_BUILD_TYPE=Release -DCABLE_MPI=OFF)

# Build
build_args=()

# Install
do_install=1

# Argument parsing adapted and stolen from http://mywiki.wooledge.org/BashFAQ/035#Complex_nonstandard_add-on_utilities
while [ ${#} -gt 0 ]; do
    case ${1} in
        -c|--clean)
            rm -rf build bin
            exit
            ;;
        -m|--mpi)
            mpi=1
            cmake_args+=(-DCABLE_MPI="ON")
            ;;
        -l|--library)
            build_args+=(--target cable_science)
            do_install=0 # Disable installation when only building the science library
            ;;
        -C|--compiler)
            compiler=${2}
            shift
            ;;
        -n|--ncpus)
            CMAKE_BUILD_PARALLEL_LEVEL=${2}
            shift
            ;;
        -h|--help)
            show_help
            exit
            ;;
        ?*)
            cmake_args+=("${1}")
            ;;
    esac
    shift
done

if hostname -f | grep gadi.nci.org.au > /dev/null; then
    : "${compiler:=intel}"

    . /etc/bashrc
    module purge
    module add cmake/3.24.2
    module add netcdf/4.6.3
    case ${compiler} in
        intel)
            module add intel-compiler/2019.5.281
            compiler_lib_install_dir=Intel
            [[ -n ${mpi} ]] && module add intel-mpi/2019.5.281
            ;;
        gnu)
            module add gcc/13.2.0
            compiler_lib_install_dir=GNU
            [[ -n ${mpi} ]] && module add openmpi/4.1.4
            ;;
        ?*)
            echo -e "\nError: compiler ${compiler} is not supported.\n"
            exit 1
            ;;
    esac

    # This is required so that the netcdf-fortran library is discoverable by
    # pkg-config:
    prepend_path PKG_CONFIG_PATH "${NETCDF_BASE}/lib/${compiler_lib_install_dir}/pkgconfig"

    if module is-loaded openmpi; then
        # This is required so that the openmpi MPI libraries are discoverable
        # via CMake's `find_package` mechanism:
        prepend_path CMAKE_PREFIX_PATH "${OPENMPI_BASE}/include/${compiler_lib_install_dir}"
    fi

elif hostname -f | grep -E '(mc16|mcmini)' > /dev/null; then
    : "${compiler:=gnu}"

    case ${compiler} in
        gnu)
            export PKG_CONFIG_PATH=/usr/local/netcdf-fortran-4.6.1-gfortran/lib/pkgconfig:${PKG_CONFIG_PATH}
            export PKG_CONFIG_PATH=${HOMEBREW_PREFIX}/lib/pkgconfig:${PKG_CONFIG_PATH}
            cmake_args+=(-DCMAKE_Fortran_COMPILER=gfortran)
            ;;
        ?*)
            echo -e "\nError: compiler ${compiler} is not supported.\n"
            exit 1
            ;;
    esac
fi

export CMAKE_BUILD_PARALLEL_LEVEL="${CMAKE_BUILD_PARALLEL_LEVEL:=${ncpus_default}}"

# Configure by default
cmake -S . -B build "${cmake_args[@]}"

# Build requested targets
cmake --build build "${build_args[@]}"

# Install if requested
if [ $do_install -eq 1 ]; then
    cmake --install build --prefix .
fi