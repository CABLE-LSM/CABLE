#!/usr/bin/env bash

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
  -h, --help    Show this screen.

Enabling debug mode:

  The release build is default. To enable debug mode, specify the CMake option
  -DCMAKE_BUILD_TYPE=Debug when invoking ${script_name}.

Enabling verbose output from Makefile builds:

  To enable more verbose output from Makefile builds, specify the CMake option
  -DCMAKE_VERBOSE_MAKEFILE=ON when invoking ${script_name}.

EOF
}

cmake_args=(-DCMAKE_BUILD_TYPE=Release -DCABLE_MPI=OFF)

# Argument parsing adapted and stolen from http://mywiki.wooledge.org/BashFAQ/035#Complex_nonstandard_add-on_utilities
while [ ${#} -gt 0 ]; do
    case ${1} in
        -c|--clean)
	    if [[ -d bin ]] ; then rm -r bin ; fi
	    if [[ -d build ]] ; then rm -r build ; fi
            exit
            ;;
        -m|--mpi)
            mpi=1
            cmake_args+=(-DCABLE_MPI="ON")
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

if hostname -f | grep gadi.nci.org.au > /dev/null ; then
    : "${compiler:=intel}"
    # Gadi at NCI, Australia

    . /etc/bashrc
    module purge
    module add cmake/3.24.2
    module add netcdf/4.9.2
    case ${compiler} in
        intel)
            module add intel-compiler-llvm/2023.0.0
            compiler_lib_install_dir=Intel
            [[ -n ${mpi} ]] && module add intel-mpi/2019.5.281

            # This is required to use ifort instead of ifx (default):
            cmake_args+=(-DCMAKE_Fortran_COMPILER="ifort")

            # This is required to build with runtime dispatch for all
            # architectures on Gadi (see
            # https://opus.nci.org.au/display/Help/Sapphire+Rapids+Compute+Nodes
            # for more information):
            cmake_args+=(-DCMAKE_Fortran_FLAGS_RELEASE_INIT="-march=broadwell -axSKYLAKE-AVX512,CASCADELAKE,SAPPHIRERAPIDS -diag-disable=15009")
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

elif hostname -f | grep -E '(mc16|mcmini)' > /dev/null ; then
    : "${compiler:=gnu}"
    # macOS, netcd-c with homebrew and local netcdf-fortran

    case ${compiler} in
        gnu)
            export PKG_CONFIG_PATH=/usr/local/netcdf-fortran-4.6.2-gfortran/lib/pkgconfig:${PKG_CONFIG_PATH}
            export PKG_CONFIG_PATH=${HOMEBREW_PREFIX}/lib/pkgconfig:${PKG_CONFIG_PATH}
            cmake_args+=(-DCMAKE_Fortran_COMPILER=gfortran)
            ;;
        ?*)
            echo -e "\nError: compiler ${compiler} is not supported.\n"
            exit 1
            ;;
    esac

elif hostname | grep biocomp > /dev/null ; then
    : "${compiler:=gnu}"
    # Linux cluster in Nancy, France, with miniconda

    eval "$(${HOME}/miniconda3/bin/conda 'shell.bash' 'hook' 2> /dev/null)"
    conda activate pystd
    case ${compiler} in
        gnu)
            export PKG_CONFIG_PATH=${HOME}/miniconda3/envs/pystd/lib/pkgconfig:${PKG_CONFIG_PATH}
            cmake_args+=(-DCMAKE_Fortran_COMPILER=gfortran)
            ;;
        ?*)
            echo -e "\nError: compiler ${compiler} is not supported yet.\n"
            exit 1
            ;;
    esac

elif uname -r | grep microsoft-standard-WSL2 > /dev/null ; then
    : "${compiler:=gnu}"
    # Ubuntu on WSL2 of Windows

    case ${compiler} in
        gnu)
            export PKG_CONFIG_PATH=/usr/lib/x86_64-linux-gnu/pkgconfig:${PKG_CONFIG_PATH}
            # pkg-config --cflags netcdf-fortran
            # only gives include path for hdf5 -> set include path by hand
            cmake_args+=(-DCMAKE_Fortran_FLAGS_RELEASE_INIT="-I/usr/include")
            cmake_args+=(-DCMAKE_Fortran_FLAGS_DEBUG_INIT="-I/usr/include")
            cmake_args+=(-DCMAKE_Fortran_COMPILER=gfortran)
            ;;
        ?*)
            echo -e "\nError: compiler ${compiler} is not supported yet.\n"
            exit 1
            ;;
    esac
fi

export CMAKE_BUILD_PARALLEL_LEVEL="${CMAKE_BUILD_PARALLEL_LEVEL:=${ncpus_default}}"

cmake -S . -B build "${cmake_args[@]}" &&\
cmake --build build &&\
cmake --install build --prefix .
