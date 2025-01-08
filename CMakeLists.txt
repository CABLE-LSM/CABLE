# Set minimum CMake version to latest version on Gadi:
cmake_minimum_required(VERSION 3.22.1)

project(
    CABLE
    LANGUAGES Fortran
)

option(CABLE_MPI "Build the MPI executable" OFF)

if(CABLE_MPI)
    # The MPI implementation is currently not supported in CABLE-POP_TRENDY:
    message(FATAL_ERROR "MPI compilation not supported.")
endif()

# third party libs
if(CABLE_MPI)
    find_package(MPI REQUIRED COMPONENTS Fortran)
endif()
find_package(PkgConfig REQUIRED)
pkg_check_modules(NETCDF REQUIRED IMPORTED_TARGET "netcdf-fortran")

set(CABLE_Intel_Fortran_FLAGS -fp-model precise -fpp -nofixed -assume byterecl -stand f08 -D__INTEL__ -D__INTEL_COMPILER__ -D__NETCDF3__)
set(CABLE_Intel_Fortran_FLAGS_DEBUG -O0 -g -traceback -fpe0 -debug extended -check all -warn all -fp-stack-check)
set(CABLE_Intel_Fortran_FLAGS_RELEASE -O3 -ip)
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    add_compile_options(
        ${CABLE_Intel_Fortran_FLAGS}
        "$<$<CONFIG:Release>:${CABLE_Intel_Fortran_FLAGS_RELEASE}>"
        "$<$<CONFIG:Debug>:${CABLE_Intel_Fortran_FLAGS_DEBUG}>"
    )
endif()

set(CABLE_GNU_Fortran_FLAGS -cpp -ffree-form -ffixed-line-length-132 -std=f2008 -D__GFORTRAN__ -D__gFortran__ -D__NETCDF3__)
set(CABLE_GNU_Fortran_FLAGS_DEBUG -O -g -pedantic-errors -Wall -W -Wno-maybe-uninitialized -fbacktrace -ffpe-trap=zero,overflow,underflow -finit-real=nan)
set(CABLE_GNU_Fortran_FLAGS_RELEASE -O3 -Wno-aggressive-loop-optimizations)
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    add_compile_options(
        ${CABLE_GNU_Fortran_FLAGS}
        "$<$<CONFIG:Release>:${CABLE_GNU_Fortran_FLAGS_RELEASE}>"
        "$<$<CONFIG:Debug>:${CABLE_GNU_Fortran_FLAGS_DEBUG}>"
    )
endif()

add_library(
    cable_common_objlib
    OBJECT
    core/biogeochem/cable_adjustJVgm.F90
    core/biogeochem/cable_optimiseJVratio.F90
    core/biogeochem/cable_phenology.F90
    core/biogeochem/casa_cable.F90
    core/biogeochem/casa_cnp.F90
    core/biogeochem/casa_inout.F90
    core/biogeochem/CASAONLY_LUC.F90
    core/biogeochem/casa_variable.F90
    core/biogeochem/POP.F90
    core/biogeochem/pop_io.F90
    core/biogeochem/POPLUC.F90
    core/biogeochem/spincasacnp.F90
    core/biogeophys/cable_air.F90
    core/biogeophys/cable_albedo.F90
    core/biogeophys/cable_c13o2_def.F90
    core/biogeophys/cable_c13o2.F90
    core/biogeophys/cable_canopy.F90
    core/biogeophys/cable_carbon.F90
    core/biogeophys/cable_cbm.F90
    core/biogeophys/cable_climate.F90
    core/biogeophys/cable_common.F90
    core/biogeophys/cable_data.F90
    core/biogeophys/cable_define_types.F90
    core/biogeophys/cable_diag.F90
    core/biogeophys/cable_radiation.F90
    core/biogeophys/cable_roughness.F90
    core/biogeophys/cable_sli_main.F90
    core/biogeophys/cable_sli_numbers.F90
    core/biogeophys/cable_sli_roots.F90
    core/biogeophys/cable_sli_solve.F90
    core/biogeophys/cable_sli_utils.F90
    core/biogeophys/cable_soilsnow.F90
    core/biogeophys/mo_c13o2_photosynthesis.f90
    core/biogeophys/mo_constants.f90
    core/biogeophys/mo_isotope.f90
    core/biogeophys/mo_isotope_luc_model.f90
    core/biogeophys/mo_isotope_pool_model.f90
    core/biogeophys/mo_kind.f90
    core/biogeophys/mo_utils.f90
    core/blaze/blaze_driver.F90
    core/blaze/blaze.F90
    core/blaze/simfire.F90
    core/utils/minpack.f90
    offline/cable_abort.F90
    offline/cable_bios_met_obs_params.F90
    offline/cable_checks.F90
    offline/cable_cru_TRENDY.F90
    offline/cable_initialise.F90
    offline/cable_input.F90
    offline/cable_iovars.F90
    offline/cable_LUC_EXPT.F90
    offline/cable_output.F90
    offline/cable_parameters.F90
    offline/cable_plume_mip.F90
    offline/cable_read.F90
    offline/cable_site.F90
    offline/cable_weathergenerator.F90
    offline/cable_write.F90
)
target_link_libraries(cable_common_objlib PRIVATE PkgConfig::NETCDF)

add_executable(
    cable
    offline/cable_driver.F90
    "$<TARGET_OBJECTS:cable_common_objlib>"
)
target_link_libraries(cable PRIVATE PkgConfig::NETCDF)
install(TARGETS cable RUNTIME)

if(CABLE_MPI)
    add_executable(
        cable-mpi
        src/offline/cable_mpicommon.F90
        src/offline/cable_mpidrv.F90
        src/offline/cable_mpimaster.F90
        src/offline/cable_mpiworker.F90
        src/science/pop/pop_mpi.F90
        "$<TARGET_OBJECTS:cable_common_objlib>"
    )
    target_link_libraries(cable-mpi PRIVATE PkgConfig::NETCDF MPI::MPI_Fortran)
    install(TARGETS cable-mpi RUNTIME)
endif()
