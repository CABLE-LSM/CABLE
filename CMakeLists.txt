# Set minimum CMake version to latest version on Gadi:
cmake_minimum_required(VERSION 3.24.2)

project(
    CABLE
    LANGUAGES Fortran
)

option(CABLE_MPI "Build the MPI executable" OFF)
option(CABLE_TESTS "Build CABLE tests" OFF)

# third party libs
if(CABLE_MPI)
    find_package(MPI REQUIRED COMPONENTS Fortran)
    find_package(PIO COMPONENTS Fortran QUIET)
    if(TARGET PIO::PIO_Fortran)
        message(STATUS "Found PIO_Fortran: ${PIO_DIR}")
    endif()
endif()
find_package(PkgConfig REQUIRED)
pkg_check_modules(NETCDF REQUIRED IMPORTED_TARGET "netcdf-fortran")

if(CABLE_TESTS)
    enable_testing()
    include(FetchContent)
    if(CABLE_MPI)
        option(FORTUNO_WITH_MPI "Fortuno: whether to build the MPI interface" ON)
        FetchContent_Declare(
            FortunoMPI
            GIT_REPOSITORY https://github.com/fortuno-repos/fortuno
            GIT_TAG main
        )
        FetchContent_MakeAvailable(FortunoMPI)
        set(fortuno_libs Fortuno::fortuno_mpi)
    else()
        option(FORTUNO_WITH_MPI "Fortuno: whether to build the MPI interface" OFF)
        FetchContent_Declare(Fortuno
            GIT_REPOSITORY https://github.com/fortuno-repos/fortuno
            GIT_TAG main
            FIND_PACKAGE_ARGS CONFIG
        )
        FetchContent_MakeAvailable(Fortuno)
        set(fortuno_libs Fortuno::fortuno_serial)
    endif ()
endif()

set(CABLE_Intel_Fortran_FLAGS -fp-model precise)
set(CABLE_Intel_Fortran_FLAGS_DEBUG -O0 -g -traceback -fpe0)
set(CABLE_Intel_Fortran_FLAGS_RELEASE -O2)
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    add_compile_options(
        ${CABLE_Intel_Fortran_FLAGS}
        "$<$<CONFIG:Release>:${CABLE_Intel_Fortran_FLAGS_RELEASE}>"
        "$<$<CONFIG:Debug>:${CABLE_Intel_Fortran_FLAGS_DEBUG}>"
    )
endif()

set(CABLE_GNU_Fortran_FLAGS -cpp -ffree-form -ffree-line-length-none)
set(CABLE_GNU_Fortran_FLAGS_DEBUG -O -g -pedantic-errors -Wall -W -Wno-maybe-uninitialized -fbacktrace -ffpe-trap=zero,overflow,underflow -finit-real=nan)
set(CABLE_GNU_Fortran_FLAGS_RELEASE -O3 -Wno-aggressive-loop-optimizations)
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    add_compile_options(
        ${CABLE_GNU_Fortran_FLAGS}
        "$<$<CONFIG:Release>:${CABLE_GNU_Fortran_FLAGS_RELEASE}>"
        "$<$<CONFIG:Debug>:${CABLE_GNU_Fortran_FLAGS_DEBUG}>"
    )
endif()

# Set the module directory to assist Spack
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include")
if(CABLE_LIBRARY)

    # Set the core sources
    set(SOURCES
        src/science/albedo/cbl_albedo.F90
        src/science/albedo/cbl_snow_albedo.F90
        src/science/albedo/cbl_soilColour_albedo.F90
        src/science/canopy/cable_canopy.F90
        src/science/canopy/cbl_dryLeaf.F90
        src/science/canopy/cbl_friction_vel.F90
        src/science/canopy/cbl_fwsoil.F90
        src/science/canopy/cbl_init_wetfac_mod.F90
        src/science/canopy/cbl_latent_heat.F90
        src/science/canopy/cbl_photosynthesis.F90
        src/science/canopy/cbl_pot_evap_snow.F90
        src/science/canopy/cbl_qsat.F90
        src/science/canopy/cbl_SurfaceWetness.F90
        src/science/canopy/cbl_wetleaf.F90
        src/science/canopy/cbl_within_canopy.F90
        src/science/canopy/cbl_zetar.F90
        src/science/casa-cnp/bgcdriver.F90
        src/science/casa-cnp/biogeochem_casa.F90
        src/science/casa-cnp/casa_cnp.F90
        src/science/casa-cnp/casa_dimension.F90
        src/science/casa-cnp/casa_feedback.F90
        src/science/casa-cnp/casa_inout.F90
        src/science/casa-cnp/casa_param.F90
        src/science/casa-cnp/casa_phenology.F90
        src/science/casa-cnp/casa_readbiome.F90
        src/science/casa-cnp/casa_rplant.F90
        src/science/casa-cnp/casa_sumcflux.F90
        src/science/casa-cnp/casa_variable.F90
        src/science/landuse/landuse3.F90
        src/science/landuse/landuse_constant.F90
        src/science/misc/cable_air.F90
        src/science/misc/cable_carbon.F90
        src/science/misc/cable_climate.F90
        src/science/pop/pop_constants.F90
        src/science/pop/pop_def.F90
        src/science/pop/pop_types.F90
        src/science/radiation/cbl_init_radiation.F90
        src/science/radiation/cbl_radiation.F90
        src/science/radiation/cbl_rhoch.F90
        src/science/radiation/cbl_sinbet.F90
        src/science/radiation/cbl_spitter.F90
        src/science/roughness/cable_roughness.F90
        src/science/roughness/roughnessHGT_effLAI_cbl.F90
        src/science/soilsnow/cbl_conductivity.F90
        src/science/soilsnow/cbl_GW.F90
        src/science/soilsnow/cbl_hyd_redistrib.F90
        src/science/soilsnow/cbl_Oldconductivity.F90
        src/science/soilsnow/cbl_remove_trans.F90
        src/science/soilsnow/cbl_smoisturev.F90
        src/science/soilsnow/cbl_snowAccum.F90
        src/science/soilsnow/cbl_snow_aging.F90
        src/science/soilsnow/cbl_snowCheck.F90
        src/science/soilsnow/cbl_snowDensity.F90
        src/science/soilsnow/cbl_snowl_adjust.F90
        src/science/soilsnow/cbl_snowMelt.F90
        src/science/soilsnow/cbl_soilfreeze.F90
        src/science/soilsnow/cbl_soilsnow_data.F90
        src/science/soilsnow/cbl_soilsnow_init_special.F90
        src/science/soilsnow/cbl_soilsnow_main.F90
        src/science/soilsnow/cbl_stempv.F90
        src/science/soilsnow/cbl_surfbv.F90
        src/science/soilsnow/cbl_thermal.F90
        src/science/soilsnow/cbl_trimb.F90
        src/params/cable_phys_constants_mod.F90
        src/params/grid_constants_cbl.F90
        src/params/cable_photo_constants_mod.F90
        src/params/cable_other_constants_mod.F90
        src/params/cable_maths_constants_mod.F90
        src/util/cable_runtime_opts_mod.F90
        src/util/cable_common.F90
        src/coupled/shared/cable_canopy_type_mod.F90
        src/coupled/shared/cable_soilsnow_type_mod.F90
        src/coupled/shared/cable_soil_type_mod.F90
        src/coupled/shared/cable_veg_type_mod.F90
        src/util/cable_climate_type_mod.F90
        src/util/masks_cbl.F90
    )

if(CABLE_LIBRARY_TARGET STREQUAL "ESM1.6")
        list(APPEND SOURCES
            src/coupled/esm/cable_iovars.F90
            src/coupled/esm/cable_surface_types.F90
            src/coupled/esm/cable_define_types.F90
            src/coupled/esm/cable_soil_params_mod.F90
            src/coupled/esm/cable_pft_params_mod.F90
            src/coupled/esm/casa_landuse.F90
            src/coupled/esm/LAI_canopy_height_cbl.F90
        )
        add_compile_options(-r8 -i8)
    
    else()
        message(STATUS "Specified library application does not exist")
    endif()
    # Can conditionally include files for other shared apps here
    # e.g. IF (CABLE_LIBRARY_TARGET STREQUAL "AM3")

    add_library(cable STATIC ${SOURCES})
    target_compile_definitions(cable PRIVATE UM_CBL)

    target_include_directories(cable PUBLIC "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")

    message(STATUS "Build the library")
    install(TARGETS cable
      EXPORT CABLEFortranTargets
      LIBRARY DESTINATION lib
      INCLUDES DESTINATION include)

    install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ TYPE INCLUDE)

else()

    add_library(
        cable_common
        STATIC
        src/science/albedo/cbl_albedo.F90
        src/science/albedo/cbl_snow_albedo.F90
        src/science/albedo/cbl_soilColour_albedo.F90
        src/science/canopy/cable_canopy.F90
        src/science/canopy/cbl_dryLeaf.F90
        src/science/canopy/cbl_friction_vel.F90
        src/science/canopy/cbl_fwsoil.F90
        src/science/canopy/cbl_init_wetfac_mod.F90
        src/science/canopy/cbl_latent_heat.F90
        src/science/canopy/cbl_photosynthesis.F90
        src/science/canopy/cbl_pot_evap_snow.F90
        src/science/canopy/cbl_qsat.F90
        src/science/canopy/cbl_SurfaceWetness.F90
        src/science/canopy/cbl_wetleaf.F90
        src/science/canopy/cbl_within_canopy.F90
        src/science/canopy/cbl_zetar.F90
        src/science/casa-cnp/bgcdriver.F90
        src/science/casa-cnp/biogeochem_casa.F90
        src/science/casa-cnp/casa_cnp.F90
        src/science/casa-cnp/casa_dimension.F90
        src/science/casa-cnp/casa_feedback.F90
        src/science/casa-cnp/casa_inout.F90
        src/science/casa-cnp/casa_param.F90
        src/science/casa-cnp/casa_phenology.F90
        src/science/casa-cnp/casa_readbiome.F90
        src/science/casa-cnp/casa_rplant.F90
        src/science/casa-cnp/casa_sumcflux.F90
        src/science/casa-cnp/casa_variable.F90
        src/science/gw_hydro/cable_gw_hydro.F90
        src/science/gw_hydro/cable_psm.F90
        src/science/landuse/landuse3.F90
        src/science/landuse/landuse_constant.F90
        src/science/misc/cable_air.F90
        src/science/misc/cable_carbon.F90
        src/science/misc/cable_climate.F90
        src/science/pop/pop_constants.F90
        src/science/pop/pop_def.F90
        src/science/pop/POP.F90
        src/science/pop/pop_io.F90
        src/science/pop/POPLUC.F90
        src/science/pop/pop_types.F90
        src/science/radiation/cbl_init_radiation.F90
        src/science/radiation/cbl_radiation.F90
        src/science/radiation/cbl_rhoch.F90
        src/science/radiation/cbl_sinbet.F90
        src/science/radiation/cbl_spitter.F90
        src/science/roughness/cable_roughness.F90
        src/science/roughness/roughnessHGT_effLAI_cbl.F90
        src/science/sli/cable_sli_main.F90
        src/science/sli/cable_sli_numbers.F90
        src/science/sli/cable_sli_roots.F90
        src/science/sli/cable_sli_solve.F90
        src/science/sli/cable_sli_utils.F90
        src/science/soilsnow/cbl_conductivity.F90
        src/science/soilsnow/cbl_GW.F90
        src/science/soilsnow/cbl_hyd_redistrib.F90
        src/science/soilsnow/cbl_Oldconductivity.F90
        src/science/soilsnow/cbl_remove_trans.F90
        src/science/soilsnow/cbl_smoisturev.F90
        src/science/soilsnow/cbl_snowAccum.F90
        src/science/soilsnow/cbl_snow_aging.F90
        src/science/soilsnow/cbl_snowCheck.F90
        src/science/soilsnow/cbl_snowDensity.F90
        src/science/soilsnow/cbl_snowl_adjust.F90
        src/science/soilsnow/cbl_snowMelt.F90
        src/science/soilsnow/cbl_soilfreeze.F90
        src/science/soilsnow/cbl_soilsnow_data.F90
        src/science/soilsnow/cbl_soilsnow_init_special.F90
        src/science/soilsnow/cbl_soilsnow_main.F90
        src/science/soilsnow/cbl_stempv.F90
        src/science/soilsnow/cbl_surfbv.F90
        src/science/soilsnow/cbl_thermal.F90
        src/science/soilsnow/cbl_trimb.F90
        src/params/cable_phys_constants_mod.F90
        src/params/grid_constants_cbl.F90
        src/params/cable_photo_constants_mod.F90
        src/params/cable_other_constants_mod.F90
        src/params/cable_maths_constants_mod.F90
        src/util/cable_runtime_opts_mod.F90
        src/util/cable_common.F90
        src/shared/casa_offline_inout.F90
        src/shared/casa_ncdf.F90
        src/offline/cable_io_decomp.F90
        src/offline/cable_iovars.F90
        src/offline/cable_surface_types.F90
        src/offline/cable_define_types.F90
        src/shared/cable_phenology.F90
        src/shared/cable_LUC_EXPT.F90
        src/offline/CASAONLY_LUC.F90
        src/offline/cable_abort.F90
        src/offline/cable_checks.F90
        src/offline/cable_cru_TRENDY.F90
        src/offline/cable_driver_common.F90
        src/offline/cable_initialise.F90
        src/offline/cable_input.F90
        src/offline/cable_metutils.F90
        src/offline/cable_mpi_stub_types.F90
        src/offline/cable_mpi.F90
        src/offline/cable_namelist_input.F90
        src/offline/cable_output.F90
        src/offline/cable_parameters.F90
        src/offline/cable_pft_params.F90
        src/offline/cable_plume_mip.F90
        src/offline/cable_read.F90
        src/offline/cable_site.F90
        src/offline/cable_serial.F90
        src/offline/cable_soil_params.F90
        src/offline/cable_weathergenerator.F90
        src/offline/cable_write.F90
        src/offline/casa_cable.F90
        src/offline/cbl_model_driver_offline.F90
        src/offline/landuse_inout.F90
        src/offline/spincasacnp.F90
        src/util/aggregator.F90
        src/util/cable_climate_type_mod.F90
        src/util/masks_cbl.F90
        src/util/cable_array_utils.F90
        src/util/cable_enum.F90
        src/util/netcdf/cable_netcdf_decomp_util.F90
        src/util/netcdf/cable_netcdf.F90
        src/util/netcdf/cable_netcdf_internal.F90
        src/util/netcdf/cable_netcdf_stub_types.F90
        src/util/netcdf/nf90/cable_netcdf_nf90.F90
    )

    target_link_libraries(cable_common PRIVATE PkgConfig::NETCDF)

    if(CABLE_MPI)
        target_compile_definitions(cable_common PRIVATE __MPI__)
        target_link_libraries(cable_common PRIVATE MPI::MPI_Fortran)
    endif()

    if(TARGET PIO::PIO_Fortran)
        target_link_libraries(cable_common PRIVATE PIO::PIO_Fortran)
        target_sources(cable_common PRIVATE src/util/netcdf/pio/cable_netcdf_pio.F90)
    else()
        target_sources(cable_common PRIVATE src/util/netcdf/pio/cable_netcdf_pio_stub.F90)
    endif()

    if(CABLE_MPI)
        add_executable(
            cable-mpi
            src/offline/cable_mpicommon.F90
            src/offline/cable_mpimaster.F90
            src/offline/cable_mpiworker.F90
            src/science/pop/pop_mpi.F90
            src/offline/cable_offline_driver.F90
        )
        target_link_libraries(cable-mpi PRIVATE cable_common MPI::MPI_Fortran)
        install(TARGETS cable-mpi RUNTIME)
    else()
        add_executable(
            cable
            src/offline/cable_mpimaster_stub.F90
            src/offline/cable_mpiworker_stub.F90
            src/offline/cable_offline_driver.F90
        )
        target_link_libraries(cable PRIVATE cable_common)
        install(TARGETS cable RUNTIME)
    endif()

    if (CABLE_TESTS)
        add_executable(
            cable-tests
            tests/cable_tests.F90
            tests/fixtures.F90
            tests/utils/file_utils.F90
            tests/test_cable_netcdf.F90
        )
        if(CABLE_MPI)
            target_sources(cable-tests PRIVATE tests/fortuno_interface_mpi.f90)
        else()
            target_sources(cable-tests PRIVATE tests/fortuno_interface_serial.f90)
        endif()
        target_link_libraries(cable-tests PRIVATE cable_common ${fortuno_libs})
        if(CABLE_MPI)
            add_test(NAME cable-tests-serial
                COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:cable-tests> ~parallel
            )
            set_tests_properties(cable-tests-serial PROPERTIES PROCESSORS 1)
            add_test(NAME cable-tests-parallel
                COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:cable-tests> parallel
            )
            set_tests_properties(cable-tests-parallel PROPERTIES PROCESSORS 4)
        else()
            add_test(NAME cable-tests-serial
                COMMAND $<TARGET_FILE:cable-tests> ~parallel
            )
        endif()
    endif()
endif()
