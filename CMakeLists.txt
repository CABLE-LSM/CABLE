# Set minimum CMake version to latest version on Gadi:
cmake_minimum_required(VERSION 3.24.2)

project(
    CABLE
    LANGUAGES Fortran
)

option(CABLE_MPI "Build the MPI executable" OFF)

# third party libs
find_package(PkgConfig REQUIRED)
pkg_check_modules(NETCDF REQUIRED IMPORTED_TARGET "netcdf-fortran")

set(CABLE_INTEL_Fortran_FLAGS -fp-model precise)
set(CABLE_INTEL_Fortran_FLAGS_DEBUG -O0 -g -traceback -fpe0)
set(CABLE_INTEL_Fortran_FLAGS_RELEASE -O2)

# CMake $<...> syntax is explained here:
# https://cmake.org/cmake/help/latest/manual/cmake-generator-expressions.7.html#introduction
set(CABLE_EXE_NAME "$<IF:$<BOOL:${CABLE_MPI}>,cable-mpi,cable>")

set(
    CABLE_SRCS_SERIAL
    src/offline/cable_driver.F90
)
set(
    CABLE_SRCS_MPI
    src/offline/cable_mpidrv.F90
    src/offline/cable_mpicommon.F90
    src/offline/cable_mpimaster.F90
    src/offline/cable_mpiworker.F90
    src/science/pop/pop_mpi.F90
)
set(
    CABLE_SRCS_COMMON
    src/science/casa-cnp/bgcdriver.F90
    src/science/casa-cnp/biogeochem_casa.F90
    src/offline/cable_abort.F90
    src/science/misc/cable_air.F90
    src/science/canopy/cable_canopy.F90
    src/science/misc/cable_carbon.F90
    src/offline/cable_checks.F90
    src/science/misc/cable_climate.F90
    src/util/cable_climate_type_mod.F90
    src/util/cable_common.F90
    src/offline/cable_cru_TRENDY.F90
    src/offline/cable_define_types.F90
    src/science/gw_hydro/cable_gw_hydro.F90
    src/offline/cable_initialise.F90
    src/offline/cable_input.F90
    src/offline/cable_iovars.F90
    src/offline/cable_LUC_EXPT.F90
    src/params/cable_maths_constants_mod.F90
    src/offline/cable_metutils.F90
    src/offline/cable_namelist_input.F90
    src/params/cable_other_constants_mod.F90
    src/offline/cable_output.F90
    src/offline/cable_parameters.F90
    src/offline/cable_pft_params.F90
    src/offline/cable_phenology.F90
    src/params/cable_photo_constants_mod.F90
    src/params/cable_phys_constants_mod.F90
    src/offline/cable_plume_mip.F90
    src/science/gw_hydro/cable_psm.F90
    src/offline/cable_read.F90
    src/science/roughness/cable_roughness.F90
    src/util/cable_runtime_opts_mod.F90
    src/offline/cable_site.F90
    src/science/sli/cable_sli_main.F90
    src/science/sli/cable_sli_numbers.F90
    src/science/sli/cable_sli_roots.F90
    src/science/sli/cable_sli_solve.F90
    src/science/sli/cable_sli_utils.F90
    src/offline/cable_soil_params.F90
    src/offline/cable_weathergenerator.F90
    src/offline/cable_write.F90
    src/offline/casa_cable.F90
    src/science/casa-cnp/casa_cnp.F90
    src/science/casa-cnp/casa_dimension.F90
    src/science/casa-cnp/casa_feedback.F90
    src/science/casa-cnp/casa_inout.F90
    src/offline/casa_ncdf.F90
    src/offline/casa_offline_inout.F90
    src/offline/CASAONLY_LUC.F90
    src/science/casa-cnp/casa_param.F90
    src/science/casa-cnp/casa_phenology.F90
    src/science/casa-cnp/casa_readbiome.F90
    src/science/casa-cnp/casa_rplant.F90
    src/science/casa-cnp/casa_sumcflux.F90
    src/science/casa-cnp/casa_variable.F90
    src/science/albedo/cbl_albedo.F90
    src/science/soilsnow/cbl_conductivity.F90
    src/science/canopy/cbl_dryLeaf.F90
    src/science/canopy/cbl_friction_vel.F90
    src/science/canopy/cbl_fwsoil.F90
    src/science/soilsnow/cbl_GW.F90
    src/science/soilsnow/cbl_hyd_redistrib.F90
    src/science/radiation/cbl_init_radiation.F90
    src/science/canopy/cbl_latent_heat.F90
    src/offline/cbl_model_driver_offline.F90
    src/science/soilsnow/cbl_Oldconductivity.F90
    src/science/canopy/cbl_photosynthesis.F90
    src/science/canopy/cbl_pot_evap_snow.F90
    src/science/canopy/cbl_qsat.F90
    src/science/radiation/cbl_radiation.F90
    src/science/soilsnow/cbl_remove_trans.F90
    src/science/radiation/cbl_rhoch.F90
    src/science/radiation/cbl_sinbet.F90
    src/science/soilsnow/cbl_smoisturev.F90
    src/science/soilsnow/cbl_snowAccum.F90
    src/science/soilsnow/cbl_snow_aging.F90
    src/science/albedo/cbl_snow_albedo.F90
    src/science/soilsnow/cbl_snowCheck.F90
    src/science/soilsnow/cbl_snowDensity.F90
    src/science/soilsnow/cbl_snowl_adjust.F90
    src/science/soilsnow/cbl_snowMelt.F90
    src/science/albedo/cbl_soilColour_albedo.F90
    src/science/soilsnow/cbl_soilfreeze.F90
    src/science/soilsnow/cbl_soilsnow_data.F90
    src/science/soilsnow/cbl_soilsnow_init_special.F90
    src/science/soilsnow/cbl_soilsnow_main.F90
    src/science/radiation/cbl_spitter.F90
    src/science/soilsnow/cbl_stempv.F90
    src/science/canopy/cbl_SurfaceWetness.F90
    src/science/soilsnow/cbl_surfbv.F90
    src/science/soilsnow/cbl_thermal.F90
    src/science/soilsnow/cbl_trimb.F90
    src/science/canopy/cbl_wetleaf.F90
    src/science/canopy/cbl_within_canopy.F90
    src/science/canopy/cbl_zetar.F90
    src/params/grid_constants_cbl.F90
    src/science/landuse/landuse3.F90
    src/science/landuse/landuse_constant.F90
    src/offline/landuse_inout.F90
    src/util/masks_cbl.F90
    src/science/pop/pop_constants.F90
    src/science/pop/pop_def.F90
    src/science/pop/pop_io.F90
    src/science/pop/POPLUC.F90
    src/science/pop/POP.F90
    src/science/pop/pop_types.F90
    src/science/roughness/roughnessHGT_effLAI_cbl.F90
    src/offline/spincasacnp.F90
)
set(
    CABLE_SRCS
    "$<IF:$<BOOL:${CABLE_MPI}>,${CABLE_SRCS_MPI},${CABLE_SRCS_SERIAL}>"
    ${CABLE_SRCS_COMMON}
)

add_executable(cable ${CABLE_SRCS})

set_target_properties(cable PROPERTIES OUTPUT_NAME ${CABLE_EXE_NAME})

if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    target_compile_options(
        cable
        PRIVATE
        ${CABLE_INTEL_Fortran_FLAGS}
        "$<$<CONFIG:Release>:${CABLE_INTEL_Fortran_FLAGS_RELEASE}>"
        "$<$<CONFIG:Debug>:${CABLE_INTEL_Fortran_FLAGS_DEBUG}>"
    )
endif()

target_link_libraries(
    cable
    PkgConfig::NETCDF
)

install(TARGETS cable RUNTIME)
