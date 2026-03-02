module cable_diagnostics_core_mod

  use cable_def_types_mod, only: met_type
  use cable_def_types_mod, only: canopy_type
  use cable_def_types_mod, only: soil_parameter_type
  use cable_def_types_mod, only: soil_snow_type
  use cable_def_types_mod, only: radiation_type
  use cable_def_types_mod, only: veg_parameter_type
  use cable_def_types_mod, only: balances_type
  use cable_def_types_mod, only: roughness_type
  use cable_def_types_mod, only: bgc_pool_type
  use cable_def_types_mod, only: mvtype, mstype

  use cable_phys_constants_mod, only: c_molar_mass

  use cable_netcdf_mod, only: CABLE_NETCDF_INT
  use cable_netcdf_mod, only: CABLE_NETCDF_FLOAT

  use aggregator_mod, only: new_aggregator

  use cable_common_module, only: cable_user
  use cable_common_module, only: gw_params

  use cable_io_vars_module, only: output, patchout
  use cable_io_vars_module, only: landpt_global
  use cable_io_vars_module, only: patch

  use cable_output_prototype_v2_mod, only: cable_output_variable_t
  use cable_output_prototype_v2_mod, only: attribute => cable_output_attribute_t
  use cable_output_prototype_v2_mod, only: DIM_PATCH => CABLE_OUTPUT_DIM_PATCH
  use cable_output_prototype_v2_mod, only: DIM_SOIL => CABLE_OUTPUT_DIM_SOIL
  use cable_output_prototype_v2_mod, only: DIM_RAD => CABLE_OUTPUT_DIM_RAD
  use cable_output_prototype_v2_mod, only: DIM_SNOW => CABLE_OUTPUT_DIM_SNOW
  use cable_output_prototype_v2_mod, only: DIM_PLANTCARBON => CABLE_OUTPUT_DIM_PLANTCARBON
  use cable_output_prototype_v2_mod, only: DIM_SOILCARBON => CABLE_OUTPUT_DIM_SOILCARBON
  use cable_output_prototype_v2_mod, only: DIM_LAND_GLOBAL => CABLE_OUTPUT_DIM_LAND_GLOBAL

  use cable_checks_module, only: ranges

  implicit none
  private

  public :: cable_diagnostics_core

contains

  function cable_diagnostics_core(met, canopy, soil, ssnow, rad, veg, bal, rough, bgc, dels) result(output_variables)
    type(met_type), intent(inout) :: met
    type(canopy_type), intent(inout) :: canopy
    type(soil_parameter_type), intent(inout) :: soil
    type(soil_snow_type), intent(inout) :: ssnow
    type(radiation_type), intent(inout) :: rad
    type(veg_parameter_type), intent(inout) :: veg
    type(balances_type), intent(inout) :: bal
    type(roughness_type), intent(inout) :: rough
    type(bgc_pool_type), intent(inout) :: bgc
    real, intent(in) :: dels

    type(cable_output_variable_t), allocatable :: output_variables(:)

    output_variables = [ &
      cable_output_variable_t( &
        field_name="fsd", &
        netcdf_name="SWdown", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%SWdown, &
        active=output%met .or. output%SWdown, &
        patchout=output%patch .or. patchout%SWdown, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(met%ofsd), &
        metadata=[ &
          attribute("units", "W/m^2"), &
          attribute("long_name", "Downward shortwave radiation") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="fld", &
        netcdf_name="LWdown", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%LWdown, &
        active=output%met .or. output%LWdown, &
        patchout=output%patch .or. patchout%LWdown, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(met%fld), &
        metadata=[ &
          attribute("units", "W/m^2"), &
          attribute("long_name", "Downward longwave radiation") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="precip", &
        netcdf_name="Rainf", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%Rainf, &
        active=output%met .or. output%Rainf, &
        patchout=output%patch .or. patchout%Rainf, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(met%precip), &
        scale=(1.0 / dels), &
        metadata=[ &
          attribute("units", "kg/m^2/s"), &
          attribute("long_name", "Rainfall+snowfall") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="precip_sn", &
        netcdf_name="Snowf", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%Snowf, &
        active=output%met .or. output%Snowf, &
        patchout=output%patch .or. patchout%Snowf, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(met%precip_sn), &
        scale=(1.0 / dels), &
        metadata=[ &
          attribute("units", "kg/m^2/s"), &
          attribute("long_name", "Snowfall") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="pmb", &
        netcdf_name="PSurf", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%PSurf, &
        active=output%met .or. output%PSurf, &
        patchout=output%patch .or. patchout%PSurf, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(met%pmb), &
        metadata=[ &
          attribute("units", "hPa"), &
          attribute("long_name", "Surface air pressure") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="tk", &
        netcdf_name="Tair", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%Tair, &
        active=output%met .or. output%Tair, &
        patchout=output%patch .or. patchout%Tair, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(met%tk), &
        metadata=[ &
          attribute("units", "K"), &
          attribute("long_name", "Surface air temperature") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="qv", &
        netcdf_name="Qair", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%Qair, &
        active=output%met .or. output%Qair, &
        patchout=output%patch .or. patchout%Qair, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(met%qv), &
        metadata=[ &
          attribute("units", "kg/kg"), &
          attribute("long_name", "Surface specific humidity") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="ua", &
        netcdf_name="Wind", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%Wind, &
        active=output%met .or. output%Wind, &
        patchout=output%patch .or. patchout%Wind, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(met%ua), &
        metadata=[ &
          attribute("units", "m/s"), &
          attribute("long_name", "Scalar surface wind speed") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="ca", &
        netcdf_name="CO2air", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%CO2air, &
        active=output%met .or. output%CO2air, &
        patchout=output%patch .or. patchout%CO2air, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        scale=1e6, & ! Convert to ppmv
        aggregator=new_aggregator(met%ca), &
        metadata=[ &
          attribute("units", "ppmv"), &
          attribute("long_name", "Surface air CO2 concentration") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="qmom", &
        netcdf_name="Qmom", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%Qmom, &
        active=output%flux .or. output%Qmom, &
        patchout=output%patch .or. patchout%Qmom, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(canopy%qmom), &
        metadata=[ &
          attribute("units", "kg/m/s2"), &
          attribute("long_name", "Surface momentum flux") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="fe", &
        netcdf_name="Qle", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%Qle, &
        active=output%flux .or. output%Qle, &
        patchout=output%patch .or. patchout%Qle, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(canopy%fe), &
        metadata=[ &
          attribute("units", "W/m^2"), &
          attribute("long_name", "Surface latent heat flux") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="fh", &
        netcdf_name="Qh", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%Qh, &
        active=output%flux .or. output%Qh, &
        patchout=output%patch .or. patchout%Qh, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(canopy%fh), &
        metadata=[ &
          attribute("units", "W/m^2"), &
          attribute("long_name", "Surface sensible heat flux") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="ga", &
        netcdf_name="Qg", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%Qg, &
        active=output%flux .or. output%Qg, &
        patchout=output%patch .or. patchout%Qg, &
        restart=.true., &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(canopy%ga), &
        metadata=[ &
          attribute("units", "W/m^2"), &
          attribute("long_name", "Surface ground heat flux") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="rnof1", &
        netcdf_name="Qs", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%Qs, &
        active=output%flux .or. output%Qs, &
        patchout=output%patch .or. patchout%Qs, &
        restart=.true., &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(ssnow%rnof1), &
        scale=(1.0 / dels), &
        metadata=[ &
          attribute("units", "kg/m^2/s"), &
          attribute("long_name", "Surface runoff") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="rnof2", &
        netcdf_name="Qsb", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%Qsb, &
        active=output%flux .or. output%Qsb, &
        patchout=output%patch .or. patchout%Qsb, &
        restart=.true., &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(ssnow%rnof2), &
        scale=(1.0 / dels), &
        metadata=[ &
          attribute("units", "kg/m^2/s"), &
          attribute("long_name", "Subsurface runoff") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="et", &
        netcdf_name="Evap", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%Evap, &
        active=output%flux .or. output%Evap, &
        patchout=output%patch .or. patchout%Evap, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(canopy%et), &
        metadata=[ &
          attribute("units", "kg/m^2/s"), &
          attribute("long_name", "Total evapotranspiration") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="epot", &
        netcdf_name="PotEvap", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%PotEvap, &
        active=output%flux .or. output%PotEvap, &
        patchout=output%patch .or. patchout%PotEvap, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(canopy%epot), &
        scale=(1.0 / dels), &
        metadata=[ &
          attribute("units", "kg/m^2/s"), &
          attribute("long_name", "Potential evaporation") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="eint", &
        netcdf_name="ECanop", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%ECanop, &
        active=output%flux .or. output%ECanop, &
        patchout=output%patch .or. patchout%ECanop, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(canopy%eint), &
        metadata=[ &
          attribute("units", "kg/m^2/s"), &
          attribute("long_name", "Wet canopy evaporation") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="tveg", &
        netcdf_name="TVeg", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%TVeg, &
        active=output%flux .or. output%TVeg, &
        patchout=output%patch .or. patchout%TVeg, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(canopy%tveg), &
        metadata=[ &
          attribute("units", "kg/m^2/s"), &
          attribute("long_name", "Vegetation transpiration") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="esoil", &
        netcdf_name="ESoil", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%ESoil, &
        active=output%flux .or. output%ESoil, &
        patchout=output%patch .or. patchout%ESoil, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(canopy%esoil), &
        metadata=[ &
          attribute("units", "kg/m^2/s"), &
          attribute("long_name", "Evaporation from soil") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="fhv", &
        netcdf_name="HVeg", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%HVeg, &
        active=output%flux .or. output%HVeg, &
        patchout=output%patch .or. patchout%HVeg, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(canopy%fhv), &
        metadata=[ &
          attribute("units", "W/m^2"), &
          attribute("long_name", "Sensible heat from vegetation") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="fhs", &
        netcdf_name="HSoil", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%HSoil, &
        active=output%flux .or. output%HSoil, &
        patchout=output%patch .or. patchout%HSoil, &
        restart=.true., &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(canopy%fhs), &
        metadata=[ &
          attribute("units", "W/m^2"), &
          attribute("long_name", "Sensible heat from soil") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="fns", &
        netcdf_name="RnetSoil", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%HSoil, &
        active=output%flux .or. output%RnetSoil, &
        patchout=output%patch .or. patchout%RnetSoil, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(canopy%fns), &
        metadata=[ &
          attribute("units", "W/m^2"), &
          attribute("long_name", "Net radiation absorbed by ground") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="wb", &
        netcdf_name="SoilMoist", &
        data_shape=[DIM_PATCH, DIM_SOIL], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%SoilMoist, &
        active=output%soil .or. output%SoilMoist, &
        patchout=output%patch .or. patchout%SoilMoist, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(ssnow%wb), &
        restart=.true., &
        metadata=[ &
          attribute("units", "m^3/m^3"), &
          attribute("long_name", "Average layer soil moisture") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="wbice", &
        netcdf_name="SoilMoistIce", &
        data_shape=[DIM_PATCH, DIM_SOIL], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%SoilMoist, &
        active=output%soil .or. output%SoilMoistIce, &
        patchout=output%patch .or. patchout%SoilMoistIce, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(ssnow%wbice), &
        restart=.true., &
        metadata=[ &
          attribute("units", "m^3/m^3"), &
          attribute("long_name", "Average layer frozen soil moisture") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="tgg", &
        netcdf_name="SoilTemp", &
        data_shape=[DIM_PATCH, DIM_SOIL], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%SoilTemp, &
        active=output%soil .or. output%SoilTemp, &
        patchout=output%patch .or. patchout%SoilTemp, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(ssnow%tgg), &
        restart=.true., &
        metadata=[ &
          attribute("units", "K"), &
          attribute("long_name", "Average layer soil temperature") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="gammzz", &
        data_shape=[DIM_PATCH, DIM_SOIL], &
        var_type=CABLE_NETCDF_FLOAT, &
        active=.false., &
        aggregation_method="mean", &
        restart=.true., &
        aggregator=new_aggregator(ssnow%gammzz), &
        metadata=[ &
          attribute("units", "J/kg/C"), &
          attribute("long_name", "Heat capacity for each soil layer") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="ssdn", &
        data_shape=[DIM_PATCH, DIM_SNOW], &
        var_type=CABLE_NETCDF_FLOAT, &
        active=.false., &
        aggregation_method="mean", &
        restart=.true., &
        aggregator=new_aggregator(ssnow%ssdn), &
        metadata=[ &
          attribute("units", "kg/m^3"), &
          attribute("long_name", "Average layer snow density") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="smass", &
        data_shape=[DIM_PATCH, DIM_SNOW], &
        var_type=CABLE_NETCDF_FLOAT, &
        active=.false., &
        aggregation_method="mean", &
        restart=.true., &
        aggregator=new_aggregator(ssnow%smass), &
        metadata=[ &
          attribute("units", "kg/m^2"), &
          attribute("long_name", "Average layer snow mass") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="tgg1", &
        netcdf_name="BaresoilT", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%BaresoilT, &
        active=output%soil .or. output%BaresoilT, &
        patchout=output%patch .or. patchout%BaresoilT, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(ssnow%tgg(:, 1)), &
        metadata=[ &
          attribute("units", "K"), &
          attribute("long_name", "Bare soil temperature") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="snowd", &
        netcdf_name="SWE", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%SWE, &
        active=output%snow .or. output%SWE, &
        patchout=output%patch .or. patchout%SWE, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(ssnow%snowd), &
        restart=.true., &
        metadata=[ &
          attribute("units", "kg/m^2"), &
          attribute("long_name", "Snow water equivalent") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="smelt", &
        netcdf_name="SnowMelt", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%SnowMelt, &
        active=output%snow .or. output%SnowMelt, &
        patchout=output%patch .or. patchout%SnowMelt, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(ssnow%smelt), &
        scale=(1.0 / dels), &
        metadata=[ &
          attribute("units", "kg/m^2/s"), &
          attribute("long_name", "Snow Melt Rate") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="tggsn", &
        data_shape=[DIM_PATCH, DIM_SNOW], &
        var_type=CABLE_NETCDF_FLOAT, &
        active=.false., &
        aggregation_method="mean", &
        restart=.true., &
        aggregator=new_aggregator(ssnow%tggsn), &
        metadata=[ &
          attribute("units", "K"), &
          attribute("long_name", "Average layer snow temperature") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="tggsn1", &
        netcdf_name="SnowT", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%SnowT, &
        active=output%snow .or. output%SnowT, &
        patchout=output%patch .or. patchout%SnowT, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(ssnow%tggsn(:, 1)), &
        metadata=[ &
          attribute("units", "K"), &
          attribute("long_name", "Snow surface temperature") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="sdepth", &
        data_shape=[DIM_PATCH, DIM_SNOW], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%SnowDepth, &
        active=.false., &
        aggregation_method="mean", &
        restart=.true., &
        aggregator=new_aggregator(ssnow%sdepth), &
        metadata=[ &
          attribute("units", "m"), &
          attribute("long_name", "Snow layer depth") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="totsdepth", &
        netcdf_name="SnowDepth", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%SnowDepth, &
        active=output%snow .or. output%SnowDepth, &
        patchout=output%patch .or. patchout%SnowDepth, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(ssnow%totsdepth), &
        metadata=[ &
          attribute("units", "m"), &
          attribute("long_name", "Snow depth") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="swnet", &
        netcdf_name="SWnet", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%SWnet, &
        active=output%radiation .or. output%SWnet, &
        patchout=output%patch .or. patchout%SWnet, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(rad%swnet), &
        metadata=[ &
          attribute("units", "W/m^2"), &
          attribute("long_name", "Net shortwave radiation absorbed by surface") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="lwnet", &
        netcdf_name="LWnet", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%LWnet, &
        active=output%radiation .or. output%LWnet, &
        patchout=output%patch .or. patchout%LWnet, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(rad%lwnet), &
        metadata=[ &
          attribute("units", "W/m^2"), &
          attribute("long_name", "Net longwave radiation absorbed by surface") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="rnet", &
        netcdf_name="Rnet", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%Rnet, &
        active=output%radiation .or. output%Rnet, &
        patchout=output%patch .or. patchout%Rnet, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(rad%rnet), &
        metadata=[ &
          attribute("units", "W/m^2"), &
          attribute("long_name", "Net radiation absorbed by surface") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="albedo_T", &
        netcdf_name="Albedo", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%Albedo, &
        active=output%radiation .or. output%Albedo, &
        patchout=output%patch .or. patchout%Albedo, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(rad%albedo_T), &
        metadata=[ &
          attribute("units", "-"), &
          attribute("long_name", "Surface albedo") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="albedo_vis", &
        netcdf_name="visAlbedo", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%visAlbedo, &
        active=output%radiation .or. output%visAlbedo, &
        patchout=output%patch .or. patchout%visAlbedo, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(rad%albedo(:, 1)), &
        metadata=[ &
          attribute("units", "-"), &
          attribute("long_name", "Surface vis albedo") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="albedo_nir", &
        netcdf_name="nirAlbedo", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%nirAlbedo, &
        active=output%radiation .or. output%nirAlbedo, &
        patchout=output%patch .or. patchout%nirAlbedo, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(rad%albedo(:, 2)), &
        metadata=[ &
          attribute("units", "-"), &
          attribute("long_name", "Surface nir albedo") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="trad", &
        netcdf_name="RadT", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%RadT, &
        active=output%radiation .or. output%RadT, &
        patchout=output%patch .or. patchout%RadT, &
        restart=.true., &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(rad%trad), &
        metadata=[ &
          attribute("units", "K"), &
          attribute("long_name", "Radiative surface temperature") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="tscrn", &
        netcdf_name="Tscrn", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%Tscrn, &
        active=output%veg .or. output%Tscrn, &
        patchout=output%patch .or. patchout%Tscrn, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(canopy%tscrn), &
        metadata=[ &
          attribute("units", "oC"), &
          attribute("long_name", "screen level air temperature") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="tscrn_max", &
        netcdf_name="Txx", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%Tscrn, &
        active=output%veg .or. output%Tex, &
        patchout=output%patch .or. patchout%Tex, &
        reduction_method="grid_cell_average", &
        aggregation_method="max", &
        aggregator=new_aggregator(canopy%tscrn), &
        metadata=[ &
          attribute("units", "oC"), &
          attribute("long_name", "max screen-level T in reporting period") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="tscrn_min", &
        netcdf_name="Tnn", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%Tscrn, &
        active=output%veg .or. output%Tex, &
        patchout=output%patch .or. patchout%Tex, &
        reduction_method="grid_cell_average", &
        aggregation_method="min", &
        aggregator=new_aggregator(canopy%tscrn), &
        metadata=[ &
          attribute("units", "oC"), &
          attribute("long_name", "min screen-level T in reporting period") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="tscrn_max_daily", &
        netcdf_name="Tmx", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        active=(output%veg .or. output%Tex) .and. output%averaging == "monthly", &
        patchout=output%patch .or. patchout%Tex, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        range=ranges%Tscrn, &
        accumulation_frequency="daily", &
        aggregator=new_aggregator(canopy%tscrn_max_daily%aggregated_data), &
        metadata=[ &
          attribute("units", "oC"), &
          attribute("long_name", "averaged daily maximum screen-level T") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="tscrn_min_daily", &
        netcdf_name="Tmn", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        active=(output%veg .or. output%Tex) .and. output%averaging == "monthly", &
        patchout=output%patch .or. patchout%Tex, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        range=ranges%Tscrn, &
        accumulation_frequency="daily", &
        aggregator=new_aggregator(canopy%tscrn_min_daily%aggregated_data), &
        metadata=[ &
          attribute("units", "oC"), &
          attribute("long_name", "averaged daily minimum screen-level T") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="qscrn", &
        netcdf_name="Qscrn", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%Qscrn, &
        active=output%veg .or. output%Qscrn, &
        patchout=output%patch .or. patchout%Qscrn, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(canopy%qscrn), &
        metadata=[ &
          attribute("units", "kg/kg"), &
          attribute("long_name", "screen level specific humidity") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="tv", &
        netcdf_name="VegT", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%VegT, &
        active=output%veg .or. output%VegT, &
        patchout=output%patch .or. patchout%VegT, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(canopy%tv), &
        metadata=[ &
          attribute("units", "K"), &
          attribute("long_name", "Average vegetation temperature") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="tvair", &
        netcdf_name="CanT", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%CanT, &
        active=output%veg .or. output%CanT, &
        patchout=output%patch .or. patchout%CanT, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(met%tvair), &
        metadata=[ &
          attribute("units", "K"), &
          attribute("long_name", "Within-canopy temperature") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="fwsoil", &
        netcdf_name="Fwsoil", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%Fwsoil, &
        active=output%veg .or. output%Fwsoil, &
        patchout=output%patch .or. patchout%Fwsoil, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(canopy%fwsoil), &
        metadata=[ &
          attribute("units", "[-]"), &
          attribute("long_name", "soil moisture modifier to stomatal conductance") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="cansto", &
        netcdf_name="CanopInt", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%CanopInt, &
        active=output%veg .or. output%CanopInt, &
        patchout=output%patch .or. patchout%CanopInt, &
        restart=.true., &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(canopy%cansto), &
        metadata=[ &
          attribute("units", "kg/m^2"), &
          attribute("long_name", "Canopy intercepted water storage") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="vlai", &
        netcdf_name="LAI", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%LAI, &
        active=output%veg .or. output%LAI, &
        patchout=output%patch .or. patchout%LAI, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(veg%vlai), &
        metadata=[ &
          attribute("units", "-"), &
          attribute("long_name", "Leaf area index") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="ebal_tot", &
        netcdf_name="Ebal", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%Ebal, &
        active=output%balances .or. output%Ebal, &
        patchout=output%patch .or. patchout%Ebal, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(bal%ebal_tot), &
        metadata=[ &
          attribute("units", "W/m^2"), &
          attribute("long_name", "Cumulative energy imbalance") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="wbal_tot", &
        netcdf_name="Wbal", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%Wbal, &
        active=output%balances .or. output%Wbal, &
        patchout=output%patch .or. patchout%Wbal, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(bal%wbal_tot), &
        metadata=[ &
          attribute("units", "kg/m^2"), &
          attribute("long_name", "Cumulative water imbalance") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="wbtot0", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        active=.false., &
        restart=.true., &
        aggregator=new_aggregator(bal%wbtot0), &
        metadata=[ &
          attribute("units", "mm"), &
          attribute("long_name", "Initial time step soil water total") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="osnowd0", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        active=.false., &
        restart=.true., &
        aggregator=new_aggregator(bal%osnowd0), &
        metadata=[ &
          attribute("units", "mm"), &
          attribute("long_name", "Initial time step snow water total") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="frday", &
        netcdf_name="LeafResp", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%AutoResp, &
        active=output%carbon .or. output%LeafResp, &
        patchout=output%patch .or. patchout%LeafResp, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(canopy%frday), &
        scale=(1.0 / c_molar_mass), &
        metadata=[ &
          attribute("units", "umol/m^2/s"), &
          attribute("long_name", "Leaf respiration") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="frs", &
        netcdf_name="HeteroResp", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%HeteroResp, &
        active=output%carbon .or. output%HeteroResp, &
        patchout=output%patch .or. patchout%HeteroResp, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(canopy%frs), &
        scale=(1.0 / c_molar_mass), &
        metadata=[ &
          attribute("units", "umol/m^2/s"), &
          attribute("long_name", "Heterotrophic respiration") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="fgpp", &
        netcdf_name="GPP", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%GPP, &
        active=output%carbon .or. output%GPP, &
        patchout=output%patch .or. patchout%GPP, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(canopy%fgpp), &
        scale=(1.0 / c_molar_mass), &
        metadata=[ &
          attribute("units", "umol/m^2/s"), &
          attribute("long_name", "Gross primary production") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="fnpp", &
        netcdf_name="NPP", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%NPP, &
        active=output%carbon .or. output%NPP, &
        patchout=output%patch .or. patchout%NPP, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(canopy%fnpp), &
        scale=(1.0 / c_molar_mass), &
        metadata=[ &
          attribute("units", "umol/m^2/s"), &
          attribute("long_name", "Net primary production") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="fra", &
        netcdf_name="AutoResp", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%AutoResp, &
        active=output%carbon .or. output%AutoResp, &
        patchout=output%patch .or. patchout%AutoResp, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(canopy%fra), &
        scale=(1.0 / c_molar_mass), &
        metadata=[ &
          attribute("units", "umol/m^2/s"), &
          attribute("long_name", "Autotrophic respiration") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="fnee", &
        netcdf_name="NEE", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%NEE, &
        active=output%flux .or. output%NEE, &
        patchout=output%patch .or. patchout%NEE, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(canopy%fnee), &
        scale=(1.0 / c_molar_mass), &
        metadata=[ &
          attribute("units", "umol/m^2/s"), &
          attribute("long_name", "Net ecosystem exchange of CO2") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="wtd", &
        netcdf_name="WatTable", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%WatTable, &
        active=output%soil .or. output%WatTable, &
        patchout=output%patch .or. patchout%WatTable, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(ssnow%wtd), &
        scale=1e-3, &
        metadata=[ &
          attribute("units", "m"), &
          attribute("long_name", "Water Table Depth") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="GWwb", &
        netcdf_name="GWMoist", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%GWwb, &
        active=output%soil .or. output%GWMoist, &
        patchout=output%patch .or. patchout%GWMoist, &
        restart=.true., &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(ssnow%GWwb), &
        metadata=[ &
          attribute("units", "mm3/mm3"), &
          attribute("long_name", "Aquifer moisture content") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="satfrac", &
        netcdf_name="SatFrac", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%SatFrac, &
        active=output%soil .or. output%SatFrac, &
        patchout=output%patch .or. patchout%SatFrac, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(ssnow%satfrac), &
        metadata=[ &
          attribute("units", "unitless"), &
          attribute("long_name", "Saturated fraction of grid cell") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="Qrecharge", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%Qrecharge, &
        active=output%soil .or. output%Qrecharge, &
        patchout=output%patch .or. patchout%Qrecharge, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(ssnow%Qrecharge), &
        metadata=[ &
          attribute("units", "mm/s"), &
          attribute("long_name", "Recharge to or from Aquifer") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="tss", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        active=.false., &
        restart=.true., &
        aggregator=new_aggregator(ssnow%tss), &
        metadata=[ &
          attribute("units", "K"), &
          attribute("long_name", "Combined soil/snow temperature") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="rtsoil", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        active=.false., &
        restart=.true., &
        aggregator=new_aggregator(ssnow%rtsoil), &
        metadata=[ &
          attribute("units", "??"), &
          attribute("long_name", "Turbulent resistance for soil") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="runoff", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        active=.false., &
        restart=.true., &
        aggregator=new_aggregator(ssnow%runoff), &
        metadata=[ &
          attribute("units", "mm/timestep"), &
          attribute("long_name", "Total runoff") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="ssdnn", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        active=.false., &
        restart=.true., &
        aggregator=new_aggregator(ssnow%ssdnn), &
        metadata=[ &
          attribute("units", "kg/m^3"), &
          attribute("long_name", "Average snow density") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="snage", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        active=.false., &
        restart=.true., &
        aggregator=new_aggregator(ssnow%snage), &
        metadata=[ &
          attribute("units", "??"), &
          attribute("long_name", "Snow age") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="osnowd", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        active=.false., &
        restart=.true., &
        aggregator=new_aggregator(ssnow%osnowd), &
        metadata=[ &
          attribute("units", "mm"), &
          attribute("long_name", "Previous time step snow depth in water equivalent") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="albsoilsn", &
        data_shape=[DIM_PATCH, DIM_RAD], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%albsoiln, &
        active=.false., &
        aggregation_method="mean", &
        restart=.true., &
        aggregator=new_aggregator(ssnow%albsoilsn), &
        metadata=[ &
          attribute("units", "-"), &
          attribute("long_name", "Combined soil/snow albedo") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="isflag", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_INT, &
        active=.false., &
        restart=.true., &
        aggregator=new_aggregator(ssnow%isflag), &
        metadata=[ &
          attribute("units", "-"), &
          attribute("long_name", "Snow layer scheme flag") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="ghflux", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        active=.false., &
        restart=.true., &
        aggregator=new_aggregator(canopy%ghflux), &
        metadata=[ &
          attribute("units", "W/m^2"), &
          attribute("long_name", "????") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="sghflux", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        active=.false., &
        restart=.true., &
        aggregator=new_aggregator(canopy%sghflux), &
        metadata=[ &
          attribute("units", "W/m^2"), &
          attribute("long_name", "????") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="dgdtg", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        active=.false., &
        restart=.true., &
        aggregator=new_aggregator(canopy%dgdtg), &
        metadata=[ &
          attribute("units", "W/m^2/K"), &
          attribute("long_name", "Derivative of ground heat flux wrt soil temperature") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="fev", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        active=.false., &
        restart=.true., &
        aggregator=new_aggregator(canopy%fev), &
        metadata=[ &
          attribute("units", "W/m^2"), &
          attribute("long_name", "Latent heat flux from vegetation") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="fes", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        active=.false., &
        restart=.true., &
        aggregator=new_aggregator(canopy%fes), &
        metadata=[ &
          attribute("units", "W/m^2"), &
          attribute("long_name", "Latent heat flux from soil") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="albedo", &
        data_shape=[DIM_PATCH, DIM_RAD], &
        var_type=CABLE_NETCDF_FLOAT, &
        active=.false., &
        restart=.true., &
        aggregator=new_aggregator(rad%albedo), &
        metadata=[ &
          attribute("units", "-"), &
          attribute("long_name", "Albedo for shortwave and NIR radiation") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="iveg", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_INT, &
        range=ranges%iveg, &
        active=output%params .or. output%iveg, &
        patchout=output%patch .or. patchout%iveg, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        restart=.true., &
        aggregator=new_aggregator(veg%iveg), &
        metadata=[ &
          attribute("units", "-"), &
          attribute("long_name", "Vegetation type") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="patchfrac", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%patchfrac, &
        active=(output%params .or. output%patchfrac) .and. (output%patch .or. patchout%patchfrac), &
        patchout=output%patch .or. patchout%patchfrac, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        restart=.true., &
        aggregator=new_aggregator(patch(:)%frac), &
        metadata=[ &
          attribute("units", "-"), &
          attribute("long_name", "Fractional cover of vegetation patches") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="isoil", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_INT, &
        range=ranges%isoil, &
        active=output%params .or. output%isoil, &
        patchout=output%patch .or. patchout%isoil, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        restart=.true., &
        aggregator=new_aggregator(soil%isoilm), &
        metadata=[ &
          attribute("units", "-"), &
          attribute("long_name", "Soil type") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="bch", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%bch, &
        active=output%params .or. output%bch, &
        patchout=output%patch .or. patchout%bch, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(soil%bch), &
        metadata=[ &
          attribute("units", "-"), &
          attribute("long_name", "Parameter b, Campbell eqn 1985") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="clay", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%clay, &
        active=output%params .or. output%clay, &
        patchout=output%patch .or. patchout%clay, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(soil%clay), &
        metadata=[ &
          attribute("units", "-"), &
          attribute("long_name", "Fraction of soil which is clay") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="sand", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%sand, &
        active=output%params .or. output%sand, &
        patchout=output%patch .or. patchout%sand, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(soil%sand), &
        metadata=[ &
          attribute("units", "-"), &
          attribute("long_name", "Fraction of soil which is sand") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="silt", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%silt, &
        active=output%params .or. output%silt, &
        patchout=output%patch .or. patchout%silt, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(soil%silt), &
        metadata=[ &
          attribute("units", "-"), &
          attribute("long_name", "Fraction of soil which is silt") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="ssat", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%ssat, &
        active=output%params .or. output%ssat, &
        patchout=output%patch .or. patchout%ssat, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(soil%ssat), &
        metadata=[ &
          attribute("units", "-"), &
          attribute("long_name", "Fraction of soil volume which is water @ saturation") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="sfc", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%sfc, &
        active=output%params .or. output%sfc, &
        patchout=output%patch .or. patchout%sfc, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(soil%sfc), &
        metadata=[ &
          attribute("units", "-"), &
          attribute("long_name", "Fraction of soil volume which is water @ field capacity") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="swilt", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%swilt, &
        active=output%params .or. output%swilt, &
        patchout=output%patch .or. patchout%swilt, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(soil%swilt), &
        metadata=[ &
          attribute("units", "-"), &
          attribute("long_name", "Fraction of soil volume which is water @ wilting point") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="hyds", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%hyds, &
        active=output%params .or. output%hyds, &
        patchout=output%patch .or. patchout%hyds, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(soil%hyds), &
        metadata=[ &
          attribute("units", "m/s"), &
          attribute("long_name", "Hydraulic conductivity @ saturation") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="sucs", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%sucs, &
        active=output%params .or. output%sucs, &
        patchout=output%patch .or. patchout%sucs, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(soil%sucs), &
        metadata=[ &
          attribute("units", "m"), &
          attribute("long_name", "Suction @ saturation") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="css", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%css, &
        active=output%params .or. output%css, &
        patchout=output%patch .or. patchout%css, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(soil%css), &
        metadata=[ &
          attribute("units", "J/kg/C"), &
          attribute("long_name", "Heat capacity of soil minerals") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="rhosoil", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%rhosoil, &
        active=output%params .or. output%rhosoil, &
        patchout=output%patch .or. patchout%rhosoil, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(soil%rhosoil), &
        metadata=[ &
          attribute("units", "kg/m^3"), &
          attribute("long_name", "Density of soil minerals") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="rs20", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%rs20, &
        active=output%params .or. output%rs20, &
        patchout=output%patch .or. patchout%rs20, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(veg%rs20), &
        metadata=[ &
          attribute("units", "-"), &
          attribute("long_name", "Soil respiration coefficient at 20C") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="albsoil", &
        data_shape=[DIM_PATCH, DIM_RAD], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%albsoil, &
        active=output%params .or. output%albsoil, &
        patchout=output%patch .or. patchout%albsoil, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        restart=.true., &
        aggregator=new_aggregator(soil%albsoil), &
        metadata=[ &
          attribute("units", "-"), &
          attribute("long_name", "Snow free shortwave soil reflectance fraction") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="hc", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%hc, &
        active=output%params .or. output%hc, &
        patchout=output%patch .or. patchout%hc, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(veg%hc), &
        metadata=[ &
          attribute("units", "m"), &
          attribute("long_name", "Height of canopy") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="canst1", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%canst1, &
        active=output%params .or. output%canst1, &
        patchout=output%patch .or. patchout%canst1, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(veg%canst1), &
        metadata=[ &
          attribute("units", "mm/LAI"), &
          attribute("long_name", "Max water intercepted by canopy") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="dleaf", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%dleaf, &
        active=output%params .or. output%dleaf, &
        patchout=output%patch .or. patchout%dleaf, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(veg%dleaf), &
        metadata=[ &
          attribute("units", "m"), &
          attribute("long_name", "Characteristic length of leaf") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="frac4", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%frac4, &
        active=output%params .or. output%frac4, &
        patchout=output%patch .or. patchout%frac4, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(veg%frac4), &
        metadata=[ &
          attribute("units", "-"), &
          attribute("long_name", "Fraction of plants which are C4") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="ejmax", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%ejmax, &
        active=output%params .or. output%ejmax, &
        patchout=output%patch .or. patchout%ejmax, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(veg%ejmax), &
        metadata=[ &
          attribute("units", "mol/m^2/s"), &
          attribute("long_name", "Max potential electron transport rate top leaf") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="vcmax", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%vcmax, &
        active=output%params .or. output%vcmax, &
        patchout=output%patch .or. patchout%vcmax, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(veg%vcmax), &
        metadata=[ &
          attribute("units", "mol/m^2/s"), &
          attribute("long_name", "Maximum RuBP carboxylation rate top leaf") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="rp20", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%rp20, &
        active=output%params .or. output%rp20, &
        patchout=output%patch .or. patchout%rp20, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(veg%rp20), &
        metadata=[ &
          attribute("units", "-"), &
          attribute("long_name", "Plant respiration coefficient at 20C") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="g0", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%g0, &
        active=output%params .or. output%g0, &
        patchout=output%patch .or. patchout%g0, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(veg%g0), &
        metadata=[ &
          attribute("units", "-"), &
          attribute("long_name", "g0 term in Medlyn Stom Cond. Param") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="g1", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%g1, &
        active=output%params .or. output%g1, &
        patchout=output%patch .or. patchout%g1, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(veg%g1), &
        metadata=[ &
          attribute("units", "-"), &
          attribute("long_name", "g1 term in Medlyn Stom Cond. Param") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="rpcoef", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%rpcoef, &
        active=output%params .or. output%rpcoef, &
        patchout=output%patch .or. patchout%rpcoef, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(veg%rpcoef), &
        metadata=[ &
          attribute("units", "1/C"), &
          attribute("long_name", "Temperature coef nonleaf plant respiration") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="shelrb", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%shelrb, &
        active=output%params .or. output%shelrb, &
        patchout=output%patch .or. patchout%shelrb, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(veg%shelrb), &
        metadata=[ &
          attribute("units", "-"), &
          attribute("long_name", "Sheltering factor") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="xfang", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%xfang, &
        active=output%params .or. output%xfang, &
        patchout=output%patch .or. patchout%xfang, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(veg%xfang), &
        metadata=[ &
          attribute("units", "-"), &
          attribute("long_name", "Leaf angle parameter") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="wai", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%wai, &
        active=output%params .or. output%wai, &
        patchout=output%patch .or. patchout%wai, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(veg%wai), &
        metadata=[ &
          attribute("units", "-"), &
          attribute("long_name", "Wood area index") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="vegcf", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%vegcf, &
        active=output%params .or. output%vegcf, &
        patchout=output%patch .or. patchout%vegcf, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(veg%vegcf), &
        metadata=[ &
          attribute("units", "-"), &
          attribute("long_name", "vegcf") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="extkn", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%extkn, &
        active=output%params .or. output%extkn, &
        patchout=output%patch .or. patchout%extkn, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(veg%extkn), &
        metadata=[ &
          attribute("units", "-"), &
          attribute("long_name", "Nitrogen extinction coef for vert. canopy profile") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="tminvj", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%tminvj, &
        active=output%params .or. output%tminvj, &
        patchout=output%patch .or. patchout%tminvj, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(veg%tminvj), &
        metadata=[ &
          attribute("units", "C"), &
          attribute("long_name", "Min temperature for the start of photosynthesis") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="tmaxvj", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%tmaxvj, &
        active=output%params .or. output%tmaxvj, &
        patchout=output%patch .or. patchout%tmaxvj, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(veg%tmaxvj), &
        metadata=[ &
          attribute("units", "C"), &
          attribute("long_name", "Max temperature for photosynthesis") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="vbeta", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%vbeta, &
        active=output%params .or. output%vbeta, &
        patchout=output%patch .or. patchout%vbeta, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(veg%vbeta), &
        metadata=[ &
          attribute("units", "-"), &
          attribute("long_name", "Stomatal sensitivity to soil water") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="xalbnir", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%xalbnir, &
        active=output%params .or. output%xalbnir, &
        patchout=output%patch .or. patchout%xalbnir, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(veg%xalbnir), &
        metadata=[ &
          attribute("units", "-"), &
          attribute("long_name", "Modifier for albedo in near ir band") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="meth", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%meth, &
        active=output%params .or. output%meth, &
        patchout=output%patch .or. patchout%meth, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(veg%meth), &
        metadata=[ &
          attribute("units", "-"), &
          attribute("long_name", "Canopy turbulence parameterisation choice") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="za_uv", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%za, &
        active=output%params .or. output%za, &
        patchout=output%patch .or. patchout%za, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(rough%za_uv), &
        metadata=[ &
          attribute("units", "m"), &
          attribute("long_name", "Reference height (lowest atm. model layer) for momentum") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="za_tq", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%za, &
        active=output%params .or. output%za, &
        patchout=output%patch .or. patchout%za, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(rough%za_tq), &
        metadata=[ &
          attribute("units", "m"), &
          attribute("long_name", "Reference height (lowest atm. model layer) for scalars") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="ratecp", &
        data_shape=[DIM_PLANTCARBON], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%ratecp, &
        active=output%params .or. output%ratecp, &
        distributed=.false., &
        parameter=.true., &
        aggregator=new_aggregator(bgc%ratecp), &
        metadata=[ &
          attribute("units", "1/year"), &
          attribute("long_name", "Plant carbon rate constant") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="ratecs", &
        data_shape=[DIM_SOILCARBON], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%ratecs, &
        active=output%params .or. output%ratecs, &
        distributed=.false., &
        parameter=.true., &
        aggregator=new_aggregator(bgc%ratecs), &
        metadata=[ &
          attribute("units", "1/year"), &
          attribute("long_name", "Soil carbon rate constant") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="cplant", &
        data_shape=[DIM_PATCH, DIM_PLANTCARBON], &
        var_type=CABLE_NETCDF_FLOAT, &
        active=.false., &
        restart=.true., &
        aggregator=new_aggregator(bgc%cplant), &
        metadata=[ &
          attribute("units", "gC/m^2"), &
          attribute("long_name", "Plant carbon stores") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="csoil", &
        data_shape=[DIM_PATCH, DIM_SOILCARBON], &
        var_type=CABLE_NETCDF_FLOAT, &
        active=.false., &
        restart=.true., &
        aggregator=new_aggregator(bgc%csoil), &
        metadata=[ &
          attribute("units", "gC/m^2"), &
          attribute("long_name", "Soil carbon stores") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="zse", &
        data_shape=[DIM_SOIL], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%zse, &
        active=output%params .or. output%zse, &
        distributed=.false., &
        parameter=.true., &
        restart=.true., &
        aggregator=new_aggregator(soil%zse), &
        metadata=[ &
          attribute("units", "m"), &
          attribute("long_name", "Depth of each soil layer") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="froot", &
        data_shape=[DIM_PATCH, DIM_SOIL], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%froot, &
        active=output%params .or. output%froot, &
        patchout=output%patch .or. patchout%froot, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(veg%froot), &
        metadata=[ &
          attribute("units", "-"), &
          attribute("long_name", "Fraction of roots in each soil layer") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="slope", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%slope, &
        active=output%params .or. output%slope, &
        patchout=output%patch .or. patchout%slope, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(soil%slope), &
        metadata=[ &
          attribute("units", "-"), &
          attribute("long_name", "Mean subgrid topographic slope") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="slope_std", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%slope_std, &
        active=output%params .or. output%slope_std, &
        patchout=output%patch .or. patchout%slope_std, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(soil%slope_std), &
        metadata=[ &
          attribute("units", "-"), &
          attribute("long_name", "Mean subgrid topographic slope_std") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="GWdz", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%GWdz, &
        active=output%params .AND. cable_user%gw_model, &
        patchout=output%patch .OR. patchout%GWdz, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(soil%GWdz), &
        metadata=[ &
          attribute("units", "m"), &
          attribute("long_name", "Mean aquifer layer thickness") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="MaxHorzDrainRate", &
        netcdf_name="Qhmax", &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%gw_default, &
        active=output%params .AND. cable_user%gw_model, &
        patchout=output%patch .OR. patchout%Qhmax, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(gw_params%MaxHorzDrainRate), &
        metadata=[ &
          attribute("units", "mm/s"), &
          attribute("long_name", "Maximum subsurface drainage") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="EfoldHorzDrainRate", &
        netcdf_name="QhmaxEfold", &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%gw_default, &
        active=output%params .AND. cable_user%gw_model, &
        patchout=output%patch .OR. patchout%QhmaxEfold, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(gw_params%EfoldHorzDrainRate), &
        metadata=[ &
          attribute("units", "m"), &
          attribute("long_name", "Maximum subsurface drainage decay rate") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="MaxSatFraction", &
        netcdf_name="SatFracmax", &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%gw_default, &
        active=output%params .AND. cable_user%gw_model, &
        patchout=output%patch .OR. patchout%SatFracmax, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(gw_params%MaxSatFraction), &
        metadata=[ &
          attribute("units", "-"), &
          attribute("long_name", "Controls max saturated fraction") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="hkrz", &
        netcdf_name="HKefold", &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%gw_default, &
        active=output%params .AND. cable_user%gw_model, &
        patchout=output%patch .OR. patchout%HKefold, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(gw_params%hkrz), &
        metadata=[ &
          attribute("units", "1/m"), &
          attribute("long_name", "Rate HK decays with depth") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="zdepth", &
        netcdf_name="HKdepth", &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%gw_default, &
        active=output%params .AND. cable_user%gw_model, &
        patchout=output%patch .OR. patchout%HKdepth, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(gw_params%zdepth), &
        metadata=[ &
          attribute("units", "m"), &
          attribute("long_name", "Depth at which HKsat(z) is HKsat(0)") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="nap", &
        data_shape=[DIM_LAND_GLOBAL], &
        var_type=CABLE_NETCDF_FLOAT, &
        active=.false., &
        restart=.true., &
        distributed=.false., &
        aggregation_method="point", &
        aggregator=new_aggregator(landpt_global(:)%nap), &
        metadata=[ &
          attribute("units", ""), &
          attribute("long_name", "") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="mvtype", &
        var_type=CABLE_NETCDF_FLOAT, &
        active=.false., &
        restart=.true., &
        distributed=.false., &
        aggregation_method="point", &
        aggregator=new_aggregator(mvtype), &
        metadata=[ &
          attribute("units", ""), &
          attribute("long_name", "Number of vegetation types") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="mstype", &
        var_type=CABLE_NETCDF_FLOAT, &
        active=.false., &
        restart=.true., &
        distributed=.false., &
        aggregation_method="point", &
        aggregator=new_aggregator(mstype), &
        metadata=[ &
          attribute("units", ""), &
          attribute("long_name", "Number of soil types") &
        ] &
      ) &
    ]

  end function cable_diagnostics_core

end module
