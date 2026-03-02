module cable_diagnostics_casa_mod
  use casavariable, only: casa_flux
  use casavariable, only: casa_pool
  use casavariable, only: casa_met

  use casaparm, only: LEAF
  use casaparm, only: WOOD
  use casaparm, only: FROOT

  use casaparm, only: MIC
  use casaparm, only: SLOW
  use casaparm, only: PASS

  use casaparm, only: METB
  use casaparm, only: STR
  use casaparm, only: CWD

  use cable_timing_utils_mod, only: seconds_per_day

  use cable_phys_constants_mod, only: c_molar_mass

  use cable_common_module, only: l_casacnp

  use cable_common_module, only: cable_user

  use cable_io_vars_module, only: output, patchout

  use cable_output_prototype_v2_mod, only: cable_output_variable_t
  use cable_output_prototype_v2_mod, only: attribute => cable_output_attribute_t
  use cable_output_prototype_v2_mod, only: DIM_PATCH => CABLE_OUTPUT_DIM_PATCH

  use cable_netcdf_mod, only: CABLE_NETCDF_FLOAT

  use aggregator_mod, only: new_aggregator

  use cable_checks_module, only: ranges

  implicit none
  private

  public :: cable_diagnostics_casa

contains

  function cable_diagnostics_casa(casaflux, casapool, casamet) result(casa_output_variables)
    type(casa_flux), intent(in) :: casaflux
    type(casa_pool), intent(in) :: casapool
    type(casa_met), intent(in) :: casamet
    type(cable_output_variable_t), allocatable :: casa_output_variables(:)

    if (.not. l_casacnp) then
      allocate(casa_output_variables(0))
      return
    end if

    casa_output_variables = [ &
      cable_output_variable_t( &
        field_name="crmplant_froot", &
        netcdf_name="RootResp", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%AutoResp, &
        active=output%carbon .or. output%casa, &
        patchout=output%patch .or. patchout%casa, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(casaflux%crmplant(:, FROOT)), &
        scale=(1.0 / seconds_per_day / c_molar_mass), &
        metadata=[ &
          attribute("units", "umol/m^2/s"), &
          attribute("long_name", "Fine Root Autotrophic respiration") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="crmplant_wood", &
        netcdf_name="StemResp", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%AutoResp, &
        active=output%carbon .or. output%casa, &
        patchout=output%patch .or. patchout%casa, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(casaflux%crmplant(:, WOOD)), &
        scale=(1.0 / seconds_per_day / c_molar_mass), &
        metadata=[ &
          attribute("units", "umol/m^2/s"), &
          attribute("long_name", "StemWood Autotrophic respiration") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="cnbp", &
        netcdf_name="NBP", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%NEE, &
        active=output%casa .or. output%NBP, &
        patchout=output%patch .or. patchout%NBP, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(casaflux%cnbp), &
        scale=(1.0 / seconds_per_day / c_molar_mass), &
        metadata=[ &
          attribute("units", "umol/m^2/s"), &
          attribute("long_name", "Net Biosphere Production") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="dCdt", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%NEE, &
        active=output%casa .or. output%dCdt, &
        patchout=output%patch .or. patchout%dCdt, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(casapool%dCdt), &
        scale=(1.0 / seconds_per_day / c_molar_mass), &
        metadata=[ &
          attribute("units", "umol/m^2/s"), &
          attribute("long_name", "Carbon accumulation rate (uptake +ve)") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="csoiltot", &
        netcdf_name="TotSoilCarb", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%TotSoilCarb, &
        active=output%casa .or. output%TotSoilCarb, &
        patchout=output%patch .or. patchout%TotSoilCarb, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(casapool%csoiltot), &
        scale=(1.0 / 1000.0), &
        metadata=[ &
          attribute("units", "kg C/m^2"), &
          attribute("long_name", "Total Soil and Litter Carbon") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="clittertot", &
        netcdf_name="TotLittCarb", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%TotLittCarb, &
        active=output%casa .or. output%TotLittCarb, &
        patchout=output%patch .or. patchout%TotLittCarb, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(casapool%clittertot), &
        scale=(1.0 / 1000.0), &
        metadata=[ &
          attribute("units", "kg C/m^2"), &
          attribute("long_name", "Total Litter Carbon") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="csoil_mic", &
        netcdf_name="SoilCarbFast", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%TotLittCarb, &
        active=output%casa .or. output%SoilCarbFast, &
        patchout=output%patch .or. patchout%SoilCarbFast, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(casapool%csoil(:, MIC)), &
        scale=(1.0 / 1000.0), &
        metadata=[ &
          attribute("units", "kg C/m^2"), &
          attribute("long_name", "Soil Carbon: Fast Turnover") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="csoil_slow", &
        netcdf_name="SoilCarbSlow", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%TotSoilCarb, &
        active=output%casa .or. output%SoilCarbSlow, &
        patchout=output%patch .or. patchout%SoilCarbSlow, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(casapool%csoil(:, SLOW)), &
        scale=(1.0 / 1000.0), &
        metadata=[ &
          attribute("units", "kg C/m^2"), &
          attribute("long_name", "Soil Carbon: Slow Turnover") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="csoil_pass", &
        netcdf_name="SoilCarbPassive", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%TotSoilCarb, &
        active=output%casa .or. output%SoilCarbPassive, &
        patchout=output%patch .or. patchout%SoilCarbPassive, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(casapool%csoil(:, PASS)), &
        scale=(1.0 / 1000.0), &
        metadata=[ &
          attribute("units", "kg C/m^2"), &
          attribute("long_name", "Soil Carbon: Passive") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="clitter_metb", &
        netcdf_name="LittCarbMetabolic", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%TotLittCarb, &
        active=output%casa .or. output%LittCarbMetabolic, &
        patchout=output%patch .or. patchout%LittCarbMetabolic, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(casapool%clitter(:, METB)), &
        scale=(1.0 / 1000.0), &
        metadata=[ &
          attribute("units", "kg C/m^2"), &
          attribute("long_name", "Litter Carbon: metabolic") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="clitter_str", &
        netcdf_name="LittCarbStructural", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%TotLittCarb, &
        active=output%casa .or. output%LittCarbStructural, &
        patchout=output%patch .or. patchout%LittCarbStructural, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(casapool%clitter(:, STR)), &
        scale=(1.0 / 1000.0), &
        metadata=[ &
          attribute("units", "kg C/m^2"), &
          attribute("long_name", "Litter Carbon: structural") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="clitter_cwd", &
        netcdf_name="LittCarbCWD", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%TotLittCarb, &
        active=output%casa .or. output%LittCarbCWD, &
        patchout=output%patch .or. patchout%LittCarbCWD, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(casapool%clitter(:, CWD)), &
        scale=(1.0 / 1000.0), &
        metadata=[ &
          attribute("units", "kg C/m^2"), &
          attribute("long_name", "Litter Carbon: CWD") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="cplant_leaf", &
        netcdf_name="PlantCarbLeaf", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%TotLittCarb, &
        active=output%casa .or. output%PlantCarbLeaf, &
        patchout=output%patch .or. patchout%PlantCarbLeaf, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(casapool%cplant(:, LEAF)), &
        scale=(1.0 / 1000.0), &
        metadata=[ &
          attribute("units", "kg C/m^2"), &
          attribute("long_name", "Plant Carbon: leaf") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="cplant_wood", &
        netcdf_name="PlantCarbWood", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%TotLittCarb, &
        active=output%casa .or. output%PlantCarbWood, &
        patchout=output%patch .or. patchout%PlantCarbWood, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(casapool%cplant(:, WOOD)), &
        scale=(1.0 / 1000.0), &
        metadata=[ &
          attribute("units", "kg C/m^2"), &
          attribute("long_name", "Plant Carbon: wood (above- and below-ground)") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="cplant_froot", &
        netcdf_name="PlantCarbFineRoot", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%TotLittCarb, &
        active=output%casa .or. output%PlantCarbFineRoot, &
        patchout=output%patch .or. patchout%PlantCarbFineRoot, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(casapool%cplant(:, FROOT)), &
        scale=(1.0 / 1000.0), &
        metadata=[ &
          attribute("units", "kg C/m^2"), &
          attribute("long_name", "Plant Carbon: Fine roots") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="cplanttot", &
        netcdf_name="TotLivBiomass", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%TotLivBiomass, &
        active=output%casa .or. output%TotLivBiomass, &
        patchout=output%patch .or. patchout%TotLivBiomass, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(casapool%cplanttot), &
        scale=(1.0 / 1000.0), &
        metadata=[ &
          attribute("units", "kg C/m^2"), &
          attribute("long_name", "Total Biomass") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="cplant_turnover_tot", &
        netcdf_name="PlantTurnover", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%NEE, &
        active=output%casa .or. output%PlantTurnover, &
        patchout=output%patch .or. patchout%PlantTurnover, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(casaflux%cplant_turnover_tot), &
        scale=(1.0 / seconds_per_day / c_molar_mass), &
        metadata=[ &
          attribute("units", "umol/m^2/s"), &
          attribute("long_name", "Total Biomass Turnover") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="Cplant_turnover_leaf", &
        netcdf_name="PlantTurnoverLeaf", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%NEE, &
        active=output%casa .or. output%PlantTurnoverLeaf, &
        patchout=output%patch .or. patchout%PlantTurnoverLeaf, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(casaflux%Cplant_turnover(:, LEAF)), &
        scale=(1.0 / seconds_per_day / c_molar_mass), &
        metadata=[ &
          attribute("units", "umol/m^2/s"), &
          attribute("long_name", "Leaf Biomass Turnover") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="Cplant_turnover_wood", &
        netcdf_name="PlantTurnoverWood", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%NEE, &
        active=output%casa .or. output%PlantTurnoverWood, &
        patchout=output%patch .or. patchout%PlantTurnoverWood, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(casaflux%Cplant_turnover(:, WOOD)), &
        scale=(1.0 / seconds_per_day / c_molar_mass), &
        metadata=[ &
          attribute("units", "umol/m^2/s"), &
          attribute("long_name", "Woody Biomass Turnover") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="Cplant_turnover_froot", &
        netcdf_name="PlantTurnoverFineRoot", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%NEE, &
        active=output%casa .or. output%PlantTurnoverFineRoot, &
        patchout=output%patch .or. patchout%PlantTurnoverFineRoot, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(casaflux%Cplant_turnover(:, FROOT)), &
        scale=(1.0 / seconds_per_day / c_molar_mass), &
        metadata=[ &
          attribute("units", "umol/m^2/s"), &
          attribute("long_name", "FineRoot Biomass Turnover") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="Cplant_turnover_disturbance", &
        netcdf_name="PlantTurnoverWoodDist", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%NEE, &
        active=output%casa .or. output%PlantTurnoverWoodDist, &
        patchout=output%patch .or. patchout%PlantTurnoverWoodDist, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(casaflux%Cplant_turnover_disturbance), &
        scale=(1.0 / seconds_per_day / c_molar_mass), &
        metadata=[ &
          attribute("units", "umol/m^2/s"), &
          attribute("long_name", "Woody Biomass Turnover (disturbance)") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="Cplant_turnover_crowding", &
        netcdf_name="PlantTurnoverWoodCrowding", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%NEE, &
        active=output%casa .or. output%PlantTurnoverWoodCrowding, &
        patchout=output%patch .or. patchout%PlantTurnoverWoodCrowding, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(casaflux%Cplant_turnover_crowding), &
        scale=(1.0 / seconds_per_day / c_molar_mass), &
        metadata=[ &
          attribute("units", "umol/m^2/s"), &
          attribute("long_name", "Woody Biomass Turnover (crowding)") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="Cplant_turnover_resource_limitation", &
        netcdf_name="PlantTurnoverWoodResourceLim", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%NEE, &
        active=output%casa .or. output%PlantTurnoverWoodResourceLim, &
        patchout=output%patch .or. patchout%PlantTurnoverWoodResourceLim, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(casaflux%Cplant_turnover_resource_limitation), &
        scale=(1.0 / seconds_per_day / c_molar_mass), &
        metadata=[ &
          attribute("units", "umol/m^2/s"), &
          attribute("long_name", "Woody Biomass Turnover (Resource Limitation)") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="areacell", &
        netcdf_name="Area", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%Area, &
        active=output%casa .or. output%Area, &
        patchout=output%patch .or. patchout%Area, &
        reduction_method="grid_cell_average", &
        aggregation_method="point", &
        aggregator=new_aggregator(casamet%areacell), &
        scale=1e-6, &
        metadata=[ &
          attribute("units", "km2"), &
          attribute("long_name", "Patch Area") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="FluxCtoLUC", &
        netcdf_name="LandUseFlux", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%NEE, &
        active=cable_user%POPLUC .or. output%LandUseFlux, &
        patchout=output%patch .or. patchout%LandUseFlux, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(casaflux%FluxCtoLUC), &
        scale=(1.0 / seconds_per_day / c_molar_mass), &
        metadata=[ &
          attribute("units", "umol/m^2/s"), &
          attribute("long_name", "Sum of wood harvest and clearing fluxes") &
        ] &
      ) &
    ]

  end function cable_diagnostics_casa

end module
