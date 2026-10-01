! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

module cable_output_bindings_mod
  !! Bindings between output catalogue names and CABLE model state.
  !!
  !! Everything that describes a variable lives in `output_catalogue.yaml`.
  !! This module supplies only what Fortran must provide: the model variable
  !! each name points at, its unit conversion, its valid range and whether it
  !! exists in the current model configuration. All model state and
  !! configuration is passed in as arguments.

  use cable_def_types_mod, only: met_type
  use cable_def_types_mod, only: canopy_type
  use cable_def_types_mod, only: soil_parameter_type
  use cable_def_types_mod, only: soil_snow_type
  use cable_def_types_mod, only: radiation_type
  use cable_def_types_mod, only: veg_parameter_type
  use cable_def_types_mod, only: balances_type
  use cable_def_types_mod, only: roughness_type
  use cable_def_types_mod, only: bgc_pool_type

  use cable_phys_constants_mod, only: c_molar_mass
  use cable_phys_constants_mod, only: HL

  use cable_common_module, only: gw_parameters_type
  use cable_checks_module, only: ranges_type
  use cable_io_vars_module, only: land_type
  use cable_io_vars_module, only: patch_type

  use aggregator_mod, only: new_aggregator

  use cable_output_mod, only: cable_output_binding_t

  implicit none
  private

  public :: cable_output_bindings

contains

  function cable_output_bindings( &
    met, canopy, soil, ssnow, rad, veg, bal, rough, bgc, ranges, gw_params, landpt_global, patch, mvtype, mstype, &
    dels, use_groundwater_model, use_popluc, calculate_soil_albedo &
  ) result(bindings)
    type(met_type), intent(inout) :: met
    type(canopy_type), intent(inout) :: canopy
    type(soil_parameter_type), intent(inout) :: soil
    type(soil_snow_type), intent(inout) :: ssnow
    type(radiation_type), intent(inout) :: rad
    type(veg_parameter_type), intent(inout) :: veg
    type(balances_type), intent(inout) :: bal
    type(roughness_type), intent(inout) :: rough
    type(bgc_pool_type), intent(inout) :: bgc
    type(ranges_type), intent(in) :: ranges
    type(gw_parameters_type), intent(inout) :: gw_params
    type(land_type), intent(inout) :: landpt_global(:)
    type(patch_type), intent(inout) :: patch(:)
    integer, intent(inout) :: mvtype
    integer, intent(inout) :: mstype
    real, intent(in) :: dels
    logical, intent(in) :: use_groundwater_model
    logical, intent(in) :: use_popluc
    logical, intent(in) :: calculate_soil_albedo
    type(cable_output_binding_t), allocatable :: bindings(:)

    ! One entry per catalogue variable, in catalogue order. Reading an entry:
    !   name=       the catalogue name this binding belongs to;
    !   aggregator= wraps the model's working variable. The aggregator keeps a
    !               pointer to it, so whatever the model stores there each
    !               time step is what gets sampled (no copy is made here);
    !   scale_by, divide_by, offset_by
    !               unit conversion: output = scale_by * value / divide_by + offset_by.
    !               Left out when no conversion is needed. Values such as dels
    !               or c_molar_mass are not known to the catalogue, which is why
    !               conversions live here and not in the YAML;
    !   range=      valid range in output units, checked if range checking is on;
    !   available=  false when the model configuration does not provide the
    !               variable (for example groundwater variables when the
    !               groundwater model is off).
    bindings = [ &
      cable_output_binding_t( &
        name="SWdown", &
        aggregator=new_aggregator(met%ofsd), &
        range=ranges%SWdown &
      ), &
      cable_output_binding_t( &
        name="LWdown", &
        aggregator=new_aggregator(met%fld), &
        range=ranges%LWdown &
      ), &
      cable_output_binding_t( &
        name="Rainf", &
        aggregator=new_aggregator(met%precip), &
        divide_by=dels, &
        range=ranges%Rainf &
      ), &
      cable_output_binding_t( &
        name="Snowf", &
        aggregator=new_aggregator(met%precip_sn), &
        divide_by=dels, &
        range=ranges%Snowf &
      ), &
      cable_output_binding_t( &
        name="PSurf", &
        aggregator=new_aggregator(met%pmb), &
        range=ranges%PSurf &
      ), &
      cable_output_binding_t( &
        name="Tair", &
        aggregator=new_aggregator(met%tk), &
        range=ranges%Tair &
      ), &
      cable_output_binding_t( &
        name="Qair", &
        aggregator=new_aggregator(met%qv), &
        range=ranges%Qair &
      ), &
      cable_output_binding_t( &
        name="Wind", &
        aggregator=new_aggregator(met%ua), &
        range=ranges%Wind &
      ), &
      cable_output_binding_t( &
        name="CO2air", &
        aggregator=new_aggregator(met%ca), &
        scale_by=1e6, &
        range=ranges%CO2air &
      ), &
      cable_output_binding_t( &
        name="Qmom", &
        aggregator=new_aggregator(canopy%qmom), &
        range=ranges%Qmom &
      ), &
      cable_output_binding_t( &
        name="Qle", &
        aggregator=new_aggregator(canopy%fe), &
        range=ranges%Qle &
      ), &
      cable_output_binding_t( &
        name="Qh", &
        aggregator=new_aggregator(canopy%fh), &
        range=ranges%Qh &
      ), &
      cable_output_binding_t( &
        name="Qg", &
        aggregator=new_aggregator(canopy%ga), &
        range=ranges%Qg &
      ), &
      cable_output_binding_t( &
        name="Qs", &
        aggregator=new_aggregator(ssnow%rnof1), &
        divide_by=dels, &
        range=ranges%Qs &
      ), &
      cable_output_binding_t( &
        name="Qsb", &
        aggregator=new_aggregator(ssnow%rnof2), &
        divide_by=dels, &
        range=ranges%Qsb &
      ), &
      cable_output_binding_t( &
        name="Evap", &
        aggregator=new_aggregator(canopy%fe), &
        divide_by=HL, &
        range=ranges%Evap &
      ), &
      cable_output_binding_t( &
        name="PotEvap", &
        aggregator=new_aggregator(canopy%epot), &
        divide_by=dels, &
        range=ranges%PotEvap &
      ), &
      cable_output_binding_t( &
        name="ECanop", &
        aggregator=new_aggregator(canopy%fevw), &
        divide_by=HL, &
        range=ranges%ECanop &
      ), &
      cable_output_binding_t( &
        name="TVeg", &
        aggregator=new_aggregator(canopy%fevc), &
        divide_by=HL, &
        range=ranges%TVeg &
      ), &
      cable_output_binding_t( &
        name="ESoil", &
        aggregator=new_aggregator(canopy%fes), &
        divide_by=HL, &
        range=ranges%ESoil &
      ), &
      cable_output_binding_t( &
        name="HVeg", &
        aggregator=new_aggregator(canopy%fhv), &
        range=ranges%HVeg &
      ), &
      cable_output_binding_t( &
        name="HSoil", &
        aggregator=new_aggregator(canopy%fhs), &
        range=ranges%HSoil &
      ), &
      cable_output_binding_t( &
        name="RnetSoil", &
        aggregator=new_aggregator(canopy%fns), &
        range=ranges%HSoil &
      ), &
      cable_output_binding_t( &
        name="SoilMoist", &
        aggregator=new_aggregator(ssnow%wb), &
        range=ranges%SoilMoist &
      ), &
      cable_output_binding_t( &
        name="SoilMoistIce", &
        aggregator=new_aggregator(ssnow%wbice), &
        range=ranges%SoilMoist &
      ), &
      cable_output_binding_t( &
        name="SoilTemp", &
        aggregator=new_aggregator(ssnow%tgg), &
        range=ranges%SoilTemp &
      ), &
      cable_output_binding_t( &
        name="gammzz", &
        aggregator=new_aggregator(ssnow%gammzz) &
      ), &
      cable_output_binding_t( &
        name="ssdn", &
        aggregator=new_aggregator(ssnow%ssdn) &
      ), &
      cable_output_binding_t( &
        name="smass", &
        aggregator=new_aggregator(ssnow%smass) &
      ), &
      cable_output_binding_t( &
        name="BaresoilT", &
        aggregator=new_aggregator(ssnow%tgg(:, 1)), &
        range=ranges%BaresoilT &
      ), &
      cable_output_binding_t( &
        name="SWE", &
        aggregator=new_aggregator(ssnow%snowd), &
        range=ranges%SWE &
      ), &
      cable_output_binding_t( &
        name="SnowMelt", &
        aggregator=new_aggregator(ssnow%smelt), &
        divide_by=dels, &
        range=ranges%SnowMelt &
      ), &
      cable_output_binding_t( &
        name="tggsn", &
        aggregator=new_aggregator(ssnow%tggsn) &
      ), &
      cable_output_binding_t( &
        name="SnowT", &
        aggregator=new_aggregator(ssnow%tggsn(:, 1)), &
        range=ranges%SnowT &
      ), &
      cable_output_binding_t( &
        name="sdepth", &
        aggregator=new_aggregator(ssnow%sdepth), &
        range=ranges%SnowDepth &
      ), &
      cable_output_binding_t( &
        name="SnowDepth", &
        aggregator=new_aggregator(ssnow%totsdepth), &
        range=ranges%SnowDepth &
      ), &
      cable_output_binding_t( &
        name="SWnet", &
        aggregator=new_aggregator(rad%swnet), &
        range=ranges%SWnet &
      ), &
      cable_output_binding_t( &
        name="LWnet", &
        aggregator=new_aggregator(rad%lwnet), &
        range=ranges%LWnet &
      ), &
      cable_output_binding_t( &
        name="Rnet", &
        aggregator=new_aggregator(rad%rnet), &
        range=ranges%Rnet &
      ), &
      cable_output_binding_t( &
        name="Albedo", &
        aggregator=new_aggregator(rad%albedo_T), &
        range=ranges%Albedo &
      ), &
      cable_output_binding_t( &
        name="visAlbedo", &
        aggregator=new_aggregator(rad%albedo(:, 1)), &
        range=ranges%visAlbedo, &
        available=calculate_soil_albedo &
      ), &
      cable_output_binding_t( &
        name="nirAlbedo", &
        aggregator=new_aggregator(rad%albedo(:, 2)), &
        range=ranges%nirAlbedo, &
        available=calculate_soil_albedo &
      ), &
      cable_output_binding_t( &
        name="RadT", &
        aggregator=new_aggregator(rad%trad), &
        range=ranges%RadT &
      ), &
      cable_output_binding_t( &
        name="Tscrn", &
        aggregator=new_aggregator(canopy%tscrn), &
        range=ranges%Tscrn &
      ), &
      cable_output_binding_t( &
        name="Txx", &
        aggregator=new_aggregator(canopy%tscrn), &
        range=ranges%Tscrn &
      ), &
      cable_output_binding_t( &
        name="Tnn", &
        aggregator=new_aggregator(canopy%tscrn), &
        range=ranges%Tscrn &
      ), &
      cable_output_binding_t( &
        name="Tmx", &
        aggregator=new_aggregator(canopy%tscrn_max_daily%aggregated_data), &
        range=ranges%Tscrn &
      ), &
      cable_output_binding_t( &
        name="Tmn", &
        aggregator=new_aggregator(canopy%tscrn_min_daily%aggregated_data), &
        range=ranges%Tscrn &
      ), &
      cable_output_binding_t( &
        name="Qscrn", &
        aggregator=new_aggregator(canopy%qscrn), &
        range=ranges%Qscrn &
      ), &
      cable_output_binding_t( &
        name="VegT", &
        aggregator=new_aggregator(canopy%tv), &
        range=ranges%VegT &
      ), &
      cable_output_binding_t( &
        name="CanT", &
        aggregator=new_aggregator(met%tvair), &
        range=ranges%CanT &
      ), &
      cable_output_binding_t( &
        name="Fwsoil", &
        aggregator=new_aggregator(canopy%fwsoil), &
        range=ranges%Fwsoil &
      ), &
      cable_output_binding_t( &
        name="CanopInt", &
        aggregator=new_aggregator(canopy%cansto), &
        range=ranges%CanopInt &
      ), &
      cable_output_binding_t( &
        name="LAI", &
        aggregator=new_aggregator(veg%vlai), &
        range=ranges%LAI &
      ), &
      cable_output_binding_t( &
        name="Ebal", &
        aggregator=new_aggregator(bal%ebal_tot), &
        range=ranges%Ebal &
      ), &
      cable_output_binding_t( &
        name="Wbal", &
        aggregator=new_aggregator(bal%wbal_tot), &
        range=ranges%Wbal &
      ), &
      cable_output_binding_t( &
        name="wbtot0", &
        aggregator=new_aggregator(bal%wbtot0) &
      ), &
      cable_output_binding_t( &
        name="osnowd0", &
        aggregator=new_aggregator(bal%osnowd0) &
      ), &
      cable_output_binding_t( &
        name="LeafResp", &
        aggregator=new_aggregator(canopy%frday), &
        divide_by=c_molar_mass, &
        range=ranges%AutoResp &
      ), &
      cable_output_binding_t( &
        name="HeteroResp", &
        aggregator=new_aggregator(canopy%frs), &
        divide_by=c_molar_mass, &
        range=ranges%HeteroResp &
      ), &
      cable_output_binding_t( &
        name="GPP", &
        aggregator=new_aggregator(canopy%fgpp), &
        divide_by=c_molar_mass, &
        range=ranges%GPP &
      ), &
      cable_output_binding_t( &
        name="NPP", &
        aggregator=new_aggregator(canopy%fnpp), &
        divide_by=c_molar_mass, &
        range=ranges%NPP &
      ), &
      cable_output_binding_t( &
        name="AutoResp", &
        aggregator=new_aggregator(canopy%fra), &
        divide_by=c_molar_mass, &
        range=ranges%AutoResp &
      ), &
      cable_output_binding_t( &
        name="NEE", &
        aggregator=new_aggregator(canopy%fnee), &
        divide_by=c_molar_mass, &
        range=ranges%NEE &
      ), &
      cable_output_binding_t( &
        name="WatTable", &
        aggregator=new_aggregator(ssnow%wtd), &
        scale_by=1e-3, &
        range=ranges%WatTable, &
        available=use_groundwater_model &
      ), &
      cable_output_binding_t( &
        name="GWMoist", &
        aggregator=new_aggregator(ssnow%GWwb), &
        range=ranges%GWwb, &
        available=use_groundwater_model &
      ), &
      cable_output_binding_t( &
        name="SatFrac", &
        aggregator=new_aggregator(ssnow%satfrac), &
        range=ranges%SatFrac, &
        available=use_groundwater_model &
      ), &
      cable_output_binding_t( &
        name="Qrecharge", &
        aggregator=new_aggregator(ssnow%Qrecharge), &
        range=ranges%Qrecharge &
      ), &
      cable_output_binding_t( &
        name="tss", &
        aggregator=new_aggregator(ssnow%tss) &
      ), &
      cable_output_binding_t( &
        name="rtsoil", &
        aggregator=new_aggregator(ssnow%rtsoil) &
      ), &
      cable_output_binding_t( &
        name="runoff", &
        aggregator=new_aggregator(ssnow%runoff) &
      ), &
      cable_output_binding_t( &
        name="ssdnn", &
        aggregator=new_aggregator(ssnow%ssdnn) &
      ), &
      cable_output_binding_t( &
        name="snage", &
        aggregator=new_aggregator(ssnow%snage) &
      ), &
      cable_output_binding_t( &
        name="osnowd", &
        aggregator=new_aggregator(ssnow%osnowd) &
      ), &
      cable_output_binding_t( &
        name="albsoilsn", &
        aggregator=new_aggregator(ssnow%albsoilsn), &
        range=ranges%albsoiln &
      ), &
      cable_output_binding_t( &
        name="isflag", &
        aggregator=new_aggregator(ssnow%isflag) &
      ), &
      cable_output_binding_t( &
        name="ghflux", &
        aggregator=new_aggregator(canopy%ghflux) &
      ), &
      cable_output_binding_t( &
        name="sghflux", &
        aggregator=new_aggregator(canopy%sghflux) &
      ), &
      cable_output_binding_t( &
        name="dgdtg", &
        aggregator=new_aggregator(canopy%dgdtg) &
      ), &
      cable_output_binding_t( &
        name="fev", &
        aggregator=new_aggregator(canopy%fev) &
      ), &
      cable_output_binding_t( &
        name="fes", &
        aggregator=new_aggregator(canopy%fes) &
      ), &
      cable_output_binding_t( &
        name="albedo", &
        aggregator=new_aggregator(rad%albedo) &
      ), &
      cable_output_binding_t( &
        name="iveg", &
        aggregator=new_aggregator(veg%iveg), &
        range=ranges%iveg &
      ), &
      cable_output_binding_t( &
        name="patchfrac", &
        aggregator=new_aggregator(patch(:)%frac), &
        range=ranges%patchfrac &
      ), &
      cable_output_binding_t( &
        name="isoil", &
        aggregator=new_aggregator(soil%isoilm), &
        range=ranges%isoil &
      ), &
      cable_output_binding_t( &
        name="bch", &
        aggregator=new_aggregator(soil%bch), &
        range=ranges%bch &
      ), &
      cable_output_binding_t( &
        name="clay", &
        aggregator=new_aggregator(soil%clay), &
        range=ranges%clay &
      ), &
      cable_output_binding_t( &
        name="sand", &
        aggregator=new_aggregator(soil%sand), &
        range=ranges%sand &
      ), &
      cable_output_binding_t( &
        name="silt", &
        aggregator=new_aggregator(soil%silt), &
        range=ranges%silt &
      ), &
      cable_output_binding_t( &
        name="ssat", &
        aggregator=new_aggregator(soil%ssat), &
        range=ranges%ssat &
      ), &
      cable_output_binding_t( &
        name="sfc", &
        aggregator=new_aggregator(soil%sfc), &
        range=ranges%sfc &
      ), &
      cable_output_binding_t( &
        name="swilt", &
        aggregator=new_aggregator(soil%swilt), &
        range=ranges%swilt &
      ), &
      cable_output_binding_t( &
        name="hyds", &
        aggregator=new_aggregator(soil%hyds), &
        range=ranges%hyds &
      ), &
      cable_output_binding_t( &
        name="sucs", &
        aggregator=new_aggregator(soil%sucs), &
        range=ranges%sucs &
      ), &
      cable_output_binding_t( &
        name="css", &
        aggregator=new_aggregator(soil%css), &
        range=ranges%css &
      ), &
      cable_output_binding_t( &
        name="rhosoil", &
        aggregator=new_aggregator(soil%rhosoil), &
        range=ranges%rhosoil &
      ), &
      cable_output_binding_t( &
        name="rs20", &
        aggregator=new_aggregator(veg%rs20), &
        range=ranges%rs20 &
      ), &
      cable_output_binding_t( &
        name="albsoil", &
        aggregator=new_aggregator(soil%albsoil), &
        range=ranges%albsoil &
      ), &
      cable_output_binding_t( &
        name="hc", &
        aggregator=new_aggregator(veg%hc), &
        range=ranges%hc &
      ), &
      cable_output_binding_t( &
        name="canst1", &
        aggregator=new_aggregator(veg%canst1), &
        range=ranges%canst1 &
      ), &
      cable_output_binding_t( &
        name="dleaf", &
        aggregator=new_aggregator(veg%dleaf), &
        range=ranges%dleaf &
      ), &
      cable_output_binding_t( &
        name="frac4", &
        aggregator=new_aggregator(veg%frac4), &
        range=ranges%frac4 &
      ), &
      cable_output_binding_t( &
        name="ejmax", &
        aggregator=new_aggregator(veg%ejmax), &
        range=ranges%ejmax &
      ), &
      cable_output_binding_t( &
        name="vcmax", &
        aggregator=new_aggregator(veg%vcmax), &
        range=ranges%vcmax &
      ), &
      cable_output_binding_t( &
        name="rp20", &
        aggregator=new_aggregator(veg%rp20), &
        range=ranges%rp20 &
      ), &
      cable_output_binding_t( &
        name="g0", &
        aggregator=new_aggregator(veg%g0), &
        range=ranges%g0 &
      ), &
      cable_output_binding_t( &
        name="g1", &
        aggregator=new_aggregator(veg%g1), &
        range=ranges%g1 &
      ), &
      cable_output_binding_t( &
        name="rpcoef", &
        aggregator=new_aggregator(veg%rpcoef), &
        range=ranges%rpcoef &
      ), &
      cable_output_binding_t( &
        name="shelrb", &
        aggregator=new_aggregator(veg%shelrb), &
        range=ranges%shelrb &
      ), &
      cable_output_binding_t( &
        name="xfang", &
        aggregator=new_aggregator(veg%xfang), &
        range=ranges%xfang &
      ), &
      cable_output_binding_t( &
        name="wai", &
        aggregator=new_aggregator(veg%wai), &
        range=ranges%wai &
      ), &
      cable_output_binding_t( &
        name="vegcf", &
        aggregator=new_aggregator(veg%vegcf), &
        range=ranges%vegcf &
      ), &
      cable_output_binding_t( &
        name="extkn", &
        aggregator=new_aggregator(veg%extkn), &
        range=ranges%extkn &
      ), &
      cable_output_binding_t( &
        name="tminvj", &
        aggregator=new_aggregator(veg%tminvj), &
        range=ranges%tminvj &
      ), &
      cable_output_binding_t( &
        name="tmaxvj", &
        aggregator=new_aggregator(veg%tmaxvj), &
        range=ranges%tmaxvj &
      ), &
      cable_output_binding_t( &
        name="vbeta", &
        aggregator=new_aggregator(veg%vbeta), &
        range=ranges%vbeta &
      ), &
      cable_output_binding_t( &
        name="xalbnir", &
        aggregator=new_aggregator(veg%xalbnir), &
        range=ranges%xalbnir &
      ), &
      cable_output_binding_t( &
        name="meth", &
        aggregator=new_aggregator(veg%meth), &
        range=ranges%meth &
      ), &
      cable_output_binding_t( &
        name="za_uv", &
        aggregator=new_aggregator(rough%za_uv), &
        range=ranges%za &
      ), &
      cable_output_binding_t( &
        name="za_tq", &
        aggregator=new_aggregator(rough%za_tq), &
        range=ranges%za &
      ), &
      cable_output_binding_t( &
        name="ratecp", &
        aggregator=new_aggregator(bgc%ratecp), &
        range=ranges%ratecp &
      ), &
      cable_output_binding_t( &
        name="ratecs", &
        aggregator=new_aggregator(bgc%ratecs), &
        range=ranges%ratecs &
      ), &
      cable_output_binding_t( &
        name="cplant", &
        aggregator=new_aggregator(bgc%cplant) &
      ), &
      cable_output_binding_t( &
        name="csoil", &
        aggregator=new_aggregator(bgc%csoil) &
      ), &
      cable_output_binding_t( &
        name="zse", &
        aggregator=new_aggregator(soil%zse), &
        range=ranges%zse &
      ), &
      cable_output_binding_t( &
        name="froot", &
        aggregator=new_aggregator(veg%froot), &
        range=ranges%froot &
      ), &
      cable_output_binding_t( &
        name="GWdz", &
        aggregator=new_aggregator(soil%GWdz), &
        range=ranges%GWdz, &
        available=use_groundwater_model &
      ), &
      cable_output_binding_t( &
        name="Qhmax", &
        aggregator=new_aggregator(gw_params%MaxHorzDrainRate), &
        range=ranges%gw_default, &
        available=use_groundwater_model &
      ), &
      cable_output_binding_t( &
        name="QhmaxEfold", &
        aggregator=new_aggregator(gw_params%EfoldHorzDrainRate), &
        range=ranges%gw_default, &
        available=use_groundwater_model &
      ), &
      cable_output_binding_t( &
        name="SatFracmax", &
        aggregator=new_aggregator(gw_params%MaxSatFraction), &
        range=ranges%gw_default, &
        available=use_groundwater_model &
      ), &
      cable_output_binding_t( &
        name="HKefold", &
        aggregator=new_aggregator(gw_params%hkrz), &
        range=ranges%gw_default, &
        available=use_groundwater_model &
      ), &
      cable_output_binding_t( &
        name="HKdepth", &
        aggregator=new_aggregator(gw_params%zdepth), &
        range=ranges%gw_default, &
        available=use_groundwater_model &
      ), &
      cable_output_binding_t( &
        name="nap", &
        aggregator=new_aggregator(landpt_global(:)%nap) &
      ), &
      cable_output_binding_t( &
        name="mvtype", &
        aggregator=new_aggregator(mvtype) &
      ), &
      cable_output_binding_t( &
        name="mstype", &
        aggregator=new_aggregator(mstype) &
      ) &
    ]
  end function cable_output_bindings

end module cable_output_bindings_mod
