! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

module cable_output_bindings_casa_mod
  !! Bindings between output catalogue names and CABLE model state.
  !!
  !! Everything that describes a variable lives in `output_catalogue.yaml`.
  !! This module supplies only what Fortran must provide: the model variable
  !! each name points at, its unit conversion, its valid range and whether it
  !! exists in the current model configuration. All model state and
  !! configuration is passed in as arguments.

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

  use cable_timing_mod, only: seconds_per_day

  use cable_phys_constants_mod, only: c_molar_mass

  use cable_checks_module, only: ranges_type

  use aggregator_mod, only: new_aggregator

  use cable_output_mod, only: cable_output_binding_t

  implicit none
  private

  public :: cable_output_bindings_casa

contains

  function cable_output_bindings_casa( &
    casaflux, casapool, casamet, ranges, casa_enabled, use_popluc &
  ) result(bindings)
    type(casa_flux), intent(in) :: casaflux
    type(casa_pool), intent(in) :: casapool
    type(casa_met), intent(in) :: casamet
    type(ranges_type), intent(in) :: ranges
    logical, intent(in) :: casa_enabled
    logical, intent(in) :: use_popluc
    type(cable_output_binding_t), allocatable :: bindings(:)

    ! Without this component none of its variables exist. They stay in the list, unavailable.
    if (.not. casa_enabled) then
      bindings = [ &
        cable_output_binding_t(name="RootResp", available=.false.), &
        cable_output_binding_t(name="StemResp", available=.false.), &
        cable_output_binding_t(name="NBP", available=.false.), &
        cable_output_binding_t(name="dCdt", available=.false.), &
        cable_output_binding_t(name="TotSoilCarb", available=.false.), &
        cable_output_binding_t(name="TotLittCarb", available=.false.), &
        cable_output_binding_t(name="SoilCarbFast", available=.false.), &
        cable_output_binding_t(name="SoilCarbSlow", available=.false.), &
        cable_output_binding_t(name="SoilCarbPassive", available=.false.), &
        cable_output_binding_t(name="LittCarbMetabolic", available=.false.), &
        cable_output_binding_t(name="LittCarbStructural", available=.false.), &
        cable_output_binding_t(name="LittCarbCWD", available=.false.), &
        cable_output_binding_t(name="PlantCarbLeaf", available=.false.), &
        cable_output_binding_t(name="PlantCarbWood", available=.false.), &
        cable_output_binding_t(name="PlantCarbFineRoot", available=.false.), &
        cable_output_binding_t(name="TotLivBiomass", available=.false.), &
        cable_output_binding_t(name="PlantTurnover", available=.false.), &
        cable_output_binding_t(name="PlantTurnoverLeaf", available=.false.), &
        cable_output_binding_t(name="PlantTurnoverWood", available=.false.), &
        cable_output_binding_t(name="PlantTurnoverFineRoot", available=.false.), &
        cable_output_binding_t(name="PlantTurnoverWoodDist", available=.false.), &
        cable_output_binding_t(name="PlantTurnoverWoodCrowding", available=.false.), &
        cable_output_binding_t(name="PlantTurnoverWoodResourceLim", available=.false.), &
        cable_output_binding_t(name="Area", available=.false.), &
        cable_output_binding_t(name="LandUseFlux", available=.false.) &
      ]
      return
    end if

    ! One entry per CASA variable, in catalogue order. Each entry has the same
    ! meaning as in cable_output_bindings_mod: the model variable (through
    ! the aggregator), its unit conversion, its valid range and its availability.
    ! Many conversions here are divide_by=(seconds_per_day * c_molar_mass), which
    ! combines two steps: seconds_per_day turns a per-day amount into a per-second
    ! one (CASA accumulates per day), and c_molar_mass, in grams per micromole of
    ! carbon, turns grams into micromoles. The result matches the catalogue's
    ! units of umol/m^2/s.
    bindings = [ &
      cable_output_binding_t( &
        name="RootResp", &
        aggregator=new_aggregator(casaflux%crmplant(:, FROOT)), &
        divide_by=(seconds_per_day * c_molar_mass), &
        range=ranges%AutoResp &
      ), &
      cable_output_binding_t( &
        name="StemResp", &
        aggregator=new_aggregator(casaflux%crmplant(:, WOOD)), &
        divide_by=(seconds_per_day * c_molar_mass), &
        range=ranges%AutoResp &
      ), &
      cable_output_binding_t( &
        name="NBP", &
        aggregator=new_aggregator(casaflux%cnbp), &
        divide_by=(seconds_per_day * c_molar_mass), &
        range=ranges%NEE &
      ), &
      cable_output_binding_t( &
        name="dCdt", &
        aggregator=new_aggregator(casapool%dCdt), &
        divide_by=(seconds_per_day * c_molar_mass), &
        range=ranges%NEE &
      ), &
      cable_output_binding_t( &
        name="TotSoilCarb", &
        aggregator=new_aggregator(casapool%csoiltot), &
        divide_by=1e3, &
        range=ranges%TotSoilCarb &
      ), &
      cable_output_binding_t( &
        name="TotLittCarb", &
        aggregator=new_aggregator(casapool%clittertot), &
        divide_by=1e3, &
        range=ranges%TotLittCarb &
      ), &
      cable_output_binding_t( &
        name="SoilCarbFast", &
        aggregator=new_aggregator(casapool%csoil(:, MIC)), &
        divide_by=1e3, &
        range=ranges%TotLittCarb &
      ), &
      cable_output_binding_t( &
        name="SoilCarbSlow", &
        aggregator=new_aggregator(casapool%csoil(:, SLOW)), &
        divide_by=1e3, &
        range=ranges%TotSoilCarb &
      ), &
      cable_output_binding_t( &
        name="SoilCarbPassive", &
        aggregator=new_aggregator(casapool%csoil(:, PASS)), &
        divide_by=1e3, &
        range=ranges%TotSoilCarb &
      ), &
      cable_output_binding_t( &
        name="LittCarbMetabolic", &
        aggregator=new_aggregator(casapool%clitter(:, METB)), &
        divide_by=1e3, &
        range=ranges%TotLittCarb &
      ), &
      cable_output_binding_t( &
        name="LittCarbStructural", &
        aggregator=new_aggregator(casapool%clitter(:, STR)), &
        divide_by=1e3, &
        range=ranges%TotLittCarb &
      ), &
      cable_output_binding_t( &
        name="LittCarbCWD", &
        aggregator=new_aggregator(casapool%clitter(:, CWD)), &
        divide_by=1e3, &
        range=ranges%TotLittCarb &
      ), &
      cable_output_binding_t( &
        name="PlantCarbLeaf", &
        aggregator=new_aggregator(casapool%cplant(:, LEAF)), &
        divide_by=1e3, &
        range=ranges%TotLittCarb &
      ), &
      cable_output_binding_t( &
        name="PlantCarbWood", &
        aggregator=new_aggregator(casapool%cplant(:, WOOD)), &
        divide_by=1e3, &
        range=ranges%TotLittCarb &
      ), &
      cable_output_binding_t( &
        name="PlantCarbFineRoot", &
        aggregator=new_aggregator(casapool%cplant(:, FROOT)), &
        divide_by=1e3, &
        range=ranges%TotLittCarb &
      ), &
      cable_output_binding_t( &
        name="TotLivBiomass", &
        aggregator=new_aggregator(casapool%cplanttot), &
        divide_by=1e3, &
        range=ranges%TotLivBiomass &
      ), &
      cable_output_binding_t( &
        name="PlantTurnover", &
        aggregator=new_aggregator(casaflux%cplant_turnover_tot), &
        divide_by=(seconds_per_day * c_molar_mass), &
        range=ranges%NEE &
      ), &
      cable_output_binding_t( &
        name="PlantTurnoverLeaf", &
        aggregator=new_aggregator(casaflux%Cplant_turnover(:, LEAF)), &
        divide_by=(seconds_per_day * c_molar_mass), &
        range=ranges%NEE &
      ), &
      cable_output_binding_t( &
        name="PlantTurnoverWood", &
        aggregator=new_aggregator(casaflux%Cplant_turnover(:, WOOD)), &
        divide_by=(seconds_per_day * c_molar_mass), &
        range=ranges%NEE &
      ), &
      cable_output_binding_t( &
        name="PlantTurnoverFineRoot", &
        aggregator=new_aggregator(casaflux%Cplant_turnover(:, FROOT)), &
        divide_by=(seconds_per_day * c_molar_mass), &
        range=ranges%NEE &
      ), &
      cable_output_binding_t( &
        name="PlantTurnoverWoodDist", &
        aggregator=new_aggregator(casaflux%Cplant_turnover_disturbance), &
        divide_by=(seconds_per_day * c_molar_mass), &
        range=ranges%NEE &
      ), &
      cable_output_binding_t( &
        name="PlantTurnoverWoodCrowding", &
        aggregator=new_aggregator(casaflux%Cplant_turnover_crowding), &
        divide_by=(seconds_per_day * c_molar_mass), &
        range=ranges%NEE &
      ), &
      cable_output_binding_t( &
        name="PlantTurnoverWoodResourceLim", &
        aggregator=new_aggregator(casaflux%Cplant_turnover_resource_limitation), &
        divide_by=(seconds_per_day * c_molar_mass), &
        range=ranges%NEE &
      ), &
      cable_output_binding_t( &
        name="Area", &
        aggregator=new_aggregator(casamet%areacell), &
        divide_by=1e6, &
        range=ranges%Area &
      ), &
      cable_output_binding_t( &
        name="LandUseFlux", &
        aggregator=new_aggregator(casaflux%FluxCtoLUC), &
        divide_by=(seconds_per_day * c_molar_mass), &
        range=ranges%NEE, &
        available=use_popluc &
      ) &
    ]
  end function cable_output_bindings_casa

end module cable_output_bindings_casa_mod
