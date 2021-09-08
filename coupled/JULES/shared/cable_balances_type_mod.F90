MODULE cable_balances_type_mod

IMPLICIT NONE

  ! Energy and water balance variables:
TYPE balances_type

  REAL, DIMENSION(:), POINTER ::                                              &
       drybal,           & ! energy balance for dry canopy
       ebal,             & ! energy balance per time step (W/m^2)
       ebal_tot,         & ! cumulative energy balance (W/m^2)
       ebal_cncheck,     & ! energy balance consistency check (W/m^2)
       ebal_tot_cncheck, & ! cumulative energy balance (W/m^2)
       ebaltr,           & ! energy balance per time step (W/m^2)
       ebal_tottr,       & ! cumulative energy balance (W/m^2)
       evap_tot,         & ! cumulative evapotranspiration (mm/dels)
       osnowd0,          & ! snow depth, first time step
       precip_tot,       & ! cumulative precipitation (mm/dels)
       rnoff_tot,        & ! cumulative runoff (mm/dels)
       wbal,             & ! water balance per time step (mm/dels)
       wbal_tot,         & ! cumulative water balance (mm/dels)
       wbtot0,           & ! total soil water (mm), first time step
       wetbal,           & ! energy balance for wet canopy
       cansto0,          & ! canopy water storage (mm)
       owbtot,           & ! total soil water (mm), first time step
       evapc_tot,        & ! cumulative evapotranspiration (mm/dels)
       evaps_tot,        & ! cumulative evapotranspiration (mm/dels)
       rnof1_tot,        & ! cumulative runoff (mm/dels)
       rnof2_tot,        & ! cumulative runoff (mm/dels)
       snowdc_tot,       & ! cumulative runoff (mm/dels)
       wbal_tot1,        & ! cumulative water balance (mm/dels)
       delwc_tot,        & ! energy balance for wet canopy
       qasrf_tot,        & ! heat advected to the snow by precip.
       qfsrf_tot,        & ! energy of snowpack phase changes
       qssrf_tot, &        ! energy of snowpack phase changes
       Radbal,                                                                &
       EbalSoil,                                                              &
       Ebalveg,                                                               &
       Radbalsum

END TYPE balances_type

!Instantiation:
TYPE(balances_type) :: bal_cbl

END MODULE cable_balances_type_mod
