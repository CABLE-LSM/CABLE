MODULE cable_met_type_mod

IMPLICIT NONE

  ! Meterological data:
TYPE met_type

  INTEGER, DIMENSION(:), POINTER ::                                           &
       year,    & ! local time year AD
       moy        ! local time month of year

  REAL, DIMENSION(:), POINTER ::                                              &
       ca,      & ! CO2 concentration (mol/mol)
       doy,     & ! local time day of year = days since 0 hr 1st Jan
       hod,     & ! local hour of day
       ofsd,    & ! downward short-wave radiation (W/m2)
       fld,     & ! downward long-wave radiation (W/m2)
       precip,  & ! rainfall (liquid+solid)(mm/dels)
       precip_sn,&! solid preipitation only (mm/dels)
       tk,      & ! surface air temperature (oK)
       tvair,   & ! within canopy air temperature (oK)
       tvrad,   & ! radiative vegetation temperature (K)
       pmb,     & ! surface air pressure (mbar)
       ua,      & ! surface wind speed (m/s)
       qv,      & ! surface specific humidity (g/g)
       qvair,   & ! within canopy specific humidity (g/g)
       da,      & ! water vap pressure deficit at ref height (Pa)
       dva,     & ! in canopy water vap pressure deficit (Pa)
       coszen,   &  ! cos(zenith angle of sun)
       Ndep,     &   ! nitrogen deposition (gN m-2 d-1)
       Pdep ! P deposition (gP m-2 d-1)

  REAL, DIMENSION(:,:), POINTER ::                                            &
       fsd  ! downward short-wave radiation (W/m2)

END TYPE met_type

!Instantiation:
TYPE(met_type) :: met_cbl

END MODULE cable_met_type_mod
