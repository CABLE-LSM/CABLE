MODULE cable_air_type_mod

IMPLICIT NONE

  ! Air variables:
TYPE air_type

  REAL, DIMENSION(:), ALLOCATABLE ::                                          &
       rho,     & ! dry air density (kg m-3)
       volm,    & ! molar volume (m3 mol-1)
       rlam,    & ! latent heat for water (j/kg)
       qsat,    & ! saturation specific humidity
       epsi,    & ! d(qsat)/dT ((kg/kg)/K)
       visc,    & ! air kinematic viscosity (m2/s)
       psyc,    & ! psychrometric constant
       dsatdk,  & ! d(es)/dT (mb/K)
       cmolar     ! conv. from m/s to mol/m2/s

END TYPE air_type

!Instantiation:
TYPE(air_type) :: air_cbl

END MODULE cable_air_type_mod
