MODULE init_cable_state_mod

IMPLICIT NONE

PRIVATE

PUBLIC :: init_cable_state 

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='init_cable_params_mod'

CONTAINS

SUBROUTINE init_cable_state( mp ) 

!subr:
USE allocate_cable_state_mod, ONLY: alloc_cable_state 
!data:
USE cable_air_type_mod,       ONLY: air_cbl
USE cable_met_type_mod,       ONLY: met_cbl
USE cable_radiation_type_mod, ONLY: rad_cbl
USE cable_roughness_type_mod, ONLY: rough_cbl
USE cable_canopy_type_mod,    ONLY: canopy_cbl
USE cable_soil_snow_type_mod, ONLY: ssnow_cbl
USE cable_bgc_pool_type_mod,  ONLY: bgc_cbl
USE cable_balances_type_mod,  ONLY: bal_cbl
USE cable_sum_flux_type_mod,  ONLY: sum_flux_cbl

IMPLICIT NONE

INTEGER :: mp

CALL alloc_cable_state ( mp, air_cbl, met_cbl, rad_cbl, rough_cbl, canopy_cbl,&
                         ssnow_cbl, bgc_cbl, bal_cbl, sum_flux_cbl )

END SUBROUTINE init_cable_state

END MODULE init_cable_state_mod
