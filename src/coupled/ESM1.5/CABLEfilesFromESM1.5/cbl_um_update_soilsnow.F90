MODULE cbl_um_update_soilsnow_mod
   
IMPLICIT NONE

CONTAINS

SUBROUTINE update_soilsnow( mp, soil, ssnow, veg_iveg )
                            

USE cable_def_types_mod,       ONLY: soil_parameter_type
USE cable_def_types_mod,       ONLY: soil_snow_type
USE cable_phys_constants_mod,  ONLY: density_ice, density_liq
USE cable_surface_types_mod,   ONLY: lakes_cable
IMPLICIT NONE

INTEGER, INTENT(IN) :: mp                ! # active land points
INTEGER, INTENT(IN) :: veg_iveg(mp)
!jhan:build fudge inOUT
TYPE(soil_snow_type), INTENT(inOUT)     :: ssnow      ! 
TYPE(soil_parameter_type), INTENT(IN) :: soil       ! soil parameters
REAL :: wbtot1(mp), wbtot2(mp) 
INTEGER :: j 

ssnow%wb_lake = 0.0
wbtot1        = 0.0
wbtot2        = 0.0

ssnow%wbliq = ssnow%wb - ssnow%wbice

WHERE( veg_iveg == lakes_cable .AND. ssnow%wb(:,1) < soil%sfc ) 
  wbtot1  = wbtot1 + REAL( ssnow%wb(:,1) ) * density_liq * soil%zse(1)
  ssnow%wb(:,1) = soil%sfc
  wbtot2  = wbtot2 + REAL( ssnow%wb(:,1) ) * density_liq * soil%zse(1)
ENDWHERE
  
ssnow%wb_lake = MAX( wbtot2 - wbtot1, 0.)

RETURN
END SUBROUTINE update_soilsnow

END MODULE cbl_um_update_soilsnow_mod


