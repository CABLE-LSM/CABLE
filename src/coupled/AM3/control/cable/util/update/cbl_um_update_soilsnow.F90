MODULE cbl_um_update_soilsnow_mod
   
IMPLICIT NONE

CONTAINS

SUBROUTINE update_soilsnow( mp, soil, ssnow, veg_iveg )
                            

USE cable_def_types_mod,       ONLY: soil_parameter_type
USE cable_def_types_mod,       ONLY: soil_snow_type
USE cable_surface_types_mod,   ONLY: lakes_cable
IMPLICIT NONE

INTEGER, INTENT(IN) :: mp                ! # active land points
INTEGER, INTENT(IN) :: veg_iveg(mp)
!jhan:build fudge inOUT
TYPE(soil_snow_type), INTENT(inOUT)     :: ssnow      ! 
TYPE(soil_parameter_type), INTENT(IN) :: soil       ! soil parameters

INTEGER :: j 

ssnow%wb_lake = 0.0
ssnow%wbtot1  = 0.0
ssnow%wbtot2  = 0.0

ssnow%wbliq = ssnow%wb - ssnow%wbice

!jhan:this was in a do loop from 1:1
! lakes: remove hard-wired number in future version
WHERE( veg_iveg == lakes_cable .AND. ssnow%wb(:,1) < soil%sfc ) 

  ssnow%wbtot1  = ssnow%wbtot1 + REAL( ssnow%wb(:,1) ) * 1000.0 *     &
                  soil%zse(1)
  ssnow%wb(:,1) = soil%sfc
  ssnow%wbtot2  = ssnow%wbtot2 + REAL( ssnow%wb(:,1) ) * 1000.0 *     &
                  soil%zse(1)
ENDWHERE
  
ssnow%wb_lake = MAX( ssnow%wbtot2 - ssnow%wbtot1, 0.)

RETURN
END SUBROUTINE update_soilsnow

END MODULE cbl_um_update_soilsnow_mod


