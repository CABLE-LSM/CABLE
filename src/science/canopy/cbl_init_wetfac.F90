MODULE cable_init_wetfac_mod
   
IMPLICIT NONE

CONTAINS

SUBROUTINE initialize_wetfac( mp, ssnow_wetfac, soil_swilt, soil_sfc,          &
                              ssnow_wb, ssnow_wbice, ssnow_snowd,              &
                              veg_iveg, met_tk, Ctfrz )

USE cable_surface_types_mod, ONLY: lakes => lakes_cable        
USE cable_def_types_mod,  ONLY : r_2
                          
IMPLICIT NONE                     

INTEGER, INTENT(IN) :: mp                ! # active land points
REAL, INTENT(INOUT) :: ssnow_wetfac(mp)

REAL, INTENT(IN) :: soil_swilt(mp)
REAL, INTENT(IN) :: soil_sfc(mp)
REAL, INTENT(IN) :: ssnow_wb(mp)
REAL, INTENT(IN) :: ssnow_wbice(mp)
REAL, INTENT(IN) :: ssnow_snowd(mp)
REAL, INTENT(IN) :: met_tk(mp)
REAL, INTENT(IN) :: Ctfrz 
INTEGER, INTENT(IN) :: veg_iveg(mp)

REAL :: wilting_pt(mp)
REAL :: wetfac_num(mp)
REAL :: wetfac_den(mp)
REAL :: ice_ratio
REAL :: ice_factor
INTEGER :: i

wilting_pt(:) = soil_swilt(:)/2.0   
wetfac_num(:) = REAL( ssnow_wb(:) ) - wilting_pt(:) 
wetfac_den(:) = REAL( soil_sfc(:) - wilting_pt(:) ) 
wetfac_den(:) = MAX( 0.0830, wetfac_den(:) )

ssnow_wetfac(:) = wetfac_num(:)  / wetfac_den(:)
ssnow_wetfac(:) = MIN( 1.0, ssnow_wetfac(:) )
ssnow_wetfac(:) = MAX( 0.0, ssnow_wetfac(:) )

DO i=1,mp
  ! Ultimately reduces surface wetness consideihg wetness locked up in ice
  IF( ssnow_wbice(i) > 0.0 ) THEN
    !ice moistrue / total  moiisture (** ?) 
    ice_ratio  = ( ssnow_wbice(i) / ssnow_wb(i) )**2                  
    !~ 1-Ice_ratio^2 
    ice_factor = 1._r_2 - MIN( 0.2_r_2, ice_ratio )
    ice_factor = REAL( MAX( 0.5_r_2, ice_factor )  )
    ssnow_wetfac(i) = ssnow_wetfac(i) * ice_factor 
  END IF

  IF( ssnow_snowd(i) > 0.1) THEN 
    ssnow_wetfac(i) = 0.9
  END IF

  IF ( veg_iveg(i) == lakes ) THEN 
    IF ( met_tk(i) >= Ctfrz + 5. ) THEN
      ssnow_wetfac(i) = 1.0
    END IF
    IF( met_tk(i) < Ctfrz + 5. ) THEN
      ssnow_wetfac(i) = 0.7
    END IF
  END IF

ENDDO

RETURN
END SUBROUTINE initialize_wetfac  

End MODULE cable_init_wetfac_mod














