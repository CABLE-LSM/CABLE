MODULE remove_trans_mod

USE cbl_ssnow_data_mod

PUBLIC  remove_trans

CONTAINS

SUBROUTINE remove_trans(dels, soil, ssnow, canopy, veg)

    USE cable_common_module, ONLY : cable_user

    ! Removes transpiration water from soil.
    REAL, INTENT(IN)                    :: dels ! integration time step (s)
    TYPE(canopy_type), INTENT(INOUT)         :: canopy
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil
    TYPE(veg_parameter_type), INTENT(INOUT)  :: veg
    INTEGER i, k

    IF (cable_user%FWSOIL_switch.NE.'Haverd2013') THEN

      ssnow%evapfbl = 0.0_r_2
      
      DO i = 1, mp
         ssnow%evapfbl(i,:) = trans_soil_water(dels, soil%swilt(i),            &
                 veg%froot(i,:), soil%zse, canopy%fevc(i), ssnow%wb(i,:))
         ssnow%wb(i,:) = ssnow%wb(i,:) - ssnow%evapfbl(i,:) / (soil%zse(:)*Cdensity_liq)
         IF ( canopy%fevc(i) > 0.0_r_2 ) THEN
            canopy%fevc(i) = SUM(ssnow%evapfbl(i,:)) * CHL / dels
         END IF
      END DO

    ELSE
       WHERE (canopy%fevc .LT. 0.0_r_2)
          canopy%fevw = canopy%fevw+canopy%fevc
          canopy%fevc = 0.0_r_2
       END WHERE
       DO k = 1,ms
          ssnow%wb(:,k) = ssnow%wb(:,k) - ssnow%evapfbl(:,k)/(soil%zse(k)*Cdensity_liq)
       ENDDO


    ENDIF

  END SUBROUTINE remove_trans

  FUNCTION trans_soil_water(dels, swilt, froot, zse, fevc, wb) RESULT(evapfbl)

    !! Calculates the amount of water removed from the soil by transpiration.
    !
    REAL, INTENT(IN)                     :: dels
       !! integration time step (s)
    REAL, INTENT(IN)                     :: swilt 
       !! wilting point (m3/m3)
    REAL(r_2), INTENT(IN)                :: fevc 
       !! transpiration (kg/m2/s)
    REAL, DIMENSION(ms), INTENT(IN)      :: froot 
       !! root fraction (-)
    REAL, DIMENSION(ms), INTENT(IN)      :: zse 
       !! soil depth (m)
    REAL(r_2), DIMENSION(ms), INTENT(IN) :: wb 
       !! water balance (m3/m3)
   
    ! Local variables
    REAL(r_2), DIMENSION(ms)   :: evapfbl
    REAL(r_2), DIMENSION(0:ms) :: diff
    REAL(r_2)                  :: xx,xxd
    INTEGER k
      
    xx = 0.; xxd = 0.; diff(:) = 0.
    IF (fevc > 0.0) THEN
      DO k = 1,ms
        ! Removing transpiration from soil:

        ! Calculate the amount (perhaps moisture/ice limited)
        ! that can be removed:
        xx = fevc * dels / CHL * froot(k) + diff(k-1)   ! kg/m2
        diff(k) = MAX( 0.0_r_2, wb(k) - 1.1 * swilt) * zse(k)*Cdensity_liq
        xxd = xx - diff(k)

        IF ( xxd > 0.0 ) THEN
          evapfbl(k) = diff(k)
          diff(k)    = xxd
        ELSE
          evapfbl(k) = xx
          diff(k)    = 0.0
        END IF

      END DO
    END IF
   
   END FUNCTION trans_soil_water

END MODULE remove_trans_mod
