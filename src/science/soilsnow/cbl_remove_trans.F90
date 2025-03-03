MODULE remove_trans_mod

USE cbl_ssnow_data_mod

PUBLIC  remove_trans

CONTAINS

SUBROUTINE remove_trans(dels, soil, ssnow, canopy, veg)

    USE cable_common_module, ONLY : redistrb, cable_user

    ! Removes transpiration water from soil.
    REAL, INTENT(IN)                    :: dels ! integration time step (s)
    TYPE(canopy_type), INTENT(INOUT)         :: canopy
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil
    TYPE(veg_parameter_type), INTENT(INOUT)  :: veg
    REAL(r_2), DIMENSION(mp,0:ms) :: diff
    REAL(r_2), DIMENSION(mp)      :: xx,xxd,evap_cur
    INTEGER k

    IF (cable_user%FWSOIL_switch.NE.'Haverd2013') THEN
       xx = 0.; xxd = 0.; diff(:,:) = 0.
       DO k = 1,ms

          ! Removing transpiration from soil:
          WHERE (canopy%fevc > 0.0 )     ! convert to mm/dels

             ! Calculate the amount (perhaps moisture/ice limited)
             ! which can be removed:
             xx = canopy%fevc * dels / CHL * veg%froot(:,k) + diff(:,k-1)   ! kg/m2
             diff(:,k) = MAX( 0.0_r_2, ssnow%wb(:,k) - soil%swilt) &      ! m3/m3
                  * soil%zse(k)*Cdensity_liq
             xxd = xx - diff(:,k)

             WHERE ( xxd .GT. 0.0 )
                ssnow%wb(:,k) = ssnow%wb(:,k) - diff(:,k) / (soil%zse(k)*Cdensity_liq)
                diff(:,k) = xxd
             ELSEWHERE
                ssnow%wb(:,k) = ssnow%wb(:,k) - xx / (soil%zse(k)*Cdensity_liq)
                diff(:,k) = 0.0
             ENDWHERE

          END WHERE

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

END MODULE remove_trans_mod
