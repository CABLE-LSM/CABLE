MODULE remove_trans_mod

USE cbl_ssnow_data_mod

PUBLIC  remove_trans

CONTAINS

SUBROUTINE remove_trans(soil, ssnow, canopy)
    !! Removes transpiration water from soil.
    !! We also attribute the negative canopy transpiration (dew) 
    !! to the wet canopy flux.

    USE cable_common_module, ONLY : cable_user

    TYPE(canopy_type), INTENT(INOUT)         :: canopy
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil
    INTEGER k


    WHERE (canopy%fevc < 0.0_r_2)
      canopy%fevw = canopy%fevw+canopy%fevc
      canopy%fevc = 0.0_r_2
    END WHERE

    DO k = 1,ms
      ssnow%wbliq(:,k) = ssnow%wbliq(:,k) - ssnow%evapfbl(:,k)/               &
                                             (soil%zse_vec(:,k)*Cdensity_liq)
      ssnow%wb(:,k)    = ssnow%wbliq(:,k) + ssnow%wbice(:,k)
    END DO

    IF (cable_user%gw_model) THEN
       ssnow%wb    = ssnow%wbliq + den_rat * ssnow%wbice
       ssnow%wmliq = ssnow%wbliq * soil%zse_vec * Cdensity_liq !mass
       ssnow%wmtot = ssnow%wmliq + ssnow%wmice  !mass

    END IF

  END SUBROUTINE remove_trans


  FUNCTION transp_soil_water(dels, swilt, froot, zse, fevc, wbliq) RESULT(evapfbl)

    !! Calculates the amount of water removed from the soil by transpiration.
    !
    REAL, INTENT(IN)                     :: dels
       !! integration time step (s)
    REAL(r_2), DIMENSION(ms), INTENT(IN) :: swilt 
       !! wilting point (m3/m3)
    REAL(r_2), INTENT(IN)                :: fevc 
       !! transpiration (kg/m2/s)
    REAL, DIMENSION(ms), INTENT(IN)      :: froot 
       !! root fraction (-)
    REAL(r_2), DIMENSION(ms), INTENT(IN) :: zse 
       !! soil depth (m)
    REAL(r_2), DIMENSION(ms), INTENT(IN) :: wbliq 
       !! liquid soil water (m3/m3)
   
    ! Local variables
    REAL(r_2), DIMENSION(ms)   :: evapfbl
    REAL(r_2), DIMENSION(0:ms) :: diff
    REAL(r_2)                  :: xx,xxd
    INTEGER k
      
    xx = 0.; xxd = 0.; diff(:) = 0.
    IF (fevc > 0.0) THEN
      DO k = 1,ms

        ! Calculate the amount (perhaps moisture/ice limited)
        ! that can be removed:
        ! xx: water demand from the transpiration and above soil layers
        ! diff(k-1): excess demand from higher soil layers
        ! diff(k): maximum water amount available for this layer (supply)
        ! xxd: demand minus supply. If the demand is larger (xxd>0), 
        ! evapfbl is limited by the supply and the excess demand is shifted 
        ! to the next layer.
        xx = fevc * dels / CHL * froot(k) + diff(k-1)   ! kg/m2
        diff(k) = MAX( 0.0_r_2, wbliq(k) - 1.1 * swilt(k)) * zse(k)*Cdensity_liq
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
   
   END FUNCTION transp_soil_water

END MODULE remove_trans_mod
