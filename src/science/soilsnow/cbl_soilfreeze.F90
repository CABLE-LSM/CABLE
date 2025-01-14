MODULE soilfreeze_mod

USE cbl_ssnow_data_mod

PUBLIC  soilfreeze

CONTAINS

SUBROUTINE soilfreeze(dels, soil, ssnow,heat_cap_lower_limit)
    USE cable_common_module
IMPLICIT NONE
    REAL, INTENT(IN)                    :: dels ! integration time step (s)
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil
    REAL(r_2), DIMENSION(mp)           :: sicefreeze
    REAL(r_2), DIMENSION(mp)           :: sicemelt
    REAL, DIMENSION(mp)           :: xx
INTEGER :: i,k
REAL :: heat_cap_lower_limit(mp,ms)
REAL :: max_arg1, max_arg2

xx = 0.
DO k = 1, ms  !loop over soil levels
  DO i = 1, mp  !loop over active tiles
   
    IF(ssnow%tgg(i,k) < CTFRZ &
       & .AND. frozen_limit*ssnow%wb(i,k) - ssnow%wbice(i,k) > .001) THEN
      
      ! Limits on water mass that can be converted to ice due frozen_limit and saturation
      sicefreeze(i) = MIN( frozen_limit * ssnow%wb(i,k) - ssnow%wbice(i,k) ,  &
        (soil%ssat(i) - ssnow%wb(i,k)) / MAX((1.0-Cdensity_ice/Cdensity_liq),1.0E-3)  )
      ! Apply limits on watermass to freeze by zero and the energy limit
      sicefreeze(i) = MIN( MAX( 0.0_r_2, sicefreeze(i) ) * soil%zse(k) * Cdensity_ice,             &
        ( CTFRZ - ssnow%tgg(i,k) ) * ssnow%gammzz(i,k) / CHLF )
        
      ssnow%wbice(i,k) = MIN( ssnow%wbice(i,k) + sicefreeze(i) / (soil%zse(k)  &
                         * Cdensity_ice), frozen_limit * ssnow%wb(i,k) )
      ssnow%wb(i,k) = ssnow%wb(i,k) + sicefreeze(i) / (soil%zse(k) * Cdensity_ice) &
                      - sicefreeze(i) / (soil%zse(k) * Cdensity_liq)
      xx(i) = soil%css(i) * soil%rhosoil(i)
      max_arg1 = heat_cap_lower_limit(i,k) 
      max_arg2 = REAL((1.0 - soil%ssat(i)) * soil%css(i) * soil%rhosoil(i) ,r_2) &
          + (ssnow%wb(i,k) - ssnow%wbice(i,k)) * REAL(Ccswat * Cdensity_liq,r_2)  &
          + ssnow%wbice(i,k) * REAL(Ccsice * Cdensity_ice,r_2)
      ssnow%gammzz(i,k) = MAX( max_arg1, max_arg2 ) * REAL( soil%zse(k),r_2 )

      IF (k == 1 .AND. ssnow%isflag(i) == 0)  &
         ssnow%gammzz(i,k) = ssnow%gammzz(i,k) + Ccgsnow * ssnow%snowd(i)
      
      ssnow%tgg(i,k) = ssnow%tgg(i,k) + REAL(sicefreeze(i))                    &
                       * CHLF / REAL(ssnow%gammzz(i,k) )

   ELSEIF( ssnow%tgg(i,k) > CTFRZ .AND. ssnow%wbice(i,k) > 0. ) THEN

      sicemelt(i) = MIN( ssnow%wbice(i,k) * soil%zse(k) * Cdensity_ice,              &
                 ( ssnow%tgg(i,k) - CTFRZ ) * ssnow%gammzz(i,k) / CHLF )

      ssnow%wbice(i,k) = MAX( 0.0_r_2, ssnow%wbice(i,k) - sicemelt(i)          &
                         / (soil%zse(k) * Cdensity_ice) )
      ssnow%wb(i,k) = ssnow%wb(i,k) - sicemelt(i) / (soil%zse(k) * Cdensity_ice) &
            + sicemelt(i) / (soil%zse(k) * Cdensity_liq)
      xx(i) = soil%css(i) * soil%rhosoil(i)
      max_arg1 = heat_cap_lower_limit(i,k) 
      max_arg2 = REAL((1.0 - soil%ssat(i)) * soil%css(i) * soil%rhosoil(i) ,r_2) &
          + (ssnow%wb(i,k) - ssnow%wbice(i,k)) * REAL(Ccswat * Cdensity_liq,r_2)  &
          + ssnow%wbice(i,k) * REAL(Ccsice * Cdensity_ice,r_2)
      
      ssnow%gammzz(i,k) = MAX( max_arg1, max_arg2 ) * REAL( soil%zse(k),r_2 )

      IF (k == 1 .AND. ssnow%isflag(i) == 0) &
         ssnow%gammzz(i,k) = ssnow%gammzz(i,k) + Ccgsnow * ssnow%snowd(i)

      ssnow%tgg(i,k) = ssnow%tgg(i,k) - REAL(sicemelt(i))                     &
                       * CHLF / REAL(ssnow%gammzz(i,k))

   ENDIF

    END DO
END DO

END SUBROUTINE soilfreeze

END MODULE soilfreeze_mod

