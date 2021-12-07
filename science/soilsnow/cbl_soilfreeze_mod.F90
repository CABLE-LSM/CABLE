MODULE soilfreeze_mod

USE cbl_ssnow_data_mod

PUBLIC  soilfreeze

CONTAINS

SUBROUTINE soilfreeze(dels, soil, ssnow,heat_cap_lower_limit)
    USE cable_common_module
    REAL, INTENT(IN)                    :: dels ! integration time step (s)
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil
    REAL(r_2), DIMENSION(mp)           :: sicefreeze
    REAL(r_2), DIMENSION(mp)           :: sicemelt
    REAL, DIMENSION(mp)           :: xx
    INTEGER k
REAL :: heat_cap_lower_limit(mp,ms)

    xx = 0.
    DO k = 1, ms

       WHERE (ssnow%tgg(:,k) < CTFRZ &
            & .AND. frozen_limit * ssnow%wb(:,k) - ssnow%wbice(:,k) > .001)

          sicefreeze = MIN( MAX( 0.0_r_2, ( frozen_limit * ssnow%wb(:,k) -      &
               ssnow%wbice(:,k) ) ) * soil%zse(k) * 1000.0,             &
               ( CTFRZ - ssnow%tgg(:,k) ) * ssnow%gammzz(:,k) / CHLF )
          ssnow%wbice(:,k) = MIN( ssnow%wbice(:,k) + sicefreeze / (soil%zse(k)  &
               * 1000.0), frozen_limit * ssnow%wb(:,k) )
          xx = soil%css * soil%rhosoil
          ssnow%gammzz(:,k) = MAX((heat_cap_lower_limit(:,k)),           &
               REAL((1.0 - soil%ssat) * soil%css * soil%rhosoil ,r_2)            &
               + (ssnow%wb(:,k) - ssnow%wbice(:,k)) * REAL(Ccswat * Cdensity_liq,r_2)   &
               + ssnow%wbice(:,k) * REAL(Ccsice * Cdensity_liq * 0.9,r_2))* &
               REAL( soil%zse(k),r_2 )

          WHERE (k == 1 .AND. ssnow%isflag == 0)
             ssnow%gammzz(:,k) = ssnow%gammzz(:,k) + Ccgsnow * ssnow%snowd
          END WHERE
          ssnow%tgg(:,k) = ssnow%tgg(:,k) + REAL(sicefreeze)                    &
               * CHLF / REAL(ssnow%gammzz(:,k) )

       ELSEWHERE( ssnow%tgg(:,k) > CTFRZ .AND. ssnow%wbice(:,k) > 0. )

          sicemelt = MIN( ssnow%wbice(:,k) * soil%zse(k) * 1000.0,              &
               ( ssnow%tgg(:,k) - CTFRZ ) * ssnow%gammzz(:,k) / CHLF )

          ssnow%wbice(:,k) = MAX( 0.0_r_2, ssnow%wbice(:,k) - sicemelt          &
               / (soil%zse(k) * 1000.0) )
          xx = soil%css * soil%rhosoil
          ssnow%gammzz(:,k) = MAX((heat_cap_lower_limit(:,k)),       &
               REAL((1.0-soil%ssat) * soil%css * soil%rhosoil,r_2)             &
               + (ssnow%wb(:,k) - ssnow%wbice(:,k)) * REAL(Ccswat*Cdensity_liq,r_2)   &
               + ssnow%wbice(:,k) * REAL(Ccsice * Cdensity_liq * 0.9,r_2))            &
               * REAL(soil%zse(k),r_2)
          WHERE (k == 1 .AND. ssnow%isflag == 0)
             ssnow%gammzz(:,k) = ssnow%gammzz(:,k) + Ccgsnow * ssnow%snowd
          END WHERE
          ssnow%tgg(:,k) = ssnow%tgg(:,k) - REAL(sicemelt)                     &
               * CHLF / REAL(ssnow%gammzz(:,k))

       END WHERE

    END DO

  END SUBROUTINE soilfreeze

END MODULE soilfreeze_mod

