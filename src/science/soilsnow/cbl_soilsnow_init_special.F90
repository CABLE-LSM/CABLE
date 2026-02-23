MODULE cbl_soil_snow_init_special_module

USE cbl_ssnow_data_mod

  IMPLICIT NONE

  PRIVATE

  PUBLIC spec_init_soil_snow

CONTAINS

SUBROUTINE spec_init_soil_snow(dels, soil, ssnow, canopy, met, bal, veg,heat_cap_lower_limit)
USE cable_common_module
REAL, INTENT(IN)                    :: dels ! integration time step (s)
TYPE(soil_parameter_type), INTENT(INOUT) :: soil
TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow
TYPE(canopy_type), INTENT(INOUT)         :: canopy
TYPE(veg_parameter_type), INTENT(INOUT)  :: veg
TYPE(met_type), INTENT(INOUT)            :: met ! all met forcing
TYPE (balances_type), INTENT(INOUT)      :: bal
INTEGER             :: k
REAL, DIMENSION(mp) :: snowmlt
REAL, DIMENSION(mp) :: totwet
REAL, DIMENSION(mp) :: weting
REAL, DIMENSION(mp) :: xx, tgg_old, tggsn_old
REAL(r_2), DIMENSION(mp) :: xxx,deltat,sinfil1,sinfil2,sinfil3
REAL                :: zsetot
INTEGER, SAVE :: ktau =0
REAL :: heat_cap_lower_limit(mp,ms)

ktau = ktau +1

IF( (cable_user%soilsnow_init_spec ) ) THEN

   IF( ktau_gl <= 1 ) THEN
      IF (cable_runtime%um) canopy%dgdtg = 0.0 ! RML added um condition
      ! after discussion with BP
      ! N.B. snmin should exceed sum of layer depths, i.e. .11 m
      ssnow%wbtot = 0.0
      DO k = 1, ms
         ssnow%wb(:,k)  = MIN( soil%ssat, MAX( REAL(ssnow%wb(:,k)), soil%swilt ) )
      END DO
      ssnow%wb(:,ms-2)  = MIN( soil%ssat, MAX( REAL(ssnow%wb(:,ms-2)),           &
           0.5 * ( soil%sfc + soil%swilt ) ) )
      ssnow%wb(:,ms-1)  = MIN( soil%ssat, MAX( REAL(ssnow%wb(:,ms-1)),           &
           0.8 * soil%sfc ) )
      ssnow%wb(:,ms)    = MIN( soil%ssat, MAX( REAL(ssnow%wb(:,ms)), soil%sfc ) )
      DO k = 1, ms
         WHERE( ssnow%tgg(:,k) <= CTFRZ .AND. ssnow%wbice(:,k) <= 0.01 )   &
              ssnow%wbice(:,k) = 0.5 * ssnow%wb(:,k)
         WHERE( ssnow%tgg(:,k) < CTFRZ)                                    &
              ssnow%wbice(:,k) = frozen_limit * ssnow%wb(:,k)
      END DO
      WHERE (soil%isoilm == 9)
         ! permanent ice: fix hard-wired number in next version
         ssnow%snowd = max_glacier_snowd
         ssnow%osnowd = max_glacier_snowd
         ssnow%tgg(:,1) = ssnow%tgg(:,1) - 1.0
         ssnow%wb(:,1) = 0.95 * soil%ssat
         ssnow%wb(:,2) = 0.95 * soil%ssat
         ssnow%wb(:,3) = 0.95 * soil%ssat
         ssnow%wb(:,4) = 0.95 * soil%ssat
         ssnow%wb(:,5) = 0.95 * soil%ssat
         ssnow%wb(:,6) = 0.95 * soil%ssat
         ssnow%wbice(:,1) = frozen_limit * ssnow%wb(:,1)
         ssnow%wbice(:,2) = frozen_limit * ssnow%wb(:,2)
         ssnow%wbice(:,3) = frozen_limit * ssnow%wb(:,3)
         ssnow%wbice(:,4) = frozen_limit * ssnow%wb(:,4)
         ssnow%wbice(:,5) = frozen_limit * ssnow%wb(:,5)
         ssnow%wbice(:,6) = frozen_limit * ssnow%wb(:,6)
      ENDWHERE
      xx=REAL(heat_cap_lower_limit(:,1))
      ssnow%gammzz(:,1) = MAX( (1.0 - soil%ssat) * soil%css * soil%rhosoil &
           & + (ssnow%wb(:,1) - ssnow%wbice(:,1) ) * Ccswat * Cdensity_liq &
           & + ssnow%wbice(:,1) * Ccsice * Cdensity_ice, xx ) * soil%zse(1)
   END IF
ENDIF  ! if(.NOT.soilsnow_init_spec )

IF (ktau <= 1)       THEN
  xx=heat_cap_lower_limit(:,1)
  ssnow%gammzz(:,1) = MAX( (1.0 - soil%ssat) * soil%css * soil%rhosoil      &
        & + (ssnow%wb(:,1) - ssnow%wbice(:,1) ) * Ccswat * Cdensity_liq           &
        & + ssnow%wbice(:,1) * Ccsice * Cdensity_ice, xx ) * soil%zse(1) +   &
        & (1. - ssnow%isflag) * Ccgsnow * ssnow%snowd
END IF

END SUBROUTINE spec_init_soil_snow

END MODULE cbl_soil_snow_init_special_module
