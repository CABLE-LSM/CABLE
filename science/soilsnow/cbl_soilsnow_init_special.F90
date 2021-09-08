MODULE cbl_soil_snow_init_special_module

  USE cable_def_types_mod, ONLY : soil_snow_type, soil_parameter_type,        &
       veg_parameter_type, canopy_type, met_type,        &
       balances_type, r_2, ms, mp

USE cable_phys_constants_mod, ONLY : CTFRZ => TFRZ
USE cable_phys_constants_mod, ONLY : CHL => HL
USE cable_phys_constants_mod, ONLY : Cdensity_liq => density_liq
USE cable_phys_constants_mod, ONLY : Ccgsnow => cgsnow
USE cable_phys_constants_mod, ONLY : Ccswat => cswat
USE cable_phys_constants_mod, ONLY : Ccsice => csice

  USE cable_common_module, ONLY: cable_user,snow_ccnsw,snmin,&
       max_ssdn,max_sconds,frozen_limit,&
       max_glacier_snowd

  IMPLICIT NONE

  PRIVATE

  PUBLIC spec_init_soil_snow
  PUBLIC spec_init_snowcheck

CONTAINS

SUBROUTINE spec_init_soil_snow(dels, soil, ssnow, canopy, met, bal, veg,heat_cap_lower_limit)
USE cable_common_module
!all subrs-implement ONLY:
USE cbl_soil_snow_subrs_module
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

IF( .NOT.cable_user%cable_runtime_coupled ) THEN

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
         ssnow%wbice(:,1) = 0.90 * ssnow%wb(:,1)
         ssnow%wbice(:,2) = 0.90 * ssnow%wb(:,2)
         ssnow%wbice(:,3) = 0.90 * ssnow%wb(:,3)
         ssnow%wbice(:,4) = 0.90 * ssnow%wb(:,4)
         ssnow%wbice(:,5) = 0.90 * ssnow%wb(:,5)
         ssnow%wbice(:,6) = 0.90 * ssnow%wb(:,6)
      ENDWHERE
      xx=REAL(heat_cap_lower_limit(:,1))
      ssnow%gammzz(:,1) = MAX( (1.0 - soil%ssat) * soil%css * soil%rhosoil &
           & + (ssnow%wb(:,1) - ssnow%wbice(:,1) ) * Ccswat * Cdensity_liq &
           & + ssnow%wbice(:,1) * Ccsice * Cdensity_liq * .9, xx ) * soil%zse(1)
   END IF
ENDIF  ! if(.NOT.cable_runtime_coupled)

IF (ktau <= 1)       THEN
  xx=heat_cap_lower_limit(:,1)
  ssnow%gammzz(:,1) = MAX( (1.0 - soil%ssat) * soil%css * soil%rhosoil      &
        & + (ssnow%wb(:,1) - ssnow%wbice(:,1) ) * Ccswat * Cdensity_liq           &
        & + ssnow%wbice(:,1) * Ccsice * Cdensity_liq * .9, xx ) * soil%zse(1) +   &
        & (1. - ssnow%isflag) * Ccgsnow * ssnow%snowd
END IF

END SUBROUTINE spec_init_soil_snow

  SUBROUTINE spec_init_snowcheck(dels, ssnow, soil, met )

    USE cable_common_module

    REAL, INTENT(IN) :: dels ! integration time step (s)

    TYPE(soil_snow_type), INTENT(INOUT) :: ssnow
    TYPE(met_type),       INTENT(INOUT) :: met ! all met forcing

    TYPE(soil_parameter_type), INTENT(INOUT) :: soil  ! soil parameters

    INTEGER :: k,j


    DO j=1,mp

       IF( ssnow%snowd(j) <= 0.0 ) THEN
          !H!ssnow%isflag(j) = 0
          !H!ssnow%ssdn(j,:) = 120.0
          !H!ssnow%ssdnn(j) = 120.0
          !H!ssnow%tggsn(j,:) = CTFRZ
          !H!ssnow%sdepth(j,1) = ssnow%snowd(j) / ssnow%ssdn(j,1)
          !H!ssnow%sdepth(j,2) = 0.
          !H!ssnow%sdepth(j,3) = 0.
          !H!ssnow%smass(j,1) = ssnow%snowd(j)
          !H!ssnow%smass(j,2) = 0.0     ! EK to fix -ve sdepth 21Dec2007
          !H!ssnow%smass(j,3) = 0.0     ! EK to fix -ve sdepth 21Dec2007
       ELSEIF( ssnow%snowd(j) < snmin * ssnow%ssdnn(j) ) THEN
          !H!IF( ssnow%isflag(j) == 1 ) THEN
          !H!   ssnow%ssdn(j,1) = ssnow%ssdnn(j)
          !H!   ssnow%tgg(j,1) = ssnow%tggsn(j,1)
          !H!ENDIF
          !H!ssnow%isflag(j) = 0
          !H!ssnow%ssdnn(j) = MIN( 400.0, MAX( 120.0, ssnow%ssdn(j,1) ) )
          !H!ssnow%tggsn(j,:) = MIN( CTFRZ,ssnow%tgg(j,1) )
          !H!ssnow%sdepth(j,1) = ssnow%snowd(j) / ssnow%ssdn(j,1)
          !H!ssnow%sdepth(j,2) = 0.0
          !H!ssnow%sdepth(j,3) = 0.0
          !H!ssnow%smass(j,1) = ssnow%snowd(j)
          !H!ssnow%smass(j,2) = 0.0
          !H!ssnow%smass(j,3) = 0.0
          !H!ssnow%ssdn(j,:) = ssnow%ssdnn(j)

          IF( .NOT.cable_user%CABLE_RUNTIME_COUPLED ) THEN
             IF( soil%isoilm(j) == 9 .AND. ktau_gl <= 2 )                       &
                                ! permanent ice: fixed hard-wired number in next version
                  ssnow%ssdnn(j) = 700.0
          ENDIF

       ELSE ! in loop: IF( ssnow%snowd(j) <= 0.0 ) THEN
          ! sufficient snow now for 3 layer snowpack

          IF( ssnow%isflag(j) == 0 ) THEN
             !H!ssnow%tggsn(j,:) = MIN( CTFRZ, ssnow%tgg(j,1) )
             !H!ssnow%ssdn(j,2) = ssnow%ssdn(j,1)
             !H!ssnow%ssdn(j,3) = ssnow%ssdn(j,1)
             IF( .NOT. cable_user%cable_runtime_coupled) THEN
                IF( soil%isoilm(j) == 9 .AND. ktau_gl <= 2 ) THEN
                   ! permanent ice: fix hard-wired number in next version
                   ssnow%ssdn(j,1)  = 450.0
                   ssnow%ssdn(j,2)  = 580.0
                   ssnow%ssdn(j,3)  = 600.0
                ENDIF
             ENDIF
             !H!ssnow%sdepth(j,1) = ssnow%t_snwlr(j)
             !H!ssnow%smass(j,1)  =  ssnow%t_snwlr(j) * ssnow%ssdn(j,1)
             !H!ssnow%smass(j,2)  = ( ssnow%snowd(j) - ssnow%smass(j,1) ) * 0.4
             !H!ssnow%smass(j,3)  = ( ssnow%snowd(j) - ssnow%smass(j,1) ) * 0.6
             !H!ssnow%sdepth(j,2) = ssnow%smass(j,2) / ssnow%ssdn(j,2)
             !H!ssnow%sdepth(j,3) = ssnow%smass(j,3) / ssnow%ssdn(j,3)
             !H!ssnow%ssdnn(j) = ( ssnow%ssdn(j,1) * ssnow%smass(j,1) +            &
             !H!     ssnow%ssdn(j,2) * ssnow%smass(j,2) +             &
             !H!     ssnow%ssdn(j,3) * ssnow%smass(j,3) )             &
             !H!     / ssnow%snowd(j)
          ENDIF
          !H!ssnow%isflag(j) = 1
       ENDIF ! END: IF( ssnow%snowd(j) <= 0.0 ) THEN
    ENDDO ! END: DO j=1,mp

  END SUBROUTINE spec_init_snowcheck
END MODULE cbl_soil_snow_init_special_module
