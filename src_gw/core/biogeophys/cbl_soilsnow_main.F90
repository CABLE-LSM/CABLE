MODULE cbl_soil_snow_main_module

  USE cable_def_types_mod, ONLY : soil_snow_type, soil_parameter_type,        &
       veg_parameter_type, canopy_type, met_type,        &
       balances_type, r_2, ms, mp

  USE cable_data_module, ONLY : issnow_type, point2constants

  USE cable_common_module, ONLY: cable_user,snow_ccnsw,snmin,&
       max_ssdn,max_sconds,frozen_limit,&
       max_glacier_snowd

  IMPLICIT NONE

  PRIVATE

  TYPE ( issnow_type ), SAVE :: C

  PUBLIC soil_snow 

CONTAINS


  ! Inputs:
  !        dt_in - time step in sec
  !        ktau_in - time step no.
  !        ga      - ground heat flux W/m^2
  !        dgdtg   -
  !        condxpr - total precip reaching the ground (liquid and solid)
  !        scondxpr - precip (solid only)
  !        fev   - transpiration (W/m2)
  !        fes   - soil evaporation (W/m2)
  !        isoil - soil type
  !        ivegt - vegetation type
  ! Output
  !        ssnow
  SUBROUTINE soil_snow(dels, soil, ssnow, canopy, met, bal, veg)
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

    CALL point2constants( C )

    ktau = ktau +1

    !jhan - make switchable
    ! appropriate for ACCESS1.0
    !max_glacier_snowd = 50000.0
    ! appropriate for ACCESS1.3
    !max_glacier_snowd = 1100.0

    zsetot = SUM(soil%zse)
    ssnow%tggav = 0.
    DO k = 1, ms
       ssnow%tggav = ssnow%tggav  + soil%zse(k)*ssnow%tgg(:,k)/zsetot
    END DO


    IF( cable_runtime%offline .OR. cable_runtime%mk3l ) THEN
       ssnow%t_snwlr = 0.05
    ENDIF

    ssnow%fwtop1 = 0.0
    ssnow%fwtop2 = 0.0
    ssnow%fwtop3 = 0.0
    ssnow%runoff = 0.0 ! initialise total runoff
    ssnow%rnof1 = 0.0 ! initialise surface runoff
    ssnow%rnof2 = 0.0 ! initialise deep drainage
    ssnow%smelt = 0.0 ! initialise snowmelt
    ssnow%dtmlt = 0.0
    ssnow%osnowd = ssnow%snowd

    IF (cable_user%soil_thermal_fix) THEN
       soil%heat_cap_lower_limit(:,:) = 0.01  !never allow /0
    ELSE
       DO k=1,ms
          soil%heat_cap_lower_limit(:,k) = soil%css(:) * soil%rhosoil(:)
       END DO
    END IF

!$ inserted block below as per MMY code -- rk4417

   IF( .NOT.cable_user%cable_runtime_coupled ) THEN

      IF( ktau_gl <= 1 ) THEN

         IF (cable_runtime%um) canopy%dgdtg = 0.0 ! RML added um condition
                                                  ! after discussion with BP
         ! N.B. snmin should exceed sum of layer depths, i.e. .11 m
         ssnow%wbtot = 0.0
         DO k = 1, ms
            ssnow%wb(:,k)  = MIN( soil%ssat, MAX( real(ssnow%wb(:,k)), soil%swilt ) )
         END DO

         ssnow%wb(:,ms-2)  = MIN( soil%ssat, MAX( real(ssnow%wb(:,ms-2)),           &
                             0.5 * ( soil%sfc + soil%swilt ) ) )
         ssnow%wb(:,ms-1)  = MIN( soil%ssat, MAX( real(ssnow%wb(:,ms-1)),           &
                             0.8 * soil%sfc ) )
         ssnow%wb(:,ms)    = MIN( soil%ssat, MAX( real(ssnow%wb(:,ms)), soil%sfc ) )

         DO k = 1, ms

            WHERE( ssnow%tgg(:,k) <= C%TFRZ .AND. ssnow%wbice(:,k) <= 0.01 )   &
               ssnow%wbice(:,k) = 0.5 * ssnow%wb(:,k)

            WHERE( ssnow%tgg(:,k) < C%TFRZ)                                    &
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

         xx=real(soil%heat_cap_lower_limit(:,1))

         ssnow%gammzz(:,1) = MAX( (1.0 - soil%ssat) * soil%css * soil%rhosoil &
              & + (ssnow%wb(:,1) - ssnow%wbice(:,1) ) * C%cswat * C%density_liq &
              & + ssnow%wbice(:,1) * C%csice * C%density_liq * .9, xx ) * soil%zse(1)
      END IF
   ENDIF  ! if(.NOT.cable_runtime_coupled)


   IF (ktau <= 1)       THEN

       xx=soil%heat_cap_lower_limit(:,1)
     ssnow%gammzz(:,1) = MAX( (1.0 - soil%ssat) * soil%css * soil%rhosoil      &
            & + (ssnow%wb(:,1) - ssnow%wbice(:,1) ) * C%cswat * C%density_liq           &
            & + ssnow%wbice(:,1) * C%csice * C%density_liq * .9, xx ) * soil%zse(1) +   &
            & (1. - ssnow%isflag) * C%cgsnow * ssnow%snowd

   END IF

!$ -------------  end of block -------------- rk4417 -----------
    
    ssnow%wbliq = ssnow%wb - ssnow%wbice

    DO k = 1, ms ! for stempv

       ! Set liquid soil water fraction (fraction of saturation value):
       ssnow%wblf(:,k) = MAX( 0.01_r_2, (ssnow%wb(:,k) - ssnow%wbice(:,k)) )    &
            & / REAL(soil%ssat,r_2)

       ! Set ice soil water fraction (fraction of saturation value):
       ssnow%wbfice(:,k) = REAL(ssnow%wbice(:,k)) / soil%ssat
    END DO

    CALL snowcheck (dels, ssnow, soil, met )

    CALL snowdensity (dels, ssnow, soil)

    CALL snow_accum (dels, canopy, met, ssnow, soil )

    CALL snow_melting (dels, snowmlt, ssnow, soil )

    ! Add snow melt to global snow melt variable:
    ssnow%smelt = snowmlt

    ! Adjust levels in the snowpack due to snow accumulation/melting,
    ! snow aging etc...
    CALL snowl_adjust(dels, ssnow, canopy )

    CALL stempv(dels, canopy, ssnow, soil)

    ssnow%tss =  (1-ssnow%isflag)*ssnow%tgg(:,1) + ssnow%isflag*ssnow%tggsn(:,1)

    CALL snow_melting (dels, snowmlt, ssnow, soil )

    ! Add new snow melt to global snow melt variable:
    ssnow%smelt = ssnow%smelt + snowmlt

    CALL remove_trans(dels, soil, ssnow, canopy, veg)

    CALL  soilfreeze(dels, soil, ssnow)

    totwet = canopy%precis + ssnow%smelt

    ! total available liquid including puddle
    weting = totwet + MAX(0._r_2,ssnow%pudsto - canopy%fesp/C%HL*dels)
    xxx=soil%ssat - ssnow%wb(:,1)

    sinfil1 = MIN( 0.95*xxx*soil%zse(1)*C%density_liq, weting) !soil capacity
    xxx=soil%ssat - ssnow%wb(:,2)
    sinfil2 = MIN( 0.95*xxx*soil%zse(2)*C%density_liq, weting - REAL(sinfil1)) !soil capacity
    xxx=soil%ssat - ssnow%wb(:,3)
    sinfil3 = MIN( 0.95*xxx*soil%zse(3)*C%density_liq,weting-REAL(sinfil1)-REAL(sinfil2))

    ! net water flux to the soil
    ssnow%fwtop1 = sinfil1 / dels - canopy%segg
    ssnow%fwtop2 = sinfil2 / dels
    ssnow%fwtop3 = sinfil3 / dels

    ! Puddle for the next time step
    ssnow%pudsto = MAX( 0._r_2, weting - sinfil1 - sinfil2 - sinfil3 )
    ssnow%rnof1 = MAX(0.,ssnow%pudsto - ssnow%pudsmx)
    ssnow%pudsto = ssnow%pudsto - ssnow%rnof1

    CALL surfbv(dels, met, ssnow, soil, veg, canopy )


    ! correction required for energy balance in online simulations
    IF( cable_runtime%um ) THEN

       !cls package - rewritten for flexibility
       canopy%fhs_cor = ssnow%dtmlt(:,1)*ssnow%dfh_dtg
       !canopy%fes_cor = ssnow%dtmlt(:,1)*(ssnow%dfe_ddq * ssnow%ddq_dtg)
       canopy%fes_cor = ssnow%dtmlt(:,1)*ssnow%dfe_dtg

       canopy%fhs = canopy%fhs+canopy%fhs_cor
       canopy%fes = canopy%fes+canopy%fes_cor

       !REV_CORR associated changes to other energy balance terms
       !NB canopy%fns changed not rad%flws as the correction term needs to
       !pass through the canopy in entirety, not be partially absorbed
       IF (cable_user%L_REV_CORR) THEN
          canopy%fns_cor = ssnow%dtmlt(:,1)*ssnow%dfn_dtg
          canopy%ga_cor = ssnow%dtmlt(:,1)*canopy%dgdtg

          canopy%fns = canopy%fns + canopy%fns_cor
          canopy%ga = canopy%ga + canopy%ga_cor

          canopy%fess = canopy%fess + canopy%fes_cor
       ENDIF
    ENDIF

    ! redistrb (set in cable.nml) by default==.FALSE.
    IF( redistrb )                                                              &
         CALL hydraulic_redistribution( dels, soil, ssnow, canopy, veg, met )

    ssnow%smelt = ssnow%smelt/dels

    ! Set weighted soil/snow surface temperature
    ssnow%tss=(1-ssnow%isflag)*ssnow%tgg(:,1) + ssnow%isflag*ssnow%tggsn(:,1)

    ssnow%wbliq = ssnow%wb - ssnow%wbice

    ssnow%wbtot = 0.0
    DO k = 1, ms
       ssnow%wbtot = ssnow%wbtot + REAL(ssnow%wb(:,k)*1000.0*soil%zse(k),r_2)
    END DO

  END SUBROUTINE soil_snow

END MODULE cbl_soil_snow_main_module
