MODULE cbl_soil_snow_main_module


IMPLICIT NONE

PRIVATE
   
PUBLIC soil_snow ! must be available outside this module
   


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
USE cbl_ssnow_data_mod
!called subrs
USE hydraulic_redistribution_mod, ONLY: hydraulic_redistribution
USE soilfreeze_mod,               ONLY: soilfreeze
USE remove_trans_mod,             ONLY: remove_trans
USE snowl_adjust_mod,             ONLY: snowl_adjust
USE snowCheck_mod,                ONLY: snowCheck
USE stempv_mod,                   ONLY: stempv
USE surfbv_mod,                   ONLY: surfbv
USE snow_melting_mod,             ONLY: snow_melting
USE snow_accum_mod,               ONLY: snow_accum
USE snowdensity_mod,              ONLY: snowDensity

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
    REAL, DIMENSION(mp) :: xx
    REAL(r_2), DIMENSION(mp) :: xxx
    REAL(r_2), DIMENSION(mp) :: deltat,sinfil1,sinfil2,sinfil3
    REAL                :: zsetot
    INTEGER, SAVE :: ktau =0
REAL :: wbliq(mp,ms)

    ktau = ktau +1
  !this is the value it is initialized with in cable_common anyway 
  max_glacier_snowd = 1100.0 ! for ACCESS1.3 onwards. = 50000.0 for ACCESS1.0

    zsetot = SUM(soil%zse)
    ssnow%tggav = 0.
    DO k = 1, ms
      ssnow%tggav = ssnow%tggav  + soil%zse(k)*ssnow%tgg(:,k)/zsetot
      soil%heat_cap_lower_limit(:,k) = MAX( 0.01, soil%css(:) * soil%rhosoil(:) )
    END DO

  IF( cable_runtime%offline .or. cable_runtime%mk3l ) THEN !in um_init for UM
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


    wbliq = ssnow%wb - ssnow%wbice

  !%cable_runtime_coupled special initalizations in um_init NA for ESM1.5

   xx=soil%css * soil%rhosoil
   IF (ktau <= 1)                                                              &
     ssnow%gammzz(:,1) = MAX( (1.0 - soil%ssat) * soil%css * soil%rhosoil      &
            & + (ssnow%wb(:,1) - ssnow%wbice(:,1) ) * Ccswat * Cdensity_liq           &
            & + ssnow%wbice(:,1) * Ccsice * Cdensity_liq * .9, xx ) * soil%zse(1) +   &
            & (1. - ssnow%isflag) * Ccgsnow * ssnow%snowd



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

   CALL stempv(dels, canopy, ssnow, soil, soil%heat_cap_lower_limit )

    ssnow%tss =  (1-ssnow%isflag)*ssnow%tgg(:,1) + ssnow%isflag*ssnow%tggsn(:,1)

    CALL snow_melting (dels, snowmlt, ssnow, soil )

    ! Add new snow melt to global snow melt variable:
    ssnow%smelt = ssnow%smelt + snowmlt

    CALL remove_trans(dels, soil, ssnow, canopy, veg)

   CALL  soilfreeze(dels, soil, ssnow, soil%heat_cap_lower_limit)


    totwet = canopy%precis + ssnow%smelt

    ! total available liquid including puddle
    weting = totwet + MAX(0._r_2,ssnow%pudsto - canopy%fesp/CHL*dels)
    xxx=soil%ssat - ssnow%wb(:,1)

    sinfil1 = MIN( 0.95*xxx*soil%zse(1)*Cdensity_liq, weting) !soil capacity
    xxx=soil%ssat - ssnow%wb(:,2)
    sinfil2 = MIN( 0.95*xxx*soil%zse(2)*Cdensity_liq, weting - REAL(sinfil1)) !soil capacity
    xxx=soil%ssat - ssnow%wb(:,3)
    sinfil3 = MIN( 0.95*xxx*soil%zse(3)*Cdensity_liq,weting-REAL(sinfil1)-REAL(sinfil2))

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
  canopy%fhs_cor = ssnow%dtmlt(:,1)*ssnow%dfh_dtg
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

    wbliq = ssnow%wb - ssnow%wbice

    ssnow%wbtot = 0.0
    DO k = 1, ms
       ssnow%wbtot = ssnow%wbtot + REAL(ssnow%wb(:,k)*1000.0*soil%zse(k),r_2)
    END DO


RETURN
END SUBROUTINE soil_snow

END MODULE cbl_soil_snow_main_module
