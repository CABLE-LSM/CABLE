!==============================================================================
! This source code is part of the
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
!
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) is located
! in each directory containing CABLE code.
!
! ==============================================================================
! Purpose: Calculates surface exchange fluxes through the solution of surface
!          energy balance and its interaction with plant physiology. Specific
!        representation of the transport of scalars within a canopy is included.
!
! Called from: cbm
!
! Contact: Yingping.Wang@csiro.au and Eva.Kowalczyk@csiro.au
!
! History: Revision of canopy temperature calculation (relative to v1.4b)
!          Reorganisation of code (dryLeaf, wetLeaf, photosynthesis subroutines
!          taken out of define_canopy)
!        : Martin De Kauwe and Jatin Kala added new switch to compute stomatal
!          conductance based on: Medlyn BE et al (2011) Global Change Biology 17:
!          2134-2144. The variables xleuning, xleuningz are no longer used, but
!          replaced with gs_coeff,gs_coeffz. If GS_SWITCH is set to "leuning",
!          gs_coeff=xleuning and gs_coeffz=xleuningz, but based on the new model
!          if set to "medlyn". Search for "Ticket #56"
!        : Vanessa Haverd added new fwsoil_switch  for response of stomatal conductance
!          to soil moisture, which resolves decoupling of transpiration and
!          photosynthesis at low soil moisture. Search for "Haverd2013".
!        : Vanessa Hverd added new logical switch cable_user%litter. When 'true',
!          leaf litter suppresses soil evaporation
!        : Vanessa Haverd added new switch to enable SLI alternative to default soil
!          module. When cable_user%soil_struc=='sli', then SLI is used to compute
!          coupled transfers of heat and water in the soil and snow and at the surface
!          and an in-canopy stability correction is applied.
!        : See http://www.geosci-model-dev.net/9/3111/2016 for full documentation
!          of last 3 changes.
!
!        : Feb 2017 - various changes to facilitate the cls (Ticket 137) and
!          REV_CORR packages
!
!        : August 2017 - GW and Or parameterisations includes and merged with
!          cls and REV_CORR
! ==============================================================================

MODULE cable_canopy_module

  USE cable_data_module, ONLY : icanopy_type, point2constants

  IMPLICIT NONE

  PUBLIC define_canopy
  PRIVATE

  TYPE( icanopy_type ) :: C


CONTAINS


  SUBROUTINE define_canopy(bal,rad,rough,air,met,dels,ssnow,soil,veg, canopy,climate)
    USE cable_def_types_mod
    USE cable_radiation_module
    USE cable_air_module
    USE cable_common_module
    USE cable_roughness_module
    USE cable_psm, ONLY: or_soil_evap_resistance,rtevap_max,&
         rt_Dff,update_or_soil_resis
    USE cable_gw_hydro_module, ONLY : pore_space_relative_humidity
    USE sli_main_mod, ONLY : sli_main


    TYPE (balances_type), INTENT(INOUT)  :: bal
    TYPE (radiation_type), INTENT(INOUT) :: rad
    TYPE (roughness_type), INTENT(INOUT) :: rough
    TYPE (air_type), INTENT(INOUT)       :: air
    TYPE (met_type), INTENT(INOUT)       :: met
    TYPE (soil_snow_type), INTENT(INOUT) :: ssnow
    TYPE (canopy_type), INTENT(INOUT)    :: canopy
    TYPE (climate_type), INTENT(IN)    :: climate

    TYPE (soil_parameter_type), INTENT(INOUT)   :: soil
    TYPE (veg_parameter_type), INTENT(INOUT)    :: veg
    !INTEGER, INTENT(IN) :: wlogn

    REAL, INTENT(IN)               :: dels ! integration time setp (s)
    INTEGER  ::                                                                 &
         iter,  & ! iteration #
         iterplus !

    REAL, DIMENSION(mp) ::                                                      &
         rt0,           & ! turbulent resistance
         ortsoil,       & ! turb. resist. prev t-step
         rt1usc,        & ! eq. 3.53, SCAM manual, 1997
         tstar,         & !
         zscrn,         & !
         qstar,         & !
         rsts,          & !
         qsurf,         & !
         qtgnet,        & !
         tss4,          & ! soil/snow temperature**4
         qstvair,       & ! sat spec humidity at leaf temperature
         xx,            & ! delta-type func 4 sparse canopy limit, p20 SCAM manual
         r_sc,          & !
         zscl,          & !
         pwet,          & !
         dq,            & ! sat sp
         dq_unsat,      & ! spec hum diff including rh at srf
         xx1,           & !
         sum_rad_rniso, & !
         sum_rad_gradis,& !
         rttsoil,       & ! REV_CORR working variable for sensitivity terms
         rhlitt,        & ! REV_CORR working variables for litter resistances
         relitt,        & !
         alpm1,         & ! REV_CORR working variables for Or scheme
         beta2,         & ! beta_div_alpm1 = beta2/alpm1 (goes to zero without
         beta_div_alpm    ! division when no canopy)

    ! temporary buffers to simplify equations
    REAL, DIMENSION(mp) ::                                                      &
         ftemp,z_eff,psim_arg, psim_1, psim_2, rlower_limit,                      &
         term1, term2, term3, term5
    ! arguments for potential_evap (sli)
    REAL(r_2), DIMENSION(mp) ::  Rn, rbh, rbw, Ta, rha,Ts, &
         kth, dz,lambdav, &
         Tsoil, Epot, Hpot, Gpot, &
         dEdrha, dEdTa, dEdTsoil, dGdTa, dGdTsoil
    REAL, DIMENSION(mp) :: qsat

    REAL, DIMENSION(:), POINTER ::                                              &
         cansat,        & ! max canopy intercept. (mm)
         dsx,           & ! leaf surface vpd
         fwsoil,        & ! soil water modifier of stom. cond
         tlfx,          & ! leaf temp prev. iter (K)
         tlfy             ! leaf temp (K)

    REAL(r_2), DIMENSION(mp) ::                                                 &
         gbvtop                   ! bnd layer cond. top leaf

    REAL(r_2), DIMENSION(:), POINTER ::                                         &
         ecy,           & ! lat heat fl dry big leaf
         hcy,           & ! veg. sens heat
         rny,           & ! net rad
         ghwet             ! cond for heat for a wet canopy

    REAL(r_2), DIMENSION(:,:), POINTER ::                                       &
         gbhu,          & ! forcedConvectionBndryLayerCond
         gbhf,          & ! freeConvectionBndryLayerCond
         csx              ! leaf surface CO2 concentration

    REAL  :: rt_min
    REAL, DIMENSION(mp)       :: zstar, rL, phist, csw, psihat,rt0bus

    INTEGER :: j

    INTEGER, SAVE :: call_number =0

    ! END header

    call_number = call_number + 1

    ! assign local ptrs to constants defined in cable_data_module
    CALL point2constants(C)

    IF( .NOT. cable_runtime%um)                                                 &
         canopy%cansto =  canopy%oldcansto

    ALLOCATE( cansat(mp), gbhu(mp,mf))
    ALLOCATE( dsx(mp), fwsoil(mp), tlfx(mp), tlfy(mp) )
    ALLOCATE( ecy(mp), hcy(mp), rny(mp))
    ALLOCATE( gbhf(mp,mf), csx(mp,mf))
    ALLOCATE( ghwet(mp))

    ! BATS-type canopy saturation proportional to LAI:
    cansat = veg%canst1 * canopy%vlaiw

    !---compute surface wetness factor, update cansto, through
    CALL surf_wetness_fact( cansat, canopy, ssnow,veg,met, soil, dels )

    canopy%fevw_pot = 0.0
    canopy%gswx = 1e-3     ! default stomatal conuctance
    gbhf = 1e-3     ! default free convection boundary layer conductance
    gbhu = 1e-3     ! default forced convection boundary layer conductance
    ssnow%evapfbl = 0.0
    ssnow%rex = 0.0
    ! Initialise in-canopy temperatures and humidity:
    csx = SPREAD(met%ca, 2, mf) ! initialise leaf surface CO2 concentration
    met%tvair = met%tk
    met%qvair = met%qv
    canopy%tv = met%tvair
    canopy%fwsoil = 1.0

    CALL define_air (met, air)

    CALL qsatfjh(qstvair,met%tvair-C%tfrz,met%pmb)

    met%dva = (qstvair - met%qvair) *  C%rmair/C%rmh2o * met%pmb * 100.0
    dsx = met%dva     ! init. leaf surface vpd
    dsx= MAX(dsx,0.0)
    tlfx = met%tk  ! initialise leaf temp iteration memory variable (K)
    tlfy = met%tk  ! initialise current leaf temp (K)

    ortsoil = ssnow%rtsoil
    IF (cable_user%soil_struc=='sli') THEN
       ssnow%tss = REAL(ssnow%Tsurface) + C%tfrz
    ELSE
       ssnow%tss =  REAL((1-ssnow%isflag))*ssnow%tgg(:,1) +                    &
            REAL(ssnow%isflag)*ssnow%tggsn(:,1)
    ENDIF
    tss4 = ssnow%tss**4
    canopy%fes = 0.
    canopy%fess = 0.
    canopy%fesp = 0.
    ssnow%potev = 0.
    canopy%fevw_pot = 0.

    !L_REV_CORR - initialise sensitivity/ACCESS correction terms
    !NB %fes_cor is NOT initialised to zero at this point
    canopy%fhs_cor = 0.0
    canopy%fns_cor = 0.0
    canopy%ga_cor = 0.0
    !canopy%fes_cor = 0.0

    !L_REV_CORR - new working variables
    rttsoil = 0.
    rhlitt = 0.
    relitt = 0.
    alpm1  = 0.
    beta2  = 0.

    CALL radiation( ssnow, veg, air, met, rad, canopy )

    canopy%zetar(:,1) = C%ZETA0 ! stability correction terms
    canopy%zetar(:,2) = C%ZETPOS + 1
    canopy%zetash(:,1) = C%ZETA0 ! stability correction terms
    canopy%zetash(:,2) = C%ZETPOS + 1


    DO iter = 1, NITER

       ! AERODYNAMIC PROPERTIES: friction velocity us, thence turbulent
       ! resistances rt0, rt1 (elements of dispersion matrix):
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.46:
       CALL comp_friction_vel()

       ! E.Kowalczyk 2014
       IF (cable_user%l_new_roughness_soil)                                     &
            CALL ruff_resist(veg, rough, ssnow, canopy)




       ! Turbulent aerodynamic resistance from roughness sublayer depth
       ! to reference height, x=1 if zref+disp>zruffs,
       ! 0 otherwise: thus rt1usc = 0 if zref+disp<zruffs
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.53:
       xx = 0.5 + SIGN(0.5,rough%zref_tq+rough%disp-rough%zruffs)

       ! correction  by Ian Harman to the 2nd psis term
       rt1usc = xx * ( LOG( rough%zref_tq/MAX( rough%zruffs-rough%disp,         &
            rough%z0soilsn ) )               &
            - psis( canopy%zetar(:,iter) )                                  &
            + psis( canopy%zetar(:,iter) * ( MAX( rough%zruffs-rough%disp,  &
            rough%z0soilsn ) )            &
            / rough%zref_tq ) ) / C%VONK
       rt_min = 5.

       !! vh_js !!
       IF (cable_user%soil_struc=='sli') THEN
          ! for stable conditions, update rough%rt0us & rough%rt1usa by replacing C%CSW by
          ! csw = cd/2* (U(hc)/ust)**2 according to Eqs 15 & 19 from notes by Ian Harman (9-9-2011)
          WHERE (canopy%vlaiw > C%LAI_thresh .AND. rough%hruff > rough%z0soilsn)
             rt0bus = (LOG(0.1*rough%hruff/rough%z0soilsn) - psis(canopy%zetash(:,iter)) + &
                  psis(canopy%zetash(:,iter)*rough%z0soilsn/(0.1*rough%hruff))) / &
                  C%vonk/rough%term6a

             zstar = rough%disp + 1.5*(veg%hc - rough%disp)

             psihat = LOG((zstar - rough%disp)/ (veg%hc - rough%disp)) + &
                  (veg%hc - zstar)/(zstar - rough%disp)
             rL = -(C%vonk*C%grav*(zstar - rough%disp)*(canopy%fh))/ &  ! 1/Monin-Obokov Length
                  MAX( (air%rho*C%capp*met%tk*canopy%us**3), 1.e-12)
             phist = 1 + 5.0*(zstar - rough%disp)*rL

             WHERE (canopy%zetar(:,iter) .GT. 1.e-6)! stable conditions

                csw = MIN(0.3*((LOG((veg%hc-rough%disp)/rough%z0m) + phist*psihat - &
                     psim(canopy%zetar(:,iter)*(veg%hc-rough%disp)/(rough%zref_tq-rough%disp))+ &
                     psim(canopy%zetar(:,iter)*rough%z0m/(rough%zref_tq-rough%disp)))/0.4)**2/2., 3.0)* c%csw
             ELSEWHERE
                csw = c%csw
             endwhere

             rough%term2  = EXP( 2. * CSW * canopy%rghlai * &
                  ( 1 - rough%disp / rough%hruff ) )
             rough%term3  = C%A33**2 * C%CTL * 2. * CSW * canopy%rghlai
             rough%term5  = MAX( ( 2. / 3. ) * rough%hruff / rough%disp, 1.0 )
             rough%term6 =  EXP( 3. * rough%coexp * ( rough%disp / rough%hruff -1. ) )

             rough%rt0us  = LOG(rough%disp/(0.1 * rough%hruff)) * &
                  EXP(2. * C%CSW * canopy%rghlai) * rough%disp &
                  / rough%hruff / (c%a33 ** 2 * c%ctl)

             rough%rt1usa = rough%term5 * ( rough%term2 - 1.0 ) / rough%term3
             rt0 = MAX(rt_min, rough%rt0us+rt0bus) / canopy%us
          ELSEWHERE
             rt0 = MAX(rt_min,rough%rt0us) / canopy%us
          ENDWHERE

       ELSE
          rt0 = MAX(rt_min,rough%rt0us / canopy%us)

          IF (cable_user%litter) THEN
             ! Mathews (2006), A process-based model of offine fuel moisture,
             !                 International Journal of Wildland Fire 15,155-168
             ! assuming here u=1.0 ms-1, bulk litter density 63.5 kgm-3
             canopy%kthLitt = 0.3_r_2 ! ~ 0.2992125984251969 = 0.2+0.14*0.045*1000.0/63.5
             canopy%DvLitt = 3.1415841138194147e-05_r_2 ! = 2.17e-5*exp(1.0*2.6)*exp(-0.5*(2.08+(1.0*2.38)))
          ENDIF

       ENDIF


       ! Aerodynamic resistance (sum 3 height integrals)/us
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.50:
       rough%rt1 = MAX(5.,(rough%rt1usa + rough%rt1usb + rt1usc) / canopy%us)

       DO j=1,mp

          IF(canopy%vlaiw(j) > C%LAI_THRESH) THEN
             ssnow%rtsoil(j) = rt0(j)
          ELSE
             ssnow%rtsoil(j) = rt0(j) + rough%rt1(j)
          ENDIF

       ENDDO


       ssnow%rtsoil = MAX(rt_min,ssnow%rtsoil)

       DO j=1,mp

          IF( ssnow%rtsoil(j) > 2.*ortsoil(j) .OR.                              &
               ssnow%rtsoil(j) < 0.5*ortsoil(j) ) THEN

             ssnow%rtsoil(j) = MAX(rt_min,0.5*(ssnow%rtsoil(j) + ortsoil(j)))

          ENDIF

       ENDDO

       IF (cable_user%or_evap) THEN
          CALL or_soil_evap_resistance(soil,air,met,canopy,ssnow,veg,rough)
       END IF

       ! Vegetation boundary-layer conductance (mol/m2/s)
       ! C%prandt = kinematic viscosity/molecular diffusivity
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.12. Top leaf:
       DO j=1,mp

          IF(canopy%vlaiw(j) > C%LAI_THRESH) THEN
             gbvtop(j) = air%cmolar(j)*C%APOL * air%visc(j) / C%prandt /        &
                  veg%dleaf(j) * (canopy%us(j) / MAX(rough%usuh(j),1.e-6)&
                  * veg%dleaf(j) / air%visc(j) )**0.5                    &
                  * C%prandt**(1.0/3.0) / veg%shelrb(j)
             gbvtop(j) = MAX (0.05_r_2,gbvtop(j) )      ! for testing (BP aug2010)

             ! Forced convection boundary layer conductance
             ! (see Wang & Leuning 1998, AFM):


             !vh! inserted 'min' to avoid floating underflow
             gbhu(j,1) = gbvtop(j)*(1.0-EXP(-MIN(canopy%vlaiw(j)                    &
                  *(0.5*rough%coexp(j)+rad%extkb(j) ),20.0))) /            &
                  (rad%extkb(j)+0.5*rough%coexp(j))

             gbhu(j,2) = (2.0/rough%coexp(j))*gbvtop(j)*  &
                  (1.0-EXP(-MIN(0.5*rough%coexp(j)*canopy%vlaiw(j),20.0))) &
                  - gbhu(j,1)
          ENDIF

       ENDDO

       rny = SUM(rad%rniso,2) ! init current estimate net rad
       hcy = 0.0              ! init current estimate lat heat
       ecy = rny - hcy        ! init current estimate lat heat

       sum_rad_rniso = SUM(rad%rniso,2)

       CALL dryLeaf( dels, rad, rough, air, met,                                &
            veg, canopy, soil, ssnow, dsx,                             &
            fwsoil, tlfx, tlfy, ecy, hcy,                              &
            rny, gbhu, gbhf, csx, cansat,                              &
            ghwet,  iter,climate )


       CALL wetLeaf( dels, rad, rough, air, met,                                &
            veg, canopy, cansat, tlfy,                                 &
            gbhu, gbhf, ghwet )


       ! Calculate latent heat from vegetation:
       ! Calculate sensible heat from vegetation:
       ! Calculate net rad absorbed by canopy:
       canopy%fev = REAL(canopy%fevc + canopy%fevw)
       ftemp = (1.0 - canopy%fwet) *  REAL(hcy) + canopy%fhvw
       canopy%fhv = REAL(ftemp)
       ftemp= (1.0-canopy%fwet)*REAL(rny)+canopy%fevw+canopy%fhvw
       canopy%fnv = REAL(ftemp)

       ! canopy rad. temperature calc from long-wave rad. balance
       sum_rad_gradis = SUM(rad%gradis,2)

       DO j=1,mp

          IF ( canopy%vlaiw(j) > C%LAI_THRESH .AND.                             &
               rough%hruff(j) > rough%z0soilsn(j) ) THEN

             rad%lwabv(j) = C%CAPP * C%rmair * ( tlfy(j) - met%tk(j) ) *        &
                  sum_rad_gradis(j)
             !! vh_js !!

             IF (  (rad%lwabv(j) / (2.0*(1.0-rad%transd(j))            &
                  * C%SBOLTZ*C%EMLEAF)+met%tvrad(j)**4) .GT. 0.0) THEN

                canopy%tv(j) = (rad%lwabv(j) / (2.0*(1.0-rad%transd(j))            &
                     * C%SBOLTZ*C%EMLEAF)+met%tvrad(j)**4)**0.25

             ELSE
                canopy%tv(j) = met%tvrad(j)
             ENDIF


          ELSE! sparse canopy

             canopy%tv(j) = met%tvrad(j)

          ENDIF

       ENDDO


       ! Calculate net rad to soil:
       canopy%fns = rad%qssabs + rad%transd*met%fld + (1.0-rad%transd)*C%EMLEAF* &
            C%SBOLTZ*canopy%tv**4 - C%EMSOIL*C%SBOLTZ* tss4


       ! Saturation specific humidity at soil/snow surface temperature:
       CALL qsatfjh(ssnow%qstss,ssnow%tss-C%tfrz,met%pmb)

       IF (cable_user%gw_model .OR.  cable_user%or_evap) &
            CALL pore_space_relative_humidity(ssnow,soil,veg)

       IF (cable_user%soil_struc=='default') THEN

          !REV_CORR - single location for calculating litter resistances
          !can go earlier in the code if needed
          IF (cable_user%litter) THEN
             rhlitt = REAL((1-ssnow%isflag))*veg%clitt*0.003/ &
                  canopy%kthLitt/(air%rho*C%CAPP)
             relitt = REAL((1-ssnow%isflag))*veg%clitt*0.003/ &
                  canopy%DvLitt
          ENDIF

          IF(cable_user%ssnow_POTEV== "P-M") THEN

             !--- uses %ga from previous timestep
             ssnow%potev = Penman_Monteith(canopy%ga)

          ELSE !by default assumes Humidity Deficit Method

             ! Humidity deficit
             ! INH: I think this should be - met%qvair
             dq = ssnow%qstss - met%qv
             dq_unsat = ssnow%rh_srf*ssnow%qstss - met%qv
             ssnow%potev =  Humidity_deficit_method(dq, dq_unsat,ssnow%qstss)

          ENDIF

          ! Soil latent heat:
          CALL latent_heat_flux()

          ! Calculate soil sensible heat:
          ! INH: I think this should be - met%tvair
          !canopy%fhs = air%rho*C%CAPP*(ssnow%tss - met%tk) /ssnow%rtsoil
          IF (cable_user%gw_model .OR. cable_user%or_evap) THEN
             canopy%fhs =  air%rho*C%CAPP*(ssnow%tss - met%tk) / &
                  (ssnow%rtsoil + ssnow%rt_qh_sublayer)
             !note if or_evap and litter are true then litter resistance is
             !incluyded above in ssnow%rt_qh_sublayer
          ELSEIF (cable_user%litter) THEN
             !! vh_js !! account for additional litter resistance to sensible heat transfer
             !! INH simplifying code using rhlitt
             canopy%fhs =  air%rho*C%CAPP*(ssnow%tss - met%tk) / &
                                !(ssnow%rtsoil + real((1-ssnow%isflag))*veg%clitt*0.003/canopy%kthLitt/(air%rho*C%CAPP))
                  (ssnow%rtsoil + rhlitt)
          ELSE
             canopy%fhs = air%rho*C%CAPP*(ssnow%tss - met%tvair) /ssnow%rtsoil
          ENDIF

       ELSE


          ! SLI SEB to get canopy%fhs, canopy%fess, canopy%ga
          ! (Based on old Tsoil, new canopy%tv, new canopy%fns)
          CALL sli_main(1,dels,veg,soil,ssnow,met,canopy,air,rad,1)

       ENDIF

       CALL within_canopy( gbhu, gbhf, rt0, rhlitt, relitt )

       ! Saturation specific humidity at soil/snow surface temperature:
       CALL qsatfjh(ssnow%qstss,ssnow%tss-C%tfrz,met%pmb)

       IF (cable_user%soil_struc=='default') THEN

          IF(cable_user%ssnow_POTEV== "P-M") THEN

             !--- uses %ga from previous timestep
             ssnow%potev = Penman_Monteith(canopy%ga)

          ELSE !by default assumes Humidity Deficit Method

             ! Humidity deficit
             dq = ssnow%qstss - met%qvair
             dq_unsat = ssnow%rh_srf*ssnow%qstss - met%qvair
             ssnow%potev =  Humidity_deficit_method(dq, dq_unsat,ssnow%qstss)

          ENDIF

          ! Soil latent heat:
          CALL latent_heat_flux()

          ! Soil sensible heat:
          !canopy%fhs = air%rho*C%CAPP*(ssnow%tss - met%tvair) /ssnow%rtsoil
          IF (cable_user%gw_model .OR. cable_user%or_evap) THEN
             canopy%fhs =  air%rho*C%CAPP*(ssnow%tss - met%tvair) / &
                  (ssnow%rtsoil + REAL(ssnow%rt_qh_sublayer))

          ELSEIF (cable_user%litter) THEN
             !! vh_js !! account for additional litter resistance to sensible heat transfer
             !! INH simplifying code using rhlitt
             canopy%fhs =  air%rho*C%CAPP*(ssnow%tss - met%tvair) / &
                                !(ssnow%rtsoil +  real((1-ssnow%isflag))*veg%clitt*0.003/canopy%kthLitt/(air%rho*C%CAPP))
                  (ssnow%rtsoil + rhlitt)
          ELSE
             canopy%fhs = air%rho*C%CAPP*(ssnow%tss - met%tvair) /ssnow%rtsoil
          ENDIF

          !! Ticket #90 ssnow%cls factor should be retained: required for energy balance
          !! INH: %cls factor included in %fes already - do not include here
          canopy%ga = canopy%fns-canopy%fhs-canopy%fes !*ssnow%cls

       ELSE

          ! SLI SEB to get canopy%fhs, canopy%fess, canopy%ga
          ! (Based on old Tsoil, new canopy%tv, new canopy%fns)
          CALL sli_main(1,dels,veg,soil,ssnow,met,canopy,air,rad,1)

       ENDIF

       ! Set total latent heat:
       canopy%fe = canopy%fev + canopy%fes

       ! Set total sensible heat:
       canopy%fh = canopy%fhv + canopy%fhs

       !---diagnostic purposes
       DO j=1,mp

          IF (ssnow%potev(j) .GE. 0.) THEN
             ssnow%potev(j) = MAX(0.00001,ssnow%potev(j))
          ELSE
             ssnow%potev(j) = MIN(-0.0002,ssnow%potev(j))
          ENDIF

          IF (canopy%fevw_pot(j) .GE. 0.) THEN
             canopy%fevw_pot(j) = MAX(0.000001,canopy%fevw_pot(j))
          ELSE
             canopy%fevw_pot(j) = MIN(-0.002,canopy%fevw_pot(j))
          ENDIF

       ENDDO


       canopy%rnet = canopy%fnv + canopy%fns

       !INH: If PM routine corrected then match changes here
       canopy%epot = ((1.-rad%transd)*canopy%fevw_pot +                         &
            rad%transd*ssnow%potev*ssnow%cls) * dels/air%rlam



       canopy%rniso = SUM(rad%rniso,2) + rad%qssabs + rad%transd*met%fld + &
            (1.0-rad%transd)*C%EMLEAF* &
            C%SBOLTZ*met%tvrad**4 - C%EMSOIL*C%SBOLTZ*met%tvrad**4

       rlower_limit = canopy%epot * air%rlam / dels
       WHERE (rlower_limit == 0 ) rlower_limit = 1.e-7 !prevent from 0. by adding 1.e-7 (W/m2)


       canopy%wetfac_cs = MAX(0., MIN(1.0,canopy%fe / rlower_limit ))

       DO j=1,mp

          IF ( canopy%wetfac_cs(j) .LE. 0. )                                    &
               canopy%wetfac_cs(j) = MAX( 0., MIN( 1.,                            &
               MAX( canopy%fev(j)/canopy%fevw_pot(j),       &
               REAL(canopy%fes(j))/ssnow%potev(j) ) ) )

       ENDDO

       CALL update_zetar()

    END DO           ! do iter = 1, NITER



    canopy%cduv = canopy%us * canopy%us / (MAX(met%ua,C%UMIN))**2

    !---diagnostic purposes
    canopy%gswx_T = rad%fvlai(:,1)/MAX( C%LAI_THRESH, canopy%vlaiw(:) )         &
         * canopy%gswx(:,1) + rad%fvlai(:,2) / MAX(C%LAI_THRESH,     &
         canopy%vlaiw(:))*canopy%gswx(:,2)

    ! The surface conductance below is required by dust scheme; it is composed from canopy and soil conductances
    canopy%gswx_T = (1.-rad%transd)*MAX(1.e-06,canopy%gswx_T ) +  &   !contribution from  canopy conductance
         rad%transd*(.01*ssnow%wb(:,1)/soil%sfc)**2 ! + soil conductance; this part is done as in Moses
    WHERE ( soil%isoilm == 9 ) canopy%gswx_T = 1.e6   ! this is a value taken from Moses for ice points

    canopy%cdtq = canopy%cduv *( LOG( rough%zref_uv / rough%z0m) -              &
         psim( canopy%zetar(:,NITER) * rough%zref_uv/rough%zref_tq )   &
         + psim( canopy%zetar(:,NITER) * rough%z0m/rough%zref_tq ) & ! new term from Ian Harman
         ) / ( LOG( rough%zref_tq /(0.1*rough%z0m) )                   &
         - psis( canopy%zetar(:,NITER))                                  &
         + psis(canopy%zetar(:,NITER)*0.1*rough%z0m/rough%zref_tq) ) ! n

    !INH - the screen level calculations should be split off into a new subroutine -------

    ! Calculate screen temperature: 1) original method from SCAM
    ! screen temp., windspeed and relative humidity at 1.5m
    ! screen temp., windspeed and relative humidity at 2.0m
    ! cls factor included in qstar
    tstar = - canopy%fh / ( air%rho*C%CAPP*canopy%us)
    qstar = - canopy%fe / ( air%rho*air%rlam *canopy%us * ssnow%cls)
    zscrn = MAX(rough%z0m,2.0-rough%disp)
    ftemp = ( LOG(rough%zref_tq/zscrn)- psis(canopy%zetar(:,iterplus)) +       &
         psis(canopy%zetar(:,iterplus) * zscrn / rough%zref_tq) ) /C%VONK

    ! Calculate screen temperature:
    canopy%tscrn = met%tk - C%tfrz - tstar * ftemp

    ! Calculate radiative/skin temperature;
    ! at this stage old soil temperature is used
    ! calculation of screen temepratures for LAI > 0.1 . Method by Ian Harman

    term1=0.
    term2=0.
    term5=0.
    term3 = 0. ! Work around for Intel compiler problem with nested whres
    r_sc = 0.  ! sum of resistance from ground to screen level
    zscl = MAX(rough%z0soilsn,2.0)

    ! assume screen temp of bareground if all these conditions are not met
    DO j=1,mp

       IF ( canopy%vlaiw(j) > C%LAI_THRESH .AND. rough%hruff(j) > 0.01) THEN

          IF ( rough%disp(j)  > 0.0 ) THEN

             term1(j) = EXP(2*C%CSW*canopy%rghlai(j)*(1-zscl(j)/rough%hruff(j)))
             term2(j) = EXP(2*C%CSW*canopy%rghlai(j) *                          &
                  (1-rough%disp(j)/rough%hruff(j)))
             term5(j) = MAX(2./3.*rough%hruff(j)/rough%disp(j), 1.)

          ENDIF

          term3(j) = C%A33**2*C%CTL*2*C%CSW*canopy%rghlai(j)

          IF( zscl(j) < rough%disp(j) ) THEN

             !Ticket #154
             !r_sc(j) = term5(j) * LOG(zscl(j)/rough%z0soilsn(j)) *              &
             !     ( EXP(2*C%CSW*canopy%rghlai(j)) - term1(j) ) / term3(j)
             r_sc(j) = term5(j) * LOG(zscl(j)/rough%z0soilsn(j)) *              &
                  ( EXP(2*C%CSW*canopy%rghlai(j)) - term2(j) ) / term3(j)
             r_sc(j) = r_sc(j) + term5(j) * LOG(rough%disp(j)/rough%z0soilsn(j)) *  &
                  ( EXP(2*C%CSW*canopy%rghlai(j)) - term1(j) ) / term3(j)

          ELSEIF( rough%disp(j) <= zscl(j) .AND.                                &
               zscl(j) < rough%hruff(j) ) THEN

             r_sc(j) = rough%rt0us(j) + term5(j) * ( term2(j) - term1(j) ) /    &
                  term3(j)

          ELSEIF( rough%hruff(j) <= zscl(j) .AND.                               &
               zscl(j) <  rough%zruffs(j) ) THEN

             r_sc(j) = rough%rt0us(j) + rough%rt1usa(j) + term5(j) *            &
                  ( zscl(j) - rough%hruff(j) ) /                           &
                  ( C%A33**2 * C%CTL * rough%hruff(j) )


          ELSEIF( zscl(j) >= rough%zruffs(j) ) THEN
             !Ticket #67 - Modify order of operations to avoid potential error
             r_sc(j) = rough%rt0us(j) + rough%rt1usa(j) + rough%rt1usb(j) +     &
                  ( LOG( (zscl(j) - rough%disp(j)) /                       &
                  MAX( rough%zruffs(j)-rough%disp(j),                      &
                  rough%z0soilsn(j) ) ) - psis( (zscl(j)-rough%disp(j))    &
                                !Ticket #67 - change order of operations to avoid /0
                                !        / (rough%zref_tq(j)/canopy%zetar(j,iterplus) ) )        &
                  * canopy%zetar(j,iterplus)/rough%zref_tq(j) )            &
                  + psis( (rough%zruffs(j) - rough%disp(j) )               &
                                !        / (rough%zref_tq(j)/canopy%zetar(j,iterplus ) ) ) )     &
                  * canopy%zetar(j,iterplus)/rough%zref_tq(j) ) )          &
                  / C%VONK

          ENDIF

          !extensions for litter and Or evaporation model
          IF (cable_user%litter) THEN
             canopy%tscrn(j) = ssnow%tss(j) + (met%tk(j) - ssnow%tss(j)) *     &
                  MIN(1., ( (r_sc(j)+rhlitt(j)*canopy%us(j))  / MAX( 1.,          &
                  rough%rt0us(j) + rough%rt1usa(j) + rough%rt1usb(j)              &
                  + rt1usc(j) + rhlitt(j)*canopy%us(j) )) ) - C%tfrz
          ELSEIF (cable_user%or_evap .OR. cable_user%gw_model) THEN
             canopy%tscrn(j) = ssnow%tss(j) + (met%tk(j) - ssnow%tss(j)) *     &
                  MIN(1., ( (ssnow%rt_qh_sublayer(j)*canopy%us(j) + r_sc(j) ) /   &
                  MAX( 1., rough%rt0us(j) + rough%rt1usa(j) + rough%rt1usb(j)     &
                  + rt1usc(j) + ssnow%rt_qh_sublayer(j)*canopy%us(j) )) ) - C%tfrz
          ELSE
             canopy%tscrn(j) = ssnow%tss(j) + (met%tk(j) - ssnow%tss(j)) *      &
                  MIN(1., (r_sc(j) / MAX( 1.,                            &
                  rough%rt0us(j) + rough%rt1usa(j) + rough%rt1usb(j)   &
                  + rt1usc(j))) )  - C%tfrz
          ENDIF

       ENDIF

    ENDDO


    !screen level humdity - this is only approximate --------------------------
    CALL qsatfjh(rsts,canopy%tscrn,met%pmb)

    qtgnet = rsts * ssnow%wetfac - met%qv

    DO j=1,mp

       IF (qtgnet(j) .GT. 0. ) THEN
          qsurf(j) = rsts(j) * ssnow%wetfac(j)
       ELSE
          qsurf(j) = 0.1*rsts(j)*ssnow%wetfac(j) + 0.9*met%qv(j)
       ENDIF

       canopy%qscrn(j) = met%qv(j) - qstar(j) * ftemp(j)

       IF( canopy%vlaiw(j) >C%LAI_THRESH .AND. rough%hruff(j) > 0.01) THEN

          !extensions for litter and Or model
          IF (cable_user%litter) THEN
             canopy%qscrn(j) = qsurf(j) + (met%qv(j) - qsurf(j)) *             &
                  MIN(1., ( ( r_sc(j)+relitt(j)*canopy%us(j) ) / MAX( 1.,         &
                  rough%rt0us(j) + rough%rt1usa(j) + rough%rt1usb(j)            &
                  + rt1usc(j) + relitt(j)*canopy%us(j) )) )

          ELSEIF (cable_user%or_evap .OR. cable_user%gw_model) THEN
             !using alpm1 as a dumy variable
             alpm1(j) = REAL(&
                  ssnow%satfrac(j)/(REAL(ssnow%rtsoil(j),r_2)+&
                  ssnow%rtevap_sat(j)) &
                  + (1.0-ssnow%satfrac(j))/(REAL(ssnow%rtsoil(j),r_2)+ ssnow%rtevap_unsat(j)) &
                  )

             canopy%qscrn(j) = qsurf(j) + (met%qv(j) - qsurf(j)) *             &
                  MIN(1., ( (r_sc(j) + canopy%us(j)/alpm1(j) ) / MAX( 1.,         &
                  rough%rt0us(j) + rough%rt1usa(j) + rough%rt1usb(j)              &
                  + rt1usc(j) + canopy%us(j)/alpm1(j) )) )

          ELSE
             canopy%qscrn(j) = qsurf(j) + (met%qv(j) - qsurf(j)) *             &
                  MIN(1., (r_sc(j) / MAX( 1.,                                     &
                  rough%rt0us(j) + rough%rt1usa(j) + rough%rt1usb(j)              &
                  + rt1usc(j))) )
          ENDIF

       ENDIF

    ENDDO

    ! Calculate dewfall: from negative lh wet canopy + neg. lh dry canopy:
    canopy%dewmm = - (MIN(0.0,canopy%fevw) + MIN(0.0_r_2,canopy%fevc)) * &
         dels * 1.0e3 / (C%RHOW*air%rlam)

    ! Add dewfall to canopy water storage:
    canopy%cansto = canopy%cansto + canopy%dewmm

    ! Modify canopy water storage for evaporation:
    canopy%cansto = MAX(canopy%cansto-MAX(0.0,REAL(canopy%fevw))*dels &
         *1.0e3/(C%RHOW*air%rlam), 0.0)

    ! Calculate canopy water storage excess:
    canopy%spill=MAX(0.0, canopy%cansto-cansat)

    ! Move excess canopy water to throughfall:
    ! %through is /dels in UM app. (unpacked in hyd driver) for STASH output
    canopy%through = canopy%through + canopy%spill

    ! Initialise 'throughfall to soil' as 'throughfall from canopy';
    ! snow may absorb
    canopy%precis = MAX(0.,canopy%through)

    ! Update canopy storage term:
    canopy%cansto=canopy%cansto - canopy%spill

    ! Calculate the total change in canopy water store (mm/dels):
    canopy%delwc = canopy%cansto-canopy%oldcansto

    ! calculate dgdtg, derivative of ghflux 3 instances
    ! d(canopy%fns)/d(ssnow%tgg)
    ! d(canopy%fhs)/d(ssnow%tgg)
    ! d(canopy%fes)/d(dq)
    !IF (cable_user%soil_struc=='default') THEN
    ssnow%dfn_dtg = (-1.)*4.*C%EMSOIL*C%SBOLTZ*tss4/ssnow%tss

    !INH: REV_CORR revised sensitivity terms working variable
    rttsoil = ssnow%rtsoil
    IF (cable_user%L_REV_CORR) THEN
       WHERE (canopy%vlaiw > C%LAI_THRESH)
          !if %vlaiw<=%LAI_THRESH then %rt1 already added to %rtsoil
          rttsoil = rttsoil + rough%rt1
       ENDWHERE
    ENDIF

    IF (cable_user%gw_model .OR. cable_user%or_evap) THEN

       ssnow%dfh_dtg = air%rho*C%CAPP/(ssnow%rtsoil+ REAL(ssnow%rt_qh_sublayer))

       !! INH simplifying code for legibility
       !ssnow%dfe_ddq = real(ssnow%satfrac)*air%rho*air%rlam*ssnow%cls/ &
       !     (ssnow%rtsoil+ real(ssnow%rtevap_sat))  +
       !     (1.0-real(ssnow%satfrac))*real(ssnow%rh_srf)*&
       !      air%rho*air%rlam*ssnow%cls/ (ssnow%rtsoil+
       !      real(ssnow%rtevap_unsat) )
       ssnow%dfe_ddq = REAL(ssnow%satfrac)/(ssnow%rtsoil+ REAL(ssnow%rtevap_sat))  &
            + (1.0-REAL(ssnow%satfrac))*REAL(ssnow%rh_srf)                   &
            / (ssnow%rtsoil+ REAL(ssnow%rtevap_unsat) )

       !mrd561 fixes.  Do same thing as INH but has been tested.
       IF (cable_user%L_REV_CORR) THEN
          alpm1  = REAL(ssnow%satfrac/(REAL(ssnow%rtsoil,r_2)+ ssnow%rtevap_sat) +     &
               (1.0-ssnow%satfrac) / (REAL(ssnow%rtsoil,r_2)+ ssnow%rtevap_unsat ) )
          beta2 = REAL(ssnow%satfrac/(REAL(ssnow%rtsoil,r_2)+ ssnow%rtevap_sat) +     &
               (1.0-ssnow%satfrac) * ssnow%rh_srf                  &
               / (REAL(ssnow%rtsoil,r_2)+ ssnow%rtevap_unsat ) )
          WHERE (canopy%vlaiw > C%LAI_THRESH)
             alpm1 = alpm1 + 1._r_2/REAL(rough%rt1,r_2)
             beta_div_alpm  = beta2 / alpm1  !might need limit here
             rttsoil = ssnow%rtsoil + rough%rt1
          ELSEWHERE!if there is no canopy then qa should not change
             beta_div_alpm=0.0  !do not divide by aplm1 prevent issues
             rttsoil = ssnow%rtsoil
          ENDWHERE
          ssnow%dfh_dtg = air%rho*C%CAPP/(rttsoil +               &
               REAL(ssnow%rt_qh_sublayer))
          ssnow%dfe_ddq = REAL(ssnow%satfrac*(1.0-REAL(beta_div_alpm,r_2)) /        &
               (REAL(ssnow%rtsoil,r_2)+ ssnow%rtevap_sat) +           &
               (1.0-ssnow%satfrac)* (ssnow%rh_srf - REAL(beta_div_alpm,r_2)) /    &
               (REAL(ssnow%rtsoil,r_2)+ ssnow%rtevap_unsat ) )

       ELSE
          ssnow%dfh_dtg = air%rho*C%CAPP/(ssnow%rtsoil+ REAL(ssnow%rt_qh_sublayer))

          ssnow%dfe_ddq = REAL(ssnow%satfrac)/(ssnow%rtsoil+ REAL(ssnow%rtevap_sat))  &
               + (1.0-REAL(ssnow%satfrac))*REAL(ssnow%rh_srf)                   &
               / (ssnow%rtsoil+ REAL(ssnow%rtevap_unsat) )
       ENDIF

       !cls applies for both REV_CORR false and true
       ssnow%dfe_ddq = ssnow%dfe_ddq*air%rho*air%rlam*ssnow%cls

       !REV_CORR: factor %wetfac needed for potev>0. and gw_model &/or snow cover
       !NB %wetfac=1. if or_evap
       IF (cable_user%L_REV_CORR) THEN
          WHERE (ssnow%potev >= 0.)
             ssnow%dfe_ddq = ssnow%dfe_ddq*ssnow%wetfac
          ENDWHERE
       ENDIF


    ELSEIF (cable_user%litter) THEN
       !!vh_js!! INH simplifying code for legibility and REV_CORR
       !ssnow%dfh_dtg = air%rho*C%CAPP/(ssnow%rtsoil+ &
       !     real((1-ssnow%isflag))*veg%clitt*0.003/canopy%kthLitt/(air%rho*C%CAPP))
       !ssnow%dfe_ddq = ssnow%wetfac*air%rho*air%rlam*ssnow%cls/ &
       !     (ssnow%rtsoil+ real((1-ssnow%isflag))*veg%clitt*0.003/canopy%DvLitt)

       !recalculated - probably not needed
       rhlitt = REAL((1-ssnow%isflag))*veg%clitt*0.003/canopy%kthLitt/(air%rho*C%CAPP)
       relitt = REAL((1-ssnow%isflag))*veg%clitt*0.003/canopy%DvLitt

       !incorporates REV_CORR changes
       ssnow%dfh_dtg = air%rho*C%CAPP/(rttsoil+rhlitt)
       ssnow%dfe_ddq = ssnow%wetfac*air%rho*air%rlam*ssnow%cls/(rttsoil+relitt)

       !REV_CORR: factor ssnow%wetfac is not applied if dew/frost i.e. potev<0
       IF (cable_user%L_REV_CORR) THEN
          WHERE (ssnow%potev < 0.)
             ssnow%dfe_ddq = air%rho*air%rlam*ssnow%cls/(rttsoil+relitt)
          ENDWHERE
       ENDIF

    ELSE
       !ssnow%dfh_dtg = air%rho*C%CAPP/ssnow%rtsoil
       !ssnow%dfe_ddq = ssnow%wetfac*air%rho*air%rlam*ssnow%cls/ssnow%rtsoil

       !incorporates REV_CORR changes
       ssnow%dfh_dtg = air%rho*C%CAPP/rttsoil
       ssnow%dfe_ddq = ssnow%wetfac*air%rho*air%rlam*ssnow%cls/rttsoil

       !REV_CORR: factor ssnow%wetfac is not applied if dew/frost i.e. potev<0
       IF (cable_user%L_REV_CORR) THEN
          WHERE (ssnow%potev < 0.)
             ssnow%dfe_ddq = air%rho*air%rlam*ssnow%cls/(rttsoil+relitt)
          ENDWHERE
       ENDIF

    ENDIF

    ssnow%ddq_dtg = (C%rmh2o/C%rmair) /met%pmb * C%TETENA*C%TETENB * C%TETENC   &
         / ( ( C%TETENC + ssnow%tss-C%tfrz )**2 )*EXP( C%TETENB *       &
         ( ssnow%tss-C%tfrz ) / ( C%TETENC + ssnow%tss-C%tfrz ) )

    !canopy%dgdtg = ssnow%dfn_dtg - ssnow%dfh_dtg - ssnow%dfe_ddq *    &
    !     ssnow%ddq_dtg

    !INH: REV_CORR Rewritten for flexibility
    ssnow%dfe_dtg = ssnow%dfe_ddq * ssnow%ddq_dtg
    canopy%dgdtg = ssnow%dfn_dtg - ssnow%dfh_dtg - ssnow%dfe_dtg

    !ENDIF

    bal%drybal = REAL(ecy+hcy) - SUM(rad%rniso,2)                               &
         + C%CAPP*C%rmair*(tlfy-met%tk)*SUM(rad%gradis,2)  ! YP nov2009

    bal%wetbal = canopy%fevw + canopy%fhvw - SUM(rad%rniso,2) * canopy%fwet      &
         + C%CAPP*C%rmair * (tlfy-met%tk) * SUM(rad%gradis,2) *          &
         canopy%fwet  ! YP nov2009

    DEALLOCATE(cansat,gbhu)
    DEALLOCATE(dsx, fwsoil, tlfx, tlfy)
    DEALLOCATE(ecy, hcy, rny)
    DEALLOCATE(gbhf, csx)
    DEALLOCATE(ghwet)

    RETURN

  CONTAINS

    ! ------------------------------------------------------------------------------

    SUBROUTINE comp_friction_vel()
      USE cable_def_types_mod, ONLY : mp
      REAL, DIMENSION(mp)  :: lower_limit, rescale

      !INH: Ticket #138 %us is defined based on U(rough%zref_uv)
      ! but zetar based on rough%zref_tq - changes to ensure consistency
      !NB no RSL incorporated here

      !psim_1 = psim(canopy%zetar(:,iter))
      psim_1 = psim(canopy%zetar(:,iter)*rough%zref_uv/rough%zref_tq)

      rescale = C%VONK * MAX(met%ua,C%UMIN)
      z_eff = rough%zref_uv / rough%z0m

      !psim_arg = canopy%zetar(:,iter) / z_eff
      psim_arg = canopy%zetar(:,iter) * rough%z0m / rough%zref_tq

      !---fix for compiler limitation. bitwise reproducable whilst we
      !---we know it to 11th decimal. psim_arg typically of a few
      !psim_arg = nint(psim_arg * 1.e11)*1.e-11

      psim_2 = psim( psim_arg )

      lower_limit = rescale / ( LOG(z_eff) - psim_1 + psim_2 )

      canopy%us = MIN(MAX(1.e-6, lower_limit ), 10.0 )

    END SUBROUTINE comp_friction_vel

    ! ------------------------------------------------------------------------------

    FUNCTION Penman_Monteith( ground_H_flux ) RESULT(ssnowpotev)
      USE cable_def_types_mod, ONLY : mp
      REAL, INTENT(IN), DIMENSION(mp)  :: ground_H_flux
      REAL, DIMENSION(MP)  ::                                                     &
           ssnowpotev,      & ! returned result of function
           sss,             & ! var for Penman-Monteith soil evap
           cc1,             & ! var for Penman-Monteith soil evap
           cc2,             & ! var for Penman-Monteith soil evap
           qsatfvar           !
      INTEGER :: j

      ! Penman-Monteith formula
      sss=air%dsatdk
      cc1=sss/(sss+air%psyc )
      cc2=air%psyc /(sss+air%psyc )

      CALL qsatfjh(qsatfvar,met%tvair-C%tfrz,met%pmb)

      !INH 10-1-2017 - this P-M implementation is incorrect over snow.
      !variable ssnowpotev is actually the latent heat flux associated with
      !potential evaporation.
      !Needs to be addressed/simplified at a later date - involves changes
      !to HDM method and latent_heat_flux() and elsewhere

      IF (cable_user%litter) THEN
         !! vh_js !!
         ssnowpotev = cc1 * (canopy%fns - ground_H_flux) + &
              cc2 * air%rho * air%rlam*(qsatfvar - met%qvair)/ &
              (ssnow%rtsoil+ REAL((1-ssnow%isflag))*veg%clitt*0.003/canopy%DvLitt)
      ELSE
         ssnowpotev = cc1 * (canopy%fns - ground_H_flux) + &
              cc2 * air%rho * air%rlam*(qsatfvar  - met%qvair)/ssnow%rtsoil
      ENDIF

    END FUNCTION Penman_Monteith


    ! ------------------------------------------------------------------------------
    ! method alternative to P-M formula above
    FUNCTION humidity_deficit_method(dq,dqu,qstss ) RESULT(ssnowpotev)

      USE cable_def_types_mod, ONLY : mp

      REAL, DIMENSION(mp) ::                                                      &
           ssnowpotev,    & !
           dq,            & ! sat spec hum diff.
           dqu,           & ! sat spec hum diff.
           qstss             !dummy var for compilation

      INTEGER :: j
      REAL, DIMENSION(mp) :: q_air

      q_air = qstss - dq

      DO j=1,mp
         !if(ssnow%snowd(j) > 1.0) dq(j) = max( -0.1e-3, dq(j))
         IF( ssnow%snowd(j)>1.0 .OR. ssnow%tgg(j,1).EQ.C%tfrz)   THEN
            dq(j) = MAX( -0.1e-3, dq(j))
            dqu(j) = MAX( -0.1e-3, dqu(j))
         END IF

         IF (dq(j) .LE. 0.0 .AND. dqu(j) .LT. dq(j)) THEN
            dqu(j) = dq(j)
         END IF

         IF (dq(j) .GE. 0.0 .AND. dqu(j) .LT. 0.0) THEN
            dqu(j) = 0.0
         ENDIF
      ENDDO

      IF (cable_user%or_evap .OR. cable_user%gw_model) THEN

         IF (cable_user%or_evap) THEN
            DO j=1,mp

               IF (veg%iveg(j) .LT. 16 .AND. ssnow%snowd(j) .LT. 1e-7) THEN

                  IF (dq(j) .LE. 0.0) THEN
                     ssnow%rtevap_sat(j) = MIN(rtevap_max,canopy%sublayer_dz(j)/rt_Dff)
                  END IF

                  IF (dqu(j) .LE. 0.0) THEN
                     ssnow%rtevap_unsat(j) = MIN(rtevap_max,canopy%sublayer_dz(j)/rt_Dff)
                  END IF

               END IF

            END DO

         END IF

         ssnowpotev = air%rho * air%rlam * ( &
              REAL(ssnow%satfrac) * dq /(ssnow%rtsoil + REAL(ssnow%rtevap_sat)) + &
              (1.0 - REAL(ssnow%satfrac))* dqu/( &
              ssnow%rtsoil + REAL(ssnow%rtevap_unsat)) )

      ELSEIF (cable_user%litter) THEN
         !! vh_js !!
         ssnowpotev =air%rho * air%rlam * dq /(ssnow%rtsoil + &
              REAL((1-ssnow%isflag))*veg%clitt*0.003/canopy%DvLitt)
      ELSE
         ssnowpotev =air%rho * air%rlam * dq /ssnow%rtsoil
      ENDIF

    END FUNCTION Humidity_deficit_method

    ! ------------------------------------------------------------------------------

    SUBROUTINE Latent_heat_flux()

      USE cable_common_module
      USE cable_def_types_mod, ONLY : mp

      REAL, DIMENSION(mp) ::                                                      &
           frescale,  flower_limit, fupper_limit

      INTEGER :: j

      !Ticket 137 - adjustments made so that one of four cases occurs
      !             i) evaporation from/dew onto surfaces with no snow, T>Tfrz
      !            ii) sublimation from/frost onto surfaces with snow cover
      !           iii) evaporation of liquid water if no snow but frozen soil
      !            iv) deposition of frost onto frozen soils if no snow cover
      !
      !IMPORTANTLY the value of %cls set here is used to control whether
      !water fluxes are from the snow pack or soil column in _soilsnow

      ! Soil latent heat:
      WHERE (ssnow%potev < 0. ) ssnow%wetfac(:) = 1.0
      canopy%fess= ssnow%wetfac * ssnow%potev
      !removed below, set wetfac to 1.0 is potev < 0.0
      !WHERE (ssnow%potev < 0. ) canopy%fess = ssnow%potev

      ! Reduce soil evap due to presence of puddle
      pwet = MAX(0.,MIN(0.2,ssnow%pudsto/MAX(1.,ssnow%pudsmx)))
      canopy%fess = canopy%fess * (1.-pwet)

      frescale = soil%zse(1) * 1000. * air%rlam / dels

      DO j=1,mp

         IF(ssnow%snowd(j) < 0.1 .AND. canopy%fess(j) .GT. 0. ) THEN

            IF (.NOT.cable_user%l_new_reduce_soilevp) THEN
               flower_limit(j) = REAL(ssnow%wb(j,1))-soil%swilt(j)/2.0
            ELSE
               ! E.Kowalczyk 2014 - reduces the soil evaporation
               flower_limit(j) = REAL(ssnow%wb(j,1))-soil%swilt(j)
            ENDIF
            fupper_limit(j) = MAX( 0.,                                        &
                 flower_limit(j) * frescale(j)                       &
                 - ssnow%evapfbl(j,1)*air%rlam(j)/dels)

            canopy%fess(j) = MIN(canopy%fess(j), REAL(fupper_limit(j),r_2))

            !Ticket 137 - case iii)
            !evaporation from frozen soils needs to respect the assumption that
            !ice fraction of soil moisture cannot exceed frozen_limit=0.85
            !see soilsnow: if frozen_limit changes need to be consistent
            fupper_limit(j) = REAL(ssnow%wb(j,1)-ssnow%wbice(j,1)/0.85)*frescale(j)
            fupper_limit(j) = MAX(REAL(fupper_limit(j),r_2),0.)

            canopy%fess(j) = MIN(canopy%fess(j), REAL(fupper_limit(j),r_2))

         END IF

         ssnow%cls(j)=1.

         !Ticket 137 - case ii) deposition of frost onto snow
         ! case of sublimation of snow overwrites later
         IF (ssnow%snowd(j) >=0.1 ) THEN
            ssnow%cls(j) = 1.1335
            canopy%fess(j) = ssnow%cls(j)*ssnow%potev(j)
         ENDIF

         !Ticket 137 - case iv) deposition of frost onto frozen soil, no snow
         IF (ssnow%snowd(j) < 0.1 .AND. ssnow%potev(j) < 0. .AND. &
              ssnow%tss(j)<C%TFRZ) THEN
            ssnow%cls(j)=1.1335
            canopy%fess(j) = ssnow%cls(j)*ssnow%potev(j)
         ENDIF

         !Ticket 137 - case ii) sublimation of snow
         IF (ssnow%snowd(j) >= 0.1 .AND. ssnow%potev(j) > 0.) THEN

            ssnow%cls(j) = 1.1335

            !INH - if changes to PM routine then matching changes here
            canopy%fess(j) = MIN( (ssnow%wetfac(j)*ssnow%potev(j))*ssnow%cls(j), &
                 ssnow%snowd(j)/dels*air%rlam(j)*ssnow%cls(j))

         ENDIF

      ENDDO

      ! Evaporation from soil puddle
      canopy%fesp = MIN(ssnow%pudsto/dels*air%rlam,MAX(pwet*ssnow%potev,0.))
      canopy%fes = canopy%fess + canopy%fesp

    END SUBROUTINE latent_heat_flux

    ! -----------------------------------------------------------------------------

    SUBROUTINE within_canopy( gbhu, gbhf, rt0, rhlitt, relitt )

      USE cable_def_types_mod, ONLY : mp, r_2

      REAL(r_2), INTENT(IN), DIMENSION(:,:) ::                                    &
           gbhu,    &  ! forcedConvectionBndryLayerCond
           gbhf        ! freeConvectionBndryLayerCond

      REAL, INTENT(IN), DIMENSION(mp) ::                                       &
           rt0,     &  ! resistance from ground to canopy air
           rhlitt,  &  ! additional litter resistance for heat (=0 by default)
           relitt      ! additional litter resistance for water

      REAL, DIMENSION(mp) ::                                                      &
           rrsw,             & ! recipr. stomatal resistance for water
           rrbw,             & ! recipr. leaf boundary layer resistance for water
           dmah,             & ! A_{H} in eq. 3.41 in SCAM, CSIRO tech report 132
           dmbh,             & ! B_{H} in eq. 3.41 in SCAM, CSIRO tech report 132
           dmch,             & ! C_{H} in eq. 3.41 in SCAM, CSIRO tech report 132
           dmae,             & ! A_{E} in eq. 3.41 in SCAM, CSIRO tech report 132
           dmbe,             & ! B_{E} in eq. 3.41 in SCAM, CSIRO tech report 132
           dmce                ! C_{E} in eq. 3.41 in SCAM, CSIRO tech report 132

      REAL  :: lower_limit, upper_limit
      REAL, DIMENSION(mp) :: fix_eqn,fix_eqn2

      INTEGER :: j

      !INH: rhlitt=relitt=0. if litter resistance not active but case included
      !dmah through to dmce are not A_{H} through C_{E} as per Eqn 3.40
      !in SCAM documentation but rt0*((1+esp)/rs + 1/rb)*A_{H} etc.
      !
      !changes from v1.4 for %cls package, litter and Or hydrology

      rrbw = SUM(gbhu+gbhf,2)/air%cmolar  ! MJT

      ! leaf stomatal resistance for water
      rrsw = SUM(canopy%gswx,2)/air%cmolar ! MJT

      IF (cable_user%or_evap) THEN
         fix_eqn(:) = rt0(:)*(REAL(ssnow%satfrac(:))/(rt0(:)+REAL(ssnow%rtevap_sat(:))) + &
              (1-REAL(ssnow%satfrac(:)))/(rt0(:)+REAL(ssnow%rtevap_unsat(:))))
         !lakes/ice rtevap=0 and wetfac is .ne. 1
         fix_eqn(:) = ssnow%wetfac(:) * fix_eqn(:)*ssnow%cls(:)   !INH correction. & M.Dekker +d wetfac

         fix_eqn2(:) = rt0(:) / (rt0(:) + REAL(ssnow%rt_qh_sublayer) )

      ELSE  !with INH corrections for litter and cls

         fix_eqn(:) = ssnow%cls(:)*rt0(:)/(rt0(:)+relitt(:))
         WHERE (ssnow%potev>0.) fix_eqn(:)=fix_eqn(:)*ssnow%wetfac(:)
         fix_eqn2(:) = rt0(:)/(rt0(:)+rhlitt(:))

      END IF

      DO j=1,mp

         IF(veg%meth(j) > 0 .AND. canopy%vlaiw(j) > C%LAI_THRESH .AND.              &
              rough%hruff(j) > rough%z0soilsn(j) ) THEN

            !   use the dispersion matrix (DM) to find the air temperature
            !   and specific humidity
            !   (Raupach, Finkele and Zhang 1997, pp 17)
            ! leaf boundary layer resistance for water
            ! A_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132
            dmah(j) = (rt0(j)+fix_eqn2(j)*rough%rt1(j))*((1.+air%epsi(j))*rrsw(j) + rrbw(j))  &
                 + air%epsi(j) * (rt0(j)*rough%rt1(j))*(rrbw(j)*rrsw(j))

            ! B_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132
            dmbh(j) = (-air%rlam(j)/C%CAPP)*(rt0(j)*rough%rt1(j))*(rrbw(j)*rrsw(j))

            ! C_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132
            dmch(j) = ((1.+air%epsi(j))*rrsw(j) + rrbw(j))*rt0(j)*rough%rt1(j)*   &
                 (canopy%fhv(j) + canopy%fhs(j))/(air%rho(j)*C%CAPP)

            ! A_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
            dmae(j) = (-air%epsi(j)*C%CAPP/air%rlam(j))*(rt0(j)*rough%rt1(j)) *   &
                 (rrbw(j)*rrsw(j))

            ! B_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
            dmbe(j) = ( rt0(j) + fix_eqn(j) * rough%rt1(j) ) *               &
                 ( (1.+air%epsi(j) ) * rrsw(j) + rrbw(j) ) +                 &
                 ( rt0(j) * rough%rt1(j) ) * ( rrbw(j) * rrsw(j) )

            ! C_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
            ! INH: includes modifications for %cls
            dmce(j) = ((1.+air%epsi(j))*rrsw(j) + rrbw(j))*rt0(j)*rough%rt1(j)*   &
                 (canopy%fev(j) + canopy%fes(j)/ssnow%cls(j)) /                   &
                 (air%rho(j)*air%rlam(j))

            ! Within canopy air temperature:
            met%tvair(j) = met%tk(j) + ( dmbe(j) * dmch(j) - dmbh(j) * dmce(j) )  &
                 / (dmah(j)*dmbe(j)-dmae(j)*dmbh(j)+1.0e-12)

            !---set limits for comparisson
            lower_limit =  MIN( ssnow%tss(j), met%tk(j)) - 5.0
            upper_limit =  MAX( ssnow%tss(j), met%tk(j)) + 5.0

            !--- tvair within these limits
            met%tvair(j) = MAX(met%tvair(j) , lower_limit)
            met%tvair(j) = MIN(met%tvair(j) , upper_limit)

            ! recalculate using canopy within temperature
            met%qvair(j) = met%qv(j) + (dmah(j)*dmce(j)-dmae(j)*dmch(j)) /        &
                 ( dmah(j)*dmbe(j)-dmae(j)*dmbh(j)+1.0e-12)
            met%qvair(j) = MAX(0.0,met%qvair(j))

            !---set limits for comparisson

            lower_limit =  MIN( ssnow%qstss(j), met%qv(j))
            upper_limit =  MAX( ssnow%qstss(j), met%qv(j))

            !--- qvair within these limits
            met%qvair(j) =  MAX(met%qvair(j),lower_limit)
            met%qvair(j) =  MIN(met%qvair(j), upper_limit)

            ! Saturated specific humidity in canopy:
            CALL qsatfjh2(qstvair(j),met%tvair(j)-C%tfrz,met%pmb(j))

            ! Saturated vapour pressure deficit in canopy:
            met%dva(j) = ( qstvair(j) - met%qvair(j) ) *  C%rmair/C%RMH2o         &
                 * met%pmb(j) * 100.
         ENDIF

      ENDDO

    END SUBROUTINE within_canopy

    ! -----------------------------------------------------------------------------

    SUBROUTINE update_zetar()

      INTEGER :: j

      ! monin-obukhov stability parameter zetar=zref/l
      ! recompute zetar for the next iteration, except on last iteration
      IF (iter < NITER) THEN ! dont compute zetar on the last iter

         iterplus = MAX(iter+1,2)
         canopy%zetar(:,iterplus) = -( C%VONK * C%GRAV * rough%zref_tq *              &
              ( canopy%fh + 0.07 * canopy%fe ) ) /          &
              ( air%rho * C%CAPP * met%tk * canopy%us**3 )

         ! stability parameter at shear height: needed for Harman in-canopy stability correction
         IF (cable_user%soil_struc=='sli') THEN
            WHERE (canopy%vlaiw > C%LAI_THRESH .AND. rough%hruff > rough%z0soilsn)
               canopy%zetash(:,iterplus) = -(C%vonk*C%grav*(0.1*rough%hruff)*(canopy%fhs+0.07*REAL(canopy%fes)))/ &
                    MAX( (air%rho*C%capp*met%tk*(canopy%us*rough%term6a)**3), 1.e-12)
            ELSEWHERE
               canopy%zetash(:,iterplus) = canopy%zetash(:,iter)
            ENDWHERE
         ENDIF

         ! case NITER=2: final zetar=C%ZETmul*zetar(2) (compute only when iter=1)
         IF (NITER == 2) THEN

            canopy%zetar(:,2) = C%ZETmul * canopy%zetar(:,2)

            ! stability parameter at shear height: needed for Harman in-canopy stability correction
            IF (cable_user%soil_struc=='sli') THEN
               canopy%zetash(:,2) = C%ZETmul * canopy%zetash(:,2)
            ENDIF


            DO j=1,mp
               IF ( (met%fsd(j,1)+met%fsd(j,2))  ==  0.0 ) &
                    canopy%zetar(j,2) = 0.5 * canopy%zetar(j,2)
            ENDDO

         END IF

         ! constrain zeta to C%ZETPOS and C%ZETNEG (set in param0)

         ! zetar too +
         canopy%zetar(:,iterplus) = MIN(C%ZETPOS,canopy%zetar(:,iterplus))
         !jhan: to get past rigorous build - however (:,i) cant be compared
         !if ( canopy%zetash(:,iterplus) .NE. C%ZETPOS ) &
         IF (cable_user%soil_struc=='sli') &
              canopy%zetash(:,iterplus) = MIN(C%ZETPOS,canopy%zetash(:,iterplus))

         ! zetar too -
         canopy%zetar(:,iterplus) = MAX(C%ZETNEG,canopy%zetar(:,iterplus))
         !if ( canopy%zetash(:,iterplus) .NE. C%ZETNEG ) &
         IF (cable_user%soil_struc=='sli') &
              canopy%zetash(:,iterplus) = MAX(C%ZETNEG,canopy%zetash(:,iterplus))

      END IF ! (iter < NITER)

    END SUBROUTINE update_zetar


    FUNCTION qsatf(j,tair,pmb) RESULT(r)
      ! MRR, 1987
      ! AT TEMP tair (DEG C) AND PRESSURE pmb (MB), GET SATURATION SPECIFIC
      ! HUMIDITY (KG/KG) FROM TETEN FORMULA

      REAL, INTENT(IN) ::                                                         &
           tair,         & ! air temperature (C)
           pmb             ! pressure PMB (mb)

      INTEGER, INTENT(IN) :: j

      REAL           :: r    ! result; sat sp humidity

      r = (C%RMH2o/C%rmair) * (C%TETENA*EXP(C%TETENB*tair/(C%TETENC+tair))) / pmb

    END FUNCTION qsatf

    ! -----------------------------------------------------------------------------

    SUBROUTINE qsatfjh(var,tair,pmb)
      USE cable_def_types_mod, ONLY : mp
      REAL, INTENT(IN), DIMENSION(mp) ::                                          &
           tair,                        & ! air temperature (C)
           pmb                            ! pressure PMB (mb)

      REAL, INTENT(OUT), DIMENSION(mp) ::                                         &
           var                            ! result; sat sp humidity

      INTEGER :: j

      DO j=1,mp

         var(j) = (C%RMH2o/C%rmair) * (C%TETENA*EXP(C%TETENB*tair(j)/(C%TETENC+tair(j))))    &
              / pmb(j)
      ENDDO

    END SUBROUTINE qsatfjh

    ! -----------------------------------------------------------------------------

    SUBROUTINE qsatfjh2(var,tair,pmb)

      REAL, INTENT(IN) ::                                                         &
           tair,         & ! air temperature (C)
           pmb             ! pressure PMB (mb)

      REAL, INTENT(OUT) ::                                                        &
           var             ! result; sat sp humidity

      var = (C%RMH2o/C%rmair) * (C%TETENA*EXP(C%TETENB*tair/(C%TETENC+tair))) / pmb

    END SUBROUTINE qsatfjh2

    ! -----------------------------------------------------------------------------

    FUNCTION psim(zeta) RESULT(r)
      USE cable_def_types_mod, ONLY : mp
      ! mrr, 16-sep-92 (from function psi: mrr, edinburgh 1977)
      ! computes integrated stability function psim(z/l) (z/l=zeta)
      ! for momentum, using the businger-dyer form for unstable cases
      ! and the Beljaars and Holtslag (1991) form for stable cases.


      REAL, INTENT(IN), DIMENSION(mp) ::  zeta       !

      ! function result
      REAL, DIMENSION(mp) :: r

      REAL, DIMENSION(mp) ::                                                      &
           x,       & !
           z,       & !
           stable,  & !
           unstable   !

      REAL, PARAMETER ::                                                          &
           gu = 16.0,  & !
           gs = 5.0

      REAL, PARAMETER ::                                                          &
           a = 1.0,       & !
           b = 0.667,     & !
           xc = 5.0,       & !
           d = 0.35         !

      z = 0.5 + SIGN(0.5,zeta)    ! z=1 in stable, 0 in unstable

      ! Beljaars and Holtslag (1991) for stable
      stable = -a*zeta - b*(zeta - xc/d)*EXP( -d*zeta) - b*xc/d
      x      = (1.0 + gu*ABS(zeta))**0.25
      unstable = ALOG((1.0+x*x)*(1.0+x)**2/8) - 2.0*ATAN(x) + C%PI_C*0.5
      r = z*stable + (1.0-z)*unstable
    END FUNCTION psim

    ! -----------------------------------------------------------------------------

    ELEMENTAL FUNCTION psis(zeta) RESULT(r)

      ! mrr, 16-sep-92 (from function psi: mrr, edinburgh 1977)
      ! computes integrated stability function psis(z/l) (z/l=zeta)
      ! for scalars, using the businger-dyer form for unstable cases
      ! and the webb form for stable cases. see paulson (1970).

      REAL, INTENT(IN)     :: zeta

      REAL, PARAMETER      ::                                                     &
           gu = 16.0,        & !
           gs = 5.0,         & !
           a = 1.0,          & !
           b = 0.667,        & !
           c = 5.0,          & !
           d = 0.35

      REAL                 ::                                                     &
           r,                & !
           stzeta,           & !
           ustzeta,          & !
           z,                & !
           y,                & !
           stable,           & !
           unstable

      z      = 0.5 + SIGN(0.5,zeta)    ! z=1 in stable, 0 in unstable

      ! Beljaars and Holtslag (1991) for stable
      stzeta = MAX(0.,zeta)
      stable = -(1.+2./3.*a*stzeta)**(3./2.) -  &
           b*(stzeta-c/d)*EXP(-d*stzeta) - b*c/d + 1.
      y      = (1.0 + gu*ABS(zeta))**0.5
      unstable = 2.0 * alog((1+y)*0.5)
      r   = z*stable + (1.0-z)*unstable

    END FUNCTION psis

    ! -----------------------------------------------------------------------------

    ELEMENTAL FUNCTION rplant(rpconst, rpcoef, tair) RESULT(z)
      REAL, INTENT(IN)     :: rpconst
      REAL, INTENT(IN)     :: rpcoef
      REAL, INTENT(IN)     :: tair
      REAL                 :: z
      z = rpconst * EXP(rpcoef * tair)
    END FUNCTION rplant

    ! -----------------------------------------------------------------------------

    SUBROUTINE wetLeaf( dels, rad, rough, air, met, veg, canopy, cansat, tlfy,     &
         gbhu, gbhf, ghwet )

      USE cable_def_types_mod

      TYPE (radiation_type), INTENT(INOUT) :: rad
      TYPE (roughness_type), INTENT(INOUT) :: rough
      TYPE (air_type),       INTENT(INOUT) :: air
      TYPE (met_type),       INTENT(INOUT) :: met
      TYPE (canopy_type),    INTENT(INOUT) :: canopy

      TYPE (veg_parameter_type), INTENT(INOUT)    :: veg

      REAL,INTENT(IN), DIMENSION(:) ::                                            &
           tlfy,          & ! leaf temp (K) - assC%UMINg the temperature of
                                ! wet leaf is equal that of dry leaf ="tlfy"
           cansat           ! max canopy intercept. (mm)

      REAL(r_2), INTENT(IN), DIMENSION(:,:) ::                                    &
           gbhu,          & ! forcedConvectionBndryLayerCond
           gbhf             ! freeConvectionBndryLayerCond

      REAL(r_2), INTENT(OUT), DIMENSION(:) ::                                     &
           ghwet            ! cond for heat for a wet canopy

      REAL, INTENT(IN)     :: dels ! integration time step (s)

      ! local variables
      REAL, DIMENSION(mp) ::                                                      &
           ccfevw,        & ! limitation term for
           gwwet,         & ! cond for water for a wet canopy
           ghrwet           ! wet canopy cond: heat & thermal rad

      !i sums, terms of convenience/readability
      REAL, DIMENSION(mp) ::                                                      &
           sum_gbh, sum_rad_rniso, sum_rad_gradis, xx1

      INTEGER :: j

      ! END header

      ghwet = 1.0e-3
      gwwet = 1.0e-3
      ghrwet= 1.0e-3
      canopy%fevw = 0.0
      canopy%fhvw = 0.0
      sum_gbh = SUM((gbhu+gbhf),2)
      sum_rad_rniso = SUM(rad%rniso,2)
      sum_rad_gradis = SUM(rad%gradis,2)

      DO j=1,mp

         IF(canopy%vlaiw(j) > C%LAI_THRESH) THEN

            ! VEG SENSIBLE & LATENT HEAT FLUXES fevw, fhvw (W/m2) for a wet canopy
            ! calculate total thermal resistance, rthv in s/m
            ghwet(j) = 2.0   * sum_gbh(j)
            gwwet(j) = 1.075 * sum_gbh(j)
            ghrwet(j) = sum_rad_gradis(j) + ghwet(j)

            ! Calculate fraction of canopy which is wet:
            canopy%fwet(j) = MAX( 0.0, MIN( 1.0,                                  &
                 0.8 * canopy%cansto(j) / MAX( cansat(j), 0.01 ) ) )

            ! Calculate lat heat from wet canopy, may be neg. if dew on wet canopy
            ! to avoid excessive evaporation:
            ccfevw(j) = MIN(canopy%cansto(j) * air%rlam(j) / dels, &
                 2.0 / (1440.0 / (dels/60.0)) * air%rlam(j) )

            canopy%fevw(j) = MIN( canopy%fwet(j) * ( air%dsatdk(j) *              &
                 ( sum_rad_rniso(j)- C%CAPP*C%rmair*( met%tvair(j)     &
                 - met%tk(j) ) * sum_rad_gradis(j) )                   &
                 + C%CAPP * C%rmair * met%dva(j) * ghrwet(j) )         &
                 / ( air%dsatdk(j)+air%psyc(j)*ghrwet(j) / gwwet(j) )  &
                 , ccfevw(j) )

            canopy%fevw_pot(j) = ( air%dsatdk(j)* (sum_rad_rniso(j) -             &
                 C%CAPP * C%rmair * ( met%tvair(j) - met%tk(j) )  &
                 *sum_rad_gradis(j) )                             &
                 + C%CAPP * C%rmair * met%dva(j) * ghrwet(j))     &
                 / (air%dsatdk(j)+air%psyc(j)*ghrwet(j)/gwwet(j) )

            ! calculate sens heat from wet canopy:
            canopy%fhvw(j) = canopy%fwet(j) * ( sum_rad_rniso(j) -C%CAPP * C%rmair&
                 * ( tlfy(j) - met%tk(j) ) * sum_rad_gradis(j) )      &
                 - canopy%fevw(j)

         ENDIF

      ENDDO

    END SUBROUTINE wetLeaf

    ! -----------------------------------------------------------------------------

  END SUBROUTINE define_canopy

  ! -----------------------------------------------------------------------------

  SUBROUTINE Surf_wetness_fact( cansat, canopy, ssnow,veg, met, soil, dels )

    USE cable_common_module
    USE cable_def_types_mod
    USE cable_gw_hydro_module, ONLY : calc_srf_wet_fraction

    TYPE (veg_parameter_type), INTENT(INOUT)    :: veg
    TYPE (soil_snow_type), INTENT(inout):: ssnow
    TYPE (soil_parameter_type), INTENT(inout)   :: soil
    TYPE (canopy_type), INTENT(INOUT)   :: canopy
    TYPE (met_type), INTENT(INOUT)   :: met

    REAL, INTENT(IN) :: dels ! integration time setp (s)

    REAL,INTENT(IN), DIMENSION(:) :: cansat ! max canopy intercept. (mm)

    !local variables
    REAL, DIMENSION(mp)  :: lower_limit, upper_limit,ftemp

    INTEGER :: j

    ! Rainfall variable is limited so canopy interception is limited,
    ! used to stabilise latent fluxes.
    ! to avoid excessive direct canopy evaporation (EK nov2007, snow scheme)
    upper_limit = 4.0 * MIN(dels,1800.0) / (60.0 * 1440.0 )
    ftemp =MIN(met%precip-met%precip_sn, upper_limit )
    ! Calculate canopy intercepted rainfall, equal to zero if temp < 0C:
    lower_limit = cansat - canopy%cansto
    upper_limit = MAX(lower_limit, 0.0)
    canopy%wcint = MERGE( MIN( upper_limit, ftemp ), 0.0,                       &
         ftemp > 0.0  .AND. met%tk > C%tfrz)  !EAK, 09/10

    ! Define canopy throughfall (100% of precip if temp < 0C, see above):
    canopy%through = met%precip_sn + MIN( met%precip - met%precip_sn ,          &
         MAX( 0.0, met%precip - met%precip_sn - canopy%wcint) )

    ! Add canopy interception to canopy storage term:
    canopy%cansto = canopy%cansto + canopy%wcint

    ! Calculate fraction of canopy which is wet:
    canopy%fwet   = MAX( 0.0, MIN( 0.9, 0.8 * canopy%cansto /                   &
         MAX( cansat, 0.01 ) ) )

    !calc the surface wetness for soil evap in this routine
    !include the default wetfac when or_evap and gw_model are not used
    CALL calc_srf_wet_fraction(ssnow,soil,met,veg)

  END SUBROUTINE Surf_wetness_fact

  ! -----------------------------------------------------------------------------
  SUBROUTINE dryLeaf( dels, rad, rough, air, met,                                &
       veg, canopy, soil, ssnow, dsx,                             &
       fwsoil, tlfx,  tlfy,  ecy, hcy,                            &
       rny, gbhu, gbhf, csx,                                      &
       cansat, ghwet, iter,climate )

    USE cable_def_types_mod
    USE cable_common_module

    TYPE (radiation_type), INTENT(INOUT) :: rad
    TYPE (roughness_type), INTENT(INOUT) :: rough
    TYPE (air_type),       INTENT(INOUT) :: air
    TYPE (met_type),       INTENT(INOUT) :: met
    TYPE (canopy_type),    INTENT(INOUT) :: canopy
    TYPE (soil_snow_type), INTENT(INOUT) :: ssnow

    TYPE (veg_parameter_type),  INTENT(INOUT)   :: veg
    TYPE (soil_parameter_type), INTENT(inout)   :: soil
    TYPE (climate_type), INTENT(IN)    :: climate

    REAL, INTENT(INOUT), DIMENSION(:) ::                                        &
         dsx,        & ! leaf surface vpd
         fwsoil,     & ! soil water modifier of stom. cond
         tlfx,       & ! leaf temp prev. iter (K)
         tlfy          ! leaf temp (K)

    REAL(R_2),INTENT(INOUT), DIMENSION(:) ::                                    &
         ecy,        & ! lat heat fl dry big leaf
         hcy,        & ! veg. sens heat
         rny         !& !

    REAL(R_2),INTENT(INOUT), DIMENSION(:,:) ::                                  &
         gbhu,       & ! forcedConvectionBndryLayerCond
         gbhf,       & ! freeConvectionBndryLayerCond
         csx           ! leaf surface CO2 concentration

    REAL,INTENT(IN), DIMENSION(:) :: cansat

    REAL(r_2), INTENT(OUT), DIMENSION(:) ::                                     &
         ghwet  ! cond for heat for a wet canopy

    INTEGER,INTENT(IN) :: iter

    REAL, INTENT(IN)     :: dels ! integration time step (s)

    !local variables
    REAL, PARAMETER  ::                                                         &
         co2cp3 = 0.0,  & ! CO2 compensation pt C3
         jtomol = 4.6e-6  ! Convert from J to Mol for light

    REAL, DIMENSION(mp) ::                                                      &
         conkct,        & ! Michaelis Menton const.
         conkot,        & ! Michaelis Menton const.
         cx1,           & ! "d_{3}" in Wang and Leuning,
         cx2,           & !     1998, appendix E
         tdiff,         & ! leaf air temp diff.
         tlfxx,         & ! leaf temp of current iteration (K)
         abs_deltlf,    & ! ABS(deltlf)
         deltlf,        & ! deltlfy of prev iter.
         deltlfy,       & ! del temp successive iter.
         gras,          & ! Grashof coeff
         evapfb,        & !
         sum_rad_rniso, & !
         sum_rad_gradis,& !
         gwwet,         & ! cond for water for a wet canopy
         ghrwet,        & ! wet canopy cond: heat & thermal rad
         sum_gbh,       & !
         ccfevw,        & ! limitation term for
                                ! wet canopy evaporation rate
         temp             !

    REAL(r_2), DIMENSION(mp)  ::                                                &
         ecx,        & ! lat. hflux big leaf
         ecx_t,      & ! lat. hflux big leaf
         hcx,        & ! sens heat fl big leaf prev iteration
         rnx,        & ! net rad prev timestep
         fwsoil_coef   !

    REAL, DIMENSION(mp,ms)  :: oldevapfbl

    REAL, DIMENSION(mp,mf)  ::                                                  &
         gw,         & ! cond for water for a dry canopy
         gh,         & ! cond for heat for a dry canopy
         ghr,        & ! dry canopy cond for heat & thermal rad
         anx,        & ! net photos. prev iteration
         an_y,       & ! net photosynthesis soln
         rdx,        & ! daytime leaf resp rate, prev iteration
         rdy,        & ! daytime leaf resp rate
         ejmax2,     & ! jmax of big leaf
         ejmxt3,     & ! jmax big leaf C3 plants
         vcmxt3,     & ! vcmax big leaf C3
         vcmxt4,     & ! vcmax big leaf C4
         vx3,        & ! carboxylation C3 plants
         vx4,        & ! carboxylation C4 plants
                                ! Ticket #56, xleuning is no longer used, we replace it with
                                ! gs_coeff,
                                ! which is computed differently based on the new GS_SWITCH. If GS_SWITCH
                                ! is "leuning", it's the same, if "medlyn", then the new Medlyn model
                                ! xleuning,   & ! leuning stomatal coeff
         gs_coeff,   & ! stom coeff, Ticket #56
         psycst,     & ! modified pych. constant
         frac42,     & ! 2D frac4
         temp2

    REAL, DIMENSION(:,:), POINTER :: gswmin ! min stomatal conductance

    REAL, DIMENSION(mp,2) ::  gsw_term, lower_limit2  ! local temp var

    REAL :: trans_mmol, conv, e_test
    REAL, PARAMETER :: KG_2_G = 1000.0
    REAL, PARAMETER :: G_WATER_TO_MOL = 1.0 / 18.01528
    REAL, PARAMETER :: MOL_2_MMOL = 1000.0
    REAL, PARAMETER :: MB_TO_PA = 100.

    INTEGER :: i, j, k, kk  ! iteration count
    REAL :: vpd, g1, ktot, fw, refill ! Ticket #56
#define VanessasCanopy
#ifdef VanessasCanopy
    REAL, DIMENSION(mp,mf)  ::                                                  &
         xleuning    ! leuning stomatal coeff
#endif

    REAL :: medlyn_lim  !INH 2018: should be a parameter in long-term
    ! END header

    ALLOCATE( gswmin(mp,mf ))

    ! Soil water limitation on stomatal conductance:
    IF( iter ==1) THEN
       IF ((cable_user%soil_struc=='default').and.(cable_user%FWSOIL_SWITCH.ne.'Haverd2013')) THEN
          IF(cable_user%FWSOIL_SWITCH == 'standard') THEN
             CALL fwsoil_calc_std( fwsoil, soil, ssnow, veg)
          ELSEIf (cable_user%FWSOIL_SWITCH == 'non-linear extrapolation') THEN
             !EAK, 09/10 - replace linear approx by polynomial fitting
             CALL fwsoil_calc_non_linear(fwsoil, soil, ssnow, veg)
          ELSEIF(cable_user%FWSOIL_SWITCH == 'Lai and Ktaul 2000') THEN
             CALL fwsoil_calc_Lai_Ktaul(fwsoil, soil, ssnow, veg)
          ! Water stress is inferred from the hydraulics approach instead.
          ELSE IF(cable_user%FWSOIL_SWITCH == 'hydraulics') THEN
             fwsoil = 1.0
          ELSE
             STOP 'fwsoil_switch failed.'
          ENDIF
          canopy%fwsoil = fwsoil
       ELSEIF ((cable_user%soil_struc=='sli').OR.(cable_user%FWSOIL_SWITCH=='Haverd2013')) THEN
          fwsoil = canopy%fwsoil
       ENDIF



    ENDIF

    ! weight min stomatal conductance by C3 an C4 plant fractions
    frac42 = SPREAD(veg%frac4, 2, mf) ! frac C4 plants
    gsw_term = SPREAD(veg%gswmin,2,mf)
    lower_limit2 = rad%scalex * gsw_term
    gswmin = MAX(1.e-6,lower_limit2)


    gw = 1.0e-3 ! default values of conductance
    gh = 1.0e-3
    ghr= 1.0e-3
    rdx = 0.0
    anx = 0.0
    rnx = SUM(rad%rniso,2)
    abs_deltlf = 999.0


    gras = 1.0e-6
    an_y= 0.0
    hcx = 0.0              ! init sens heat iteration memory variable
    hcy = 0.0
    rdy = 0.0
    ecx = SUM(rad%rniso,2) ! init lat heat iteration memory variable
    tlfxx = tlfx
    psycst(:,:) = SPREAD(air%psyc,2,mf)
    canopy%fevc = 0.0
    ssnow%evapfbl = 0.0

    ghwet = 1.0e-3
    gwwet = 1.0e-3
    ghrwet= 1.0e-3
    canopy%fevw = 0.0
    canopy%fhvw = 0.0
    sum_gbh = SUM((gbhu+gbhf),2)
    sum_rad_rniso = SUM(rad%rniso,2)
    sum_rad_gradis = SUM(rad%gradis,2)

    DO kk=1,mp

       IF(canopy%vlaiw(kk) <= C%LAI_THRESH) THEN
          rnx(kk) = 0.0 ! intialise
          ecx(kk) = 0.0 ! intialise
          ecy(kk) = ecx(kk) ! store initial values
          abs_deltlf(kk)=0.0
          rny(kk) = rnx(kk) ! store initial values
          ! calculate total thermal resistance, rthv in s/m
       END IF

    ENDDO

    deltlfy = abs_deltlf
    k = 0


    !kdcorbin, 08/10 - doing all points all the time
    DO WHILE (k < C%MAXITER)
       k = k + 1
       DO i=1,mp



          IF (canopy%vlaiw(i) > C%LAI_THRESH .AND. abs_deltlf(i) > 0.1) THEN

             ghwet(i) = 2.0   * sum_gbh(i)
             gwwet(i) = 1.075 * sum_gbh(i)
             ghrwet(i) = sum_rad_gradis(i) + ghwet(i)

             ! Calculate fraction of canopy which is wet:
             canopy%fwet(i) = MAX( 0.0, MIN( 1.0, 0.8 * canopy%cansto(i)/ MAX(  &
                  cansat(i),0.01 ) ) )

             ! Calculate lat heat from wet canopy, may be neg.
             ! if dew on wet canopy to avoid excessive evaporation:
             ccfevw(i) = MIN(canopy%cansto(i) * air%rlam(i) / dels,             &
                  2.0 / (1440.0 / (dels/60.0)) * air%rlam(i) )

             ! Grashof number (Leuning et al, 1995) eq E4:
             gras(i) = MAX(1.0e-6,                                              &
                  1.595E8* ABS( tlfx(i)-met%tvair(i))* (veg%dleaf(i)**3.0) )

             ! See Appendix E in (Leuning et al, 1995):
             gbhf(i,1) = rad%fvlai(i,1) * air%cmolar(i) * 0.5*C%dheat           &
                  * ( gras(i)**0.25 ) / veg%dleaf(i)
             gbhf(i,2) = rad%fvlai(i,2) * air%cmolar(i) * 0.5 * C%dheat         &
                  * ( gras(i)**0.25 ) / veg%dleaf(i)
             gbhf(i,:) = MAX( 1.e-6_r_2, gbhf(i,:) )

             ! Conductance for heat:
             gh(i,:) = 2.0 * (gbhu(i,:) + gbhf(i,:))

             ! Conductance for heat and longwave radiation:
             ghr(i,:) = rad%gradis(i,:)+gh(i,:)

             ! Leuning 2002 (P C & E) equation for temperature response
             ! used for Vcmax for C3 plants:
             temp(i) =  xvcmxt3(tlfx(i)) * veg%vcmax(i) * (1.0-veg%frac4(i))

             vcmxt3(i,1) = rad%scalex(i,1) * temp(i)
             vcmxt3(i,2) = rad%scalex(i,2) * temp(i)

             ! Temperature response Vcmax, C4 plants (Collatz et al 1989):
             temp(i) = xvcmxt4(tlfx(i)-C%tfrz) * veg%vcmax(i) * veg%frac4(i)
             vcmxt4(i,1) = rad%scalex(i,1) * temp(i)
             vcmxt4(i,2) = rad%scalex(i,2) * temp(i)

             ! Leuning 2002 (P C & E) equation for temperature response
             ! used for Jmax for C3 plants:
             temp(i) = xejmxt3(tlfx(i)) * veg%ejmax(i) * (1.0-veg%frac4(i))
             ejmxt3(i,1) = rad%scalex(i,1) * temp(i)
             ejmxt3(i,2) = rad%scalex(i,2) * temp(i)

             ! Difference between leaf temperature and reference temperature:
             tdiff(i) = tlfx(i) - C%TREFK

             ! Michaelis menten constant of Rubisco for CO2:
             conkct(i) = veg%conkc0(i) * EXP( ( veg%ekc(i) / (C%rgas*C%trefk) ) &
                  * ( 1.0 - C%trefk/tlfx(i) ) )

             ! Michaelis menten constant of Rubisco for oxygen:
             conkot(i) = veg%conko0(i) * EXP( ( veg%eko(i) / (C%rgas*C%trefk) ) &
                  * ( 1.0 - C%trefk/tlfx(i) ) )

             ! Store leaf temperature
             tlfxx(i) = tlfx(i)

             ! "d_{3}" in Wang and Leuning, 1998, appendix E:
             cx1(i) = conkct(i) * (1.0+0.21/conkot(i))
             cx2(i) = 2.0 * C%gam0 * ( 1.0 + C%gam1 * tdiff(i)                  &
                  + C%gam2 * tdiff(i) * tdiff(i) )

             ! All equations below in appendix E in Wang and Leuning 1998 are
             ! for calculating anx, csx and gswx for Rubisco limited,
             ! RuBP limited, sink limited
             temp2(i,1) = rad%qcan(i,1,1) * jtomol * (1.0-veg%frac4(i))
             temp2(i,2) = rad%qcan(i,2,1) * jtomol * (1.0-veg%frac4(i))
             vx3(i,1)  = ej3x(temp2(i,1),veg%alpha(i),veg%convex(i),ejmxt3(i,1))
             vx3(i,2)  = ej3x(temp2(i,2),veg%alpha(i),veg%convex(i),ejmxt3(i,2))
             temp2(i,1) = rad%qcan(i,1,1) * jtomol * veg%frac4(i)
             temp2(i,2) = rad%qcan(i,2,1) * jtomol * veg%frac4(i)
             vx4(i,1)  = ej4x(temp2(i,1),veg%alpha(i),veg%convex(i),vcmxt4(i,1))
             vx4(i,2)  = ej4x(temp2(i,2),veg%alpha(i),veg%convex(i),vcmxt4(i,2))

             rdx(i,1) = (veg%cfrd(i)*Vcmxt3(i,1) + veg%cfrd(i)*vcmxt4(i,1))
             rdx(i,2) = (veg%cfrd(i)*vcmxt3(i,2) + veg%cfrd(i)*vcmxt4(i,2))

             !Vanessa - the trunk does not contain xleauning as of Ticket#56 inclusion
             !as well as other inconsistencies here that need further investigation. In the
             !interests of getting this into the trunk ASAP just isolate this code for now
             !default side of this condition is to use trunk version

             !#ifdef VanessasCanopy


             IF (cable_user%CALL_climate) THEN

                ! Atkins et al. 2015, Table S4,
                ! modified by saling factor to reduce leaf respiration to
                ! expected proportion of GPP
                !Broad-leaved trees: Rdark a25 =
                !1.2818 + (0.0116 × Vcmax,a25) – (0.0334 × TWQ)
                !C3 herbs/grasses: Rdark,a25 =
                !1.6737 + (0.0116 × Vcmax,a25) – (0.0334 × TWQ)
                !Needle-leaved trees: Rdark,a25 =
                !1.2877 + (0.0116 × Vcmax,a25) – (0.0334 × TWQ)
                !Shrubs: Rdark,a25 = 1.5758 + (0.0116 × Vcmax,a25) – (0.0334 × TWQ)

                IF (veg%iveg(i).EQ.2 .OR. veg%iveg(i).EQ. 4  ) THEN ! broadleaf forest

                   rdx(i,1) = 0.60*(1.2818e-6+0.0116*veg%vcmax(i)- &
                        0.0334*climate%qtemp_max_last_year(i)*1e-6)
                   rdx(i,2) = rdx(i,1)

                ELSEIF (veg%iveg(i).EQ.1 .OR. veg%iveg(i).EQ. 3  ) THEN ! needleleaf forest
                   rdx(i,1) = 1.0*(1.2877e-6+0.0116*veg%vcmax(i)- &
                        0.0334*climate%qtemp_max_last_year(i)*1e-6)
                   rdx(i,2) = rdx(i,1)

                ELSEIF (veg%iveg(i).EQ.6 .OR. veg%iveg(i).EQ.8 .OR. &
                     veg%iveg(i).EQ. 9  ) THEN ! C3 grass, tundra, crop
                   rdx(i,1) = 0.60*(1.6737e-6+0.0116*veg%vcmax(i)- &
                        0.0334*climate%qtemp_max_last_year(i)*1e-6)
                   rdx(i,2) = rdx(i,1)

                ELSE  ! shrubs and other (C4 grass and crop)
                   rdx(i,1) = 0.60*(1.5758e-6+0.0116*veg%vcmax(i)- &
                        0.0334*climate%qtemp_max_last_year(i)*1e-6)
                   rdx(i,2) = rdx(i,1)
                ENDIF


                ! modify for leaf area and instanteous temperature response (Rd25 -> Rd)
                rdx(i,1) = rdx(i,1) * xrdt(tlfx(i)) * rad%scalex(i,1)
                rdx(i,2) = rdx(i,2) * xrdt(tlfx(i)) * rad%scalex(i,2)



                ! reduction of daytime leaf dark-respiration to account for
                !photo-inhibition
                !Mercado, L. M., Huntingford, C., Gash, J. H. C., Cox, P. M.,
                ! and Jogireddy, V.:
                ! Improving the representation of radiation
                !interception and photosynthesis for climate model applications,
                !Tellus B, 59, 553-565, 2007.
                ! Equation 3
                ! (Brooks and Farquhar, 1985, as implemented by Lloyd et al., 1995).
                ! Rc = Rd 0 < Io < 10 μmol quantam−2s−1
                ! Rc = [0.5 − 0.05 ln(Io)] Rd Io > 10μmol quantam−2s−1

                IF (jtomol*1.0e6*rad%qcan(i,1,1).GT.10.0) &
                     rdx(i,1) = rdx(i,1) * &
                     (0.5 - 0.05*LOG(jtomol*1.0e6*rad%qcan(i,1,1)))

                IF (jtomol*1.0e6*rad%qcan(i,1,2).GT.10.0) &
                     rdx(i,2) = rdx(i,2) * &
                     (0.5 - 0.05*LOG(jtomol*1.0e6*rad%qcan(i,1,2)))

!!$                xleuning(i,1) = ( fwsoil(i) / ( csx(i,1) - co2cp3 ) )              &
!!$                     * ( veg%a1gs(i) / ( 1.0 + dsx(i)/veg%d0gs(i)))
!!$                xleuning(i,2) = ( fwsoil(i) / ( csx(i,2) - co2cp3 ) )              &
!!$                     * ( veg%a1gs(i) / ( 1.0 + dsx(i)/veg%d0gs(i)))

             ELSE !cable_user%call_climate

!!$!Vanessa:note there is no xleuning to go into photosynthesis etc anymore
!!$             gs_coeff = xleuning

                !#else
                rdx(i,1) = (veg%cfrd(i)*vcmxt3(i,1) + veg%cfrd(i)*vcmxt4(i,1))
                rdx(i,2) = (veg%cfrd(i)*vcmxt3(i,2) + veg%cfrd(i)*vcmxt4(i,2))

             ENDIF !cable_user%call_climate

             ! Ticket #56 added switch for Belinda Medlyn's model
             IF (cable_user%GS_SWITCH == 'leuning') THEN
                gs_coeff(i,1) = ( fwsoil(i) / ( csx(i,1) - co2cp3 ) )              &
                     * ( veg%a1gs(i) / ( 1.0 + dsx(i)/veg%d0gs(i)))

                gs_coeff(i,2) = ( fwsoil(i) / ( csx(i,2) - co2cp3 ) )              &
                     * ( veg%a1gs(i) / ( 1.0 + dsx(i)/veg%d0gs(i)))

             ! Medlyn BE et al (2011) Global Change Biology 17: 2134-2144.
             ELSEIF(cable_user%GS_SWITCH == 'medlyn' .AND. &
                    cable_user%FWSOIL_SWITCH /= 'hydraulics') THEN

                gswmin = veg%g0(i)

                IF (dsx(i) < 50.0) THEN
                   vpd  = 0.05 ! kPa
                ELSE
                   vpd = dsx(i) * 1E-03 ! Pa -> kPa
                END IF

                g1 = veg%g1(i)

                gs_coeff(i,1) = (1.0 + (g1 * fwsoil(i)) / SQRT(vpd)) / csx(i,1)
                gs_coeff(i,2) = (1.0 + (g1 * fwsoil(i)) / SQRT(vpd)) / csx(i,2)

                !INH 2018: enforce gs_coeff to vary proportionally to fwsoil in dry soil conditions
                ! required to avoid transpiration without soil water extraction
                medlyn_lim = 0.05
                IF (fwsoil(i) <= medlyn_lim) THEN
                   gs_coeff(i,1) = (fwsoil(i) / medlyn_lim + (g1 * fwsoil(i)) / SQRT(vpd)) / csx(i,1)
                   gs_coeff(i,2) = (fwsoil(i) / medlyn_lim + (g1 * fwsoil(i)) / SQRT(vpd)) / csx(i,2)
                END IF

             ELSE IF (cable_user%GS_SWITCH == 'medlyn' .AND. &
                     cable_user%FWSOIL_SWITCH == 'hydraulics') THEN

                ! need to use gmin ...
                gswmin = 0.0

                CALL calc_hydr_conduc(canopy, ssnow, rad, veg, veg%kp_sat(i), i)

                ! Sensitivity of stomata to leaf water potential [0-1]
                fw = f_tuzet(canopy%psi_leaf_prev(i), veg%sf(i), veg%psi_f(i))

                ! convert to conductance to CO2
                gs_coeff(i,1) = (veg%g1(i) / csx(i,1) * fw) / C%RGSWC
                gs_coeff(i,2) = (veg%g1(i) / csx(i,2) * fw) / C%RGSWC

             ELSE
               print*, cable_user%GS_SWITCH, cable_user%FWSOIL_SWITCH, veg%iveg(i)
               STOP 'gs_model_switch failed.'
            ENDIF ! IF (cable_user%GS_SWITCH == 'leuning') THEN
            !#endif

          ENDIF !IF (canopy%vlaiw(i) > C%LAI_THRESH .AND. abs_deltlf(i) > 0.1)

       ENDDO !i=1,mp

       CALL photosynthesis( csx(:,:),                                           &
            SPREAD( cx1(:), 2, mf ),                            &
            SPREAD( cx2(:), 2, mf ),                            &
            gswmin(:,:), rdx(:,:), vcmxt3(:,:),                 &
            vcmxt4(:,:), vx3(:,:), vx4(:,:),                    &
                                ! Ticket #56, xleuning replaced with gs_coeff here
            gs_coeff(:,:), rad%fvlai(:,:),&
            SPREAD( abs_deltlf, 2, mf ),                        &
            anx(:,:), fwsoil(:) )

       DO i=1,mp

          IF (canopy%vlaiw(i) > C%LAI_THRESH .AND. abs_deltlf(i) > 0.1 ) THEN

             DO kk=1,mf

                IF(rad%fvlai(i,kk)>C%LAI_THRESH) THEN

                   csx(i,kk) = met%ca(i) - C%RGBWC*anx(i,kk) / (                &
                        gbhu(i,kk) + gbhf(i,kk) )
                   csx(i,kk) = MAX( 1.0e-4_r_2, csx(i,kk) )


                   ! Ticket #56, xleuning replaced with gs_coeff here
                   canopy%gswx(i,kk) = MAX( 1.e-3, gswmin(i,kk)*fwsoil(i) +     &
                        MAX( 0.0, C%RGSWC * gs_coeff(i,kk) *     &
                        anx(i,kk) ) )


                   !Recalculate conductance for water:
                   gw(i,kk) = 1.0 / ( 1.0 / canopy%gswx(i,kk) +                 &
                        1.0 / ( 1.075 * ( gbhu(i,kk) + gbhf(i,kk) ) ) )



                   gw(i,kk) = MAX( gw(i,kk), 0.00001 )

                   ! Modified psychrometric constant
                   ! (Monteith and Unsworth, 1990)
                   psycst(i,kk) = air%psyc(i) * REAL( ghr(i,kk) / gw(i,kk) )

                ENDIF

             ENDDO

             ecx(i) = ( air%dsatdk(i) * ( rad%rniso(i,1) - C%capp * C%rmair     &
                  * ( met%tvair(i) - met%tk(i) ) * rad%gradis(i,1) )        &
                  + C%capp * C%rmair * met%dva(i) * ghr(i,1) )              &
                  / ( air%dsatdk(i) + psycst(i,1) ) + ( air%dsatdk(i)       &
                  * ( rad%rniso(i,2) - C%capp * C%rmair * ( met%tvair(i) -  &
                  met%tk(i) ) * rad%gradis(i,2) ) + C%capp * C%rmair *      &
                  met%dva(i) * ghr(i,2) ) /                                 &
                  ( air%dsatdk(i) + psycst(i,2) )

             IF (cable_user%fwsoil_switch=='Haverd2013') THEN
                ! avoid root-water extraction when fwsoil is zero
                IF (fwsoil(i).LT.1e-6) THEN
                   anx(i,:) =  - rdx(i,:)
                   ecx(i) = 0.0
                ENDIF

                canopy%fevc(i) = ecx(i)*(1.0-canopy%fwet(i))

                CALL getrex_1d(ssnow%wbliq(i,:),&
                     ssnow%rex(i,:), &
                     canopy%fwsoil(i), &
                     REAL(veg%froot(i,:),r_2),&
                     soil%ssat_vec(i,:), &
                     soil%swilt_vec(i,:), &
                     MAX(REAL(canopy%fevc(i)/air%rlam(i)/1000_r_2,r_2),0.0_r_2), &
                     REAL(veg%gamma(i),r_2), &
                     REAL(soil%zse,r_2), REAL(dels,r_2), REAL(veg%zr(i),r_2))

                fwsoil(i) = canopy%fwsoil(i)
                ssnow%evapfbl(i,:) = ssnow%rex(i,:)*dels*1000_r_2 ! mm water &
                !(root water extraction) per time step

             ! PH: mgk576, 13/10/17
             ! This is over the combined direct & diffuse leaves due to the
             ! way the loops fall above
             ELSEIF (cable_user%FWSOIL_SWITCH == 'hydraulics') THEN

                ! Transpiration: W m-2 -> kg m-2 s-1 -> mmol m-2 s-1
                IF (ecx(i) > 0.0) THEN
                   conv = KG_2_G * G_WATER_TO_MOL * MOL_2_MMOL
                   trans_mmol = ecx(i) / air%rlam(i) * conv
                ELSE
                   trans_mmol = 0.0
                END IF

                ! Calculate the leaf water potential.
                CALL calc_psi_leaf(canopy, trans_mmol, dels, veg%Cl(i), i)

                ! Flux from stem to leaf (mmol s-1) = change in leaf storage,
                ! plus transpiration
                CALL calc_flux_to_leaf(canopy, trans_mmol, dels, veg%Cl(i), i)

                ! Update stem water potential
                CALL update_stem_wp(canopy, ssnow, dels, veg%Cs(i), i)

                ! Flux from the soil to the stem = change in storage +
                ! flux_to_leaf
                CALL calc_flux_to_stem(canopy, dels, veg%Cs(i), i)

                ! Update psi_stem
                canopy%psi_stem_prev(i) = canopy%psi_stem(i)
                CALL update_stem_wp_again(canopy, ssnow, dels, veg%Cs(i), &
                                          veg%Cl(i), trans_mmol, i)

                ! Flux from the soil to the stem = change in storage +
                ! flux_to_leaf
                CALL calc_flux_to_stem_again(canopy, dels, veg%Cs(i), &
                                             trans_mmol, i)

                !print*, canopy%psi_soil_prev(i), ssnow%weighted_psi_soil(i), &
                !         canopy%psi_stem_prev(i) , canopy%psi_stem(i), &
                !         canopy%psi_leaf_prev(i), canopy%psi_leaf(i)

                ! store current water potentials for next time step
                canopy%psi_leaf_prev(i) = canopy%psi_leaf(i)
                canopy%psi_soil_prev(i) = ssnow%weighted_psi_soil(i)
                canopy%psi_stem_prev(i) = canopy%psi_stem(i)

                ! Force overnight refilling - we need the real time
                ! as this won't work with 30 min and hourly
                !IF ( (met%hod(i) >= 12 .AND. met%hod(i) < 13) .AND. &
                !IF ( (met%hod(i) >= 6 .AND. met%hod(i) < 7) .AND. &
                !      canopy%psi_stem(i) > -4.0) THEN
                !
                !    !print*, canopy%psi_stem(i)
                !    refill = abs(canopy%psi_stem(i) - &
                !                 ssnow%weighted_psi_soil(i)) * 0.7
                !    canopy%psi_stem(i) = canopy%psi_stem(i) + refill
                !
                !    ! Ensure we can't refill above psi_soil
                !    canopy%psi_stem(i) = min(canopy%psi_stem(i), &
                !                             ssnow%weighted_psi_soil(i))
                !    print*, "****", canopy%psi_stem_prev(i), canopy%psi_stem(i)
                !    canopy%psi_stem_prev(i) = canopy%psi_stem(i)
                !    !print*, canopy%psi_stem(i)
                !    !print*, " "
                !
                !ENDIF

                !IF (ecx(i) > 0.0 .AND. canopy%fwet(i) < 1.0) THEN
                !    evapfb(i) = ( 1.0 - canopy%fwet(i)) * REAL( ecx(i) ) *dels      &
                !         / air%rlam(i)
                !
                !    DO kk = 1,ms
                !
                !       ssnow%evapfbl(i,kk) = MIN( evapfb(i) * ssnow%fraction_uptake(i,k),      &
                !            MAX( 0.0, REAL( ssnow%wb(i,kk) ) -     &
                !            1.1 * soil%swilt(i) ) *                &
                !            soil%zse(kk) * 1000.0 )
                !
                !   ENDDO
                !
                !      IF (cable_user%soil_struc=='default') THEN
                !       canopy%fevc(i) = SUM(ssnow%evapfbl(i,:))*air%rlam(i)/dels
                !         ecx(i) = canopy%fevc(i) / (1.0-canopy%fwet(i))
                !   ENDIF
                !
                ! ENDIF

                ! We shouldn't be using this, but I need to figure out how to
                ! correctly do the above...
                canopy%fevc(i) = ecx(i)*(1.0-canopy%fwet(i))

                CALL getrex_1d(ssnow%wbliq(i,:),&
                     ssnow%rex(i,:), &
                     canopy%fwsoil(i), &
                     REAL(veg%froot(i,:),r_2),&
                     soil%ssat_vec(i,:), &
                     soil%swilt_vec(i,:), &
                     MAX(REAL(canopy%fevc(i)/air%rlam(i)/1000_r_2,r_2),0.0_r_2), &
                     REAL(veg%gamma(i),r_2), &
                     REAL(soil%zse,r_2), REAL(dels,r_2), REAL(veg%zr(i),r_2))

                canopy%fwsoil(i) = 1.0
                ssnow%evapfbl(i,:) = ssnow%rex(i,:)*dels*1000_r_2 ! mm water &

             ELSE



                IF (ecx(i) > 0.0 .AND. canopy%fwet(i) < 1.0) THEN
                   evapfb(i) = ( 1.0 - canopy%fwet(i)) * REAL( ecx(i) ) *dels      &
                        / air%rlam(i)

                   DO kk = 1,ms

                      ssnow%evapfbl(i,kk) = MIN( evapfb(i) * veg%froot(i,kk),      &
                           MAX( 0.0, REAL( ssnow%wb(i,kk) ) -     &
                           1.1 * soil%swilt(i) ) *                &
                           soil%zse(kk) * 1000.0 )

                   ENDDO
                   IF (cable_user%soil_struc=='default' .AND. cable_user%fwsoil_switch.NE.'hydraulics') THEN
                   !IF (cable_user%soil_struc=='default') THEN
                      canopy%fevc(i) = SUM(ssnow%evapfbl(i,:))*air%rlam(i)/dels
                      ecx(i) = canopy%fevc(i) / (1.0-canopy%fwet(i))
                   !ELSEIF (cable_user%soil_struc=='sli') THEN
                   ELSEIF (cable_user%soil_struc=='default' .AND. cable_user%fwsoil_switch.EQ.'hydraulics') THEN
                      canopy%fevc(i) = SUM(ssnow%evapfbl(i,:))*air%rlam(i)/dels
                      ecx(i) = canopy%fevc(i) / (1.0-canopy%fwet(i))
                   ELSEIF (cable_user%soil_struc=='sli') THEN
                      canopy%fevc(i) = ecx(i)*(1.0-canopy%fwet(i))
                   ENDIF

                ENDIF

             ENDIF
             ! Update canopy sensible heat flux:
             hcx(i) = (SUM(rad%rniso(i,:))-ecx(i)                               &
                  - C%capp*C%rmair*(met%tvair(i)-met%tk(i))                       &
                  * SUM(rad%gradis(i,:)))                                         &
                  * SUM(gh(i,:))/ SUM(ghr(i,:))

             ! Update leaf temperature:
             tlfx(i)=met%tvair(i)+REAL(hcx(i))/(C%capp*C%rmair*SUM(gh(i,:)))

             ! Update net radiation for canopy:
             rnx(i) = SUM( rad%rniso(i,:)) -                                    &
                  C%CAPP * C%rmair *( tlfx(i)-met%tk(i) ) *                 &
                  SUM( rad%gradis(i,:) )

             ! Update leaf surface vapour pressure deficit:
             dsx(i) = met%dva(i) + air%dsatdk(i) * (tlfx(i)-met%tvair(i))

             dsx(i)=  MAX(dsx(i),0.0)

             ! Store change in leaf temperature between successive iterations:
             deltlf(i) = tlfxx(i)-tlfx(i)
             abs_deltlf(i) = ABS(deltlf(i))

          ENDIF !lai/abs_deltlf

       ENDDO !i=1,mp
       ! Where leaf temp change b/w iterations is significant, and
       ! difference is smaller than the previous iteration, store results:
       DO i=1,mp

          IF ( abs_deltlf(i) < ABS( deltlfy(i) ) ) THEN

             deltlfy(i) = deltlf(i)
             tlfy(i) = tlfx(i)
             rny(i) = rnx(i)
             hcy(i) = hcx(i)
             ecy(i) = ecx(i)
             rdy(i,1) = rdx(i,1)
             rdy(i,2) = rdx(i,2)
             an_y(i,1) = anx(i,1)
             an_y(i,2) = anx(i,2)

             ! save last values calculated for ssnow%evapfbl
             oldevapfbl(i,:) = ssnow%evapfbl(i,:)

          ENDIF

          IF( abs_deltlf(i) > 0.1 )                                             &

                                ! after 4 iterations, take mean of current & previous estimates
                                ! as the next estimate of leaf temperature, to avoid oscillation
               tlfx(i) = ( 0.5 * ( MAX( 0, k-5 ) / ( k - 4.9999 ) ) ) *tlfxx(i) + &
               ( 1.0 - ( 0.5 * ( MAX( 0, k-5 ) / ( k - 4.9999 ) ) ) )   &
               * tlfx(i)

          IF(k==1) THEN

             ! take the first iterated estimates as the defaults
             tlfy(i) = tlfx(i)
             rny(i) = rnx(i)
             hcy(i) = hcx(i)
             ecy(i) = ecx(i)
             rdy(i,:) = rdx(i,:)
             an_y(i,:) = anx(i,:)
             ! save last values calculated for ssnow%evapfbl
             oldevapfbl(i,:) = ssnow%evapfbl(i,:)

          END IF

       END DO !over mp

    END DO  ! DO WHILE (ANY(abs_deltlf > 0.1) .AND.  k < C%MAXITER)


    ! dry canopy flux
    canopy%fevc = (1.0-canopy%fwet) * ecy

    IF (cable_user%fwsoil_switch.NE.'Haverd2013' .AND. cable_user%fwsoil_switch.NE.'hydraulics') THEN

       ! Recalculate ssnow%evapfbl as ecy may not be updated with the ecx
       ! calculated in the last iteration.
       ! DO NOT use simple scaling as there are times that ssnow%evapfbl is zero.
       ! ** ssnow%evapfbl(i,:) = ssnow%evapfbl(i,:) * ecy(i) / ecx(i) **
       DO i = 1, mp

          IF( ecy(i) > 0.0 .AND. canopy%fwet(i) < 1.0 ) THEN

             IF( ABS( ecy(i) - ecx(i) ) > 1.0e-6 ) THEN

                IF( ABS( canopy%fevc(i) - ( SUM( oldevapfbl(i,:)) * air%rlam(i)    &
                     /dels ) ) > 1.0e-4 ) THEN

                   PRINT *, 'Error! oldevapfbl not right.', ktau_gl, i
                   PRINT *, 'ecx, ecy = ', ecx(i), ecy(i)
                   PRINT *, 'or in mm = ', ecx(i) * ( 1.0 - canopy%fwet(i) )       &
                        / air%rlam(i) * dels,                   &
                        ecy(i) * ( 1.0 - canopy%fwet(i) ) /     &
                        air%rlam(i) * dels

                   PRINT *,'fevc = ', canopy%fevc(i), SUM( oldevapfbl(i,:) ) *     &
                        air%rlam(i) / dels
                   PRINT *, 'fwet = ', canopy%fwet(i)
                   PRINT *, 'oldevapfbl = ', oldevapfbl(i,:)
                   PRINT *, 'ssnow%evapfbl before rescaling: ',                    &
                        ssnow%evapfbl(i,:)
                   ! STOP

                ELSE

                   ssnow%evapfbl(i,:) = oldevapfbl(i,:)

                END IF

             END IF

          END IF

       END DO

    ENDIF

    canopy%frday = 12.0 * SUM(rdy, 2)
    !! vh !! inserted min to avoid -ve values of GPP
    canopy%fpn = MIN(-12.0 * SUM(an_y, 2), canopy%frday)
    canopy%evapfbl = ssnow%evapfbl


    DEALLOCATE( gswmin )

  END SUBROUTINE dryLeaf
  ! -----------------------------------------------------------------------------


  ! Ticket #56, xleuningz repalced with gs_coeffz
  SUBROUTINE photosynthesis( csxz, cx1z, cx2z, gswminz,                          &
       rdxz, vcmxt3z, vcmxt4z, vx3z,                       &
       vx4z, gs_coeffz, vlaiz, deltlfz, anxz, fwsoilz )
    USE cable_def_types_mod, ONLY : mp, mf, r_2

    REAL(r_2), DIMENSION(mp,mf), INTENT(IN) :: csxz

    REAL, DIMENSION(mp,mf), INTENT(IN) ::                                       &
         cx1z,       & !
         cx2z,       & !
         gswminz,    & !
         rdxz,       & !
         vcmxt3z,    & !
         vcmxt4z,    & !
         vx4z,       & !
         vx3z,       & !
         gs_coeffz,  & ! Ticket #56, xleuningz repalced with gs_coeffz
         vlaiz,      & !
         deltlfz       !

    REAL, DIMENSION(mp,mf), INTENT(INOUT) :: anxz

    ! local variables
    REAL(r_2), DIMENSION(mp,mf) ::                                              &
         coef0z,coef1z,coef2z, ciz,delcxz,                                        &
         anrubiscoz,anrubpz,ansinkz

    REAL, DIMENSION(mp) :: fwsoilz

    REAL, PARAMETER  :: effc4 = 4000.0  ! Vc=effc4*Ci*Vcmax (see
    ! Bonan,LSM version 1.0, p106)

    INTEGER :: i,j

    DO i=1,mp

       IF (SUM(vlaiz(i,:)) .GT. C%LAI_THRESH) THEN

          DO j=1,mf

             IF( vlaiz(i,j) .GT. C%LAI_THRESH .AND. deltlfz(i,j) .GT. 0.1) THEN

                ! Rubisco limited:
                coef2z(i,j) = gswminz(i,j)*fwsoilz(i) / C%RGSWC + gs_coeffz(i,j) * &
                     ( vcmxt3z(i,j) - ( rdxz(i,j)-vcmxt4z(i,j) ) )

                coef1z(i,j) = (1.0-csxz(i,j)*gs_coeffz(i,j)) *                  &
                     (vcmxt3z(i,j)+vcmxt4z(i,j)-rdxz(i,j))             &
                     + (gswminz(i,j)*fwsoilz(i)/C%RGSWC)*(cx1z(i,j)-csxz(i,j)) &
                     - gs_coeffz(i,j)*(vcmxt3z(i,j)*cx2z(i,j)/2.0      &
                     + cx1z(i,j)*(rdxz(i,j)-vcmxt4z(i,j) ) )


                coef0z(i,j) = -(1.0-csxz(i,j)*gs_coeffz(i,j)) *                 &
                     (vcmxt3z(i,j)*cx2z(i,j)/2.0                       &
                     + cx1z(i,j)*( rdxz(i,j)-vcmxt4z(i,j ) ) )         &
                     -( gswminz(i,j)*fwsoilz(i)/C%RGSWC ) * cx1z(i,j)*csxz(i,j)


                ! kdcorbin,09/10 - new calculations
                IF( ABS(coef2z(i,j)) .GT. 1.0e-9 .AND. &
                     ABS(coef1z(i,j)) .LT. 1.0e-9) THEN

                   ! no solution, give it a huge number as
                   ! quadratic below cannot handle zero denominator
                   ciz(i,j) = 99999.0

                   anrubiscoz(i,j) = 99999.0 ! OR do ciz=0 and calc anrubiscoz

                ENDIF

                ! solve linearly
                IF( ABS( coef2z(i,j) ) < 1.e-9 .AND.                            &
                     ABS( coef1z(i,j) ) >= 1e-9 ) THEN

                   ! same reason as above
                   ciz(i,j) = -1.0 * coef0z(i,j) / coef1z(i,j)

                   ciz(i,j) = MAX( 0.0_r_2, ciz(i,j) )

                   anrubiscoz(i,j) = vcmxt3z(i,j)*(ciz(i,j)-cx2z(i,j) / 2.0 ) / &
                        ( ciz(i,j) + cx1z(i,j)) + vcmxt4z(i,j) -   &
                        rdxz(i,j)

                ENDIF

                ! solve quadratic (only take the more positive solution)
                IF( ABS( coef2z(i,j) ) >= 1.e-9 ) THEN

                   delcxz(i,j) = coef1z(i,j)**2 -4.0 * coef0z(i,j)              &
                        * coef2z(i,j)

                   ciz(i,j) = ( -coef1z(i,j) + SQRT( MAX( 0.0_r_2 ,             &
                        delcxz(i,j) ) ) ) / ( 2.0*coef2z(i,j) )

                   ciz(i,j) = MAX( 0.0_r_2, ciz(i,j) )   ! must be positive, why?

                   anrubiscoz(i,j) = vcmxt3z(i,j) * ( ciz(i,j) - cx2z(i,j)      &
                        / 2.0)  / ( ciz(i,j) + cx1z(i,j) ) +       &
                        vcmxt4z(i,j) - rdxz(i,j)

                ENDIF

                ! RuBP limited:
                coef2z(i,j) = gswminz(i,j)*fwsoilz(i) / C%RGSWC + gs_coeffz(i,j) &
                     * ( vx3z(i,j) - ( rdxz(i,j) - vx4z(i,j) ) )

                coef1z(i,j) = ( 1.0 - csxz(i,j) * gs_coeffz(i,j) ) *            &
                     ( vx3z(i,j) + vx4z(i,j) - rdxz(i,j) )             &
                     + ( gswminz(i,j)*fwsoilz(i) / C%RGSWC ) *          &
                     ( cx2z(i,j) - csxz(i,j) ) - gs_coeffz(i,j)        &
                     * ( vx3z(i,j) * cx2z(i,j) / 2.0 + cx2z(i,j) *     &
                     ( rdxz(i,j) - vx4z(i,j) ) )

                coef0z(i,j) = -(1.0-csxz(i,j)*gs_coeffz(i,j)) *   &
                     (vx3z(i,j)*cx2z(i,j)/2.0                          &
                     + cx2z(i,j)*(rdxz(i,j)-vx4z(i,j)))                &
                     - (gswminz(i,j)*fwsoilz(i)/C%RGSWC)*cx2z(i,j)*csxz(i,j)


                !Ticket #117 - initialize at all times
                ciz(i,j) = 99999.0
                anrubpz(i,j)  = 99999.0

                ! solve linearly
                IF( ABS( coef2z(i,j) ) < 1.e-9 .AND.                            &
                     ABS( coef1z(i,j) ) >= 1.e-9) THEN

                   ciz(i,j) = -1.0 * coef0z(i,j) / coef1z(i,j)

                   ciz(i,j) = MAX(0.0_r_2,ciz(i,j))

                   anrubpz(i,j) = vx3z(i,j)*(ciz(i,j)-cx2z(i,j)/2.0) /          &
                        (ciz(i,j)+cx2z(i,j)) +vx4z(i,j)-rdxz(i,j)

                ENDIF

                ! solve quadratic (only take the more positive solution)
                IF ( ABS( coef2z(i,j)) >= 1.e-9 ) THEN

                   delcxz(i,j) = coef1z(i,j)**2 -4.0*coef0z(i,j)*coef2z(i,j)

                   ciz(i,j) = (-coef1z(i,j)+SQRT(MAX(0.0_r_2,delcxz(i,j))))     &
                        /(2.0*coef2z(i,j))

                   ciz(i,j) = MAX(0.0_r_2,ciz(i,j))

                   anrubpz(i,j)  = vx3z(i,j)*(ciz(i,j)-cx2z(i,j)/2.0) /         &
                        (ciz(i,j)+cx2z(i,j)) +vx4z(i,j)-rdxz(i,j)

                ENDIF

                ! Sink limited:
                coef2z(i,j) = gs_coeffz(i,j)

                coef1z(i,j) = gswminz(i,j)*fwsoilz(i)/C%RGSWC + gs_coeffz(i,j)   &
                     * (rdxz(i,j) - 0.5*vcmxt3z(i,j))                  &
                     + effc4 * vcmxt4z(i,j) - gs_coeffz(i,j)           &
                     * csxz(i,j) * effc4 * vcmxt4z(i,j)

                coef0z(i,j) = -( gswminz(i,j)*fwsoilz(i)/C%RGSWC )*csxz(i,j)*effc4 &
                     * vcmxt4z(i,j) + ( rdxz(i,j)                      &
                     - 0.5 * vcmxt3z(i,j)) * gswminz(i,j)*fwsoilz(i)/C%RGSWC

                ! no solution, give it a huge number
                IF( ABS( coef2z(i,j) ) < 1.0e-9 .AND.                           &
                     ABS( coef1z(i,j)) < 1.0e-9 ) THEN

                   ciz(i,j) = 99999.0
                   ansinkz(i,j)  = 99999.0

                ENDIF

                ! solve linearly
                IF( ABS( coef2z(i,j) ) < 1.e-9 .AND.                            &
                     ABS( coef1z(i,j) ) >= 1.e-9 ) THEN

                   ciz(i,j) = -1.0 * coef0z(i,j) / coef1z(i,j)
                   ansinkz(i,j)  = ciz(i,j)

                ENDIF

                ! solve quadratic (only take the more positive solution)
                IF( ABS( coef2z(i,j) ) >= 1.e-9 ) THEN

                   delcxz(i,j) = coef1z(i,j)**2 -4.0*coef0z(i,j)*coef2z(i,j)

                   ciz(i,j) = (-coef1z(i,j)+SQRT (MAX(0.0_r_2,delcxz(i,j)) ) )  &
                        / ( 2.0 * coef2z(i,j) )

                   ansinkz(i,j) = ciz(i,j)

                ENDIF

                ! minimal of three limited rates
                anxz(i,j) = MIN(anrubiscoz(i,j),anrubpz(i,j),ansinkz(i,j))


             ENDIF

          ENDDO

       ENDIF

    ENDDO



  END SUBROUTINE photosynthesis

  ! ------------------------------------------------------------------------------

  FUNCTION ej3x(parx,alpha,convex,x) RESULT(z)

    REAL, INTENT(IN)     :: parx
    REAL, INTENT(IN)     :: alpha
    REAL, INTENT(IN)     :: convex
    REAL, INTENT(IN)     :: x
    REAL                 :: z

    z = MAX(0.0,                                                                &
         0.25*((alpha*parx+x-SQRT((alpha*parx+x)**2 -                      &
         4.0*convex*alpha*parx*x)) /(2.0*convex)) )
  END FUNCTION ej3x

  ! ------------------------------------------------------------------------------

  FUNCTION ej4x(parx,alpha,convex,x) RESULT(z)

    REAL, INTENT(IN)     :: parx
    REAL, INTENT(IN)     :: alpha
    REAL, INTENT(IN)     :: convex
    REAL, INTENT(IN)     :: x
    REAL                 :: z

    z = MAX(0.0,                                                                &
         (alpha*parx+x-SQRT((alpha*parx+x)**2 -                           &
         4.0*convex*alpha*parx*x))/(2.0*convex))

  END FUNCTION ej4x

  ! ------------------------------------------------------------------------------

  ! Explicit array dimensions as temporary work around for NEC inlining problem
  FUNCTION xvcmxt4(x) RESULT(z)

    REAL, PARAMETER      :: q10c4 = 2.0
    REAL, INTENT(IN) :: x
    REAL :: z

    z = q10c4 ** (0.1 * x - 2.5) /                                              &
         ((1.0 + EXP(0.3 * (13.0 - x))) * (1.0 + EXP(0.3 * (x - 36.0))))

  END FUNCTION xvcmxt4

  ! ------------------------------------------------------------------------------

  FUNCTION xvcmxt3(x) RESULT(z)

    !  leuning 2002 (p c & e) equation for temperature response
    !  used for vcmax for c3 plants
    REAL, INTENT(IN) :: x
    REAL :: xvcnum,xvcden,z

    REAL, PARAMETER  :: EHaVc  = 73637.0  ! J/mol (Leuning 2002)
    REAL, PARAMETER  :: EHdVc  = 149252.0 ! J/mol (Leuning 2002)
    REAL, PARAMETER  :: EntropVc = 486.0  ! J/mol/K (Leuning 2002)
    REAL, PARAMETER  :: xVccoef = 1.17461 ! derived parameter
    ! xVccoef=1.0+exp((EntropJx*C%TREFK-EHdJx)/(Rconst*C%TREFK))

    xvcnum=xvccoef*EXP( ( ehavc / ( C%rgas*C%TREFK ) )* ( 1.-C%TREFK/x ) )
    xvcden=1.0+EXP( ( entropvc*x-ehdvc ) / ( C%rgas*x ) )
    z = MAX( 0.0,xvcnum / xvcden )

  END FUNCTION xvcmxt3

  ! ------------------------------------------------------------------------------
  REAL FUNCTION xrdt(x)

    !  Atkins et al. (Eq 1, New Phytologist (2015) 206: 614–636)
    !variable Q10 temperature of dark respiration
    ! Originally from Tjoelker et al. 2001

    REAL, INTENT(IN) :: x


    xrdt = (3.09 - 0.043*((x-273.15)+25.)/2.0)**((x-273.15 -25.0)/10.0)

  END FUNCTION xrdt

  ! ------------------------------------------------------------------------------

  FUNCTION xejmxt3(x) RESULT(z)

    !  leuning 2002 (p c & e) equation for temperature response
    !  used for jmax for c3 plants

    REAL, INTENT(IN) :: x
    REAL :: xjxnum,xjxden,z

    REAL, PARAMETER  :: EHaJx  = 50300.0  ! J/mol (Leuning 2002)
    REAL, PARAMETER  :: EHdJx  = 152044.0 ! J/mol (Leuning 2002)
    REAL, PARAMETER  :: EntropJx = 495.0  ! J/mol/K (Leuning 2002)
    REAL, PARAMETER  :: xjxcoef = 1.16715 ! derived parameter

    xjxnum = xjxcoef*EXP( ( ehajx / ( C%rgas*C%TREFK ) ) * ( 1.-C%TREFK / x ) )
    xjxden=1.0+EXP( ( entropjx*x-ehdjx) / ( C%rgas*x ) )
    z = MAX(0.0, xjxnum/xjxden)

  END FUNCTION xejmxt3

  ! ------------------------------------------------------------------------------

  SUBROUTINE fwsoil_calc_std(fwsoil, soil, ssnow, veg)
    USE cable_def_types_mod
    USE cable_common_module, ONLY : cable_user
    TYPE (soil_snow_type), INTENT(INOUT):: ssnow
    TYPE (soil_parameter_type), INTENT(INOUT)   :: soil
    TYPE (veg_parameter_type), INTENT(INOUT)    :: veg
    REAL, INTENT(OUT), DIMENSION(:):: fwsoil ! soil water modifier of stom. cond
    REAL, DIMENSION(mp) :: rwater ! soil water availability

    !note even though swilt_vec is defined in default model it is r_2
    !and even using real(_vec) gives results different from trunk (rounding
    !errors)

    IF (.NOT.cable_user%gw_model) THEN

       rwater = MAX(1.0e-9,                                                    &
            SUM(veg%froot * MAX(1.0e-9,MIN(1.0, REAL(ssnow%wb) -                   &
            SPREAD(soil%swilt, 2, ms))),2) /(soil%sfc-soil%swilt))

    ELSE
       rwater = MAX(1.0e-9,                                                    &
            SUM(veg%froot * MAX(1.0e-9,MIN(1.0, REAL((ssnow%wbliq -                 &
            soil%swilt_vec)/(soil%sfc_vec-soil%swilt_vec)) )),2) )

    ENDIF

    ! Remove vbeta #56
    IF(cable_user%GS_SWITCH == 'medlyn') THEN
       fwsoil = MAX(1.0e-4,MIN(1.0, rwater))
    ELSE
       fwsoil = MAX(1.0e-9,MIN(1.0, veg%vbeta * rwater))
    ENDIF

  END SUBROUTINE fwsoil_calc_std

  ! ------------------------------------------------------------------------------

  SUBROUTINE fwsoil_calc_non_linear(fwsoil, soil, ssnow, veg)
    USE cable_def_types_mod
    TYPE (soil_snow_type), INTENT(INOUT):: ssnow
    TYPE (soil_parameter_type), INTENT(INOUT)   :: soil
    TYPE (veg_parameter_type), INTENT(INOUT)    :: veg
    REAL, INTENT(OUT), DIMENSION(:):: fwsoil ! soil water modifier of stom. cond
    REAL, DIMENSION(mp) :: rwater ! soil water availability
    REAL, DIMENSION(mp,3)          :: xi, ti, si
    INTEGER :: j

    rwater = MAX(1.0e-9,                                                    &
         SUM(veg%froot * MAX(0.0,MIN(1.0, REAL(ssnow%wb) -                   &
         SPREAD(soil%swilt, 2, ms))),2) /(soil%sfc-soil%swilt))

    fwsoil = 1.

    rwater = soil%swilt + rwater * (soil%sfc-soil%swilt)

    xi(:,1) = soil%swilt
    xi(:,2) = soil%swilt + (soil%sfc - soil%swilt)/2.0
    xi(:,3) = soil%sfc

    ti(:,1) = 0.
    ti(:,2) = 0.9
    ti(:,3) = 1.0

    si(:,1) = (rwater - xi(:,2)) / ( xi(:,1) - xi(:,2)) *                       &
         (rwater - xi(:,3)) / ( xi(:,1) - xi(:,3))

    si(:,2) = (rwater - xi(:,1)) / ( xi(:,2) - xi(:,1)) *                       &
         (rwater - xi(:,3)) / ( xi(:,2) - xi(:,3))

    si(:,3) = (rwater - xi(:,1)) / ( xi(:,3) - xi(:,1)) *                       &
         (rwater - xi(:,2)) / ( xi(:,3) - xi(:,2))

    DO j=1,mp
       IF (rwater(j) < soil%sfc(j) - 0.02)                                      &
            fwsoil(j) = MAX(0.,MIN(1., ti(j,1)*si(j,1) +                          &
            ti(j,2)*si(j,2) + ti(j,3)*si(j,3)))

    ENDDO

  END SUBROUTINE fwsoil_calc_non_linear

  ! ------------------------------------------------------------------------------

  ! ypw 19/may/2010 soil water uptake efficiency (see Lai and Ktaul 2000)
  SUBROUTINE fwsoil_calc_Lai_Ktaul(fwsoil, soil, ssnow, veg)
    USE cable_def_types_mod
    TYPE (soil_snow_type), INTENT(INOUT):: ssnow
    TYPE (soil_parameter_type), INTENT(INOUT)   :: soil
    TYPE (veg_parameter_type), INTENT(INOUT)    :: veg
    REAL, INTENT(OUT), DIMENSION(:):: fwsoil ! soil water modifier of stom. cond
    INTEGER   :: ns
    REAL, PARAMETER ::rootgamma = 0.01   ! (19may2010)
    REAL, DIMENSION(mp)  :: dummy, normFac
    !--- local level dependent rwater
    REAL, DIMENSION(mp,ms)  :: frwater

    fwsoil(:) = 0.0
    normFac(:) = 0.0

    DO ns=1,ms

       dummy(:) = rootgamma/MAX(1.0e-3_r_2,ssnow%wb(:,ns)-soil%swilt(:))

       frwater(:,ns) = MAX(1.0e-4_r_2,((ssnow%wb(:,ns)-soil%swilt(:))/soil%ssat(:)) &
            ** dummy)

       fwsoil(:) = MIN(1.0,MAX(fwsoil(:),frwater(:,ns)))

       normFac(:) = normFac(:) + frwater(:,ns) * veg%froot(:,ns)

    ENDDO

  END SUBROUTINE fwsoil_calc_Lai_Ktaul

  ! ------------------------------------------------------------------------------
  SUBROUTINE fwsoil_calc_sli(fwsoil, soil, ssnow, veg)
    USE cable_def_types_mod
    TYPE (soil_snow_type), INTENT(INOUT):: ssnow
    TYPE (soil_parameter_type), INTENT(INOUT)   :: soil
    TYPE (veg_parameter_type), INTENT(INOUT)    :: veg
    REAL, INTENT(OUT), DIMENSION(:):: fwsoil ! soil water modifier of stom. cond
    REAL, DIMENSION(mp,ms):: tmp2d1, tmp2d2, delta_root, alpha2a_root, alpha2_root
    ! Lai and Katul formulation for root efficiency function  vh 17/07/09
    alpha2a_root = MAX(ssnow%wb-soil%swilt_vec, 0.001_r_2)/(soil%ssat_vec)
    tmp2d1 = ssnow%wb -soil%swilt_vec
    tmp2d2 = SPREAD(veg%gamma,2,ms)/tmp2d1*LOG(alpha2a_root)
    WHERE ((tmp2d1>0.001) .AND. (tmp2d2 > -10.0))
       alpha2_root = EXP(tmp2d2)
    ELSEWHERE
       alpha2_root = 0.0
    ENDWHERE

    WHERE (veg%froot>0.0)
       delta_root = 1.0
    ELSEWHERE
       delta_root = 0.0
    ENDWHERE

    fwsoil  = MAXVAL(alpha2_root*delta_root, 2)
    fwsoil  = MAX(0.0, fwsoil)

  END SUBROUTINE fwsoil_calc_sli

  !*********************************************************************************************************************

  SUBROUTINE getrex_1d(theta, rex, fws, Fs, thetaS, thetaw, Etrans, gamma, dx, dt, zr)

    ! root extraction : Haverd et al. 2013
    USE cable_def_types_mod, ONLY: r_2

    IMPLICIT NONE

    REAL(r_2), DIMENSION(:), INTENT(IN)    :: theta      ! volumetric soil moisture
    REAL(r_2), DIMENSION(:), INTENT(INOUT)   :: rex    ! water extraction per layer
    REAL(r_2),               INTENT(INOUT) :: fws    ! stomatal limitation factor
    REAL(r_2), DIMENSION(:), INTENT(IN)    :: Fs     ! root length density
    REAL(r_2), DIMENSION(:), INTENT(IN)    :: thetaS ! saturation soil moisture
    REAL(r_2), DIMENSION(:), INTENT(IN)    :: thetaw ! soil moisture at wiliting point
    REAL(r_2),               INTENT(IN)    :: Etrans ! total transpiration
    REAL(r_2),               INTENT(IN)    :: gamma  ! skew of Li & Katul alpha2 function
    REAL(r_2), DIMENSION(:), INTENT(IN)    :: dx     ! layer thicknesses (m)
    REAL(r_2),               INTENT(IN)    :: dt
    REAL(r_2),               INTENT(IN)    :: zr

    ! Gets rate of water extraction compatible with CABLE stomatal conductance model
    ! theta(:) - soil moisture(m3 m-3)
    ! rex(:)   - rate of water extraction by roots from layers (cm/h).
    REAL(r_2), DIMENSION(1:SIZE(theta)) ::  lthetar, alpha_root, delta_root, layer_depth
    REAL(r_2)                       :: trex, e3, one, zero
    INTEGER :: k

    e3 = 0.001
    one = 1.0;
    zero = 0.0;

    layer_depth(1) = 0.0
    DO k=2,SIZE(theta)
       layer_depth(k) = SUM(dx(1:k-1))
    ENDDO

    !theta(:)   = S(:)*thetaS(:)
    lthetar(:) = LOG(MAX(theta(:)-thetaw(:),e3)/thetaS(:))

    WHERE ((theta(:)-thetaw(:)) > e3)
       alpha_root(:) = EXP( gamma/MAX(theta(:)-thetaw(:), e3) * lthetar(:) )
    ELSEWHERE
       alpha_root(:) = zero
    endwhere

    WHERE (Fs(:) > zero .AND. layer_depth < zr )  ! where there are roots and we are aobe max rooting depth
       delta_root(:) = one
    ELSEWHERE
       delta_root(:) = zero
    endwhere

    rex(:) = alpha_root(:)*Fs(:)

    trex = SUM(rex(:))
    IF (trex > zero) THEN
       rex(:) = rex(:)/trex
    ELSE
       rex(:) = zero
    ENDIF
    rex(:) = Etrans*rex(:)


    ! reduce extraction efficiency where total extraction depletes soil moisture below wilting point
    WHERE (((rex*dt) > (theta(:)-thetaw(:))*dx(:)) .AND. ((rex*dt) > zero))
       alpha_root = alpha_root*(theta(:)-thetaw(:))*dx(:)/(1.1_r_2*rex*dt)
    endwhere
    rex(:) = alpha_root(:)*Fs(:)

    trex = SUM(rex(:))
    IF (trex > zero) THEN
       rex(:) = rex(:)/trex
    ELSE
       rex(:) = zero
    ENDIF
    rex(:) = Etrans*rex(:)

    ! check that the water available in each layer exceeds the extraction
    !if (any((rex*dt) > (theta(:)-0.01_r_2)*dx(:))) then
    IF (ANY(((rex*dt) > MAX((theta(:)-thetaw(:)),zero)*dx(:)) .AND. (Etrans > zero))) THEN
       fws = zero
       ! distribute extraction according to available water
       !rex(:) = (theta(:)-0.01_r_2)*dx(:)
       rex(:) = MAX((theta(:)-thetaw(:))*dx(:),zero)
       trex = SUM(rex(:))
       IF (trex > zero) THEN
          rex(:) = rex(:)/trex
       ELSE
          rex(:) = zero
       ENDIF
       rex(:) = Etrans*rex(:)
    ELSE
       fws    = MAXVAL(alpha_root(2:)*delta_root(2:))
    ENDIF

  END SUBROUTINE getrex_1d
  !*********************************************************************************************************************

  ! ----------------------------------------------------------------------------
  FUNCTION f_tuzet(psi_leaf, sf, psi_f) RESULT(fw)
     ! Empirical logistic function to describe the sensitivity of stomata
     ! to leaf water potential.
     !
     ! Sigmoid function assumes that stomata are insensitive to psi_leaf at
     ! values close to zero and that stomata rapidly close with decreasing
     ! psi_leaf.
     !
     ! Reference:
     ! ----------
     ! * Tuzet et al. (2003) A coupled model of stomatal conductance,
     !   photosynthesis and transpiration. Plant, Cell and Environment 26,
     !   1097–1116
     !
     ! Martin De Kauwe, 3rd June, 2019

     IMPLICIT NONE

     REAL             :: num, den, fw
     REAL, INTENT(IN) :: psi_leaf, sf, psi_f

     !print*, "tuzet", sf, psi_f, psi_leaf
     num = 1.0 + EXP(sf * psi_f)
     den = 1.0 + EXP(sf * (psi_f - psi_leaf))
     fw = num / den
     !fw = MAX(1.0e-9, MIN(1.0, num / den))

  END FUNCTION f_tuzet
  ! ----------------------------------------------------------------------------

  ! ----------------------------------------------------------------------------
  SUBROUTINE calc_hydr_conduc(canopy, ssnow, rad, veg, kp_sat, i)
     ! Calculate conductance terms (root to stem & stem to leaf)
     !
     ! Martin De Kauwe, 3rd June, 2019

     USE cable_common_module
     USE cable_def_types_mod

     IMPLICIT NONE

     TYPE (canopy_type), INTENT(INOUT)          :: canopy
     TYPE (soil_snow_type), INTENT(INOUT)       :: ssnow
     TYPE (radiation_type), INTENT(INOUT)       :: rad
     TYPE (veg_parameter_type), INTENT(INOUT)   :: veg

     INTEGER, INTENT(IN) :: i ! patch
     REAL, INTENT(IN)    :: kp_sat
     REAL                :: ksoil, kroot2stem, kplant

     ! Convert total below ground resistance to leaf-specific resistance.
     ! Belowground resistance is calculated on a ground area basis;
     ! multiplying by LAI converts to leaf area. his assumes that each canopy
     ! layer is connected to each soil layer, so that the roots in each soil
     ! layer supply water to each canopy layer, and that the fraction of roots
     ! supplying each canopy layer is the same as the leaf area in that layer.
     IF (canopy%vlaiw(i) > 0.0) THEN
        ssnow%tot_bg_resist(i) = ssnow%tot_bg_resist(i) * canopy%vlaiw(i)
     END IF

     !print*, "LAI check", canopy%vlaiw(i), canopy%vlaiw

     ! soil-to-root hydraulic conductance (mmol m-2 leaf area s-1 MPa-1)

     ! weird issue with patchfrac, not sure if this is an initialisiation call
     ! issue because canopy is called before soilsnow? Check. Add this so it
     ! runs
     !IF (ssnow%tot_bg_resist(i) < 0.0000000001) THEN
     !     ssnow%tot_bg_resist(i) = 1E9
     !END IF


     ksoil = 1.0 / ssnow%tot_bg_resist(i)

     ! Plant hydraulic conductance (mmol m-2 s-1 MPa-1). NB. depends on stem
     ! water potential from the previous timestep.
     kplant = kp_sat * fsig_hydr(canopy%psi_stem_prev(i), veg%X_hyd(i), &
                                 veg%p50(i), veg%s50(i))

     ! Conductance from root surface to the stem water pool (assumed to be
     ! halfway to the leaves)
     kroot2stem = 2.0 * kplant

     ! Conductance from soil to stem water store (mmol m-2 s-1 MPa-1)
     canopy%ksoil2stem(i) = 1.0 / (1.0 / ksoil + 1.0 / kroot2stem)

     ! Conductance from stem water store to leaf (mmol m-2 s-1 MPa-1)
     canopy%kstem2leaf(i) = 2.0 * kplant

  END SUBROUTINE calc_hydr_conduc
  ! ----------------------------------------------------------------------------

  ! ----------------------------------------------------------------------------
  FUNCTION fsig_hydr(psi_stem_prev, X_hyd, p50, s50) RESULT(relk)
     ! Calculate relative plant conductance as a function of xylem pressure
     !
     ! Martin De Kauwe, 3rd June, 2019

     IMPLICIT NONE

     REAL             :: PX, V, p, relk, PX50
     REAL, INTENT(IN) :: psi_stem_prev, X_hyd, p50, s50

     ! xylem pressure
     PX = ABS(psi_stem_prev)

     ! the xylem pressure (P) x% of the conductivity is lost
     PX50 = ABS(p50)

     V = (X_hyd - 100.) * LOG(1.0 - X_hyd / 100.)
     p = (PX / PX50)**((PX50 * s50) / V)

     ! relative conductance (K/Kmax) as a funcion of xylem pressure
     relk = (1. - X_hyd / 100.)**p
     !relk = max(1.0e-9, min(1.0, relk))

  END FUNCTION fsig_hydr
  ! ----------------------------------------------------------------------------

  ! ----------------------------------------------------------------------------
  SUBROUTINE calc_flux_to_leaf(canopy, transpiration, dels, Cl, i)
     ! Flux from stem to leaf = change in leaf storage, plus transpiration
     !
     ! Martin De Kauwe, 3rd June, 2019

     USE cable_common_module
     USE cable_def_types_mod

     IMPLICIT NONE

     TYPE (canopy_type), INTENT(INOUT)    :: canopy

     REAL, INTENT(IN)    :: transpiration, Cl
     REAL, INTENT(IN)    :: dels ! integration time step (s)
     INTEGER, INTENT(IN) :: i

     ! is there conductance in the trunk?
     IF (canopy%kstem2leaf(i) * canopy%vlaiw(i) > 1E-09) THEN

        ! sapflow rate from stem to leaf within the time step
        canopy%flx_to_leaf(i) = (canopy%psi_leaf(i) - &
                                 canopy%psi_leaf_prev(i)) * &
                                 Cl / dels + canopy%vlaiw(i) * transpiration

     ! no conductance in the trunk
     ELSE
        canopy%flx_to_leaf(i) = 0.0
     ENDIF

  END SUBROUTINE calc_flux_to_leaf
  ! ----------------------------------------------------------------------------

  ! ----------------------------------------------------------------------------
  SUBROUTINE calc_flux_to_stem(canopy, dels, Cs, i)
     ! Calculate the flux from the root to the stem, i.e. the root water
     ! uptake (mmol s-1) = change in stem storage plus flux_to_leaf
     !
     ! Martin De Kauwe, 3rd June, 2019

     USE cable_common_module
     USE cable_def_types_mod

     IMPLICIT NONE

     TYPE (canopy_type), INTENT(INOUT)    :: canopy

     INTEGER, INTENT(IN) :: i
     REAL, INTENT(IN)    :: dels ! integration time step (s)
     REAL, INTENT(IN)    :: Cs

     ! plant can take up water
     IF (canopy%ksoil2stem(i) * canopy%vlaiw(i) > 1E-09) THEN
        canopy%flx_to_stem(i) = (canopy%psi_stem(i) - &
                                 canopy%psi_stem_prev(i)) * &
                                 Cs / dels + canopy%flx_to_leaf(i)

     ! plant cannot take up water, change of psi_stem is solely due to
     ! flux_to_leaf (J_rl)
     ELSE
        canopy%flx_to_stem(i) = 0.0
     ENDIF

  END SUBROUTINE calc_flux_to_stem
  ! ----------------------------------------------------------------------------

  ! ----------------------------------------------------------------------------
  SUBROUTINE calc_flux_to_stem_again(canopy, dels, Cs, transpiration, i)
     ! Calculate the flux from the root to the stem, i.e. the root water
     ! uptake (mmol s-1) = change in stem storage plus flux_to_leaf
     !
     ! Martin De Kauwe, 3rd June, 2019

     USE cable_common_module
     USE cable_def_types_mod

     IMPLICIT NONE

     TYPE (canopy_type), INTENT(INOUT)    :: canopy

     INTEGER, INTENT(IN) :: i
     REAL, INTENT(IN)    :: dels ! integration time step (s)
     REAL, INTENT(IN)    :: Cs
     REAL, INTENT(IN)    :: transpiration

     ! plant can take up water
     IF (canopy%ksoil2stem(i) * canopy%vlaiw(i) > 1E-09) THEN
        canopy%flx_to_stem(i) = (canopy%psi_stem(i) - &
                                  canopy%psi_stem_prev(i) * &
                                  Cs / dels + (transpiration * canopy%vlaiw(i)))

     ! plant cannot take up water, change of psi_stem is solely due to
     ! flux_to_leaf (J_rl)
     ELSE
        canopy%flx_to_stem(i) = 0.0
     ENDIF

  END SUBROUTINE calc_flux_to_stem_again
  ! ----------------------------------------------------------------------------


  ! ----------------------------------------------------------------------------
  SUBROUTINE update_stem_wp(canopy, ssnow, dels, Cs, i)
     ! Calculate the flux from the stem to the leaf = change in leaf storage
     ! plus transpiration
     !
     !
     ! Reference:
     ! ==========
     ! * Xu, X., Medvigy, D., Powers, J. S., Becknell, J. M. and Guan, K.
     !   (2016), Diversity in plant hydraulic traits explains seasonal and
     !    inter-annual variations of vegetation dynamics in seasonally dry
     !    tropical forests. New Phytol, 212: 80–95. doi:10.1111/nph.14009.
     !
     ! Can write the dynamic equation as: dpsi_leaf_dt = b + a*psi_leaf
     ! Then it follows (Xu et al. 2016, Appendix, and Code).”
     !
     ! Martin De Kauwe, 3rd June, 2019

     USE cable_common_module
     USE cable_def_types_mod

     IMPLICIT NONE

     TYPE (canopy_type), INTENT(INOUT)    :: canopy
     TYPE (soil_snow_type), INTENT(INOUT) :: ssnow

     REAL                :: ap, bp
     REAL, INTENT(IN)    :: dels ! integration time step (s)
     REAL, INTENT(IN)    :: Cs
     INTEGER, INTENT(IN) :: i

     ! plant can take up water
     IF (canopy%ksoil2stem(i) * canopy%vlaiw(i) > 1E-09) THEN
        ap = - canopy%vlaiw(i) * canopy%ksoil2stem(i) / Cs
        bp = (canopy%psi_soil_prev(i) - canopy%flx_to_leaf(i)) / Cs
        canopy%psi_stem(i) = ((ap * canopy%psi_stem_prev(i) + bp) * &
                             EXP(ap * dels) - bp) / ap

     ! plant cannot take up water, change of psi_stem is solely due to
     ! flux_to_leaf (J_rl)
     ELSE
        canopy%psi_stem(i) = canopy%psi_stem_prev(i) - &
                                 canopy%flx_to_leaf(i) * dels / Cs

     ENDIF

     !if ( canopy%psi_stem(i) < -100.) THEN
      !  print*, canopy%ksoil2stem(i) * canopy%vlaiw(i), canopy%psi_stem(i), canopy%flx_to_leaf(i)
     !endif

  END SUBROUTINE update_stem_wp
  ! ----------------------------------------------------------------------------

  ! ----------------------------------------------------------------------------
  SUBROUTINE update_stem_wp_again(canopy, ssnow, dels, Cs, Cl, transpiration, i)
     ! Calculate the flux from the stem to the leaf = change in leaf storage
     ! plus transpiration
     !
     !  This is a simplified equation based on Xu et al., using the water
     !  potentials from the previous timestep
     !
     ! Reference:
     ! ==========
     ! * Xu, X., Medvigy, D., Powers, J. S., Becknell, J. M. and Guan, K.
     !   (2016), Diversity in plant hydraulic traits explains seasonal and
     !    inter-annual variations of vegetation dynamics in seasonally dry
     !    tropical forests. New Phytol, 212: 80–95. doi:10.1111/nph.14009.
     !
     ! Can write the dynamic equation as: dpsi_leaf_dt = b + a*psi_leaf
     ! Then it follows (Xu et al. 2016, Appendix, and Code).”
     !
     ! Martin De Kauwe, 3rd June, 2019

     USE cable_common_module
     USE cable_def_types_mod

     IMPLICIT NONE

     TYPE (canopy_type), INTENT(INOUT)    :: canopy
     TYPE (soil_snow_type), INTENT(INOUT) :: ssnow

     REAL                :: ap, bp
     REAL, INTENT(IN)    :: dels ! integration time step (s)
     REAL, INTENT(IN)    :: Cs, Cl, transpiration
     INTEGER, INTENT(IN) :: i

     ! plant can take up water
     IF (canopy%ksoil2stem(i) * canopy%vlaiw(i) > 1E-09) THEN
        ap = - canopy%vlaiw(i) * canopy%ksoil2stem(i) / (Cs + Cl)
        bp = (canopy%psi_soil_prev(i) - &
              (canopy%vlaiw(i) * transpiration)) / (Cs + Cl)
        canopy%psi_stem(i) = ((ap * canopy%psi_stem_prev(i) + bp) * &
                             EXP(ap * dels) - bp) / ap

     ! plant cannot take up water, change of psi_stem is solely due to
     ! flux_to_leaf (J_rl)
     ELSE
        canopy%psi_stem(i) = canopy%psi_stem_prev(i) - &
                              (canopy%vlaiw(i) * transpiration) * dels  / &
                              (Cs + Cl)

     ENDIF
     !print*, "***", canopy%psi_stem(i)

  END SUBROUTINE update_stem_wp_again
  ! ----------------------------------------------------------------------------


  ! ----------------------------------------------------------------------------
  SUBROUTINE calc_psi_leaf(canopy, transpiration, dels, Cl, i)
     ! Calculate the leaf water potential (MPa)
     !
     !
     ! Reference:
     ! ==========
     ! * Xu, X., Medvigy, D., Powers, J. S., Becknell, J. M. and Guan, K.
     !   (2016), Diversity in plant hydraulic traits explains seasonal and
     !    inter-annual variations of vegetation dynamics in seasonally dry
     !    tropical forests. New Phytol, 212: 80–95. doi:10.1111/nph.14009.
     !
     ! Can write the dynamic equation as: dpsi_leaf_dt = b + a*psi_leaf
     ! Then it follows (Xu et al. 2016, Appendix, and Code).”
     !
     ! Martin De Kauwe, 3rd June, 2019

     USE cable_common_module
     USE cable_def_types_mod

     IMPLICIT NONE

     TYPE (canopy_type), INTENT(INOUT)    :: canopy

     INTEGER, INTENT(IN) :: i
     REAL, INTENT(IN)    :: transpiration
     REAL                :: ap, bp
     REAL, INTENT(IN)    :: dels ! integration time step (s)
     REAL, INTENT(IN)    :: Cl

     ! is there is conductance in the trunk?
     IF (canopy%kstem2leaf(i) * canopy%vlaiw(i) > 1E-09) THEN

        ap = - canopy%vlaiw(i) * canopy%kstem2leaf(i) / Cl
        bp = (canopy%psi_stem_prev(i) * canopy%vlaiw(i) * &
              canopy%kstem2leaf(i) - &
              canopy%vlaiw(i) * transpiration) / Cl
        canopy%psi_leaf(i) = ((ap * canopy%psi_leaf_prev(i) + bp) *  &
                                EXP(ap * dels) - bp) / ap

     ! No conductance in the trunk, delta psi_leaf is due to transpiration
     ELSE
        canopy%psi_leaf(i) = (canopy%psi_leaf_prev(i) - &
                              canopy%vlaiw(i) * transpiration * dels) / Cl
     ENDIF

     !ap = -(canopy%vlaiw(i) * canopy%kstem2leaf(i) / Cl)
     !bp = (canopy%vlaiw(i) * canopy%kstem2leaf(i) * canopy%psi_stem_prev(i) - &
      !      canopy%vlaiw(i) * transpiration) / Cl

     !canopy%psi_leaf(i) = ((ap * canopy%psi_leaf_prev(i) + bp) *  &
      !                        EXP(ap * dels) - bp) / ap


  END SUBROUTINE calc_psi_leaf
  ! ----------------------------------------------------------------------------

END MODULE cable_canopy_module
