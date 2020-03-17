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
!          replaced with gs_coeff, gs_coeffz. If GS_SWITCH is set to "leuning",
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
! ==============================================================================

MODULE cable_canopy_module

  use cable_data_module, only : icanopy_type, point2constants
  use cable_io_vars_module, only: wlogn
  use cable_common_module, only: cable_user

  implicit none

  public :: define_canopy, xvcmxt3, xejmxt3, ej3x, xrdt, xgmesT, &
         xvcmxt3_acclim, xejmxt3_acclim

  private

  type( icanopy_type ) :: C

CONTAINS


  SUBROUTINE define_canopy(bal, rad, rough, air, met, dels, ssnow, soil, veg, canopy, climate)
    
    USE cable_def_types_mod
    USE cable_radiation_module
    USE cable_air_module
    USE cable_common_module
    USE cable_roughness_module
    USE sli_utils, ONLY : potential_evap

    implicit none

    TYPE(balances_type),       INTENT(INOUT) :: bal
    TYPE(radiation_type),      INTENT(INOUT) :: rad
    TYPE(roughness_type),      INTENT(INOUT) :: rough
    TYPE(air_type),            INTENT(INOUT) :: air
    TYPE(met_type),            INTENT(INOUT) :: met
    REAL,                      INTENT(IN)    :: dels ! integration time setp (s)
    TYPE(soil_snow_type),      INTENT(INOUT) :: ssnow
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil
    TYPE(veg_parameter_type),  INTENT(INOUT) :: veg
    TYPE(canopy_type),         INTENT(INOUT) :: canopy
    TYPE(climate_type),        INTENT(IN)    :: climate
    ! INTEGER,                  INTENT(IN)    :: wlogn

    INTEGER :: &
         iter,  & ! iteration #
         iterplus !

    REAL, DIMENSION(mp) :: &
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
         dq,            & ! sat sp
         sum_rad_rniso, & !
         sum_rad_gradis   !

    ! temporary buffers to simplify equations
    REAL, DIMENSION(mp) :: &
         ftemp, z_eff, rlower_limit, &
         term1, term2, term3, term5
    ! arguments for potential_evap (sli)
    REAL(r_2), DIMENSION(mp) ::  Rn, Ta

    REAL, DIMENSION(:), POINTER :: &
         cansat => null(),        & ! max canopy intercept. (mm)
         dsx => null(),           & ! leaf surface vpd
         fwsoil => null(),        & ! soil water modifier of stom. cond
         tlfx => null(),          & ! leaf temp prev. iter (K)
         tlfy => null()             ! leaf temp (K)

    REAL(r_2), DIMENSION(mp) :: &
         gbvtop                   ! bnd layer cond. top leaf

    REAL(r_2), DIMENSION(:), POINTER :: &
         ecy => null(),           & ! lat heat fl dry big leaf
         hcy => null(),           & ! veg. sens heat
         rny => null(),           & ! net rad
         ghwet => null()             ! cond for heat for a wet canopy

    REAL(r_2), DIMENSION(:,:), POINTER :: &
         gbhu => null(),          & ! forcedConvectionBndryLayerCond
         gbhf => null(),          & ! freeConvectionBndryLayerCond
         csx => null()              ! leaf surface CO2 concentration
    REAL(r_2), DIMENSION(:,:), POINTER ::  gmes => null()  ! mesophyll conductance      
    !mesophyll conductance extinction coeff't (Sun et al. 2014 SI Eq S7)      
    !REAL :: gmax0 ! max mesophyll conductacne at canopy top  
    REAL :: rt_min
    REAL, DIMENSION(mp) :: zstar, rL, phist, csw, psihat,rt0bus

    INTEGER :: j

    INTEGER, SAVE :: call_number =0

    ! END header

    call_number = call_number + 1

    ! assign local ptrs to constants defined in cable_data_module
    CALL point2constants(C)

    IF (.NOT. cable_runtime%um) canopy%cansto = canopy%oldcansto

    ALLOCATE(cansat(mp), gbhu(mp,mf))
    ALLOCATE(dsx(mp), fwsoil(mp), tlfx(mp), tlfy(mp))
    ALLOCATE(ecy(mp), hcy(mp), rny(mp))
    ALLOCATE(gbhf(mp,mf), csx(mp,mf))
    ALLOCATE(ghwet(mp))
    ALLOCATE(gmes(mp,mf))

    gmes = 0.0_r_2
    
    ! BATS-type canopy saturation proportional to LAI:
    cansat = veg%canst1 * canopy%vlaiw

    !---compute surface wetness factor, update cansto, through
    CALL surf_wetness_fact(cansat, canopy, ssnow, veg, met, soil, dels)

    canopy%fevw_pot = 0.0
    canopy%gswx     = 1.0e-3     ! default stomatal conuctance
    gbhf            = 1.0e-3_r_2 ! default free convection boundary layer conductance
    gbhu            = 1.0e-3_r_2 ! default forced convection boundary layer conductance
    ssnow%evapfbl   = 0.0
    ssnow%rex       = 0.0_r_2
    ! Initialise in-canopy temperatures and humidity:
    csx = SPREAD(real(met%ca,r_2), 2, mf) ! initialise leaf surface CO2 concentration
    met%tvair = met%tk
    met%qvair = met%qv
    canopy%tv = met%tvair

    ! Initialise sunlit and shaded Ac and Aj diagnostics
    canopy%A_sh           = 0.0_r_2
    canopy%A_sl           = 0.0_r_2
    canopy%A_shC          = 0.0_r_2
    canopy%A_slC          = 0.0_r_2
    canopy%A_shJ          = 0.0_r_2
    canopy%A_slJ          = 0.0_r_2
    canopy%GPP_sh         = 0.0_r_2
    canopy%GPP_sl         = 0.0_r_2
    canopy%fevc_sh        = 0.0_r_2
    canopy%fevc_sl        = 0.0_r_2
    canopy%eta_A_cs       = 0.0_r_2
    canopy%eta_GPP_cs     = 0.0_r_2
    canopy%eta_fevc_cs    = 0.0_r_2
    canopy%eta_A_cs_sh    = 0.0_r_2
    canopy%eta_A_cs_sl    = 0.0_r_2
    canopy%eta_fevc_cs_sh = 0.0_r_2
    canopy%eta_fevc_cs_sl = 0.0_r_2

    canopy%dAdcs   = 0.0_r_2
    canopy%cs      = 0.0_r_2
    canopy%cs_sl   = 0.0_r_2
    canopy%cs_sh   = 0.0_r_2
    canopy%tlf     = 0.0_r_2
    canopy%dlf     = 0.0_r_2
    canopy%evapfbl = 0.0_r_2
    
    ! 13C
    canopy%An        = 0.0_r_2
    canopy%Rd        = 0.0_r_2
    canopy%isc3      = .true.
    canopy%vcmax     = 0.0_r_2
    canopy%gammastar = 0.0_r_2
    canopy%gsc       = 0.0_r_2
    canopy%gbc       = 0.0_r_2
    canopy%gac       = 0.0_r_2
    canopy%ci        = 0.0_r_2

    CALL define_air(met, air)

    CALL qsatfjh(qstvair, met%tvair-C%tfrz, met%pmb)

    if (cable_user%perturb_dva_by_T) then
       CALL qsatfjh(qstvair, met%tvair-C%tfrz + cable_user%dva_T_perturbation, met%pmb)        
    endif
    met%dva = (qstvair - met%qvair) *  C%rmair/C%rmh2o * met%pmb * 100.0
    dsx     = met%dva     ! init. leaf surface vpd
    dsx     = max(dsx, 0.0)
    tlfx    = met%tk  ! initialise leaf temp iteration memory variable (K)
    tlfy    = met%tk  ! initialise current leaf temp (K)

    ortsoil = ssnow%rtsoil
    IF (cable_user%soil_struc=='default') then
       ssnow%tss = real((1-ssnow%isflag))*ssnow%tgg(:,1) + real(ssnow%isflag)*ssnow%tggsn(:,1)
    elseif (cable_user%soil_struc=='sli') then
       ssnow%tss = real(ssnow%Tsurface) + C%tfrz
    endif
    tss4 = ssnow%tss**4
    canopy%fes  = 0.0_r_2
    canopy%fess = 0.0_r_2
    canopy%fesp = 0.0_r_2
    ssnow%potev = 0.0
    canopy%fevw_pot = 0.0

    CALL radiation(ssnow, veg, air, met, rad, canopy)

    canopy%zetar(:,1)  = C%ZETA0 ! stability correction terms
    canopy%zetar(:,2)  = C%ZETPOS + 1.0
    canopy%zetash(:,1) = C%ZETA0 ! stability correction terms
    canopy%zetash(:,2) = C%ZETPOS + 1.0

    ! !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! ! SPECIAL for global offine met: set screen t and humidity as
    ! ! met inputs
    ! canopy%tscrn=   met%tk
    ! canopy%qscrn  = met%qv
    ! !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    DO iter = 1, NITER
       ! AERODYNAMIC PROPERTIES: friction velocity us, thence turbulent
       ! resistances rt0, rt1 (elements of dispersion matrix):
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.46:
       CALL comp_friction_vel()

       ! E.Kowalczyk 2014
       IF (cable_user%l_new_roughness_soil) &
            CALL ruff_resist(veg, rough, ssnow, canopy) 
       
       ! Turbulent aerodynamic resistance from roughness sublayer depth
       ! to reference height, x=1 if zref+disp>zruffs,
       ! 0 otherwise: thus rt1usc = 0 if zref+disp<zruffs
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.53:
       xx = 0.5 + sign(0.5,rough%zref_tq+rough%disp-rough%zruffs)

       ! correction  by Ian Harman to the 2nd psis term
       rt1usc = xx * ( LOG( rough%zref_tq / MAX(rough%zruffs-rough%disp, &
            rough%z0soilsn) ) &
            - psis(canopy%zetar(:,iter)) &
            + psis(canopy%zetar(:,iter) * MAX(rough%zruffs-rough%disp, &
            rough%z0soilsn) / rough%zref_tq) ) / C%VONK
       rt_min = 5.0

       !! vh_js !!
       IF (cable_user%soil_struc=='sli') THEN
          ! for stable conditions, update rough%rt0us & rough%rt1usa by replacing C%CSW by
          ! csw = cd/2* (U(hc)/ust)**2 according to Eqs 15 & 19 from notes by Ian Harman (9-9-2011)
          WHERE ((canopy%vlaiw > C%LAI_thresh) .and. (rough%hruff > rough%z0soilsn))
             rt0bus = (LOG(0.1*rough%hruff/rough%z0soilsn) - psis(canopy%zetash(:,iter)) + &
                  psis(canopy%zetash(:,iter)*rough%z0soilsn/(0.1*rough%hruff))) / &
                  C%vonk/rough%term6a

             zstar = rough%disp + 1.5*(veg%hc - rough%disp)

             psihat = log((zstar - rough%disp)/ (veg%hc - rough%disp)) + &
                  (veg%hc - zstar)/(zstar - rough%disp)
             rL = -(C%vonk*C%grav*(zstar - rough%disp)*(canopy%fh))/ &  ! 1/Monin-Obokov Length
                  max( (air%rho*C%capp*met%tk*canopy%us**3), 10.e-12)
             phist = 1.0 + 5.0*(zstar - rough%disp)*rL

             where (canopy%zetar(:,iter) .gt. 1.0e-6)! stable conditions

                csw = min(0.3*((log((veg%hc-rough%disp)/rough%z0m) + phist*psihat - &
                     psim(canopy%zetar(:,iter)*(veg%hc-rough%disp)/(rough%zref_tq-rough%disp))+ &
                     psim(canopy%zetar(:,iter)*rough%z0m/(rough%zref_tq-rough%disp)))/0.4)**2/2.0, 3.0)* c%csw
             elsewhere
                csw = c%csw
             endwhere

             rough%term2 = EXP( 2.0 * CSW * canopy%rghlai * &
                  ( 1.0 - rough%disp / rough%hruff ) )
             rough%term3 = C%A33**2 * C%CTL * 2.0 * CSW * canopy%rghlai
             rough%term5 = MAX( (2.0/3.0) * rough%hruff / rough%disp, 1.0 )
             rough%term6 =  EXP( 3.0 * rough%coexp * ( rough%disp / rough%hruff -1.0 ) )

             rough%rt0us  = log(rough%disp/(0.1 * rough%hruff)) * &
                  EXP(2.0 * C%CSW * canopy%rghlai) * rough%disp &
                  / rough%hruff / (c%a33**2 * c%ctl)

             rough%rt1usa = rough%term5 * ( rough%term2 - 1.0 ) / rough%term3
             rt0 = max(rt_min, rough%rt0us+rt0bus) / canopy%us
          ELSEWHERE
             rt0 = max(rt_min, rough%rt0us) / canopy%us
          ENDWHERE

       ELSE
          
          rt0 = max(rt_min,rough%rt0us / canopy%us)

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
       rough%rt1 = MAX(5.0,(rough%rt1usa + rough%rt1usb + rt1usc) / canopy%us)
       
       DO j=1, mp
          IF(canopy%vlaiw(j) > C%LAI_THRESH) THEN
             ssnow%rtsoil(j) = rt0(j)
          ELSE
             ssnow%rtsoil(j) = rt0(j) + rough%rt1(j)
          ENDif
       ENDDO

       ssnow%rtsoil = max(rt_min,ssnow%rtsoil)

       DO j=1, mp
          IF ( (ssnow%rtsoil(j) > 2.*ortsoil(j)) .OR. &
               (ssnow%rtsoil(j) < 0.5*ortsoil(j)) ) THEN
             ssnow%rtsoil(j) = MAX(rt_min, 0.5*(ssnow%rtsoil(j) + ortsoil(j)))
          ENDIF
       ENDDO

       ! Vegetation boundary-layer conductance (mol/m2/s)
       ! C%prandt = kinematic viscosity/molecular diffusivity
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.12. Top leaf:
       DO j=1, mp

          IF (canopy%vlaiw(j) > C%LAI_THRESH) THEN
             gbvtop(j) = real( air%cmolar(j)*C%APOL * air%visc(j) / C%prandt / &
                  veg%dleaf(j) * (canopy%us(j) / MAX(rough%usuh(j),1.e-6) &
                  * veg%dleaf(j) / air%visc(j) )**0.5 &
                  * C%prandt**(1.0/3.0) / veg%shelrb(j), r_2 )
             gbvtop(j) = MAX(0.05_r_2, gbvtop(j)) ! for testing (BP aug2010)

             ! Forced convection boundary layer conductance
             ! (see Wang & Leuning 1998, AFM):
             !vh! inserted 'min' to avoid floating underflow
             gbhu(j,1) = gbvtop(j)*(1.0_r_2 - EXP(-real( min(canopy%vlaiw(j) * &
                  (0.5*rough%coexp(j)+rad%extkb(j)), 20.0), r_2))) / &
                  real(rad%extkb(j)+0.5*rough%coexp(j), r_2)
             gbhu(j,2) = real((2.0/rough%coexp(j))*gbvtop(j), r_2) * &
                  (1.0_r_2-EXP(-real( min(0.5*rough%coexp(j)*canopy%vlaiw(j), 20.0), r_2))) &
                  - gbhu(j,1)

             if (cable_user%amphistomatous) then
                gbhu(j,1) = gbhu(j,1) * 2.0_r_2
                gbhu(j,2) = gbhu(j,2) * 2.0_r_2
             endif

             if (cable_user%explicit_gm) then
                ! JK: gmmax now comes from param file
                !gmax0 = 2.47/10.1 ! molm-2s-1 
                !if (veg%iveg(j).eq.1) then
                !   gmax0 = 1.21/10.1
                !elseif (veg%iveg(j).eq.2) then
                !   gmax0 = 1.36/10.1
                !elseif (veg%iveg(j).eq.3) then
                !   gmax0 = 2.14/10.1
                !elseif (veg%iveg(j).eq.7) then
                !   gmax0 = 0.3
                !endif
                gmes(j,1) = real(veg%gmmax(j) * rad%scalex(j,1), r_2) * xgmesT(tlfx(j))
                gmes(j,2) = real(veg%gmmax(j) * rad%scalex(j,2), r_2) * xgmesT(tlfx(j))
             endif              !
          ENDIF ! canopy%vlaiw(j) > C%LAI_THRESH

       ENDDO ! j=1, mp

       rny = sum(real(rad%rniso,r_2), 2) ! init current estimate net rad
       hcy = 0.0_r_2          ! init current estimate sensible heat
       ecy = rny - hcy        ! init current estimate latent heat

       sum_rad_rniso = sum(rad%rniso,2)
       CALL dryLeaf( dels, rad, rough, air, met, &
            veg, canopy, soil, ssnow, dsx,       &
            fwsoil, tlfx, tlfy, ecy, hcy,        &
            rny, gbhu, gbhf, csx, cansat,        &
            ghwet, iter, climate, gmes )

       CALL wetLeaf( dels, rad, rough, air, met, &
            veg, canopy, cansat, tlfy,           &
            gbhu, gbhf, ghwet )


       ! Calculate latent heat from vegetation:
       ! Calculate sensible heat from vegetation:
       ! Calculate net rad absorbed by canopy:
       canopy%fev = real(canopy%fevc) + canopy%fevw
       canopy%fhv = (1.0-canopy%fwet) * real(hcy) + canopy%fhvw
       canopy%fnv = (1.0-canopy%fwet) * real(rny) + canopy%fevw + canopy%fhvw
       !MCMCMC
       ! canopy rad. temperature calc from long-wave rad. balance
       sum_rad_gradis = SUM(rad%gradis,2)

       DO j=1,mp

          IF ( canopy%vlaiw(j) > C%LAI_THRESH .AND. &
               rough%hruff(j) > rough%z0soilsn(j) ) THEN

             rad%lwabv(j) = C%CAPP * C%rmair * ( tlfy(j) - met%tk(j) ) * &
                  sum_rad_gradis(j)
             !! vh_js !!

             if (  (rad%lwabv(j) / (2.0*(1.0-rad%transd(j)) &
                  * C%SBOLTZ*C%EMLEAF)+met%tk(j)**4) .gt. 0.0) then

                canopy%tv(j) = (rad%lwabv(j) / (2.0*(1.0-rad%transd(j)) &
                     * C%SBOLTZ*C%EMLEAF)+met%tk(j)**4)**0.25

             else
                canopy%tv(j) = met%tk(j)
             endif


          ELSE! sparse canopy

             canopy%tv(j) = met%tk(j)

          ENDIF

       ENDDO


       ! Calculate net rad to soil:
       canopy%fns = rad%qssabs + rad%transd*met%fld + (1.0-rad%transd)*C%EMLEAF* &
            C%SBOLTZ*canopy%tv**4 - C%EMSOIL*C%SBOLTZ* tss4


       ! Saturation specific humidity at soil/snow surface temperature:
       call qsatfjh(ssnow%qstss,ssnow%tss-C%tfrz,met%pmb)


       If (cable_user%soil_struc=='default') THEN

          IF(cable_user%ssnow_POTEV== "P-M") THEN

             !--- uses %ga from previous timestep
             ssnow%potev = Penman_Monteith(canopy%ga)

          ELSE !by default assumes Humidity Deficit Method

             ! Humidity deficit
             dq = ssnow%qstss - met%qv
             ssnow%potev =  Humidity_deficit_method(dq, ssnow%qstss)

          ENDIF

          ! Soil latent heat:
          CALL latent_heat_flux()

          ! Calculate soil sensible heat:
          !canopy%fhs = air%rho*C%CAPP*(ssnow%tss - met%tk) /ssnow%rtsoil
          IF (cable_user%litter) THEN
             !! vh_js !! account for additional litter resistance to sensible heat transfer
             canopy%fhs =  air%rho*C%CAPP*(ssnow%tss - met%tk) / &
                  (ssnow%rtsoil + real((1-ssnow%isflag)) * &
                  real(veg%clitt)*0.003/canopy%kthLitt/(air%rho*C%CAPP))
          ELSE
             canopy%fhs = air%rho*C%CAPP*(ssnow%tss - met%tk) /ssnow%rtsoil
          ENDIF

       ELSEIF (cable_user%soil_struc=='sli') THEN
          ! SLI SEB to get canopy%fhs, canopy%fess, canopy%ga
          ! (Based on old Tsoil, new canopy%tv, new canopy%fns)
          CALL sli_main(1,dels,veg,soil,ssnow,met,canopy,air,rad,1)
       ENDIF

!!$!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$	 ! SPECIAL for global offine met: set screen t and humidity as met inputs
!!$
!!$       tstar = - (canopy%fhv+canopy%fhs) / ( air%rho*C%capp*canopy%us)
!!$       qstar = - (canopy%fev+canopy%fes) / ( air%rho*air%rlam *canopy%us)
!!$       zscrn = max(rough%z0m,1.8-rough%disp)
!!$       denom = ( log(rough%zref_tq/zscrn)- psis(canopy%zetar(:,iter)) + &
!!$            psis(canopy%zetar(:,iter) * zscrn / rough%zref_tq) ) /C%vonk
!!$
!!$       where (canopy%zetar(:,iter) > 0.7)
!!$          zeta2=canopy%zetar(:,iter) * zscrn / rough%zref_tq
!!$          denom =alpha1* ((canopy%zetar(:,iter)**beta1* &
!!$               (1.0+gamma1*canopy%zetar(:,iter)**(1.0-beta1))) &
!!$               - (zeta2**beta1*(1.0+gamma1*zeta2**(1.0-beta1)))) /C%vonk
!!$       endwhere
!!$
!!$       where (abs( tstar*denom ).lt. 5.0)
!!$          met%tk = canopy%tscrn + tstar*denom
!!$          met%qv = max(canopy%qscrn + qstar*denom, 0.0)
!!$       endwhere
!!$!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       CALL within_canopy( gbhu, gbhf )

       ! Saturation specific humidity at soil/snow surface temperature:
       call qsatfjh(ssnow%qstss,ssnow%tss-C%tfrz,met%pmb)

       IF (cable_user%soil_struc=='default') THEN

          IF(cable_user%ssnow_POTEV== "P-M") THEN

             !--- uses %ga from previous timestep
             ssnow%potev = Penman_Monteith(canopy%ga)

          ELSE !by default assumes Humidity Deficit Method

             ! Humidity deficit
             dq = ssnow%qstss - met%qvair
             ssnow%potev =  Humidity_deficit_method(dq, ssnow%qstss)

          ENDIF

          ! Soil latent heat:
          CALL latent_heat_flux()

          ! Soil sensible heat:
          !canopy%fhs = air%rho*C%CAPP*(ssnow%tss - met%tvair) /ssnow%rtsoil
          IF (cable_user%litter) THEN
             !! vh_js !! account for additional litter resistance to sensible heat transfer
             canopy%fhs =  air%rho*C%CAPP*(ssnow%tss - met%tvair) / &
                  (ssnow%rtsoil +  real((1-ssnow%isflag))*veg%clitt*0.003/canopy%kthLitt/(air%rho*C%CAPP))
          else
             canopy%fhs = air%rho*C%CAPP*(ssnow%tss - met%tvair) /ssnow%rtsoil
          ENDIF

          !! Ticket #90 ssnow%cls factor should be retained: required for energy balance
          canopy%ga = canopy%fns-canopy%fhs-canopy%fes !*ssnow%cls

       ELSEIF (cable_user%soil_struc=='sli') THEN
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
             ssnow%potev(j) = max(0.00001,ssnow%potev(j))
          ELSE
             ssnow%potev(j) = min(-0.0002,ssnow%potev(j))
          ENDIF

          IF (canopy%fevw_pot(j) .ge. 0.) then
             canopy%fevw_pot(j) = max(0.000001,canopy%fevw_pot(j))
          ELSE
             canopy%fevw_pot(j) = min(-0.002,canopy%fevw_pot(j))
          ENDIF

       ENDDO

       canopy%rnet = canopy%fnv + canopy%fns

       canopy%epot = ((1.-rad%transd)*canopy%fevw_pot + &
            rad%transd*ssnow%potev*ssnow%cls) * dels/air%rlam

       canopy%rniso = sum(rad%rniso,2) + rad%qssabs + rad%transd*met%fld + &
            (1.0-rad%transd)*C%EMLEAF* &
            C%SBOLTZ*met%tk**4 - C%EMSOIL*C%SBOLTZ*met%tk**4

       rlower_limit = canopy%epot * air%rlam / dels
       where (rlower_limit == 0 ) rlower_limit = 1.e-7 !prevent from 0. by adding 1.e-7 (W/m2)

       canopy%wetfac_cs = max(0., min(1.0,canopy%fe / rlower_limit ))

       DO j=1,mp
          IF ( canopy%wetfac_cs(j) .LE. 0. ) &
               canopy%wetfac_cs(j) = MAX( 0., MIN( 1., &
               MAX( canopy%fev(j)/canopy%fevw_pot(j), &
               real(canopy%fes(j))/ssnow%potev(j) ) ) )
       ENDDO

       CALL update_zetar()

    END DO           ! do iter = 1, NITER

    canopy%cduv = canopy%us * canopy%us / (max(met%ua,C%UMIN))**2

    !---diagnostic purposes
    canopy%gswx_T = rad%fvlai(:,1)/MAX( C%LAI_THRESH, canopy%vlaiw(:) ) &
         * canopy%gswx(:,1) + rad%fvlai(:,2) / MAX(C%LAI_THRESH, &
         canopy%vlaiw(:))*canopy%gswx(:,2)

    ! The surface conductance below is required by dust scheme; it is composed from canopy and soil conductances

    !vh ! this line is causing a floating overflow error
    !canopy%gswx_T = (1.-rad%transd)*max(1.e-06,canopy%gswx_T ) + &   !contribution from  canopy conductance
     !    rad%transd*(.01*ssnow%wb(:,1)/soil%sfc)**2 ! + soil conductance; this part is done as in Moses
    where ( soil%isoilm == 9 ) canopy%gswx_T = 1.e6   ! this is a value taken from Moses for ice points

    canopy%cdtq = canopy%cduv *( LOG( rough%zref_uv / rough%z0m) - &
         psim( canopy%zetar(:,NITER) * rough%zref_uv/rough%zref_tq ) &
         + psim( canopy%zetar(:,NITER) * rough%z0m/rough%zref_tq ) & ! new term from Ian Harman
         ) / ( LOG( rough%zref_tq /(0.1*rough%z0m) ) &
         - psis( canopy%zetar(:,NITER)) &
         + psis(canopy%zetar(:,NITER)*0.1*rough%z0m/rough%zref_tq) ) ! n


    ! Calculate screen temperature: 1) original method from SCAM
    ! screen temp., windspeed and relative humidity at 1.5m
    ! screen temp., windspeed and relative humidity at 2.0m
    tstar = - canopy%fh / ( air%rho*C%CAPP*canopy%us)
    qstar = - canopy%fe / ( air%rho*air%rlam *canopy%us)
    zscrn = MAX(rough%z0m,2.0-rough%disp)
    ftemp = ( LOG(rough%zref_tq/zscrn)- psis(canopy%zetar(:,iterplus)) + &
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
    r_sc = 0.
    zscl = MAX(rough%z0soilsn,2.0)

    ! assume screen temp of bareground if all these conditions are not met
    DO j=1,mp

       IF ( canopy%vlaiw(j) > C%LAI_THRESH .and. rough%hruff(j) > 0.01) THEN

          IF ( rough%disp(j)  > 0.0 ) then

             term1(j) = EXP(2*C%CSW*canopy%rghlai(j)*(1-zscl(j)/rough%hruff(j)))
             term2(j) = EXP(2*C%CSW*canopy%rghlai(j) * &
                  (1-rough%disp(j)/rough%hruff(j)))
             term5(j) = MAX(2./3.*rough%hruff(j)/rough%disp(j), 1.)

          ENDIF

          term3(j) = C%A33**2*C%CTL*2*C%CSW*canopy%rghlai(j)

          IF( zscl(j) < rough%disp(j) ) THEN

             r_sc(j) = term5(j) * LOG(zscl(j)/rough%z0soilsn(j)) * &
                  ( EXP(2*C%CSW*canopy%rghlai(j)) - term1(j) ) / term3(j)

          ELSEIF( rough%disp(j) <= zscl(j) .AND. &
               zscl(j) < rough%hruff(j) ) THEN

             r_sc(j) = rough%rt0us(j) + term5(j) * ( term2(j) - term1(j) ) / &
                  term3(j)

          ELSEIF( rough%hruff(j) <= zscl(j) .AND. &
               zscl(j) <  rough%zruffs(j) ) THEN

             r_sc(j) = rough%rt0us(j) + rough%rt1usa(j) + term5(j) * &
                  ( zscl(j) - rough%hruff(j) ) / &
                  ( C%A33**2 * C%CTL * rough%hruff(j) )


          ELSEIF( zscl(j) >= rough%zruffs(j) ) THEN

             r_sc(j) = rough%rt0us(j) + rough%rt1usa(j) + rough%rt1usb(j) + &
                  ( LOG( (zscl(j) - rough%disp(j)) / &
                  MAX( rough%zruffs(j)-rough%disp(j), &
                  rough%z0soilsn(j) ) ) - psis( (zscl(j)-rough%disp(j)) &
                  / (rough%zref_tq(j)/canopy%zetar(j,iterplus) ) ) &
                  + psis( (rough%zruffs(j) - rough%disp(j) ) &
                  / (rough%zref_tq(j)/canopy%zetar(j,iterplus ) ) ) ) &
                  / C%VONK

          ENDIF

          canopy%tscrn(j) = ssnow%tss(j) + (met%tk(j) - ssnow%tss(j)) * &
               MIN(1.,r_sc(j) / MAX( 1., &
               rough%rt0us(j) + rough%rt1usa(j) + rough%rt1usb(j) &
               + rt1usc(j))) - C%tfrz
       ENDIF

    ENDDO

    CALL qsatfjh(rsts,canopy%tscrn,met%pmb)

    qtgnet = rsts * ssnow%wetfac - met%qv

    DO j=1,mp

       IF (qtgnet(j) .GT. 0. ) THEN
          qsurf(j) = rsts(j) * ssnow%wetfac(j)
       ELSE
          qsurf(j) = 0.1*rsts(j)*ssnow%wetfac(j) + 0.9*met%qv(j)
       ENDIF

       canopy%qscrn(j) = met%qv(j) - qstar(j) * ftemp(j)

       IF( canopy%vlaiw(j) >C%LAI_THRESH .and. rough%hruff(j) > 0.01) &

            canopy%qscrn(j) = qsurf(j) + (met%qv(j) - qsurf(j)) * MIN( 1., &
            r_sc(j) / MAX( 1., rough%rt0us(j) + &
            rough%rt1usa(j) + rough%rt1usb(j) + rt1usc(j) ) )

    ENDDO


    ! Calculate dewfall: from negative lh wet canopy + neg. lh dry canopy:
    canopy%dewmm = - (min(0.0,canopy%fevw) + min(0.0_r_2,canopy%fevc)) * &
         dels * 1.0e3 / (C%RHOW*air%rlam)

    ! Add dewfall to canopy water storage:
    canopy%cansto = canopy%cansto + canopy%dewmm

    ! Modify canopy water storage for evaporation:
    canopy%cansto = MAX(canopy%cansto-MAX(0.0,REAL(canopy%fevw))*dels &
         *1.0e3/(C%RHOW*air%rlam), 0.0)

    ! Calculate canopy water storage excess:
    canopy%spill=max(0.0, canopy%cansto-cansat)

    ! Move excess canopy water to throughfall:
    ! %through is /dels in UM app. (unpacked in hyd driver) for STASH output
    canopy%through = canopy%through + canopy%spill

    ! Initialise 'throughfall to soil' as 'throughfall from canopy';
    ! snow may absorb
    canopy%precis = max(0.,canopy%through)

    ! Update canopy storage term:
    canopy%cansto=canopy%cansto - canopy%spill

    ! Calculate the total change in canopy water store (mm/dels):
    canopy%delwc = canopy%cansto-canopy%oldcansto

    ! calculate dgdtg, derivative of ghflux 3 instances
    ! d(canopy%fns)/d(ssnow%tgg)
    ! d(canopy%fhs)/d(ssnow%tgg)
    ! d(canopy%fes)/d(dq)
    IF (cable_user%soil_struc=='default') THEN
       ssnow%dfn_dtg = (-1.)*4.*C%EMSOIL*C%SBOLTZ*tss4/ssnow%tss

       IF (cable_user%litter) THEN
          !!vh_js!!
          ssnow%dfh_dtg = air%rho*C%CAPP/(ssnow%rtsoil+ &
               real((1-ssnow%isflag))*veg%clitt*0.003/canopy%kthLitt/(air%rho*C%CAPP))
          ssnow%dfe_ddq = ssnow%wetfac*air%rho*air%rlam*ssnow%cls/ &
               (ssnow%rtsoil+ real((1-ssnow%isflag))*veg%clitt*0.003/canopy%DvLitt)

       ELSE
          ssnow%dfh_dtg = air%rho*C%CAPP/ssnow%rtsoil
          ssnow%dfe_ddq = ssnow%wetfac*air%rho*air%rlam*ssnow%cls/ssnow%rtsoil
       ENDIF

       ssnow%ddq_dtg = (C%rmh2o/C%rmair) /met%pmb * C%TETENA*C%TETENB * C%TETENC &
            / ( ( C%TETENC + ssnow%tss-C%tfrz )**2 )*EXP( C%TETENB * &
            ( ssnow%tss-C%tfrz ) / ( C%TETENC + ssnow%tss-C%tfrz ) )
       canopy%dgdtg = ssnow%dfn_dtg - ssnow%dfh_dtg - ssnow%dfe_ddq * &
            ssnow%ddq_dtg
    ENDIF

    bal%drybal = REAL(ecy+hcy) - SUM(rad%rniso,2) &
         + C%CAPP*C%rmair*(tlfy-met%tk)*SUM(rad%gradis,2)  ! YP nov2009

    bal%wetbal = canopy%fevw + canopy%fhvw - SUM(rad%rniso,2) * canopy%fwet &
         + C%CAPP*C%rmair * (tlfy-met%tk) * SUM(rad%gradis,2) * &
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

      use cable_def_types_mod, only : mp

      implicit none

      real, dimension(mp) :: lower_limit, rescale
      real, dimension(mp) :: psim_1, psim_2, z_eff, psim_arg

      psim_1 = psim(canopy%zetar(:,iter))

      rescale  = C%vonk * max(met%ua,C%umin)
      z_eff    = rough%zref_uv / rough%z0m
      psim_arg = canopy%zetar(:,iter) / z_eff
      !---fix for compiler limitation. bitwise reproducable whilst we
      !---we know it to 11th decimal. psim_arg typically of a few
      !psim_arg = nint(psim_arg * 1.e11)*1.e-11
      psim_2   = psim(psim_arg)

      lower_limit = rescale / ( log(z_eff) - psim_1 + psim_2 )
      canopy%us   = max(1.0e-6, lower_limit )

    END SUBROUTINE comp_friction_vel

    
    ! ------------------------------------------------------------------------------

    
    FUNCTION Penman_Monteith(ground_H_flux) RESULT(ssnowpotev)
      
      USE cable_def_types_mod, only : mp

      implicit none
      
      real, dimension(mp), intent(in)  :: ground_H_flux
      real, dimension(mp)              :: ssnowpotev ! returned result of function

      real, dimension(mp) :: &
           sss,             & ! var for Penman-Monteith soil evap
           cc1,             & ! var for Penman-Monteith soil evap
           cc2,             & ! var for Penman-Monteith soil evap
           qsatfvar           !

      ! Penman-Monteith formula
      sss = air%dsatdk
      cc1 = sss / (sss + air%psyc)
      cc2 = air%psyc / (sss + air%psyc)

      call qsatfjh(qsatfvar, met%tvair-c%tfrz, met%pmb)

      if (cable_user%litter) then
         !! vh_js !!
         ssnowpotev = cc1 * (canopy%fns - ground_H_flux) + &
              cc2 * air%rho * air%rlam*(qsatfvar - met%qvair) / &
              (ssnow%rtsoil + real((1-ssnow%isflag))*veg%clitt*0.003/canopy%DvLitt)
      else
         ssnowpotev = cc1 * (canopy%fns - ground_h_flux) + &
              cc2 * air%rho * air%rlam*(qsatfvar - met%qvair)/ssnow%rtsoil
      endif

    END FUNCTION Penman_Monteith


    ! ------------------------------------------------------------------------------

    
    ! method alternative to P-M formula above
    FUNCTION humidity_deficit_method(dq, qstss) RESULT(ssnowpotev)

      use cable_def_types_mod, only: mp

      implicit none

      real, dimension(mp), intent(in) :: &
           dq,   & ! sat spec hum diff.
           qstss   ! dummy var for compilation
      real, dimension(mp)             :: ssnowpotev

      real, dimension(mp) :: idq
      integer :: j

      idq = dq
      do j=1, mp
         if ((ssnow%snowd(j)>1.0) .or. (ssnow%tgg(j,1).eq.c%tfrz)) &
              idq(j) = max(-0.1e-3, idq(j))
      enddo

      if (cable_user%litter) then
         !! vh_js !!
         ssnowpotev = air%rho * air%rlam * idq / &
              (ssnow%rtsoil + real((1-ssnow%isflag))*veg%clitt*0.003/canopy%DvLitt)
      else
         ssnowpotev = air%rho * air%rlam * idq / ssnow%rtsoil
      endif

    END FUNCTION humidity_deficit_method

    
    ! ------------------------------------------------------------------------------

    
    SUBROUTINE Latent_heat_flux()

      use cable_common_module
      use cable_def_types_mod, only : mp

      implicit none
      
      real(r_2), dimension(mp) :: frescale, flower_limit, fupper_limit, pwet

      integer :: j

      ! Soil latent heat:
      canopy%fess = real(ssnow%wetfac * ssnow%potev, r_2)
      where (ssnow%potev < 0.0) canopy%fess = real(ssnow%potev, r_2)

      ! Reduce soil evap due to presence of puddle
      pwet = max(0.0_r_2, min(0.2_r_2, &
           real(ssnow%pudsto,r_2) / max(1.0_r_2, real(ssnow%pudsmx,r_2)) ) )
      canopy%fess = canopy%fess * (1.0_r_2-pwet)

      frescale = real(soil%zse(1) * 1000. * air%rlam / dels, r_2)

      do j=1, mp

         if ((ssnow%snowd(j) < 0.1) .and. (canopy%fess(j) .gt. 0.0_r_2)) then
            if (.not. cable_user%l_new_reduce_soilevp) then
               flower_limit(j) = ssnow%wb(j,1) - real(soil%swilt(j),r_2)/2.0_r_2
            else
               ! E.Kowalczyk 2014 - reduces the soil evaporation
               flower_limit(j) = ssnow%wb(j,1) - real(soil%swilt(j),r_2)
            endif
            
            fupper_limit(j) = max(0.0_r_2, &
                 flower_limit(j) * frescale(j) &
                 - real(ssnow%evapfbl(j,1)*air%rlam(j)/dels, r_2))
            canopy%fess(j) = min(canopy%fess(j), fupper_limit(j))

            fupper_limit(j) = (ssnow%wb(j,1) - real(ssnow%wbice(j,1),r_2)) * frescale(j)
            canopy%fess(j) = min(canopy%fess(j), fupper_limit(j))
         end if

         ssnow%cls(j) = 1.

         if ((ssnow%snowd(j) >= 0.1) .and. (ssnow%potev(j) > 0.0)) then
            ssnow%cls(j) = 1.1335
            canopy%fess(j) = real(min( (ssnow%wetfac(j)*ssnow%potev(j))*ssnow%cls(j), &
                 ssnow%snowd(j)/dels*air%rlam(j)*ssnow%cls(j) ), r_2)
         endif

      enddo

      ! Evaporation from soil puddle
      canopy%fesp = min(real(ssnow%pudsto/dels*air%rlam, r_2), &
           max(pwet*real(ssnow%potev,r_2), 0.0_r_2))
      canopy%fes  = canopy%fess + canopy%fesp

    END SUBROUTINE latent_heat_flux

    
    ! -----------------------------------------------------------------------------

    
    SUBROUTINE within_canopy(gbhu, gbhf)

      USE cable_def_types_mod, only : mp, r_2

      REAL(r_2), INTENT(IN), DIMENSION(:,:) :: &
           gbhu, & ! forcedConvectionBndryLayerCond
           gbhf    ! freeConvectionBndryLayerCond

      REAL, DIMENSION(mp) :: &
           rrsw, & ! recipr. stomatal resistance for water
           rrbw, & ! recipr. leaf boundary layer resistance for water
           dmah, & ! A_{H} in eq. 3.41 in SCAM, CSIRO tech report 132
           dmbh, & ! B_{H} in eq. 3.41 in SCAM, CSIRO tech report 132
           dmch, & ! C_{H} in eq. 3.41 in SCAM, CSIRO tech report 132
           dmae, & ! A_{E} in eq. 3.41 in SCAM, CSIRO tech report 132
           dmbe, & ! B_{E} in eq. 3.41 in SCAM, CSIRO tech report 132
           dmce    ! C_{E} in eq. 3.41 in SCAM, CSIRO tech report 132
      REAL    :: lower_limit, upper_limit
      INTEGER :: j

      rrbw = real(sum(gbhu+gbhf,2)) / air%cmolar  ! MJT

      ! leaf stomatal resistance for water
      rrsw = sum(canopy%gswx,2) / air%cmolar ! MJT

      do j=1, mp

         if ( (veg%meth(j) > 0.0) .and. (canopy%vlaiw(j) > C%lai_thresh) .and. &
              (rough%hruff(j) > rough%z0soilsn(j)) ) then

            !   use the dispersion matrix (DM) to find the air temperature
            !   and specific humidity
            !   (Raupach, Finkele and Zhang 1997, pp 17)
            ! leaf boundary layer resistance for water
            ! A_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132
            dmah(j) = (rt0(j)+rough%rt1(j))*((1.+air%epsi(j))*rrsw(j) + rrbw(j)) &
                 + air%epsi(j) * (rt0(j)*rough%rt1(j))*(rrbw(j)*rrsw(j))

            ! B_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132
            dmbh(j) = (-air%rlam(j)/C%CAPP)*(rt0(j)*rough%rt1(j))*(rrbw(j)*rrsw(j))

            ! C_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132
            dmch(j) = ((1.+air%epsi(j))*rrsw(j) + rrbw(j))*rt0(j)*rough%rt1(j)* &
                 (canopy%fhv(j) + canopy%fhs(j))/(air%rho(j)*C%CAPP)

            ! A_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
            dmae(j) = (-air%epsi(j)*C%CAPP/air%rlam(j))*(rt0(j)*rough%rt1(j)) * &
                 (rrbw(j)*rrsw(j))

            ! B_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
            dmbe(j) = ( rt0(j) + ssnow%wetfac(j) * rough%rt1(j) ) * &
                 ( (1.+air%epsi(j) ) * rrsw(j) + rrbw(j) ) + &
                 ( rt0(j) * rough%rt1(j) ) * ( rrbw(j) * rrsw(j) )

            ! C_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
            dmce(j) = ((1.+air%epsi(j))*rrsw(j) + rrbw(j))*rt0(j)*rough%rt1(j)* &
                 (canopy%fev(j) + real(canopy%fes(j)))/(air%rho(j)*air%rlam(j))

            ! Within canopy air temperature:
            met%tvair(j) = met%tk(j) + ( dmbe(j) * dmch(j) - dmbh(j) * dmce(j) ) &
                 / (dmah(j)*dmbe(j)-dmae(j)*dmbh(j)+1.0e-12)


            !---set limits for comparisson
            lower_limit =  min(ssnow%tss(j), met%tk(j)) - 5.0
            upper_limit =  max(ssnow%tss(j), met%tk(j)) + 5.0

            !--- tvair within these limits
            met%tvair(j) = max(met%tvair(j), lower_limit)
            met%tvair(j) = min(met%tvair(j), upper_limit)

            ! recalculate using canopy within temperature
            met%qvair(j) = met%qv(j) + (dmah(j)*dmce(j)-dmae(j)*dmch(j)) / &
                 ( dmah(j)*dmbe(j)-dmae(j)*dmbh(j)+1.0e-12)
            met%qvair(j) = max(0.0, met%qvair(j))

            ! isothermal test vh
            if (cable_user%within_canopy_isothermal) then
               met%qvair(j) = met%qv(j)
               met%tvair(j) = met%tk(j)
            endif

            !---set limits for comparisson
            lower_limit =  min(ssnow%qstss(j), met%qv(j))
            upper_limit =  max(ssnow%qstss(j), met%qv(j))

            !--- qvair within these limits
            met%qvair(j) = max(met%qvair(j), lower_limit)
            met%qvair(j) = min(met%qvair(j), upper_limit)

           ! Saturated specific humidity in canopy:
            call qsatfjh2(qstvair(j), met%tvair(j)-C%tfrz, met%pmb(j))

            ! Saturated vapour pressure deficit in canopy:
            met%dva(j) = (qstvair(j) - met%qvair(j)) *  C%rmair/C%RMH2O &
                 * met%pmb(j) * 100.
         endif

      enddo

    END SUBROUTINE within_canopy

    
    ! -----------------------------------------------------------------------------

    SUBROUTINE update_zetar()

      INTEGER :: j

      ! monin-obukhov stability parameter zetar=zref/l
      ! recompute zetar for the next iteration, except on last iteration
      IF (iter < NITER) THEN ! dont compute zetar on the last iter

         iterplus = MAX(iter+1,2)
         canopy%zetar(:,iterplus) = -( C%VONK * C%GRAV * rough%zref_tq * &
              ( canopy%fh + 0.07 * canopy%fe ) ) / &
              ( air%rho * C%CAPP * met%tk * canopy%us**3 )
         ! stability parameter at shear height: needed for Harman in-canopy stability correction
         IF (cable_user%soil_struc=='sli') THEN
            WHERE (canopy%vlaiw > C%LAI_THRESH .and. rough%hruff > rough%z0soilsn)
               canopy%zetash(:,iterplus) = -(C%vonk*C%grav*(0.1*rough%hruff)*(canopy%fhs+0.07*real(canopy%fes)))/ &
                    max( (air%rho*C%capp*met%tk*(canopy%us*rough%term6a)**3), 1.e-12)
            ELSEWHERE
               canopy%zetash(:,iterplus) = canopy%zetash(:,iter)
            ENDWHERE
         ENDIF

         ! case NITER=2: final zetar=C%ZETmul*zetar(2) (compute only when iter=1)
         IF (NITER == 2) THEN

            canopy%zetar(:,2) = C%ZETmul * canopy%zetar(:,2)

            ! stability parameter at shear height: needed for Harman in-canopy stability correction
            IF (cable_user%soil_struc=='sli') then
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
         if (cable_user%soil_struc=='sli') &
              canopy%zetash(:,iterplus) = MIN(C%ZETPOS,canopy%zetash(:,iterplus))

         ! zetar too -
         canopy%zetar(:,iterplus) = MAX(C%ZETNEG,canopy%zetar(:,iterplus))
         if (cable_user%soil_struc=='sli') &
              canopy%zetash(:,iterplus) = MAX(C%ZETNEG,canopy%zetash(:,iterplus))

      END IF ! (iter < NITER)

    END SUBROUTINE update_zetar

    
    ! -----------------------------------------------------------------------------

    
    FUNCTION qsatf(j,tair,pmb) RESULT(r)
      ! MRR, 1987
      ! AT TEMP tair (DEG C) AND PRESSURE pmb (MB), GET SATURATION SPECIFIC
      ! HUMIDITY (KG/KG) FROM TETEN FORMULA

      REAL, INTENT(IN) :: &
           tair,         & ! air temperature (C)
           pmb             ! pressure PMB (mb)

      INTEGER, INTENT(IN) :: j

      REAL           :: r    ! result; sat sp humidity

      r = (C%RMH2o/C%rmair) * (C%TETENA*EXP(C%TETENB*tair/(C%TETENC+tair))) / pmb
    END FUNCTION qsatf

    
    ! -----------------------------------------------------------------------------

    
    SUBROUTINE qsatfjh(var, tair, pmb)
      
      USE cable_def_types_mod, only: mp

      implicit none
      
      REAL, INTENT(OUT), DIMENSION(mp) :: var ! result; sat sp humidity
      REAL, INTENT(IN),  DIMENSION(mp) :: &
           tair, & ! air temperature (C)
           pmb     ! pressure PMB (mb)

      INTEGER :: j

      DO j=1,mp
         var(j) = (C%RMH2O/C%rmair) * (C%TETENA*EXP(C%TETENB*tair(j)/(C%TETENC+tair(j)))) &
              / pmb(j)
      ENDDO

    END SUBROUTINE qsatfjh

    
    ! -----------------------------------------------------------------------------

    
    SUBROUTINE qsatfjh2(var, tair, pmb)

      implicit none

      REAL, INTENT(OUT) :: var ! result; sat sp humidity
      REAL, INTENT(IN)  :: &
           tair, & ! air temperature (C)
           pmb     ! pressure PMB (mb)

      var = (C%RMH2O/C%rmair) * (C%TETENA*EXP(C%TETENB*tair/(C%TETENC+tair))) / pmb

    END SUBROUTINE qsatfjh2

    
    ! -----------------------------------------------------------------------------

    
    FUNCTION psim(zeta) RESULT(r)
      ! mrr, 16-sep-92 (from function psi: mrr, edinburgh 1977)
      ! computes integrated stability function psim(z/l) (z/l=zeta)
      ! for momentum, using the businger-dyer form for unstable cases
      ! and the Beljaars and Holtslag (1991) form for stable cases.
      
      use cable_def_types_mod, only : mp
      use mo_constants, only: pi => pi_dp

      implicit none

      real, dimension(mp), intent(in) :: zeta !
      real, dimension(mp)             :: r    ! function result

      real, dimension(mp) :: &
           x,       & !
           z,       & !
           stable,  & !
           unstable   !

      real, parameter :: &
           gu = 16.0 ! ,  & !
           ! gs = 5.0
      real, parameter :: &
           a  = 1.0,   & !
           b  = 0.667, & !
           xc = 5.0,   & !
           d  = 0.35     !

      z = 0.5 + sign(0.5,zeta) ! z=1 in stable, 0 in unstable

      ! Beljaars and Holtslag (1991) for stable
      stable   = -a*zeta - b*(zeta-xc/d)*exp(-d*zeta) - b*xc/d
      x        = (1.0 + gu*abs(zeta))**0.25
      unstable = log((1.0+x*x) * (1.0+x)**2 / 8.) - 2.0*atan(x) + pi*0.5
      r        = z*stable + (1.0-z)*unstable

    END FUNCTION psim

    
    ! -----------------------------------------------------------------------------

    ELEMENTAL FUNCTION psis(zeta) RESULT(r)
      ! mrr, 16-sep-92 (from function psi: mrr, edinburgh 1977)
      ! computes integrated stability function psis(z/l) (z/l=zeta)
      ! for scalars, using the businger-dyer form for unstable cases
      ! and the webb form for stable cases. see paulson (1970).

      implicit none

      real, intent(in) :: zeta
      real             :: r

      REAL, PARAMETER :: &
           gu = 16.0,  & !
           ! gs = 5.0,   & !
           a = 1.0,    & !
           b = 0.667,  & !
           c = 5.0,    & !
           d = 0.35

      REAL :: &
           stzeta,   & !
           z,        & !
           y,        & !
           stable,   & !
           unstable

      z = 0.5 + sign(0.5,zeta) ! z=1 in stable, 0 in unstable

      ! Beljaars and Holtslag (1991) for stable
      stzeta = max(0., zeta)
      stable = -(1.+2./3.*a*stzeta)**(3./2.) - &
           b*(stzeta-c/d)*exp(-d*stzeta) - b*c/d + 1.
      y      = (1.0 + gu*abs(zeta))**0.5
      unstable = 2.0 * alog((1.0+y)*0.5)
      r      = z*stable + (1.0-z)*unstable

    END FUNCTION psis

    
    ! -----------------------------------------------------------------------------

    
    ELEMENTAL FUNCTION rplant(rpconst, rpcoef, tair) result(z)

      implicit none
      
      real, intent(in) :: rpconst
      real, intent(in) :: rpcoef
      real, intent(in) :: tair
      real             :: z
      
      z = rpconst * exp(rpcoef * tair)
      
    END FUNCTION rplant
    

    ! -----------------------------------------------------------------------------

    SUBROUTINE wetLeaf( dels, rad, rough, air, met, veg, canopy, cansat, tlfy, &
         gbhu, gbhf, ghwet )

      USE cable_def_types_mod

      TYPE (radiation_type), INTENT(INOUT) :: rad
      TYPE (roughness_type), INTENT(INOUT) :: rough
      TYPE (air_type),       INTENT(INOUT) :: air
      TYPE (met_type),       INTENT(INOUT) :: met
      TYPE (canopy_type),    INTENT(INOUT) :: canopy

      TYPE (veg_parameter_type), INTENT(INOUT)    :: veg

      REAL,INTENT(IN), DIMENSION(:) :: &
           tlfy,          & ! leaf temp (K) - assC%UMINg the temperature of
                                ! wet leaf is equal that of dry leaf ="tlfy"
           cansat           ! max canopy intercept. (mm)

      REAL(r_2), INTENT(IN), DIMENSION(:,:) :: &
           gbhu,          & ! forcedConvectionBndryLayerCond
           gbhf             ! freeConvectionBndryLayerCond

      REAL(r_2), INTENT(OUT), DIMENSION(:) :: &
           ghwet            ! cond for heat for a wet canopy

      REAL, INTENT(IN)     :: dels ! integration time step (s)

      ! local variables
      REAL, DIMENSION(mp) :: &
           ccfevw,        & ! limitation term for
           gwwet,         & ! cond for water for a wet canopy
           ghrwet           ! wet canopy cond: heat & thermal rad

      !i sums, terms of convenience/readability
      REAL, DIMENSION(mp) :: &
           sum_gbh, sum_rad_rniso, sum_rad_gradis

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
            canopy%fwet(j) = MAX( 0.0, MIN( 1.0, &
                 0.8 * canopy%cansto(j) / MAX( cansat(j), 0.01 ) ) )


            ! Calculate lat heat from wet canopy, may be neg. if dew on wet canopy
            ! to avoid excessive evaporation:
            ccfevw(j) = MIN(canopy%cansto(j) * air%rlam(j) / dels, &
                 2.0 / (1440.0 / (dels/60.0)) * air%rlam(j) )

            canopy%fevw(j) = MIN( canopy%fwet(j) * ( air%dsatdk(j) * &
                 ( sum_rad_rniso(j)- C%CAPP*C%rmair*( met%tvair(j) &
                 - met%tk(j) ) * sum_rad_gradis(j) ) &
                 + C%CAPP * C%rmair * met%dva(j) * ghrwet(j) ) &
                 / ( air%dsatdk(j)+air%psyc(j)*ghrwet(j) / gwwet(j) ) &
                 , ccfevw(j) )

            canopy%fevw_pot(j) = ( air%dsatdk(j)* (sum_rad_rniso(j) - &
                 C%CAPP * C%rmair * ( met%tvair(j) - met%tk(j) ) &
                 *sum_rad_gradis(j) ) &
                 + C%CAPP * C%rmair * met%dva(j) * ghrwet(j)) &
                 / (air%dsatdk(j)+air%psyc(j)*ghrwet(j)/gwwet(j) )

            ! calculate sens heat from wet canopy:
            canopy%fhvw(j) = canopy%fwet(j) * ( sum_rad_rniso(j) -C%CAPP * C%rmair &
                 * ( tlfy(j) - met%tk(j) ) * sum_rad_gradis(j) ) &
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

    TYPE (veg_parameter_type), INTENT(INOUT)    :: veg
    TYPE (soil_snow_type), intent(inout):: ssnow
    TYPE (soil_parameter_type), intent(inout)   :: soil
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
    upper_limit = max(lower_limit, 0.0)
    canopy%wcint = MERGE( MIN( upper_limit, ftemp ), 0.0, &
         ftemp > 0.0  .AND. met%tk > C%tfrz)  !EAK, 09/10



    ! Define canopy throughfall (100% of precip if temp < 0C, see above):
    canopy%through = met%precip_sn + MIN( met%precip - met%precip_sn , &
         MAX( 0.0, met%precip - met%precip_sn - canopy%wcint) )


    ! Add canopy interception to canopy storage term:
    canopy%cansto = canopy%cansto + canopy%wcint

    ! Calculate fraction of canopy which is wet:
    canopy%fwet   = MAX( 0.0, MIN( 0.9, 0.8 * canopy%cansto / &
         MAX( cansat, 0.01 ) ) )

    ssnow%wetfac = MAX( 1.e-6, MIN( 1.0, &
         ( REAL (ssnow%wb(:,1) ) - soil%swilt/ 2.0 ) &
         / ( soil%sfc - soil%swilt/2.0 ) ) )

    DO j=1,mp

       IF ( ssnow%wbice(j,1) > real(tiny(1.0),r_2) ) &
            ssnow%wetfac(j) = ssnow%wetfac(j) * real(MAX( 0.5_r_2, 1._r_2 - MIN( 0.2_r_2, &
            ( ssnow%wbice(j,1) / ssnow%wb(j,1) )**2 ) ) )

       IF( ssnow%snowd(j) > 0.1) ssnow%wetfac(j) = 0.9

       IF ( veg%iveg(j) == 16 .and. met%tk(j) >= C%tfrz + 5. ) &
            ssnow%wetfac(j) = 1.0 ! lakes: hard-wired number to be removed

       IF( veg%iveg(j) == 16 .and. met%tk(j) < C%tfrz + 5. ) &
            ssnow%wetfac(j) = 0.7 ! lakes: hard-wired number to be removed

    ENDDO


    ! owetfac introduced to reduce sharp changes in dry regions,
    ! especially in offline runs in which there may be discrepancies b/n
    ! timing of precip and temperature change (EAK apr2009)
    ssnow%wetfac = 0.5*(ssnow%wetfac + ssnow%owetfac)

  END SUBROUTINE Surf_wetness_fact

  ! -----------------------------------------------------------------------------

  SUBROUTINE dryLeaf( dels, rad, rough, air, met, &
       veg, canopy, soil, ssnow, dsx, &
       fwsoil, tlfx,  tlfy,  ecy, hcy, &
       rny, gbhu, gbhf, csx, &
       cansat, ghwet, iter, climate, gmes )

    use cable_def_types_mod
    use cable_common_module

    implicit none

    real,                      intent(in)    :: dels ! integration time step (s)
    type(radiation_type),      intent(inout) :: rad
    type(roughness_type),      intent(inout) :: rough
    type(air_type),            intent(inout) :: air
    type(met_type),            intent(inout) :: met
    type(veg_parameter_type),  intent(inout) :: veg
    type(canopy_type),         intent(inout) :: canopy
    type(soil_parameter_type), intent(inout) :: soil
    type(soil_snow_type),      intent(inout) :: ssnow
    real,      dimension(:),   intent(inout) :: &
         dsx,        & ! leaf surface vpd
         fwsoil,     & ! soil water modifier of stom. cond
         tlfx,       & ! leaf temp prev. iter (K)
         tlfy          ! leaf temp (K)
    real(r_2), dimension(:),   intent(inout) :: &
         ecy,        & ! lat heat fl dry big leaf
         hcy,        & ! veg. sens heat
         rny         !& !
    real(r_2), dimension(:,:), intent(inout) :: &
         gbhu,       & ! forcedConvectionBndryLayerCond
         gbhf,       & ! freeConvectionBndryLayerCond
         csx           ! leaf surface CO2 concentration
    real,      dimension(:),   intent(in)    :: cansat
    real(r_2), dimension(:),   intent(out)   :: ghwet  ! cond for heat for a wet canopy
    integer,                   intent(in)    :: iter
    type(climate_type),        intent(in)    :: climate
    real(r_2), dimension(:,:), intent(inout) :: gmes             ! mesophyll conductance

    !local variables
    real, parameter  :: &
         jtomol = 4.6e-6  ! Convert from J to Mol for light

    real, dimension(mp) :: &
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
         temp,          & !
         temp_sun,      & !
         temp_shade       !

    real(r_2), dimension(mp)  :: &
         ecx,        & ! lat. hflux big leaf
         hcx,        & ! sens heat fl big leaf prev iteration
         rnx           ! net rad prev timestep

    real, dimension(mp,ms)  :: oldevapfbl

    real, dimension(mp,mf)  :: &
         gw,         & ! cond for water for a dry canopy
         gh,         & ! cond for heat for a dry canopy
         ghr,        & ! dry canopy cond for heat & thermal rad
         anx,        & ! net photos. prev iteration
         an_y,       & ! net photosynthesis soln
         rdx3,       & ! daytime leaf resp rate, prev iteration
         rdx4,       & ! daytime leaf resp rate, prev iteration
         rdx,        & ! daytime leaf resp rate, prev iteration
         rdy,        & ! daytime leaf resp rate
         ejmxt3,     & ! jmax big leaf C3 plants
         vcmxt3,     & ! vcmax big leaf C3
         vcmxt4,     & ! vcmax big leaf C4
         vx3,        & ! carboxylation C3 plants
         vx4,        & ! carboxylation C4 plants
         co2cp,      & ! CO2 compensation point (needed for Leuing stomatal conductance)
         ! Ticket #56, xleuning is no longer used, we replace it with gs_coeff,
         ! which is computed differently based on the new GS_SWITCH. If GS_SWITCH
         ! is "leuning", it's the same, if "medlyn", then the new Medlyn model
         ! xleuning,   & ! leuning stomatal coeff
         gs_coeff,   & ! stom coeff, Ticket #56
         psycst,     & ! modified pych. constant
         frac42,     & ! 2D frac4
         temp2,      &
         anrubiscox, & ! net photosynthesis (rubisco limited)
         anrubpx,    & ! net photosynthesis (rubp limited)
         ansinkx,    & ! net photosynthesis (sink limited)
         anrubiscoy, & ! net photosynthesis (rubisco limited)
         anrubpy       ! net photosynthesis (rubp limited)

   real(r_2), dimension(mp,mf)  ::  dAnrubiscox, & ! CO2 elasticity of net photosynthesis (rubisco limited)
         dAnrubpx, &   !  (rubp limited)
         dAnsinkx, &  ! (sink limited)
         dAnx, & !(actual rate)
         dAnrubiscoy, & ! CO2 elasticity of net photosynthesis (rubisco limited)
         dAnrubpy, &   !  (rubp limited)
         dAnsinky, &  ! (sink limited)
         dAn_y, & !(actual rate)
         eta_y, &
         eta_x

    real ::  gam0,    &
             conkc0,  &
             conko0,  &
             egam,    &
             ekc,     &
             eko

    real, dimension(:,:), pointer :: gswmin => null() ! min stomatal conductance

    real, dimension(mp,2) ::  gsw_term, lower_limit2  ! local temp var

    integer :: i, k, kk  ! iteration count
    real :: vpd, g1 ! Ticket #56

    ! END header

    allocate(gswmin(mp,mf))

    ! Soil water limitation on stomatal conductance:
    if (iter==1) then
       if ((cable_user%soil_struc=='default') .and. (cable_user%fwsoil_switch.ne.'Haverd2013')) then
          if (cable_user%fwsoil_switch == 'standard') then
             call fwsoil_calc_std(fwsoil, soil, ssnow, veg)
          elseif (cable_user%fwsoil_switch == 'non-linear extrapolation') then
             !EAK, 09/10 - replace linear approx by polynomial fitting
             call fwsoil_calc_non_linear(fwsoil, soil, ssnow, veg)
          elseif (cable_user%fwsoil_switch == 'Lai and Ktaul 2000') then
             call fwsoil_calc_Lai_Ktaul(fwsoil, soil, ssnow, veg)
          else
             stop 'fwsoil_switch failed.'
          endif
          canopy%fwsoil = fwsoil
       elseif ((cable_user%soil_struc=='sli') .or. (cable_user%fwsoil_switch=='Haverd2013')) then
          fwsoil = canopy%fwsoil
       endif
    endif

    ! weight min stomatal conductance by C3 an C4 plant fractions
    frac42       = spread(veg%frac4, 2, mf) ! frac C4 plants
    gsw_term     = spread(veg%gswmin, 2, mf)
    lower_limit2 = rad%scalex * gsw_term
    gswmin       = max(1.e-6, lower_limit2)

    gw          = 1.0e-3 ! default values of conductance
    gh          = 1.0e-3
    ghr         = 1.0e-3
    rdx         = 0.0
    anx         = 0.0
    anrubiscox  = 0.0
    anrubpx     = 0.0
    ansinkx     = 0.0
    dAnrubiscox = 0.0
    dAnrubpx    = 0.0
    dAnx        = 0.0
    eta_x       = 0.0
    rnx         = SUM(rad%rniso,2)
    abs_deltlf  = 999.0

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
             canopy%fwet(i) = MAX( 0.0, MIN( 1.0, 0.8 * canopy%cansto(i)/ MAX( &
                  cansat(i),0.01 ) ) )

             ! Calculate lat heat from wet canopy, may be neg.
             ! if dew on wet canopy to avoid excessive evaporation:
             ccfevw(i) = MIN(canopy%cansto(i) * air%rlam(i) / dels, &
                  2.0 / (1440.0 / (dels/60.0)) * air%rlam(i) )

             ! Grashof number (Leuning et al, 1995) eq E4:
             gras(i) = MAX(1.0e-6, &
                  1.595E8* ABS( tlfx(i)-met%tvair(i))* (veg%dleaf(i)**3.0) )

             ! See Appendix E in (Leuning et al, 1995):
             gbhf(i,1) = rad%fvlai(i,1) * air%cmolar(i) * 0.5*C%dheat &
                  * ( gras(i)**0.25 ) / veg%dleaf(i)
             gbhf(i,2) = rad%fvlai(i,2) * air%cmolar(i) * 0.5 * C%dheat &
                  * ( gras(i)**0.25 ) / veg%dleaf(i)
             gbhf(i,:) = MAX( 1.e-6_r_2, gbhf(i,:) )

             ! Conductance for heat:
             gh(i,:) = 2.0 * (gbhu(i,:) + gbhf(i,:))

             if (cable_user%amphistomatous) then
                gbhf(i,:) = 2.0*gbhf(i,:)
                gh(i,:) =  (gbhu(i,:) + gbhf(i,:))
             endif

             ! Conductance for heat and longwave radiation:
             ghr(i,:) = rad%gradis(i,:)+gh(i,:)


             
             ! Choose Ci-based (implicit gm) or Cc-based (explicit gm) parameters
             if (cable_user%explicit_gm) then
                gam0    = C%gam0cc
                conkc0  = C%conkc0cc
                conko0  = C%conko0cc
                egam    = C%egamcc
                ekc     = C%ekccc
                eko     = C%ekocc
             else
                gam0   = C%gam0
                conkc0 = C%conkc0
                conko0 = C%conko0
                egam   = C%egam
                ekc    = C%ekc
                eko    = C%eko
             endif
             ! JK: veg%vcmax_sun and veg%vcmax_shade are Cc-based if cable_user%explicit_gm = TRUE
             !     and Ci-based otherwise. If cable_user%coordinate_photosyn = FALSE,
             !     veg%vcmax_sun = veg%vcmax_shade = veg%vcmaxcc if cable_user%explicit_gm = TRUE
             !     and veg%vcmax_sun = veg%vcmax_shade = veg%vcmax otherwise.
             !     See also Subroutine 'casa_feedback' in casa_cable.F90. Same applies to ejmax.
             
             
             ! Leuning 2002 (PCE) equation for temperature response
             ! used for Vcmax for C3 plants:
             if (.not.cable_user%acclimate_photosyn) then
                if (cable_user%perturb_biochem_by_T) then
                   temp_sun(i)   = xvcmxt3(tlfx(i)+cable_user%Ta_perturbation) * &
                        veg%vcmax_sun(i) * (1.0-veg%frac4(i))
                   temp_shade(i) = xvcmxt3(tlfx(i)+cable_user%Ta_perturbation) * &
                        veg%vcmax_shade(i) * (1.0-veg%frac4(i))
                else
                   temp_sun(i)   = xvcmxt3(tlfx(i)) * veg%vcmax_sun(i) * (1.0-veg%frac4(i))
                   temp_shade(i) = xvcmxt3(tlfx(i)) * veg%vcmax_shade(i) * (1.0-veg%frac4(i))
                endif
             else
                if (cable_user%perturb_biochem_by_T) then
                   call xvcmxt3_acclim(tlfx(i)+cable_user%Ta_perturbation, climate%mtemp(i) , temp(i))
                else
                   call xvcmxt3_acclim(tlfx(i), climate%mtemp(i) , temp(i))
                endif
                temp_sun(i)   = temp(i) * veg%vcmax_sun(i) * (1.0-veg%frac4(i))
                temp_shade(i) = temp(i) * veg%vcmax_shade(i) * (1.0-veg%frac4(i))
             endif
             vcmxt3(i,1) = rad%scalex(i,1) * temp_sun(i)
             vcmxt3(i,2) = rad%scalex(i,2) * temp_shade(i)

             ! Temperature response Vcmax, C4 plants (Collatz et al 1989):
             temp_sun(i)   = xvcmxt4(tlfx(i)-C%tfrz) * veg%vcmax_sun(i) * veg%frac4(i)
             temp_shade(i) = xvcmxt4(tlfx(i)-C%tfrz) * veg%vcmax_shade(i) * veg%frac4(i)
             vcmxt4(i,1) = rad%scalex(i,1) * temp_sun(i)
             vcmxt4(i,2) = rad%scalex(i,2) * temp_shade(i)

             
             ! Leuning 2002 (PCE) equation for temperature response
             ! used for Jmax for C3 plants:
             if (.not.cable_user%acclimate_photosyn) then
                if (cable_user%perturb_biochem_by_T) then
                   temp_sun(i)   = xejmxt3(tlfx(i)+cable_user%Ta_perturbation) * &
                        veg%ejmax_sun(i) * (1.0-veg%frac4(i))
                   temp_shade(i) = xejmxt3(tlfx(i)+cable_user%Ta_perturbation) * &
                        veg%ejmax_shade(i) * (1.0-veg%frac4(i))
                else
                   temp_sun(i)   = xejmxt3(tlfx(i)) * &
                        veg%ejmax_sun(i) * (1.0-veg%frac4(i))
                   temp_shade(i) = xejmxt3(tlfx(i)) * &
                        veg%ejmax_shade(i) * (1.0-veg%frac4(i))
                endif
             else
                if (cable_user%perturb_biochem_by_T) then
                   call xejmxt3_acclim(tlfx(i)+cable_user%Ta_perturbation, &
                        climate%mtemp(i), climate%mtemp_max20(i), temp(i))
                else
                   call xejmxt3_acclim(tlfx(i), climate%mtemp(i), climate%mtemp_max20(i), temp(i))
                endif
                temp_sun(i)   = temp(i) * veg%ejmax_sun(i) * (1.0-veg%frac4(i))
                temp_shade(i) = temp(i) * veg%ejmax_shade(i) * (1.0-veg%frac4(i))
             endif
             ejmxt3(i,1) = rad%scalex(i,1) * temp_sun(i)
             ejmxt3(i,2) = rad%scalex(i,2) * temp_shade(i)
             

             ! Difference between leaf temperature and reference temperature:
             tdiff(i) = tlfx(i) - C%TREFK
             
             ! Michaelis menten constant of Rubisco for CO2:
             conkct(i) = conkc0 * EXP( ( ekc / (C%rgas*C%trefk) ) &
                                         * ( 1.0 - C%trefk/tlfx(i) ) )

             ! Michaelis menten constant of Rubisco for oxygen:
             conkot(i) = conko0 * EXP( ( eko / (C%rgas*C%trefk) ) &
                                         * ( 1.0 - C%trefk/tlfx(i) ) )

             ! Store leaf temperature
             tlfxx(i) = tlfx(i)

             ! "d_{3}" in Wang and Leuning, 1998, appendix E:
             cx1(i) = conkct(i) * (1.0+0.21/conkot(i))
             !cx2(i) = 2.0 * C%gam0 * ( 1.0 + C%gam1 * tdiff(i) &
             !                             + C%gam2 * tdiff(i) * tdiff(i) )
             cx2(i) = 2.0 * gam0 * EXP( ( egam / (C%rgas*C%trefk) ) &
                  * ( 1.0 - C%trefk/tlfx(i) ) )

             if (cable_user%perturb_biochem_by_T) then
                cx2(i) = 2.0 * gam0 * EXP( ( egam / (C%rgas*C%trefk) ) &
                     * ( 1.0 - C%trefk/(tlfx(i)+cable_user%Ta_perturbation) ) )

                
                conkot(i) = conko0 * EXP( ( eko / (C%rgas*C%trefk) ) &
                     * ( 1.0 - C%trefk/(tlfx(i)+cable_user%Ta_perturbation) ) )

                conkct(i) = conkc0 * EXP( ( ekc / (C%rgas*C%trefk) ) &
                     * ( 1.0 - C%trefk/(tlfx(i)+cable_user%Ta_perturbation) ) )

                cx1(i) = conkct(i) * (1.0+0.21/conkot(i))
             endif

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


              if (cable_user%CALL_climate) then
!!$                 if (veg%iveg(i).eq.1) then
!!$
!!$
!!$                    temp2(i,1) = rad%qcan(i,1,1) * jtomol * (1.0-veg%frac4(i))
!!$                    temp2(i,2) = rad%qcan(i,2,1) * jtomol * (1.0-veg%frac4(i))

!!$
!!$                    vx3(i,1)  = ej3x(temp2(i,1),climate%frec(i)*veg%alpha(i), &
!!$                         veg%convex(i),ejmxt3(i,1))
!!$                    vx3(i,2)  = ej3x(temp2(i,2),climate%frec(i)*veg%alpha(i), &
!!$                         veg%convex(i),ejmxt3(i,2))
!!$
!!$
!!$                    temp(i) =  xvcmxt3(tlfx(i)) * veg%vcmax(i) * (1.0-veg%frac4(i))
!!$
!!$                    vcmxt3(i,1) = rad%scalex(i,1) * temp(i) * climate%frec(i)
!!$                    vcmxt3(i,2) = rad%scalex(i,2) * temp(i) * climate%frec(i)
!!$                 endif


                 ! Atkin et al. 2015, Table S4,
                 ! modified by scaling factor to reduce leaf respiration to
                 ! expected proportion of GPP
                 !Broad-leaved trees: Rdark a25 =
                 !1.2818 + (0.0116 * Vcmax,a25) + (0.0334 * TWQ)
                 !C3 herbs/grasses: Rdark,a25 =
                 !1.6737 + (0.0116 * Vcmax,a25) + (0.0334 * TWQ)
                 !Needle-leaved trees: Rdark,a25 =
                 !1.2877 + (0.0116 * Vcmax,a25) + (0.0334 * TWQ)
                 !Shrubs: Rdark,a25 = 1.5758 + (0.0116 * Vcmax,a25) + (0.0334 * TWQ)

             if (veg%iveg(i).eq.2  ) then ! evergreen broadleaf forest

                rdx(i,1) = 1.0*(1.2818e-6+0.0116*veg%vcmax(i)- &
                     0.0334*climate%qtemp_max_last_year(i)*1e-6)

                rdx(i,2) = rdx(i,1)

             elseif ( veg%iveg(i).eq. 4  ) then ! decid broadleaf forest

                rdx(i,1) = 1.0*(1.2818e-6+0.0116*veg%vcmax(i)- &
                     0.0334*climate%qtemp_max_last_year(i)*1e-6)

                rdx(i,2) = rdx(i,1)

             elseif (veg%iveg(i).eq.1   ) then ! evergreen needleleaf forest
                 rdx(i,1) = 1.0*(1.2877e-6+0.0116*veg%vcmax(i)- &
                      0.0334*climate%qtemp_max_last_year(i)*1e-6)

                 rdx(i,2) = rdx(i,1)


              elseif ( veg%iveg(i).eq. 3  ) then ! decid needleleaf forest
                 rdx(i,1) = 1.0*(1.2877e-6+0.0116*veg%vcmax(i)- &
                      0.0334*climate%qtemp_max_last_year(i)*1e-6)

                 rdx(i,2) = rdx(i,1)

              elseif (veg%iveg(i).eq.6 .or. veg%iveg(i).eq.8 .or. &
                 veg%iveg(i).eq. 9  ) then ! C3 grass, tundra, crop
                 rdx(i,1) = 0.8*(1.6737e-6+0.0116*veg%vcmax(i)- &
                      0.0334*climate%qtemp_max_last_year(i)*1e-6)

                 rdx(i,2) = rdx(i,1)

              else  ! shrubs and other (C4 grass and crop)
                 rdx(i,1) = 0.7*(1.5758e-6+0.0116*veg%vcmax(i)- &
                      0.0334*climate%qtemp_max_last_year(i)*1e-6)
                 rdx(i,2) = rdx(i,1)
              endif
              veg%cfrd(i) = rdx(i,1) / veg%vcmax(i)
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
                 ! Rc = Rd 0 < Io < 10 umol quanta m-2 s-1
                 ! Rc = [0.5 * 0.05 ln(Io)] Rd Io > 10 umol quanta m-2 s-1

                 if (jtomol*1.0e6*rad%qcan(i,1,1).gt.10.0) &
                      rdx(i,1) = rdx(i,1) * &
                      (0.5 - 0.05*log(jtomol*1.0e6*rad%qcan(i,1,1)))

                 if (jtomol*1.0e6*rad%qcan(i,1,2).gt.10.0) &
                      rdx(i,2) = rdx(i,2) * &
                      (0.5 - 0.05*log(jtomol*1.0e6*rad%qcan(i,1,2)))

                 ! special for YP photosynthesis
                 rdx3(i,1) = rdx(i,1);
                 rdx3(i,2) = rdx(i,2);
                 rdx4(i,1) = rdx(i,1);
                 rdx4(i,2) = rdx(i,2);
             else !cable_user%call_climate

                 rdx(i,1) = (veg%cfrd(i)*vcmxt3(i,1) + veg%cfrd(i)*vcmxt4(i,1))
                 rdx(i,2) = (veg%cfrd(i)*vcmxt3(i,2) + veg%cfrd(i)*vcmxt4(i,2))
                 ! special for YP photosynthesis
                 rdx3(i,1) = veg%cfrd(i)*vcmxt3(i,1) 
                 rdx3(i,2) = veg%cfrd(i)*vcmxt3(i,2) 
                 rdx4(i,1) = veg%cfrd(i)*vcmxt4(i,1) 
                 rdx4(i,2) = veg%cfrd(i)*vcmxt4(i,2) 

             endif !cable_user%call_climate

             ! apply fwsoil to gmes
             ! if (any(gmes(i,:) < real(tiny(1.0),r_2) .and. gmes(i,:) > 0._r_2)) print*, 'CA03 ', gmes(i,:)
             ! if (fwsoil(i) < tiny(1.0)) print*, 'CA04 ', i, fwsoil(i)
             ! if (veg%gmmax(i) < tiny(1.0)) print*, 'CA05 ', veg%gmmax(i)
             gmes(i,1) = MAX(gmes(i,1) * real(fwsoil(i),r_2), real(0.15 * veg%gmmax(i),r_2))
             gmes(i,2) = MAX(gmes(i,2) * real(fwsoil(i),r_2), real(0.15 * veg%gmmax(i),r_2))

            ! Ticket #56 added switch for Belinda Medlyn's model
             IF (cable_user%GS_SWITCH == 'leuning') THEN

                ! vh: added calculation of CO2 compensation point
                if ((vcmxt3(i,1)+vcmxt4(i,1)).gt.1.0e-8) then
                   co2cp(i,1) = (cx2(i)/2.0 +  conkct(i) * (1.0 + 0.21/conkot(i))* &
                        rdx(i,1)/(vcmxt3(i,1)+vcmxt4(i,1)))/ &
                        (1.0 - rdx(i,1)/(vcmxt3(i,1)+vcmxt4(i,1)))
                else
                   co2cp(i,1) = 0.0
                endif
                
                if ((vcmxt3(i,2)+vcmxt4(i,2)).gt.1.0e-8) then
                   co2cp(i,2) = (cx2(i)/2.0 +  conkct(i) * (1.0 + 0.21/conkot(i))* &
                        rdx(i,2)/(vcmxt3(i,2)+vcmxt4(i,2)))/ &
                        (1.0 - rdx(i,2)/(vcmxt3(i,2)+vcmxt4(i,2)))
                else
                   co2cp(i,2) = 0.0
                endif


                ! vh set to zero for now as analytic solutions for Anet, gsc
                ! have not been developd with cs-co2cp in the denominator
                co2cp(i,1) = 0.0
                co2cp(i,2) = 0.0
                
               
                ! vh re-use g0 param here (not gswmin:just to reduce number of adjustable params)
                gswmin(i,1) = veg%g0(i)*rad%scalex(i,1)
                gswmin(i,2) = veg%g0(i)*rad%scalex(i,2)
                vpd =  dsx(i)

                ! vh re-write so that a1 and d0 are not correlated
                gs_coeff(i,1) = ( fwsoil(i) / ( csx(i,1) - co2cp(i,1) ) ) &
                     * ( 1.0 / ( 1.0/veg%a1gs(i) + (vpd/veg%d0gs(i))))

                gs_coeff(i,2) = ( fwsoil(i) / ( csx(i,2) - co2cp(i,2) ) ) &
                     * ( 1.0 / ( 1.0/veg%a1gs(i) + (vpd/veg%d0gs(i))))
                
!!$                gs_coeff(i,1) = ( fwsoil(i) / ( csx(i,1) - co2cp3 ) ) &
!!$                          * ( veg%a1gs(i) / ( 1.0 + dsx(i)/veg%d0gs(i)))
!!$
!!$                gs_coeff(i,2) = ( fwsoil(i) / ( csx(i,2) - co2cp3 ) ) &
!!$                          * ( veg%a1gs(i) / ( 1.0 + dsx(i)/veg%d0gs(i)))

            ! Medlyn BE et al (2011) Global Change Biology 17: 2134-2144.
            ELSEIF(cable_user%GS_SWITCH == 'medlyn') THEN
                 gswmin(i,1) = veg%g0(i)
                 gswmin(i,2) = veg%g0(i)

                IF (dsx(i) < 50.0) THEN
                    vpd  = 0.05 ! kPa
                ELSE
                    vpd = dsx(i) * 1E-03 ! Pa -> kPa
                    !vpd = met%dva(i) * 1E-03 ! Pa -> kPa
                END IF

                g1 = veg%g1(i)

                !gs_coeff(i,1) = (1.0* fwsoil(i) + (g1 * fwsoil(i)) / SQRT(vpd)) / csx(i,1)
                !gs_coeff(i,2) = (1.0* fwsoil(i) + (g1 * fwsoil(i)) / SQRT(vpd)) / csx(i,2)

                if (fwsoil(i) .LE. 0.05) then
                   gs_coeff(i,1) = (1.0* fwsoil(i) + (g1 * fwsoil(i)) / SQRT(vpd)) / csx(i,1)
                   gs_coeff(i,2) = (1.0* fwsoil(i) + (g1 * fwsoil(i)) / SQRT(vpd)) / csx(i,2)
                else
                   gs_coeff(i,1) = (1.0 + (g1 * fwsoil(i)) / SQRT(vpd)) / csx(i,1)
                   gs_coeff(i,2) = (1.0 + (g1 * fwsoil(i)) / SQRT(vpd)) / csx(i,2)
                endif

            ELSE
                STOP 'gs_model_switch failed.'
            ENDIF ! IF (cable_user%GS_SWITCH == 'leuning') THEN
!#endif

          ENDIF !IF (canopy%vlaiw(i) > C%LAI_THRESH .AND. abs_deltlf(i) > 0.1)

       ENDDO !i=1,mp


       if (cable_user%explicit_gm) then
          CALL photosynthesis_gm( csx(:,:), &
               spread(cx1(:),2,mf), &
               spread(cx2(:),2,mf), &
               gswmin(:,:), rdx(:,:), vcmxt3(:,:), &
               vcmxt4(:,:), vx3(:,:), vx4(:,:), &
               ! Ticket #56, xleuning replaced with gs_coeff here
               gs_coeff(:,:), rad%fvlai(:,:), &
               spread(abs_deltlf,2,mf), &
               anx(:,:), fwsoil(:), gmes(:,:), met, &
               anrubiscox(:,:), anrubpx(:,:), ansinkx(:,:), eta_x(:,:), dAnx(:,:) )
       ELSE
          CALL photosynthesis( csx(:,:), &
               spread(cx1(:),2,mf), &
               spread(cx2(:),2,mf), &
               gswmin(:,:), rdx(:,:), vcmxt3(:,:), &
               vcmxt4(:,:), vx3(:,:), vx4(:,:), &
               ! Ticket #56, xleuning replaced with gs_coeff here
               gs_coeff(:,:), rad%fvlai(:,:), &
               spread(abs_deltlf,2,mf), &
               anx(:,:), fwsoil(:), met, &
               anrubiscox(:,:), anrubpx(:,:), ansinkx(:,:), eta_x(:,:), dAnx(:,:) )
          ! CALL photosynthesis_orig( csx(:,:), &
          !      SPREAD( cx1(:), 2, mf ), &
          !      SPREAD( cx2(:), 2, mf ), &
          !      gswmin(:,:), rdx(:,:), vcmxt3(:,:), &
          !      vcmxt4(:,:), vx3(:,:), vx4(:,:), &
          !      ! Ticket #56, xleuning replaced with gs_coeff here
          !      gs_coeff(:,:), rad%fvlai(:,:), &
          !      SPREAD( abs_deltlf, 2, mf ), &
          !      anx(:,:), fwsoil(:), met , anrubiscox(:,:), anrubpx(:,:) &
          !      , dAnrubiscox(:,:), dAnrubpx(:,:))
       ENDIF


       DO i=1,mp

          IF (canopy%vlaiw(i) > C%LAI_THRESH .AND. abs_deltlf(i) > 0.1 ) Then

             DO kk=1, mf

                IF(rad%fvlai(i,kk)>C%LAI_THRESH) THEN

                   csx(i,kk) = met%ca(i) - C%RGBWC*anx(i,kk) / ( &
                        gbhu(i,kk) + gbhf(i,kk) )
                   csx(i,kk) = MAX( 1.0e-4_r_2, csx(i,kk) )

                   ! Ticket #56, xleuning replaced with gs_coeff here
                   if (cable_user%g0_switch == 'default') then
                      canopy%gswx(i,kk) = MAX( 1.e-3, gswmin(i,kk)*fwsoil(i) + &
                           MAX( 0.0, C%RGSWC * gs_coeff(i,kk) * &
                           anx(i,kk) ) )
                   elseif (cable_user%g0_switch == 'maximum') then
                      ! set gsw to maximum of g0*fwsoil and humidity-dependent term,
                      ! according to third formulation suggested by Lombardozzi et al.,
                      ! GMD 10, 321-331, 2017, and applied in CLM
                      canopy%gswx(i,kk) = MAX( 1.e-3, max(gswmin(i,kk)*fwsoil(i), &
                           MAX( 0.0, C%RGSWC * gs_coeff(i,kk) * &
                           anx(i,kk) )) )
                                              
                   endif
                   !Recalculate conductance for water:
                   gw(i,kk) = 1.0 / ( 1.0 / canopy%gswx(i,kk) + &
                        1.0 / ( 1.075 * ( gbhu(i,kk) + gbhf(i,kk) ) ) )

                   gw(i,kk) = MAX( gw(i,kk), 0.00001 )

                   ! Modified psychrometric constant
                   ! (Monteith and Unsworth, 1990)
                   psycst(i,kk) = air%psyc(i) * REAL( ghr(i,kk) / gw(i,kk) )

                ENDIF

             ENDDO

             ecx(i) = ( air%dsatdk(i) * ( rad%rniso(i,1) - C%capp * C%rmair &
                  * ( met%tvair(i) - met%tk(i) ) * rad%gradis(i,1) ) &
                  + C%capp * C%rmair * met%dva(i) * ghr(i,1) ) &
                  / ( air%dsatdk(i) + psycst(i,1) ) + ( air%dsatdk(i) &
                  * ( rad%rniso(i,2) - C%capp * C%rmair * ( met%tvair(i) - &
                  met%tk(i) ) * rad%gradis(i,2) ) + C%capp * C%rmair * &
                  met%dva(i) * ghr(i,2) ) / &
                  ( air%dsatdk(i) + psycst(i,2) )

             IF (cable_user%fwsoil_switch=='Haverd2013') then
                 ! avoid root-water extraction when fwsoil is zero
                 if (fwsoil(i).lt.1e-6) then
                   anx(i,:) =  - rdx(i,:)
                   ecx(i) = 0.0
                 endif

                canopy%fevc(i) = ecx(i)*(1.0-canopy%fwet(i))

                call getrex_1d(real(ssnow%wb(i,:)-ssnow%wbice(i,:),r_2), ssnow%rex(i,:), &
                     canopy%fwsoil(i), &
                     real(veg%froot(i,:),r_2), SPREAD(real(soil%ssat(i),r_2),1,ms), &
                      SPREAD(real(soil%swilt(i),r_2),1,ms), &
                      max(real(canopy%fevc(i)/air%rlam(i)/1000_r_2,r_2),0.0_r_2), &
                     veg%gamma(i), &
                     real(soil%zse,r_2), real(dels,r_2), veg%zr(i))

                fwsoil(i) = canopy%fwsoil(i)

                where (ssnow%rex(i,:) > real(tiny(1.0),r_2)) &
                     ssnow%evapfbl(i,:) = real(ssnow%rex(i,:))*dels*1000. ! mm water &
                !(root water extraction) per time step

                if (cable_user%Cumberland_soil) then
                   canopy%fwsoil(i) = max(canopy%fwsoil(i), 0.6)
                   fwsoil(i) = canopy%fwsoil(i)
                endif

             ELSE

                IF (ecx(i) > 0.0 .AND. canopy%fwet(i) < 1.0) Then
                   evapfb(i) = ( 1.0 - canopy%fwet(i)) * REAL( ecx(i) ) *dels &
                        / air%rlam(i)

                   DO kk = 1,ms

                      ssnow%evapfbl(i,kk) = MIN( evapfb(i) * veg%froot(i,kk), &
                           MAX( 0.0, REAL( ssnow%wb(i,kk) ) - &
                           1.1 * soil%swilt(i) ) * &
                           soil%zse(kk) * 1000.0 )

                   ENDDO
                   IF (cable_user%soil_struc=='default') then
                      canopy%fevc(i) = SUM(ssnow%evapfbl(i,:))*air%rlam(i)/dels
                      ecx(i) = canopy%fevc(i) / (1.0-canopy%fwet(i))
                   ELSEIF (cable_user%soil_struc=='sli') then
                      canopy%fevc(i) = ecx(i)*(1.0-canopy%fwet(i))
                   ENDIF

                ENDIF

             ENDIF
             ! Update canopy sensible heat flux:
             hcx(i) = (SUM(rad%rniso(i,:))-ecx(i) &
                  - C%capp*C%rmair*(met%tvair(i)-met%tk(i)) &
                  * SUM(rad%gradis(i,:))) &
                  * SUM(gh(i,:))/ SUM(ghr(i,:))

             ! Update leaf temperature:
             tlfx(i)=met%tvair(i)+REAL(hcx(i))/(C%capp*C%rmair*SUM(gh(i,:)))

             ! Update net radiation for canopy:
             rnx(i) = SUM( rad%rniso(i,:)) - &
                  C%CAPP * C%rmair *( tlfx(i)-met%tk(i) ) * &
                  SUM( rad%gradis(i,:) )

             ! Update leaf surface vapour pressure deficit:
             dsx(i) = met%dva(i) + air%dsatdk(i) * (tlfx(i)-met%tvair(i))
             if (cable_user%perturb_dva_by_T) then
                dsx(i) = met%dva(i) + air%dsatdk(i) * (tlfx(i) -met%tvair(i) + cable_user%dva_T_perturbation )
             endif

             dsx(i)=  max(dsx(i),0.0)

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
             anrubiscoy(i,1) = anrubiscox(i,1)
             anrubiscoy(i,2) = anrubiscox(i,2)
             anrubpy(i,1) = anrubpx(i,1)
             anrubpy(i,2) = anrubpx(i,2)
             dAn_y(i,1) = dAnx(i,1)
             dAn_y(i,2) = dAnx(i,2)
             eta_y(i,1) = eta_x(i,1)
             eta_y(i,2) = eta_x(i,2)
             dAnrubiscoy(i,1) = dAnrubiscox(i,1)
             dAnrubiscoy(i,2) = dAnrubiscox(i,2)
             dAnrubpy(i,1) = dAnrubpx(i,1)
             dAnrubpy(i,2) = dAnrubpx(i,2)
             dAnsinky(i,1) = dAnsinkx(i,1)
             dAnsinky(i,2) = dAnsinkx(i,2)

             ! save last values calculated for ssnow%evapfbl
             oldevapfbl(i,:) = ssnow%evapfbl(i,:)

          ENDIF

          IF ( abs_deltlf(i) > 0.1 ) &
               ! after 4 iterations, take mean of current & previous estimates
               ! as the next estimate of leaf temperature, to avoid oscillation
               tlfx(i) = ( 0.5 * ( MAX( 0, k-5 ) / ( k - 4.9999 ) ) ) *tlfxx(i) + &
               ( 1.0 - ( 0.5 * ( MAX( 0, k-5 ) / ( k - 4.9999 ) ) ) ) &
               * tlfx(i)

          IF (k==1) THEN
             ! take the first iterated estimates as the defaults
             tlfy(i) = tlfx(i)
             rny(i) = rnx(i)
             hcy(i) = hcx(i)
             ecy(i) = ecx(i)
             rdy(i,:) = rdx(i,:)
             an_y(i,:) = anx(i,:)
             anrubiscoy(i,:) = anrubiscox(i,:)
             anrubpy(i,:) = anrubpx(i,:)
             dAn_y(i,:) = dAnx(i,:)
             eta_y(i,:) = eta_x(i,:)
             dAnrubiscoy(i,:) = dAnrubiscox(i,:)
             dAnrubpy(i,:) = dAnrubpx(i,:)
             dAnsinky(i,:) = dAnsinkx(i,:)
             ! save last values calculated for ssnow%evapfbl
             oldevapfbl(i,:) = ssnow%evapfbl(i,:)
          END IF

       END DO !over mp

    END DO  ! DO WHILE (ANY(abs_deltlf > 0.1) .AND.  k < C%MAXITER)

!!$  write(3333,"(200e16.6)") tlfx
!!$  write(3334,"(200e16.6)") met%tvair
!!$  write(3335,"(200e16.6)") met%tk
!!$  write(3336,"(200e16.6)") canopy%gswx(:,1)
!!$  write(3337,"(200e16.6)") canopy%gswx(:,2)
!!$  write(3338,"(200e16.6)") gbhu(:,1)
!!$  write(3339,"(200e16.6)") gbhu(:,2)
!!$  write(3340,"(200e16.6)") gbhf(:,1)
!!$  write(3341,"(200e16.6)") gbhf(:,2)
!!$  write(3342,"(200e16.6)") rad%qcan(:,1,1)
!!$  write(3343,"(200e16.6)") rad%qcan(:,2,1)
!!$  write(3344,"(200e16.6)") gw(:,1)
!!$  write(3345,"(200e16.6)") gw(:,2)
!!$  write(3346,"(200e16.6)") met%dva
!!$  write(3347,"(200e16.6)") met%ua
!!$  write(3348,"(200e16.6)") rad%fvlai(:,1)
!!$  write(3349,"(200e16.6)") rad%fvlai(:,2)
!!$  write(3350,"(200e16.6)") dsx
!!$  write(3351,"(200e16.6)") ecx
!!$  write(3352,"(200e16.6)") met%hod
!!$  write(3353,"(200e16.6)") an_y(:,1)
!!$  write(3354,"(200e16.6)") an_y(:,2)
!!$  write(3355,"(200e16.6)") veg%vcmax_sun
!!$  write(3356,"(200e16.6)") veg%ejmax_sun
!!$  write(3357,"(200e16.6)") gmes(:,1)
!!$  write(3358,"(200e16.6)") gmes(:,2)
!!$  write(3359,"(200e16.6)") met%ca
!!$  write(3360,"(200e16.6)") csx(:,1)
!!$
!!$  write(3361,"(200e16.6)") csx(:,1) - C%RGBWC*anx(:,1) / ( &
!!$                      max( canopy%gswx(:,1), 0.01) )
!!$
!!$write(3362,"(200e16.6)") csx(:,1) - C%RGBWC*anx(:,1) / ( &
!!$                       max(canopy%gswx(:,1),0.01) )- anx(:,1) / ( &
!!$                       max(gmes(:,1),0.01) )
!!$ write(3363,"(200e16.6)") canopy%fwsoil
!!$ write(3364,"(200e16.6)")  gs_coeff(:,1)
!!$ write(3365,"(200e16.6)")  gs_coeff(:,2)
!!$ write(3366,"(200e16.6)") veg%vcmax_shade
!!$ write(3367,"(200e16.6)") veg%ejmax_shade
!!$ write(3378,"(200e16.6)") anrubiscoy(:,2)
!!$ write(3369,"(200e16.6)") anrubpy(:,2)
    ! dry canopy flux
    canopy%fevc = (1.0-canopy%fwet) * ecy

    IF (cable_user%fwsoil_switch.ne.'Haverd2013') then

       ! Recalculate ssnow%evapfbl as ecy may not be updated with the ecx
       ! calculated in the last iteration.
       ! DO NOT use simple scaling as there are times that ssnow%evapfbl is zero.
       ! ** ssnow%evapfbl(i,:) = ssnow%evapfbl(i,:) * ecy(i) / ecx(i) **
       DO i = 1, mp

          IF( ecy(i) > 0.0 .AND. canopy%fwet(i) < 1.0 ) THEN

             IF( ABS( ecy(i) - ecx(i) ) > 1.0e-6 ) THEN

                IF( ABS( canopy%fevc(i) - ( SUM( oldevapfbl(i,:)) * air%rlam(i) &
                     /dels ) ) > 1.0e-4 ) THEN

                   PRINT *, 'Error! oldevapfbl not right.', ktau_gl, i
                   PRINT *, 'ecx, ecy = ', ecx(i), ecy(i)
                   PRINT *, 'or in mm = ', ecx(i) * ( 1.0 - canopy%fwet(i) ) &
                        / air%rlam(i) * dels, &
                        ecy(i) * ( 1.0 - canopy%fwet(i) ) / &
                        air%rlam(i) * dels

                   PRINT *,'fevc = ', canopy%fevc(i), SUM( oldevapfbl(i,:) ) * &
                        air%rlam(i) / dels
                   PRINT *, 'fwet = ', canopy%fwet(i)
                   PRINT *, 'oldevapfbl = ', oldevapfbl(i,:)
                   PRINT *, 'ssnow%evapfbl before rescaling: ', &
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
    canopy%fpn = min(-12.0 * SUM(an_y, 2), canopy%frday)

    ! additional diagnostic variables for assessing contributions of rubisco and rubp limited photosynthesis to
    ! net photosynthesis in sunlit and shaded leaves.
    canopy%A_sh = an_y(:,2)
    canopy%A_sl = an_y(:,1)

    canopy%GPP_sh = an_y(:,2) + rdy(:,2)
    canopy%GPP_sl = an_y(:,1) + rdy(:,1)

    canopy%A_shC = anrubiscoy(:,2)
    canopy%A_shJ = anrubpy(:,2)

    where (anrubiscoy(:,2) .gt. an_y(:,2)) canopy%A_shC = 0.0_r_2
    where (anrubpy(:,2) .gt. an_y(:,2))    canopy%A_shJ = 0.0_r_2

    canopy%A_slC = anrubiscoy(:,1)
    canopy%A_slJ = anrubpy(:,1)

    where (anrubiscoy(:,1) .gt. an_y(:,1)) canopy%A_slC = 0.0
    where (anrubpy(:,1)    .gt. an_y(:,1)) canopy%A_slJ = 0.0

    canopy%eta_A_cs = canopy%A_sh * min(eta_y(:,2),5.0_r_2) + canopy%A_sl * min(eta_y(:,1),5.0_r_2)
    where ((canopy%A_sl .gt. 0.0_r_2) .and. (canopy%A_sh .gt. 0.0_r_2))
       canopy%eta_GPP_cs = canopy%GPP_sh  * min(eta_y(:,2),5.0_r_2) * canopy%GPP_sh/canopy%A_sh + &
            canopy%GPP_sl * min(eta_y(:,1),5.0_r_2)* canopy%GPP_sl/canopy%A_sl
    elsewhere ((canopy%A_sl .le. 0.0_r_2) .and. (canopy%A_sh .gt. 0.0_r_2))
        canopy%eta_GPP_cs = canopy%GPP_sh * min(eta_y(:,2),5.0_r_2) * canopy%GPP_sh/canopy%A_sh
    elsewhere
        canopy%eta_GPP_cs = canopy%eta_A_cs
    end where

    where (canopy%A_sl .gt. 0.0_r_2)
       canopy%eta_A_cs_sl    =  min(eta_y(:,1),5.0_r_2)
       canopy%eta_fevc_cs_sl = (min(eta_y(:,1),5.0_r_2) - 1.0_r_2) * &
            MAX(0.0_r_2, gs_coeff(:,1)*an_y(:,1)) / canopy%gswx(:,1)
    elsewhere
       canopy%eta_A_cs_sl    = 0.0_r_2
       canopy%eta_fevc_cs_sl = 0.0_r_2
    endwhere

    where (canopy%A_sh .gt. 0.0_r_2)
       canopy%eta_A_cs_sh    =  min(eta_y(:,2),5.0_r_2)
       canopy%eta_fevc_cs_sh = (min(eta_y(:,2),5.0_r_2) - 1.0_r_2) * &
            MAX(0.0_r_2, gs_coeff(:,2)*an_y(:,2)) / canopy%gswx(:,2)
    elsewhere
       canopy%eta_A_cs_sh    = 0.0_r_2
       canopy%eta_fevc_cs_sh = 0.0_r_2
    endwhere

    where ((rad%fvlai(:,1)+rad%fvlai(:,2)).gt.0.01)
       canopy%eta_fevc_cs =  (canopy%eta_fevc_cs_sh*canopy%fevc* rad%fvlai(:,2)/(rad%fvlai(:,1)+rad%fvlai(:,2)) + &
         canopy%eta_fevc_cs_sl*canopy%fevc* rad%fvlai(:,1)/(rad%fvlai(:,1)+rad%fvlai(:,2)))
    elsewhere
       canopy%eta_fevc_cs = 0.0_r_2
    endwhere

    canopy%dAdcs = canopy%A_sh*dAn_y(:,2)   + canopy%A_sl*dAn_y(:,1)
    canopy%cs    = canopy%A_sh*csx(:,2)*1.0e6_r_2 + canopy%A_sl*csx(:,1)*1.0e6_r_2
    canopy%cs_sl = csx(:,1)*1.0e6_r_2
    canopy%cs_sh = csx(:,2)*1.0e6_r_2
    canopy%tlf   = tlfy
    canopy%dlf   = dsx

    canopy%evapfbl = ssnow%evapfbl

    ! 13C
    canopy%An        = real(an_y,r_2)
    canopy%Rd        = real(rdy,r_2)
    canopy%isc3      = (1.0-veg%frac4) > epsilon(1.0)
    canopy%vcmax     = real(merge(vcmxt3, vcmxt4, spread(canopy%isc3,2,mf)), r_2)
    canopy%gammastar = real(spread(cx2*0.5,2,mf), r_2)
    canopy%gsc       = real(canopy%gswx / C%rgswc, r_2)
    canopy%gbc       = real((gbhu + gbhf) / C%rgbwc, r_2)
    canopy%gac       = huge(1.0_r_2)
    ! replace 1.0_r_2/canopy%gac by tiny(1.0_r_2) to avoid underflow
    canopy%ci        = real(spread(met%ca,2,mf),r_2) - &
         canopy%An / (1._r_2/(tiny(1.0_r_2)+1.0_r_2/canopy%gbc+1.0_r_2/canopy%gsc))

    deallocate( gswmin )

  END SUBROUTINE dryLeaf
  
  ! -----------------------------------------------------------------------------


   ! Ticket #56, xleuningz repalced with gs_coeffz
   SUBROUTINE photosynthesis_orig( csxz, cx1z, cx2z, gswminz, &
        rdxz, vcmxt3z, vcmxt4z, vx3z, &
        vx4z, gs_coeffz, vlaiz, deltlfz, anxz, fwsoilz, met, &
        anrubiscoz, anrubpz , dAmc, dAme )
     
     use cable_def_types_mod, only : mp, mf, r_2, met_type

     implicit none

     real(r_2), dimension(mp,mf), intent(in)    :: csxz
     type (met_type),             intent(in)    :: met
     real,      dimension(mp,mf), intent(in)    :: &
          cx1z, & !
          cx2z, & !
          gswminz, & !
          rdxz, & !
          vcmxt3z, & !
          vcmxt4z, & !
          vx4z, & !
          vx3z, & !
          gs_coeffz, & ! Ticket #56, xleuningz repalced with gs_coeffz
          vlaiz, & !
          deltlfz       !

     real,      dimension(mp,mf), intent(inout) :: anxz, anrubiscoz, anrubpz
     real(r_2), dimension(mp,mf), intent(out)   :: dAmc, dAme
     
     ! local variables
     real(r_2), dimension(mp,mf) :: &
          coef0z, coef1z, coef2z, ciz, delcxz
     real(r_2) :: gamma, beta, gammast, g0, X, Rd, cs
     real(r_2) :: Am
     real, dimension(mp,mf) :: ansinkz
     real, dimension(mp) :: fwsoilz
     real, parameter  :: effc4 = 4000.0  ! Vc=effc4*Ci*Vcmax (see Bonan,LSM version 1.0, p106)
     integer :: i, j

     !MC - fwsoil not defined. Stopped single/double precision rewrite.
     !     subroutine not used at the moment
     stop 2

     dAmc = 0.0_r_2
     dAme = 0.0_r_2

     DO i=1,mp

        IF (sum(vlaiz(i,:)) .GT. C%LAI_THRESH) THEN

           DO j=1, mf

              IF ((vlaiz(i,j) .GT. C%LAI_THRESH) .AND. (deltlfz(i,j) .GT. 0.1)) THEN

                 ! Rubisco limited:
                 coef2z(i,j) = gswminz(i,j)*fwsoilz(i) / C%RGSWC + gs_coeffz(i,j) * &
                      ( vcmxt3z(i,j) - ( rdxz(i,j)-vcmxt4z(i,j) ) )

                 coef1z(i,j) = (1.0_r_2-csxz(i,j)*gs_coeffz(i,j)) * &
                      (vcmxt3z(i,j)+vcmxt4z(i,j)-rdxz(i,j)) &
                      + (gswminz(i,j)*fwsoilz(i)/C%RGSWC)*(cx1z(i,j)-csxz(i,j)) &
                      - gs_coeffz(i,j)*(vcmxt3z(i,j)*cx2z(i,j)/2.0_r_2 &
                      + cx1z(i,j)*(rdxz(i,j)-vcmxt4z(i,j) ) )

                 coef0z(i,j) = -(1.0_r_2-csxz(i,j)*gs_coeffz(i,j)) * &
                      (vcmxt3z(i,j)*cx2z(i,j)/2.0_r_2 &
                      + cx1z(i,j)*( rdxz(i,j)-vcmxt4z(i,j ) ) ) &
                      -( gswminz(i,j)*fwsoilz(i)/C%RGSWC ) * cx1z(i,j)*csxz(i,j)

                ! kdcorbin,09/10 - new calculations
                if ( (abs(coef2z(i,j)) .gt. 1.0e-9_r_2) .and. &
                     (abs(coef1z(i,j)) .lt. 1.0e-9_r_2) ) then
                   ! no solution, give it a huge number as
                   ! quadratic below cannot handle zero denominator
                   ciz(i,j)        = 99999.0_r_2
                   anrubiscoz(i,j) = real(ciz(i,j)) ! OR do ciz=0 and calc anrubiscoz
                endif

                ! solve linearly
                if ( (abs(coef2z(i,j)) < 1.0e-9_r_2) .and. &
                     (abs(coef1z(i,j)) >= 1.0e-9_r_2) ) then
                   ! same reason as above
                   ciz(i,j) = -1.0_r_2 * coef0z(i,j) / coef1z(i,j)
                   ciz(i,j) = max(0.0_r_2, ciz(i,j))
                   ciz(i,j) = vcmxt3z(i,j)*(ciz(i,j)-cx2z(i,j) / 2.0_r_2) / &
                        (ciz(i,j) + cx1z(i,j)) + vcmxt4z(i,j) - rdxz(i,j)
                   anrubiscoz(i,j) = real(ciz(i,j)) ! OR do ciz=0 and calc anrubiscoz
                endif

                ! solve quadratic (only take the more positive solution)
                if (abs(coef2z(i,j)) >= 1.0e-9_r_2) then

                   delcxz(i,j) = coef1z(i,j)**2 - 4.0_r_2 * coef0z(i,j) * coef2z(i,j)

                   ciz(i,j) = ( -coef1z(i,j) + sqrt( max(0.0_r_2, delcxz(i,j)) ) ) / &
                        (2.0_r_2*coef2z(i,j))
                   ciz(i,j) = max(0.0_r_2, ciz(i,j)) ! must be positive (concentration always +ve)

                   ciz(i,j) = vcmxt3z(i,j) * (ciz(i,j) - cx2z(i,j) / 2.0_r_2) / &
                        (ciz(i,j) + cx1z(i,j)) + vcmxt4z(i,j) - rdxz(i,j)
                   anrubiscoz(i,j) = real(ciz(i,j)) ! OR do ciz=0 and calc anrubiscoz
                endif


                ! RuBP limited:
                coef2z(i,j) = gswminz(i,j)*fwsoilz(i) / C%RGSWC + gs_coeffz(i,j) &
                              * ( vx3z(i,j) - ( rdxz(i,j) - vx4z(i,j) ) )

                coef1z(i,j) = ( 1.0 - csxz(i,j) * gs_coeffz(i,j) ) * &
                              ( vx3z(i,j) + vx4z(i,j) - rdxz(i,j) ) &
                              + ( gswminz(i,j)*fwsoilz(i) / C%RGSWC ) * &
                              ( cx2z(i,j) - csxz(i,j) ) - gs_coeffz(i,j) &
                              * ( vx3z(i,j) * cx2z(i,j) / 2.0 + cx2z(i,j) * &
                              ( rdxz(i,j) - vx4z(i,j) ) )

                coef0z(i,j) = -(1.0-csxz(i,j)*gs_coeffz(i,j)) * &
                              (vx3z(i,j)*cx2z(i,j)/2.0 &
                              + cx2z(i,j)*(rdxz(i,j)-vx4z(i,j))) &
                          - (gswminz(i,j)*fwsoilz(i)/C%RGSWC)*cx2z(i,j)*csxz(i,j)


                !kdcorbin, 09/10 - new calculations
                ! no solution, give it a huge number
                IF( ABS( coef2z(i,j) ) < 1.0e-9 .AND. &
                     ABS( coef1z(i,j) ) < 1.0e-9 ) THEN

                   ciz(i,j)     = 99999.0_r_2
                   anrubpz(i,j) = real(ciz(i,j))

                ENDIF

                ! solve linearly
                IF( ABS( coef2z(i,j) ) < 1.e-9 .AND. &
                     ABS( coef1z(i,j) ) >= 1.e-9) THEN

                   ciz(i,j) = -1.0 * coef0z(i,j) / coef1z(i,j)

                   ciz(i,j) = MAX(0.0_r_2,ciz(i,j))

                   ciz(i,j) = vx3z(i,j)*(ciz(i,j)-cx2z(i,j)/2.0) / &
                        (ciz(i,j)+cx2z(i,j)) +vx4z(i,j)-rdxz(i,j)
                   anrubpz(i,j) = real(ciz(i,j))

                ENDIF

                ! solve quadratic (only take the more positive solution)
                IF ( ABS( coef2z(i,j)) >= 1.e-9 ) THEN

                   delcxz(i,j) = coef1z(i,j)**2 -4.0*coef0z(i,j)*coef2z(i,j)

                   ciz(i,j) = (-coef1z(i,j)+SQRT(MAX(0.0_r_2,delcxz(i,j)))) &
                        /(2.0*coef2z(i,j))

                   ciz(i,j) = MAX(0.0_r_2,ciz(i,j))

                   ciz(i,j) = vx3z(i,j)*(ciz(i,j)-cx2z(i,j)/2.0) / &
                        (ciz(i,j)+cx2z(i,j)) +vx4z(i,j)-rdxz(i,j)
                   anrubpz(i,j) = real(ciz(i,j))

                   gamma =   vx3z(i,j)
                   beta = cx2z(i,j)
                   X = gs_coeffz(i,j)
                   g0 = gswminz(i,j)*fwsoilz(i) / C%RGSWC
                   cs = csxz(i,j)
                   gammast = cx2z(i,j)/2.0
                   Rd = rdxz(i,j)
                   CALL fAndAn1(cs, g0, X*cs, gamma, beta, gammast, Rd, &
                        Am, dAme(i,j))

                ENDIF

                ! Sink limited:
                coef2z(i,j) = gs_coeffz(i,j)

                coef1z(i,j) = gswminz(i,j)*fwsoilz(i)/C%RGSWC + gs_coeffz(i,j) &
                              * (rdxz(i,j) - 0.5*vcmxt3z(i,j)) &
                              + effc4 * vcmxt4z(i,j) - gs_coeffz(i,j) &
                              * csxz(i,j) * effc4 * vcmxt4z(i,j)

                coef0z(i,j) = -( gswminz(i,j)*fwsoilz(i)/C%RGSWC )*csxz(i,j)*effc4 &
                              * vcmxt4z(i,j) + ( rdxz(i,j) &
                            - 0.5 * vcmxt3z(i,j)) * gswminz(i,j)*fwsoilz(i)/C%RGSWC

                ! no solution, give it a huge number
                IF( ABS( coef2z(i,j) ) < 1.0e-9 .AND. &
                     ABS( coef1z(i,j)) < 1.0e-9 ) THEN

                   ciz(i,j)     = 99999.0_r_2
                   ansinkz(i,j) = real(ciz(i,j))

                ENDIF

                ! solve linearly
                IF( ABS( coef2z(i,j) ) < 1.e-9 .AND. &
                     ABS( coef1z(i,j) ) >= 1.e-9 ) THEN

                   ciz(i,j) = -1.0 * coef0z(i,j) / coef1z(i,j)
                   ansinkz(i,j) = real(ciz(i,j))

                ENDIF

                ! solve quadratic (only take the more positive solution)
                IF( ABS( coef2z(i,j) ) >= 1.e-9 ) THEN

                   delcxz(i,j) = coef1z(i,j)**2 -4.0*coef0z(i,j)*coef2z(i,j)

                   ciz(i,j) = (-coef1z(i,j)+SQRT (MAX(0.0_r_2,delcxz(i,j)) ) ) &
                        / ( 2.0 * coef2z(i,j) )

                   ansinkz(i,j) = real(ciz(i,j))

                ENDIF

                ! minimal of three limited rates


             ENDIF
             anxz(i,j) = MIN(anrubiscoz(i,j),anrubpz(i,j),ansinkz(i,j))
             ! if (i==1 .and. j==2 .and. floor(met%hod(1))==12) then
             !    write(3372,"(200e16.6)") anrubiscoz(i,j), &
             !         anrubpz(i,j), &
             !         ansinkz(i,j),  anxz(i,j)
             ! elseif (i.eq.4 .and.j==2.and. floor(met%hod(1))==12) then
             !    write(3373,"(200e16.6)") anrubiscoz(i,j), &
             !         anrubpz(i,j), &
             !         ansinkz(i,j),  anxz(i,j)
             ! endif

          ENDDO

       ENDIF

    ENDDO
    
  END SUBROUTINE photosynthesis_orig

  
  ! ------------------------------------------------------------------------------

  
  ! vh adaptation of photosynthesis calculation to
  ! account for finite mesophyll conductance (gm) (cable_user%explicit_gm = TRUE)
  SUBROUTINE photosynthesis_gm( csxz, cx1z, cx2z, gswminz, &
       rdxz, vcmxt3z, vcmxt4z, vx3z, &
       vx4z, gs_coeffz, vlaiz, deltlfz, anxz, fwsoilz, &
       gmes,met, anrubiscoz, anrubpz, ansinkz, eta, dA )
    
    use cable_def_types_mod, only : mp, mf, r_2, met_type
    use cable_common_module, only: cable_user

    implicit none
    
    real(r_2), dimension(mp,mf), intent(in)    :: csxz, gmes
    type(met_type),              intent(in)    :: met
    real,      dimension(mp,mf), intent(in)    :: &
         cx1z,       & !
         cx2z,       & !
         rdxz,       & !
         vcmxt3z,    & !
         vcmxt4z,    & !
         vx4z,       & !
         vx3z,       & !
         gs_coeffz,  & ! Ticket #56, xleuningz repalced with gs_coeffz
         vlaiz,      & !
         deltlfz
    real,      dimension(mp,mf), intent(inout) :: gswminz
    real,      dimension(mp,mf), intent(inout) :: anxz, anrubiscoz, anrubpz, ansinkz
    real(r_2), dimension(mp,mf), intent(out)   :: eta, dA

    ! local variables
    real(r_2), dimension(mp,mf) :: dAmp, dAme, dAmc, eta_p, eta_e, eta_c
    real, dimension(mp) :: fwsoilz
    real, parameter  :: effc4 = 4000.0  ! Vc=effc4*Ci*Vcmax (see Bonan,LSM version 1.0, p106)
    real(r_2) :: gamma, beta, gammast, gm, g0, X, Rd, cs
    real(r_2) :: cc
    real(r_2) :: Am, gsm
    integer :: i, j
    ! for minimum of 3 rates and corresponding elasticities
    real,      dimension(3) :: tmp3
    real(r_2), dimension(3) :: dtmp3
    integer,   dimension(1) :: ii

    do i=1, mp
       
       do j=1, mf
          
          if (vlaiz(i,j) .gt. C%lai_thresh) then
             
             if (deltlfz(i,j) .gt. 0.1) then
                
                anxz(i,j)       = -rdxz(i,j)
                anrubiscoz(i,j) = -rdxz(i,j)
                anrubpz(i,j)    = -rdxz(i,j)
                ansinkz(i,j)    = -rdxz(i,j)
                dAmc(i,j)       = 0.0_r_2
                dAme(i,j)       = 0.0_r_2
                dAmp(i,j)       = 0.0_r_2
                dA(i,j)         = 0.0_r_2
                eta_c(i,j)      = 0.0_r_2
                eta_e(i,j)      = 0.0_r_2
                eta_p(i,j)      = 0.0_r_2
                eta(i,j)        = 0.0_r_2

                ! C3, Rubisco limited, accounting for explicit mesophyll conductance
                if ((vcmxt3z(i,j) .gt. 1.0e-8) .and. (gs_coeffz(i,j) .gt. 100.)) then
                   cs      = csxz(i,j)
                   g0      = real(gswminz(i,j) * fwsoilz(i) / C%RGSWC, r_2)
                   X       = real(gs_coeffz(i,j), r_2)
                   gamma   = real(vcmxt3z(i,j), r_2)
                   beta    = real(cx1z(i,j), r_2)
                   gammast = real(cx2z(i,j) / 2.0, r_2)
                   Rd      = real(rdxz(i,j), r_2)
                   gm      = gmes(i,j)
                   
                   if (trim(cable_user%g0_switch) == 'default') then
                      ! get partial derivative of A wrt cs
                      call fAmdAm1(cs, g0, X*cs, gamma, beta, gammast, Rd, &
                           gm, Am, dAmc(i,j))
                   elseif (trim(cable_user%g0_switch) == 'maximum') then
                      ! set g0 to zero initially
                      call fAmdAm1(cs, 0.0_r_2, X*cs, gamma, beta, gammast, Rd, &
                           gm, Am, dAmc(i,j))
                      ! repeat calculation if g0 > A*X
                      if (g0 .gt. Am*X) then
                         call fAmdAm1(cs, g0, 0.0_r_2, gamma, beta, gammast, Rd, &
                              gm, Am, dAmc(i,j))
                      endif
                   endif
                   anrubiscoz(i,j) = real(Am)
                   
                   gsm = (gm * (g0+X*Am)) / (gm + (g0+X*Am))
                   cc  = cs - Am / max(gsm, 1.0e-4_r_2)
                   if (Am > 0.0_r_2) eta_c(i,j) = dAmc(i,j) * cs / Am
                endif  ! C3

                ! C4, Rubisco limited, accounting for explicit mesophyll conductance
                if (vcmxt4z(i,j) .gt. 1.0e-8) then
                   anrubiscoz(i,j) = vcmxt4z(i,j)- rdxz(i,j)
                   dAmc(i,j)  = 0.0_r_2
                   eta_c(i,j) = 0.0_r_2
                endif

                ! C3, RuBP regeneration-limited, accounting for explicit mesophyll conductance
                if ( (vcmxt3z(i,j) .gt. 0.0) .and. (gs_coeffz(i,j) .gt. 100.) .and. &
                     (vx3z(i,j) .gt. 1.0e-8) ) then
                   cs      = csxz(i,j)
                   g0      = real(gswminz(i,j) * fwsoilz(i) / C%RGSWC, r_2)
                   X       = real(gs_coeffz(i,j), r_2)
                   gamma   = real(vx3z(i,j), r_2)
                   beta    = real(cx2z(i,j), r_2)
                   gammast = real(cx2z(i,j) / 2.0, r_2)
                   Rd      = real(rdxz(i,j), r_2)
                   gm      = gmes(i,j)

                   if (trim(cable_user%g0_switch) == 'default') then
                      call fAmdAm1(cs, g0, X*cs, gamma, beta, gammast, Rd, &
                           gm, Am, dAme(i,j))
                   elseif (trim(cable_user%g0_switch) == 'maximum') then
                      call fAmdAm1(cs, 0.0_r_2, X*cs, gamma, beta, gammast, Rd, &
                           gm, Am, dAme(i,j))
                      ! repeat calculation if g0 > A*X
                      if (g0 .gt. Am*X) then
                         call fAmdAm1(cs, g0, 0.0_r_2, gamma, beta, gammast, Rd, &
                              gm, Am, dAme(i,j))
                      endif
                   endif
                   anrubpz(i,j) = real(Am)
                   
                   if (Am > 0.0_r_2) eta_e(i,j) = dAme(i,j) * cs / Am
                   gsm = (gm * (g0+X*Am)) / (gm + (g0+X*Am))
                   cc  = cs - Am / max(gsm, 1.0e-4_r_2)
                endif  ! C3

                ! C4, RuBP limited, accounting for explicit mesophyll conductance
                if (vx4z(i,j) .gt. 1.0e-8) then
                   anrubpz(i,j) = vx4z(i,j) - rdxz(i,j)
                   dame(i,j)  = 0.0_r_2
                   eta_e(i,j) = 0.0_r_2
                endif

                ! Sink limited, accounting for explicit mesophyll conductance
                if (vcmxt3z(i,j) .gt. 1.0e-10) then ! C3
                   
                   ansinkz(i,j) = 0.5 * vcmxt3z(i,j) - rdxz(i,j)
                   dAmp(i,j)  = 0.0_r_2
                   eta_p(i,j) = 0.0_r_2
                   
                elseif ((vcmxt4z(i,j) .gt. 1.0e-10) .and. (gs_coeffz(i,j) .gt. 100.)) then ! C4
                   cs         = csxz(i,j)
                   g0         = 0.0_r_2
                   X          = real(gs_coeffz(i,j), r_2)
                   gamma      = real(effc4 * vcmxt4z(i,j), r_2)
                   beta       = 0.0_r_2
                   gammast    = 0.0_r_2
                   Rd         = real(rdxz(i,j), r_2)
                   gm         = gmes(i,j)

                   if (trim(cable_user%g0_switch) == 'default') then
                      call fAndAn2(cs, g0, X*cs, gamma, beta, gammast, Rd, &
                           Am, dAmp(i,j))
                   elseif (trim(cable_user%g0_switch) == 'maximum') then
                      call fAndAn2(cs, 0.0_r_2, X*cs, gamma, beta, gammast, Rd, &
                           Am, dAmp(i,j))
                      if (g0 .gt. Am*X) then
                         call fAndAn2(cs, g0, 0.0_r_2, gamma, beta, gammast, Rd, &
                              Am, dAmp(i,j))
                      endif
                   endif
                   ! CALL fAmdAm2(cs, g0, X*cs, gamma, beta, gammast, Rd, &
                   !      gm, Am, dAmp(i,j))  ! vh ! not sure why this one is commented out: need to test
                   ansinkz(i,j) = real(Am)
                   
                   dAmp(i,j)  = 0.0_r_2
                   eta_p(i,j) = 0.0_r_2
                   if (Am > 0.0_r_2) eta_p(i,j) = dAmp(i,j) * cs / Am
                endif ! C3/C4

             endif ! deltlfz > 0.1

             ! minimal of three limited rates
             ! anxz(i,j) = MIN(anrubiscoz(i,j),anrubpz(i,j),ansinkz(i,j))
             ! if ( anxz(i,j) .eq. anrubiscoz(i,j)) then
             !    dA(i,j) = dAmc(i,j)
             !    eta(i,j) = eta_c(i,j)
             ! elseif ( anxz(i,j) .eq. anrubpz(i,j)) then
             !    dA(i,j) = dAme(i,j)
             !    eta(i,j) = eta_e(i,j)
             ! elseif ( anxz(i,j) .eq. ansinkz(i,j)) then
             !    dA(i,j) = dAmp(i,j)
             !    eta(i,j) = eta_p(i,j)
             ! endif
             tmp3      = (/ anrubiscoz(i,j), anrubpz(i,j), ansinkz(i,j) /)
             ii        = minloc(tmp3)
             anxz(i,j) = tmp3(ii(1))
             dtmp3     = (/ dAmc(i,j), dAme(i,j), dAmp(i,j) /)
             dA(i,j)   = dtmp3(ii(1))
             dtmp3     = (/ eta_c(i,j), eta_e(i,j), eta_p(i,j) /)
             eta(i,j)  = dtmp3(ii(1))

          else ! vlaiz(i,j) .gt. C%lai_thresh
             
             anxz(i,:)       = 0.0
             anrubiscoz(i,:) = 0.0
             anrubpz(i,:)    = 0.0
             ansinkz(i,:)    = 0.0
             dAmc(i,:)       = 0.0_r_2
             dAme(i,:)       = 0.0_r_2
             dAmp(i,:)       = 0.0_r_2
             dA(i,:)         = 0.0_r_2
             eta_c(i,:)      = 0.0_r_2
             eta_e(i,:)      = 0.0_r_2
             eta_p(i,:)      = 0.0_r_2
             eta(i,:)        = 0.0_r_2

          endif ! vlaiz(i,j) .gt. C%lai_thresh

       enddo ! j=1, mf
       
    enddo ! i=1, mp
    
  END SUBROUTINE photosynthesis_gm


  ! ------------------------------------------------------------------------------


  SUBROUTINE photosynthesis( csxz, cx1z, cx2z, gswminz, &
       rdxz, vcmxt3z, vcmxt4z, vx3z,                    &
       vx4z, gs_coeffz, vlaiz, deltlfz, anxz, fwsoilz,  &
       met, anrubiscoz, anrubpz, ansinkz, eta, dA )

    use cable_def_types_mod, only : mp, mf, r_2, met_type
    use cable_common_module, only: cable_user

    implicit none

    real(r_2), dimension(mp,mf), intent(in) :: csxz
    type(met_type),              intent(in) :: met
    real,      dimension(mp,mf), intent(in) :: &
         cx1z,       & !
         cx2z,       & !
         rdxz,       & !
         vcmxt3z,    & !
         vcmxt4z,    & !
         vx4z,       & !
         vx3z,       & !
         gswminz,    &
         vlaiz,      & !
         deltlfz
    real,      dimension(mp,mf), intent(inout) :: gs_coeffz ! Ticket #56, xleuningz replaced with gs_coeffz
    real,      dimension(mp,mf), intent(inout) :: anxz, anrubiscoz, anrubpz, ansinkz
    real(r_2), dimension(mp,mf), intent(out)   :: eta, dA

    ! local variables
    real(r_2), dimension(mp,mf) :: &
         coef0z, coef1z, coef2z, ciz, delcxz, &
         dAmc, dAme, dAmp, eta_c, eta_e, eta_p

    real, dimension(mp) :: fwsoilz

    real, parameter  :: effc4 = 4000.0  ! Vc=effc4*Ci*Vcmax (see Bonan,LSM version 1.0, p106)

    real(r_2) :: gamma, beta, gammast, g0, X, Rd, cs
    real(r_2) :: Am
    integer :: i, j
    ! for minimum of 3 rates and corresponding elasticities
    real,      dimension(3) :: tmp3
    real(r_2), dimension(3) :: dtmp3
    integer,   dimension(1) :: ii

    do i=1, mp
       do j=1, mf
          if (vlaiz(i,j) .gt. C%lai_thresh) then
             if (deltlfz(i,j) .gt. 0.1) then
                anxz(i,j)       = -rdxz(i,j)
                anrubiscoz(i,j) = -rdxz(i,j)
                anrubpz(i,j)    = -rdxz(i,j)
                ansinkz(i,j)    = -rdxz(i,j)
                dAmc(i,j)       = 0.0_r_2
                dAme(i,j)       = 0.0_r_2
                dAmp(i,j)       = 0.0_r_2
                dA(i,j)         = 0.0_r_2
                eta_c(i,j)      = 0.0_r_2
                eta_e(i,j)      = 0.0_r_2
                eta_p(i,j)      = 0.0_r_2
                eta(i,j)        = 0.0_r_2

                ! Rubisco limited:

                if ((vcmxt3z(i,j) .gt. 1.0e-8) .and. (gs_coeffz(i,j) .gt. 100.)) then ! C3

                   cs      = csxz(i,j)
                   g0      = real(gswminz(i,j) * fwsoilz(i) / C%RGSWC, r_2)
                   X       = real(gs_coeffz(i,j), r_2)
                   gamma   = real(vcmxt3z(i,j), r_2)
                   beta    = real(cx1z(i,j), r_2)
                   gammast = real(cx2z(i,j) / 2.0, r_2)
                   Rd      = real(rdxz(i,j), r_2)

                   if (trim(cable_user%g0_switch) == 'default') then
                      ! get An and partial derivative of A wrt cs
                      call fAndAn1(cs, g0, X*cs, gamma, beta, gammast, Rd, &
                           Am, dAmc(i,j))
                   elseif (trim(cable_user%g0_switch) == 'maximum') then
                      ! set g0 to zero initially
                      call fAndAn1(cs, 0.0_r_2, X*cs, gamma, beta, gammast, Rd, &
                           Am, dAmc(i,j))
                      ! repeat calculation if g0 > A*X
                      if (g0 .gt. Am*X) then
                         call fAndAn1(cs, g0, 0.0_r_2, gamma, beta, gammast, Rd, &
                              Am, dAmc(i,j))
                      endif
                   endif

                   anrubiscoz(i,j) = real(Am)
                   if (Am > 0.0_r_2) eta_c(i,j) = dAmc(i,j) * cs / Am

                endif  ! C3, Rubisco limited

                ! C4
                if (vcmxt4z(i,j) .gt. 1.0e-8) then
                   anrubiscoz(i,j) = vcmxt4z(i,j) - rdxz(i,j)
                   dAmc(i,j)  = 0.0_r_2
                   eta_c(i,j) = 0.0_r_2
                endif

                if ((vcmxt3z(i,j) .gt. 0.0) .and. (gs_coeffz(i,j) .gt. 100.) .and. &
                     (vx3z(i,j) .gt. 1.0e-8)) then ! C3

                   cs      = csxz(i,j)
                   g0      = real(gswminz(i,j) * fwsoilz(i) / C%RGSWC, r_2)
                   X       = real(gs_coeffz(i,j), r_2)
                   gamma   = real(vx3z(i,j), r_2)
                   beta    = real(cx2z(i,j), r_2)
                   gammast = real(cx2z(i,j) / 2.0, r_2)
                   Rd      = real(rdxz(i,j), r_2)

                   if (trim(cable_user%g0_switch) == 'default') then
                      call fAndAn1(cs, g0, X*cs, gamma, beta, gammast, Rd, &
                           Am, dAme(i,j))
                   elseif (trim(cable_user%g0_switch) == 'maximum') then
                      call fAndAn1(cs, 0.0_r_2, X*cs, gamma, beta, gammast, Rd, &
                           Am, dAme(i,j))
                      ! repeat calculation if g0 > A*X
                      if (g0 .gt. Am*X) then
                         call fAndAn1(cs, g0, 0.0_r_2, gamma, beta, gammast, Rd, &
                              Am, dAme(i,j))
                      endif
                   endif

                   anrubpz(i,j) = real(Am)
                   if (Am > 0) eta_e(i,j) = dAme(i,j) * cs / Am

                endif ! C3, RuBp

                ! C4 RubP calculation
                g0 = gswminz(i,j)*fwsoilz(i) / C%RGSWC
                IF (vx4z(i,j) .gt. 1.0e-8) then
                   anrubpz(i,j) = vx4z(i,j) - rdxz(i,j)
                   dAme(i,j)    = 0.0_r_2
                   eta_e(i,j)   = 0.0_r_2
                ENDIF

                ! Sink limited:
                if (vcmxt3z(i,j).gt.1e-10) then ! C3

                   ansinkz(i,j) = 0.5 * vcmxt3z(i,j) - rdxz(i,j)
                   dAmp(i,j)    = 0.0_r_2
                   eta_p(i,j)   = 0.0_r_2

                elseif ((vcmxt4z(i,j) .gt. 1.0e-10) .and. (gs_coeffz(i,j) .gt. 100.)) then ! C4

                   coef2z(i,j) = real(gs_coeffz(i,j), r_2)
                   cs      = csxz(i,j)
                   if (trim(cable_user%g0_switch) == 'default') then
                      g0 = real(gswminz(i,j) * fwsoilz(i) / C%RGSWC, r_2)
                   elseif (trim(cable_user%g0_switch) == 'maximum') then
                      g0 = 0.0_r_2
                   endif
                   X       = real(gs_coeffz(i,j), r_2)
                   gamma   = real(effc4 * vcmxt4z(i,j), r_2)
                   Rd      = real(rdxz(i,j), r_2)

                   coef1z(i,j) = g0 + X * Rd + gamma * (1.0_r_2 - X * cs)
                   coef0z(i,j) = g0 * (Rd - cs * gamma)

                   ! no solution, give it a huge number
                   if ((abs(coef2z(i,j)) < 1.0e-9_r_2) .and. (abs(coef1z(i,j)) < 1.0e-9_r_2)) then
                      ciz(i,j)     = 99999.0_r_2
                      ansinkz(i,j) = real(ciz(i,j))
                   endif

                   ! solve linearly
                   if ( (abs(coef2z(i,j)) < 1.0e-9_r_2) .and. &
                        (abs(coef1z(i,j)) >= 1.0e-9_r_2) ) then
                      ciz(i,j)     = -1.0_r_2 * coef0z(i,j) / coef1z(i,j)
                      ansinkz(i,j) = real(ciz(i,j))
                   endif

                   ! solve quadratic (only take the more positive solution)
                   if (abs(coef2z(i,j)) >= 1.0e-9_r_2) then
                      delcxz(i,j)  = coef1z(i,j)**2 - 4.0_r_2*coef0z(i,j)*coef2z(i,j)
                      ciz(i,j)     = (-coef1z(i,j) + sqrt(max(0.0_r_2, delcxz(i,j)))) / &
                           (2.0_r_2 * coef2z(i,j))
                      ansinkz(i,j) = real(ciz(i,j))
                   endif

                   if (trim(cable_user%g0_switch) == 'maximum') then
                      ! repeat calculation of ansinkz, if g0 > A*X
                      if (g0 .gt. ciz(i,j)*X) then
                         X = 0.0_r_2
                         coef2z(i,j) = X
                         coef1z(i,j) = g0 + X * Rd + gamma  * (1.0_r_2 - X * cs)
                         coef0z(i,j) = -g0 * cs * gamma + Rd * g0


                         if ((abs(coef2z(i,j)) < 1.0e-9_r_2) .and. (abs(coef1z(i,j)) < 1.0e-9_r_2)) then
                            ciz(i,j)     = 99999.0_r_2
                            ansinkz(i,j) = real(ciz(i,j))
                         endif

                         ! solve linearly
                         if ( (abs(coef2z(i,j)) < 1.0e-9_r_2) .and. &
                              (abs(coef1z(i,j)) >= 1.0e-9_r_2) ) then
                            ciz(i,j)     = -1.0_r_2 * coef0z(i,j) / coef1z(i,j)
                            ansinkz(i,j) = real(ciz(i,j))
                         endif

                         ! solve quadratic (only take the more positive solution)
                         if (abs(coef2z(i,j)) >= 1.0e-9_r_2) then
                            delcxz(i,j)  = coef1z(i,j)**2 - 4.0_r_2*coef0z(i,j)*coef2z(i,j)
                            ciz(i,j)     = (-coef1z(i,j) + sqrt(max(0.0_r_2, delcxz(i,j)))) / &
                                 (2.0_r_2 * coef2z(i,j))
                            ansinkz(i,j) = real(ciz(i,j))
                         endif
                      endif ! g0 > A*X   
                   endif ! g0_switch==maximum

                   dAmp(i,j)  = 0.0_r_2
                   eta_p(i,j) = 0.0_r_2

                endif ! C3 / C4

             endif ! deltlfz > 0.1

             ! minimal of three limited rates
             ! anxz(i,j) = min(anrubiscoz(i,j), anrubpz(i,j), ansinkz(i,j))
             ! if (anxz(i,j) .eq. anrubiscoz(i,j)) then
             !    dA(i,j)  = dAmc(i,j)
             !    eta(i,j) = eta_c(i,j)
             ! elseif (anxz(i,j) .eq. anrubpz(i,j)) then
             !    dA(i,j)  = dAme(i,j)
             !    eta(i,j) = eta_e(i,j)
             ! elseif (anxz(i,j) .eq. ansinkz(i,j)) then
             !    dA(i,j)  = dAmp(i,j)
             !    eta(i,j) = eta_p(i,j)
             ! endif          
             tmp3      = (/ anrubiscoz(i,j), anrubpz(i,j), ansinkz(i,j) /)
             ii        = minloc(tmp3)
             anxz(i,j) = tmp3(ii(1))
             dtmp3     = (/ dAmc(i,j), dAme(i,j), dAmp(i,j) /)
             dA(i,j)   = dtmp3(ii(1))
             dtmp3     = (/ eta_c(i,j), eta_e(i,j), eta_p(i,j) /)
             eta(i,j)  = dtmp3(ii(1))

          else ! vlaiz > c%lai_thresh

             anxz(i,:)       = 0.0
             anrubiscoz(i,:) = 0.0
             anrubpz(i,:)    = 0.0
             ansinkz(i,:)    = 0.0
             dAmc(i,:)       = 0.0_r_2
             dAme(i,:)       = 0.0_r_2
             dAmp(i,:)       = 0.0_r_2
             dA(i,:)         = 0.0_r_2
             eta_c(i,:)      = 0.0_r_2
             eta_e(i,:)      = 0.0_r_2
             eta_p(i,:)      = 0.0_r_2
             eta(i,:)        = 0.0_r_2

          endif ! vlaiz > c%lai_thresh

       enddo ! j=1, mf

    enddo ! i=1, mp

  END SUBROUTINE photosynthesis


  !---------------------------------------------------------------------------------------


  FUNCTION ej3x(parx,alpha,convex,x) RESULT(z)

    REAL, INTENT(IN)     :: parx
    REAL, INTENT(IN)     :: alpha
    REAL, INTENT(IN)     :: convex
    REAL, INTENT(IN)     :: x
    REAL                 :: z

    z = MAX(0.0, &
         0.25*((alpha*parx+x-sqrt((alpha*parx+x)**2 - &
         4.0*convex*alpha*parx*x)) /(2.0*convex)) )
  END FUNCTION ej3x

  ! ------------------------------------------------------------------------------

  FUNCTION ej4x(parx,alpha,convex,x) RESULT(z)

    REAL, INTENT(IN)     :: parx ! Q mol photon m-2 s-1
    REAL, INTENT(IN)     :: alpha ! quantum efficiency in mol C (mol photon)-1
    REAL, INTENT(IN)     :: convex ! convexity parameter
    REAL, INTENT(IN)     :: x  ! vcmax
    REAL                 :: z  ! rubsico-limited gross photosynthesis


    z = MAX(0.0, &
         (alpha*parx+x-sqrt((alpha*parx+x)**2 - &
         4.0*convex*alpha*parx*x))/(2.0*convex))

  END FUNCTION ej4x

  ! ------------------------------------------------------------------------------

  ! Explicit array dimensions as temporary work around for NEC inlining problem
  FUNCTION xvcmxt4(x) RESULT(z)

    REAL, PARAMETER      :: q10c4 = 2.0
    REAL, INTENT(IN) :: x
    REAL :: z

!!$    z = q10c4 ** (0.1 * x - 2.5) / &
!!$         ((1.0 + exp(0.3 * (13.0 - x))) * (1.0 + exp(0.3 * (x - 36.0))))

    z = q10c4 ** (0.1 * x - 2.5) / &
         ((1.0 + exp(0.3 * (13.0 - x))) * (1.0 + exp(0.2 * (x - 38.0))))


  END FUNCTION xvcmxt4

  ! ------------------------------------------------------------------------------

  FUNCTION xgmesT(x) RESULT(z)

    USE cable_def_types_mod, only: r_2
    
    ! Temperature response of mesophyll conductance
    ! parameters as in Knauer et al. 2019 (from Bernacchi et al. 2002)
    REAL, INTENT(IN) :: x
    REAL(r_2)        :: z
    
    REAL(r_2) :: xgmes, TrefK, Rgas

    REAL(r_2), PARAMETER  :: EHa    = 49.6e3_r_2  ! J/mol 
    REAL(r_2), PARAMETER  :: EHd    = 437.4e3_r_2 ! J/mol 
    REAL(r_2), PARAMETER  :: Entrop = 1.4e3_r_2   ! J/mol/K
    CALL point2constants(C)

    TrefK = real(C%TrefK, r_2)
    Rgas  = real(C%Rgas, r_2)
    
    xgmes = exp(EHa * (x - TrefK) / (TrefK * Rgas * x )) * &
            (1.0_r_2 + exp((TrefK * Entrop - EHd) / (TrefK * Rgas))) / &
            (1.0_r_2 + exp((x * Entrop - EHd) / (x * Rgas)))

    z = max(0.0_r_2, xgmes)

  END FUNCTION xgmesT

  ! ------------------------------------------------------------------------------

  FUNCTION xvcmxt3(x) RESULT(z)

    !  Leuning 2002 (PCE) equation for temperature response
    !  used for Vcmax for c3 plants
    REAL, INTENT(IN) :: x
    REAL :: xvcnum,xvcden,z

    REAL, PARAMETER  :: EHaVc  = 73637.0  ! J/mol (Leuning 2002)
    REAL, PARAMETER  :: EHdVc  = 149252.0 ! J/mol (Leuning 2002)
    REAL, PARAMETER  :: EntropVc = 486.0  ! J/mol/K (Leuning 2002)
    REAL, PARAMETER  :: xVccoef = 1.17461 ! derived parameter
    ! xVccoef=1.0+exp((EntropJx*C%TREFK-EHdJx)/(Rconst*C%TREFK))
    CALL point2constants(C)
    xvcnum=xvccoef*exp( ( ehavc / ( C%rgas*C%TREFK ) )* ( 1.-C%TREFK/x ) )
    xvcden=1.0+exp( ( entropvc*x-ehdvc ) / ( C%rgas*x ) )
    z = max( 0.0, xvcnum / xvcden )


  END FUNCTION xvcmxt3

  ! ------------------------------------------------------------------------------
  elemental pure subroutine xvcmxt3_acclim(Tk, Tgrowth , trf)

    ! acclimated temperature response of Vcmax. Kumarathunge et al., New. Phyt., 2019,
    ! Eq 7 and Table 2
    REAL, INTENT(IN) :: Tk, Tgrowth  ! instantaneous T in K, growth T in degC
    REAL, INTENT(OUT) :: trf
    REAL:: xVccoef, EHaVc, EHdVc,  EntropVc, aKK, bKK, rgas, TREFK, xvcnum, xvcden


    EHaVc = 42.6 * 1000.0 + 1.14*Tgrowth*1000.0
    EHdVc  = 200000.0
    aKK = 645.13
    bKK = -0.38
    rgas = 8.314
    TREFK = 298.15

    entropvc = (aKK + bKK * Tgrowth)
    xVccoef  = 1.0 + exp((entropvc * TREFK - EHdVc)/  ( rgas*TREFK ) )
    xvcnum = xVccoef * exp( ( EHaVc / ( rgas*TREFK ) )* ( 1.-TREFK/Tk ) )
    xvcden=1.0 + exp( ( entropvc*Tk-EHdVc ) / ( rgas*Tk ) )
    trf = max( real(0.0), xvcnum / xvcden )

  end subroutine xvcmxt3_acclim
  ! ------------------------------------------------------------------------------
 REAL FUNCTION xrdt(x)

   !  Atkin et al. (Eq 1, New Phytologist (2015) 206: 614-636)
   !variable Q10 temperature of dark respiration
   ! Originally from Tjoelker et al. 2001

    REAL, INTENT(IN) :: x


    xrdt = (3.09 - 0.043*((x-273.15)+25.)/2.0)**((x-273.15 -25.0)/10.0)

  END FUNCTION xrdt

  ! ------------------------------------------------------------------------------

  FUNCTION xejmxt3(x) RESULT(z)

    ! Leuning 2002 (PCE) equation for temperature response
    ! used for jmax for C3 plants

    REAL, INTENT(IN) :: x
    REAL :: xjxnum,xjxden,z

    REAL, PARAMETER  :: EHaJx  = 50300.0  ! J/mol (Leuning 2002)
    REAL, PARAMETER  :: EHdJx  = 152044.0 ! J/mol (Leuning 2002)
    REAL, PARAMETER  :: EntropJx = 495.0  ! J/mol/K (Leuning 2002)
    REAL, PARAMETER  :: xjxcoef = 1.16715 ! derived parameter
    CALL point2constants(C)
    xjxnum = xjxcoef*exp( ( ehajx / ( C%rgas*C%TREFK ) ) * ( 1.-C%TREFK / x ) )
    xjxden=1.0+exp( ( entropjx*x-ehdjx) / ( C%rgas*x ) )
    z = max(0.0, xjxnum/xjxden)


  END FUNCTION xejmxt3

  ! ------------------------------------------------------------------------------
  elemental pure subroutine xejmxt3_acclim(Tk, Tgrowth, Thome, trf)

    ! acclimated temperature response of Jmax. Kumarathunge et al., New. Phyt., 2019,
    ! Eq 7 and Table 2
    REAL, INTENT(IN) :: Tk, Tgrowth, Thome  ! instantaneous T in K, home and growth T in degC
    REAL, INTENT(OUT) :: trf
    REAL:: xVccoef, EHaVc, EHdVc,  EntropVc, aKK, bKK, cKK, rgas, TREFK, xvcnum, xvcden


    EHaVc = 40.71 * 1000.0
    EHdVc  = 200000.0
    aKK = 658.77
    bKK = -0.84
    cKK = -0.52
    rgas = 8.314
    TREFK = 298.15
    
    entropvc = (aKK + bKK * Thome) + cKK * (Tgrowth - Thome)
    xVccoef  = 1.0 + exp((entropvc * TREFK - EHdVc)/  ( rgas*TREFK ) )
    xvcnum = xVccoef * exp( ( EHaVc / ( rgas*TREFK ) )* ( 1.-TREFK/Tk ) )
    xvcden=1.0 + exp( ( entropvc*Tk-EHdVc ) / ( rgas*Tk ) )
    trf = max( real(0.0), xvcnum / xvcden )

  end subroutine xejmxt3_acclim

  
  ! ------------------------------------------------------------------------------
  

  SUBROUTINE fwsoil_calc_std(fwsoil, soil, ssnow, veg)
    
    use cable_def_types_mod
    use cable_common_module, only : cable_user

    implicit none
    
    real, dimension(:),        intent(out)   :: fwsoil ! soil water modifier of stom. cond
    type(soil_parameter_type), intent(inout) :: soil
    type(soil_snow_type),      intent(inout) :: ssnow
    type(veg_parameter_type),  intent(inout) :: veg
    
    real, dimension(mp) :: rwater ! soil water availability

    rwater = max(1.0e-9, &
         sum(veg%froot * max(1.0e-9, min(1.0, real(ssnow%wb) - &
         spread(soil%swilt, 2, ms))),2) / (soil%sfc-soil%swilt))

   ! Remove vbeta #56
   if (cable_user%gs_switch == 'medlyn') then
      fwsoil = max(1.0e-4, min(1.0, rwater))
   else
      fwsoil = max(1.0e-9, min(1.0, veg%vbeta * rwater))
   endif

  END SUBROUTINE fwsoil_calc_std

  
  ! ------------------------------------------------------------------------------
  

  SUBROUTINE fwsoil_calc_non_linear(fwsoil, soil, ssnow, veg)
    
    use cable_def_types_mod

    implicit none
    
    real, dimension(:),        intent(out)   :: fwsoil ! soil water modifier of stom. cond
    type(soil_parameter_type), intent(inout) :: soil
    type(soil_snow_type),      intent(inout) :: ssnow
    type(veg_parameter_type),  intent(inout) :: veg
    
    real, dimension(mp) :: rwater ! soil water availability
    real, dimension(mp,3)          :: xi, ti, si
    integer :: j

    rwater = max(1.0e-9, &
         sum(veg%froot * max(0.0, min(1.0, real(ssnow%wb) - &
         spread(soil%swilt, 2, ms))),2) / (soil%sfc-soil%swilt))

    fwsoil = 1.

    rwater = soil%swilt + rwater * (soil%sfc-soil%swilt)

    xi(:,1) = soil%swilt
    xi(:,2) = soil%swilt + (soil%sfc - soil%swilt)/2.0
    xi(:,3) = soil%sfc

    ti(:,1) = 0.
    ti(:,2) = 0.9
    ti(:,3) = 1.0

    si(:,1) = (rwater - xi(:,2)) / ( xi(:,1) - xi(:,2)) * &
         (rwater - xi(:,3)) / ( xi(:,1) - xi(:,3))

    si(:,2) = (rwater - xi(:,1)) / ( xi(:,2) - xi(:,1)) * &
         (rwater - xi(:,3)) / ( xi(:,2) - xi(:,3))

    si(:,3) = (rwater - xi(:,1)) / ( xi(:,3) - xi(:,1)) * &
         (rwater - xi(:,2)) / ( xi(:,3) - xi(:,2))

    DO j=1, mp
       IF (rwater(j) < soil%sfc(j) - 0.02) &
            fwsoil(j) = max(0., min(1., ti(j,1)*si(j,1) + &
            ti(j,2)*si(j,2) + ti(j,3)*si(j,3)))
    ENDDO

  END SUBROUTINE fwsoil_calc_non_linear

  
  ! ------------------------------------------------------------------------------
  

  ! ypw 19/may/2010 soil water uptake efficiency (see Lai and Ktaul 2000)
  SUBROUTINE fwsoil_calc_Lai_Ktaul(fwsoil, soil, ssnow, veg)
    
    use cable_def_types_mod

    implicit none
    
    real, dimension(:),        intent(out)   :: fwsoil ! soil water modifier of stom. cond
    type(soil_parameter_type), intent(inout) :: soil
    type(soil_snow_type),      intent(inout) :: ssnow
    type(veg_parameter_type),  intent(inout) :: veg
    
    integer :: ns
    real, parameter :: rootgamma = 0.01   ! (19may2010)
    real, dimension(mp) :: dummy, normfac
    !--- local level dependent rwater
    real, dimension(mp,ms) :: frwater

    fwsoil(:)  = 0.0
    normFac(:) = 0.0

    DO ns=1, ms

       dummy(:) = rootgamma / max(1.0e-3,real(ssnow%wb(:,ns))-soil%swilt(:))

       frwater(:,ns) = max( 1.0e-4, &
            ((real(ssnow%wb(:,ns))-soil%swilt(:))/soil%ssat(:))**dummy )

       fwsoil(:) = min(1.0, max(fwsoil(:), frwater(:,ns)))

       normFac(:) = normFac(:) + frwater(:,ns) * veg%froot(:,ns)

    ENDDO

  END SUBROUTINE fwsoil_calc_Lai_Ktaul

  
  ! ------------------------------------------------------------------------------
  
  SUBROUTINE fwsoil_calc_sli(fwsoil, soil, ssnow, veg)
    
    use cable_def_types_mod

    implicit none
    
    real, dimension(:),        intent(out)   :: fwsoil ! soil water modifier of stom. cond
    type(soil_parameter_type), intent(inout) :: soil
    type(soil_snow_type),      intent(inout) :: ssnow
    type(veg_parameter_type),  intent(inout) :: veg
    
    real(r_2), dimension(mp,ms) :: tmp1, tmp2, delta_root, alpha2a_root, alpha2_root
    
    ! Lai and Katul formulation for root efficiency function - vh 17/07/09
    alpha2a_root = max(ssnow%wb-soil%swilt_vec, 0.001_r_2) / (soil%ssat_vec)
    tmp1         = ssnow%wb - soil%swilt_vec
    tmp2         = spread(veg%gamma,2,ms) / tmp1 * log(alpha2a_root)
    where ((tmp1 > 0.001_r_2) .and. (tmp2 > -10.0_r_2))
       alpha2_root = exp(tmp2)
    elsewhere
       alpha2_root = 0.0_r_2
    endwhere

    WHERE (veg%froot > 0.0)
       delta_root = 1.0_r_2
    ELSEWHERE
       delta_root = 0.0_r_2
    ENDWHERE

    fwsoil  = maxval(real(alpha2_root*delta_root), dim=2)
    fwsoil  = max(0.0, fwsoil)

  END SUBROUTINE fwsoil_calc_sli

  
  ! ------------------------------------------------------------------------------

  
  SUBROUTINE getrex_1d(theta, rex, fws, Fs, thetaS, thetaw, Etrans, gamma, dx, dt, zr)

    ! root extraction : Haverd et al. 2013
    USE cable_def_types_mod, only: r_2

    IMPLICIT NONE

    REAL(r_2), DIMENSION(:), INTENT(IN)    :: theta  ! volumetric soil moisture
    REAL(r_2), DIMENSION(:), INTENT(INOUT) :: rex    ! water extraction per layer
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
    REAL(r_2), DIMENSION(1:size(theta)) :: lthetar, alpha_root, delta_root, layer_depth
    REAL(r_2)                           :: trex, e3, one, zero
    INTEGER :: k

    e3   = 0.001_r_2
    one  = 1.0_r_2
    zero = 0.0_r_2

    layer_depth(1) = 0.0_r_2
    do k=2, size(theta)
       layer_depth(k) = sum(dx(1:k-1))
    enddo

    !theta(:)   = S(:)*thetaS(:)
    lthetar(:) = log(max(theta(:)-thetaw(:),e3)/thetaS(:))

    where ((theta(:)-thetaw(:)) > e3)
       alpha_root(:) = exp( gamma/max(theta(:)-thetaw(:), e3) * lthetar(:) )
    elsewhere
       alpha_root(:) = zero
    endwhere

    where (Fs(:) > zero .and. layer_depth < zr )  ! where there are roots and we are above max rooting depth
       delta_root(:) = one
    elsewhere
       delta_root(:) = zero
    endwhere

    rex(:) = alpha_root(:)*Fs(:)

    trex = sum(rex(:))
    if (trex > zero) then
       rex(:) = rex(:)/trex
    else
       rex(:) = zero
    endif
    rex(:) = Etrans*rex(:)

    ! reduce extraction efficiency where total extraction depletes soil moisture below wilting point
    where (((rex*dt) > (theta(:)-thetaw(:))*dx(:)) .and. ((rex*dt) > zero))
      alpha_root = alpha_root*(theta(:)-thetaw(:))*dx(:)/(1.1_r_2*rex*dt)
    endwhere
    rex(:) = alpha_root(:)*Fs(:)

    trex = sum(rex(:))
    if (trex > zero) then
       rex(:) = rex(:)/trex
    else
       rex(:) = zero
    endif
    rex(:) = Etrans*rex(:)

    ! check that the water available in each layer exceeds the extraction
    !if (any((rex*dt) > (theta(:)-0.01_r_2)*dx(:))) then
    if (any(((rex*dt) > max((theta(:)-thetaw(:)),zero)*dx(:)) .and. (Etrans > zero))) then
       !MC - ask Vanessa why fwsoil=0. if any (not all) demand exceeds water in layer.
       !     Should it be in else clause of trex>zero?
       !MC - also how is fwsoil coming back after it has been zero?
       !     Should it be fws=one here and fws=zero in else clause of trex>zero?
       fws = zero
       ! distribute extraction according to available water
       ! rex(:) = (theta(:)-0.01_r_2)*dx(:)
       rex(:) = max((theta(:)-thetaw(:))*dx(:),zero)
       trex = sum(rex(:))
       if (trex > zero) then
          rex(:) = rex(:)/trex
       else
          rex(:) = zero
       endif
       rex(:) = Etrans*rex(:)
    else
       fws    = maxval(alpha_root(2:)*delta_root(2:))
    endif

  END SUBROUTINE getrex_1d
!*****************************************************************************************
SUBROUTINE cubic_root_solver(a0,a1,a2,x1,x2,x3)
USE cable_def_types_mod, only: r_2
REAL(r_2), INTENT(IN) :: a0,a1,a2
REAL(r_2), INTENT(OUT) :: x1,x2,x3
REAL(r_2) :: Q, R, theta, a, b, c
real :: pi_c = 3.1415927
a = a2
b = a1
c = a0

Q = (a**2 - 3*b)/9
R = (2*a**3 - 9*a*b + 27*c)/54

if (R**2 .lt. Q**2) then
   theta = acos(R/(Q**3)**0.5)
   x1 = -2 * (Q**0.5)*cos(theta/3) - a/3
   x2 = -2 * (Q**0.5)*cos((theta+2*pi_c)/3) - a/3
   x3 = -2 * (Q**0.5)*cos((theta-2*pi_c)/3) - a/3

else
   x1 = 9999.0
   x2 = 9999.0
   x3 = 9999.0

endif



END SUBROUTINE cubic_root_solver
!*****************************************************************************************


! functions for infinite mesophyll conductance
elemental pure subroutine fabc(Cs, g0, x, gamma, beta, Gammastar, Rd, a, b, c)
  
  use cable_def_types_mod, only: r_2

  implicit none
  
  real(r_2), intent(in)  :: Cs, g0, x, gamma, beta, Gammastar, Rd
  real(r_2), intent(out) :: a, b, c
  
  a = (1.0_r_2-x)*Cs - x*beta
  b = -g0*Cs**2 + ((1.0_r_2-x)*(Rd-gamma)-g0*beta)*Cs - x*(gamma*Gammastar+Rd*beta)
  c = -g0*(Rd-gamma)*Cs**2 - g0*(gamma*Gammastar+Rd*beta)*Cs
  
end subroutine fabc


elemental pure subroutine fabc2(Cs, g0, x, gamma, beta, Gammastar, Rd, a, b, c)
  
  use cable_def_types_mod, only: r_2

  implicit none
  
  real(r_2), intent(in)  :: Cs, g0, x, gamma, beta, Gammastar, Rd
  real(r_2), intent(out) :: a, b, c
  
  a = x
  b = (g0+gamma*(1.0_r_2-x))*Cs - x*(beta-Rd)
  c =  -gamma*g0*Cs**2 - g0*(beta-Rd)*Cs
  
end subroutine fabc2


elemental pure subroutine fAn(a, b, c, A2)
  
  use cable_def_types_mod, only: r_2

  implicit none
  
  real(r_2), intent(in) :: a,b,c
  real(r_2), intent(out) :: A2
  
  real(r_2) :: s2

  s2 = b**2 - 4.0_r_2*a*c
  A2 = (-b - sqrt(s2))/(2.0_r_2*a)
  
end subroutine fAn


elemental pure subroutine fAn2(a, b, c, A2)
  
  use cable_def_types_mod, only: r_2

  implicit none
  
  real(r_2), intent(in)  :: a, b, c
  real(r_2), intent(out) :: A2
  
  real(r_2) :: s2
  
  s2 = b**2 - 4.0_r_2*a*c
  A2 = (-b + sqrt(s2))/(2.0_r_2*a)
  
end subroutine fAn2


elemental pure subroutine fdabc(Cs, g0, x, gamma, beta, Gammastar, Rd, da, db, dc)
  
  use cable_def_types_mod, only: r_2

  implicit none
  
  real(r_2), intent(in)  :: Cs, g0, x, gamma, beta, Gammastar, Rd
  real(r_2), intent(out) :: da, db, dc

  da = 1.0_r_2-x
  db = -2.0_r_2*g0*Cs + (1.0_r_2-x)*(Rd-gamma) - g0*beta
  dc = -2.0_r_2*g0*(Rd-gamma)*Cs - g0*(gamma*Gammastar+Rd*beta)
  
end subroutine fdabc


elemental pure subroutine fdabc2(Cs, g0, x, gamma, beta, Gammastar, Rd, da, db, dc)
  
  use cable_def_types_mod, only: r_2

  implicit none
  
  real(r_2), intent(in)  :: Cs, g0, x, gamma, beta, Gammastar, Rd
  real(r_2), intent(out) :: da, db, dc
  
  da = 0.
  db = g0 + gamma*(1.0_r_2-x)
  dc = -2.0_r_2*gamma*g0*Cs - g0*(beta-Rd)
  
end subroutine fdabc2


elemental pure subroutine fdAn(a, b, c, da, db, dc, dA2)
  
  use cable_def_types_mod, only: r_2

  implicit none
  
  real(r_2), intent(in) :: a, b, c, da, db, dc
  real(r_2), intent(out) :: dA2
  
  real(r_2) :: s, p
  
  s = sqrt(b**2 - 4.0_r_2*a*c)
  p = (2.0_r_2*b*db - 4.0_r_2*c*da - 4.0_r_2*a*dc)/(2.0_r_2*s)
  dA2 = (-db - p)/(2.0_r_2*a) - (-b - s)/(2.0_r_2*a**2)*da
  
end subroutine fdAn


elemental pure subroutine fdAn2(a, b, c, da, db, dc, dA2)
  
  use cable_def_types_mod, only: r_2

  implicit none
  
  real(r_2), intent(in)  :: a, b, c, da, db, dc
  real(r_2), intent(out) :: dA2

  real(r_2) :: s, p
  
  s = sqrt(b**2 - 4.0_r_2*a*c)
  p = (2.0_r_2*b*db - 4.0_r_2*c*da - 4.0_r_2*a*dc)/(2.0_r_2*s)
  dA2 = (-db + p)/(2.0_r_2*a) - (-b + s)/(2.0_r_2*a**2)*da
  
end subroutine fdAn2


elemental pure subroutine fAndAn1(Cs, g0, x, gamma, beta, Gammastar, Rd, An, dAn)
  
  use cable_def_types_mod, only: r_2

  implicit none
  
  real(r_2), intent(in)  :: Cs, g0, x, gamma, beta, Gammastar, Rd
  real(r_2), intent(out) :: An, dAn
  
  real(r_2) :: a, b, c, da, db, dc
  
  call fabc(Cs, g0, x, gamma, beta, Gammastar, Rd,a,b,c)
  call fAn(a, b, c, An)
  call fdabc(Cs, g0, x, gamma, beta, Gammastar, Rd, da, db, dc)
  call fdAn(a, b, c, da, db, dc, dAn)
  
end subroutine fAndAn1


elemental pure subroutine fAndAn2(Cs, g0, x, gamma, beta, Gammastar, Rd, An, dAn)
  
  use cable_def_types_mod, only: r_2

  implicit none
  
  real(r_2), intent(in)  :: Cs, g0, x, gamma, beta, Gammastar, Rd
  real(r_2), intent(out) :: An, dAn
  
  real(r_2) :: a, b, c, da, db, dc
  
  call fabc2(Cs, g0, x, gamma, beta, Gammastar, Rd,a,b,c)
  call fAn2(a, b, c, An)
  call fdabc2(Cs, g0, x, gamma, beta, Gammastar, Rd, da, db, dc)
  call fdAn2(a, b, c, da, db, dc, dAn)
  
end subroutine fAndAn2


elemental pure subroutine fAndAn(Cs, g0, x, Gammastar, Rd, gammac, betac, gammae, betae, &
     flag_eps, &
     Anc, Ane, An, dAnc, dAne, dAn)

  use cable_def_types_mod, only: r_2

  implicit none
  
  real(r_2), intent(in)  :: Cs, g0, x, Gammastar, Rd, gammac, betac, gammae, betae
  logical,   intent(in)  :: flag_eps
  real(r_2), intent(out) :: Anc, Ane, An, dAnc, dAne, dAn
  
  call fAndAn1(Cs, g0, x, gammac, betac, Gammastar, Rd, Anc, dAnc)
  call fAndAn1(Cs, g0, x, gammae, betae, Gammastar, Rd, Ane, dAne)

  ! An = min(Anc,Ane)
  ! if (An==Anc) then
  !    dAn = dAnc
  ! else
  !    dAn = dAne
  ! endif
  if (Anc < Ane) then
     An  = Anc
     dAn = dAnc
  else
     An  = Ane
     dAn = dAne
  endif

  if (flag_eps) then !  return elasticity instead of derivative
     dAnc = dAnc * Cs / Anc
     dAne = dAne * Cs / Ane
     dAn  = dAn  * Cs / An
  endif
  
end subroutine fAndAn


! elemental pure subroutines for finite mesophyll conductance   
elemental pure subroutine fabcd(Cs, g0, x, gamma, beta, Gammastar, Rd, gm, a, b, c1, d)
  
  use cable_def_types_mod, only: r_2

  implicit none

  real(r_2), intent(in)  :: Cs, g0, x, gamma, beta, Gammastar, Rd, gm
  real(r_2), intent(out) :: a, b, c1, d
  
  a  = x
  b  = (gm+g0-gm*x)*Cs + x*(Rd-gamma) - gm*x*beta
  c1 = -gm*g0*Cs**2 + ((gm+g0-gm*x)*(Rd-gamma)-gm*g0*beta)*Cs - &
       gm*x*(gamma*Gammastar+Rd*beta)
  d  = -gm*g0*(Rd-gamma)*Cs**2 - gm*g0*(gamma*Gammastar+Rd*beta)*Cs
  
end subroutine fabcd


elemental pure subroutine fabcm(Cs, g0, x, gamma, beta, Gammastar, Rd, gm, a, b, c1)
  
  use cable_def_types_mod, only: r_2

  implicit none

  real(r_2), intent(in)  :: Cs, g0, x, gamma, beta, Gammastar, Rd, gm
  real(r_2), intent(out) :: a,b,c1
  
  a  = x*(gm+gamma)
  b  = ((gm+g0)*gamma + g0*gm - gm*gamma*x)*Cs - gm*x*(beta-Rd)
  c1 = -gm*g0*gamma*Cs**2 - gm*g0*(beta-Rd)*Cs
  
end subroutine fabcm


elemental pure subroutine fpq(a, b, c, d, p, q)
  
  use cable_def_types_mod, only: r_2

  implicit none

  real(r_2), intent(in) :: a,b,c,d
  real(r_2), intent(out) :: p, q
  
  p = (3.0_r_2*a*c - b**2)/(3.0_r_2*a**2)
  q = (2.0_r_2*b**3 - 9.0_r_2*a*b*c + 27.0_r_2*a**2*d)/(27.0_r_2*a**3)
  
end subroutine fpq


elemental pure subroutine fdabcd(Cs, g0, x, gamma, beta, Gammastar, Rd, gm, da, db, dc, dd)
  
  use cable_def_types_mod, only: r_2

  implicit none

  real(r_2), intent(in)  :: Cs, g0, x, gamma, beta, Gammastar, Rd, gm
  real(r_2), intent(out) :: da, db, dc, dd
  
  da = 0.0_r_2
  db = gm+g0-gm*x
  dc = -2.0_r_2*gm*g0*Cs + (gm+g0-gm*x)*(Rd-gamma)-gm*g0*beta
  dd = -2.0_r_2*gm*g0*(Rd-gamma)*Cs - gm*g0*(gamma*Gammastar+Rd*beta)
  
end subroutine fdabcd


elemental pure subroutine fdabcm(Cs, g0, x, gamma, beta, Gammastar, Rd, gm, da, db, dc)
  
  use cable_def_types_mod, only: r_2

  implicit none
  
  real(r_2), intent(in)  :: Cs, g0, x, gamma, beta, Gammastar, Rd, gm
  real(r_2), intent(out) :: da, db, dc
  
  da = 0.0_r_2
  db = (gm+g0)*gamma + g0*gm - gm*gamma*x
  dc = -2.0_r_2*gm*g0*gamma*Cs - gm*g0*(beta-Rd)
  
end subroutine fdabcm

elemental pure subroutine fdpq(a, b, c, d, da, db, dc, dd, dp, dq)
  
  use cable_def_types_mod, only: r_2

  implicit none

  real(r_2), intent(in)  :: a, b, c, d, da, db, dc, dd
  real(r_2), intent(out) :: dp, dq
  
  dp = (3.0_r_2*da*c + 3.0_r_2*a*dc - 2.0_r_2*b*db)/(3.0_r_2*a**2) - &
       2.0_r_2*(3.0_r_2*a*c - b**2)/(3.0_r_2*a**3)*da
  dq = (6.0_r_2*b**2*db - 9.0_r_2*da*b*c - 9.0_r_2*a*db*c - 9.0_r_2*a*b*dc + &
       54.0_r_2*a*da*d + 27.0_r_2*a**2*dd) / (27.0_r_2*a**3) - &
       3.0_r_2*(2.0_r_2*b**3 - 9.0_r_2*a*b*c + 27.0_r_2*a**2*d) / &
       (27.0_r_2*a**4)*da
  
end subroutine fdpq


elemental pure subroutine fAm(a, b, c1, d, p, q, Am)
  
  use cable_def_types_mod, only: r_2
  use mo_constants, only: pi => pi_dp

  implicit none

  real(r_2), intent(in)  :: a, b, c1, d, p, q
  real(r_2), intent(out) :: Am
  
  real(r_2) :: p3, pq, k
 
  p3 = -p/3.0_r_2
  pq = min(3.0_r_2*q/(2.0_r_2*p)*sqrt(1.0_r_2/p3),0.999999999999_r_2)
  k  = 1.0_r_2
  Am = 2.0_r_2*sqrt(p3)*cos(acos(pq)/3.0_r_2 - 2.0_r_2*pi*k/3.0_r_2) - b/(3.0_r_2*a)
  
end subroutine fAm


elemental pure subroutine fAm2(a, b, c1, Am)
  
  use cable_def_types_mod, only: r_2
  
  implicit none

  real(r_2), intent(in) :: a, b, c1
  real(r_2), intent(out) :: Am
  
  real(r_2) :: s2
  
  s2 = b**2 - 4.0_r_2*a*c1
  Am = (-b + sqrt(s2))/(2.0_r_2*a)

end subroutine fAm2


elemental pure subroutine fdAm(a, b, c1, d, p, q, da, db, dc, dd, dp, dq, dAm)
  
  use cable_def_types_mod, only: r_2
  use mo_constants, only: pi => pi_dp

  implicit none

  real(r_2), intent(in)  :: a, b, c1, d, p, q, da, db, dc, dd, dp, dq
  real(r_2), intent(out) :: dAm
  
  real(r_2) :: k, p3, pq
  
  p3 = -p/3.0_r_2
  pq = min(3.0_r_2*q/(2.0_r_2*p)*sqrt(1.0_r_2/p3),0.999999999999_r_2)
  k  = 1.0_r_2
  dAm = -1.0_r_2/(3.0_r_2*sqrt(p3))*cos(acos(pq)/3.0_r_2 - 2.0_r_2*pi*k/3.0_r_2)*dp + &
       2.0_r_2*sqrt(p3) * (sin(acos(pq)/3.0_r_2 - 2.0_r_2*pi*k/3.0_r_2) / &
       (3.0_r_2*sqrt(1.0_r_2-pq**2)) * &
       3.0_r_2/2.0_r_2*(dq/p*sqrt(1.0_r_2/p3) - q/p**2*dp*sqrt(1.0_r_2/p3) + &
       q/p*3.0_r_2/2.0_r_2*1.0_r_2/sqrt(1.0_r_2/p3)*1.0_r_2/p**2*dp)) - &
       db/(3.0_r_2*a) + b/(3.0_r_2*a**2)*da

end subroutine fdAm


elemental pure subroutine fdAm2(a, b, c1, da, db, dc, dAm)
  
  use cable_def_types_mod, only: r_2
  
  implicit none

  real(r_2), intent(in) :: a, b, c1, da, db, dc
  real(r_2), intent(out) :: dAm
  
  real(r_2) :: s, p

  s = sqrt(b**2 - 4.0_r_2*a*c1)
  p = (2.0_r_2*b*db - 4.0_r_2*c1*da - 4.0_r_2*a*dc)/(2.0_r_2*s)
  dAm = (-db + p)/(2.0_r_2*a) - (-b + s)/(2.0_r_2*a**2)*da

end subroutine fdAm2


elemental pure subroutine fAmdAm1(Cs, g0, x, gamma, beta, Gammastar, Rd, gm, &
     Am, dAm)
  
  use cable_def_types_mod, only: r_2

  implicit none
  
  real(r_2), intent(in)  :: Cs, g0, x, gamma, beta, Gammastar, Rd, gm
  real(r_2), intent(out) :: Am, dAm
  
  real(r_2) :: a, b, c, d, p, q, da, db, dc, dd, dp, dq
  
  call fabcd(Cs, g0, x, gamma, beta, Gammastar, Rd, gm, a, b, c, d)
  call fpq(a, b, c, d, p, q)
  call fAm(a, b, c, d, p, q, Am)
  call fdabcd(Cs, g0, x, gamma, beta, Gammastar, Rd, gm, da, db, dc, dd)
  call fdpq(a, b, c, d, da, db, dc, dd, dp, dq)
  call fdAm(a, b, c, d, p, q, da, db, dc, dd, dp, dq, dAm)

end subroutine fAmdAm1


elemental pure subroutine fAmdAm2(Cs, g0, x, gamma, beta, Gammastar, Rd, gm, &
     Am, dAm)
  
  use cable_def_types_mod, only: r_2

  implicit none
  
  real(r_2), intent(in)  :: Cs, g0, x, gamma, beta, Gammastar, Rd, gm
  real(r_2), intent(out) :: Am, dAm
  
  real(r_2) :: a, b, c1, da, db, dc
  
  call fabcm(Cs, g0, x, gamma, beta, Gammastar, Rd, gm,a,b,c1)
  call fAm2(a, b, c1, Am)
  call fdabcm(Cs, g0, x, gamma, beta, Gammastar, Rd, gm, da, db, dc)
  call fdAm2(a, b, c1, da, db, dc, dAm)
  
end subroutine fAmdAm2


elemental pure subroutine  fAmdAm(Cs, g0, x, Gammastar, Rd, gm, gammac, betac, gammae, betae, &
     flag_eps, &
     Amc, Ame, Am, dAmc, dAme, dAm)
  
  use cable_def_types_mod, only: r_2

  implicit none
  
  real(r_2), intent(in)  :: Cs, g0, x, Gammastar, &
       Rd, gm, gammac, betac,gammae, betae
  logical,   intent(in)  :: flag_eps
  real(r_2), intent(out) :: Amc, Ame, Am, dAmc, dAme, dAm
  
  call fAmdAm1(Cs, g0, x, gammac, betac, Gammastar, Rd, gm, Amc, dAmc)
  call fAmdAm1(Cs, g0, x, gammae, betae, Gammastar, Rd, gm, Ame, dAme)

  ! Am = min(Amc,Ame)
  ! if (Am==Amc) then
  !    dAm = dAmc
  ! else
  !    dAm = dAme
  ! endif
  if (Amc < Ame) then
     Am  = Amc
     dAm = dAmc
  else
     Am  = Ame
     dAm = dAme
  endif
  
  if (flag_eps)  then !  return elasticity instead of derivative
     dAmc = dAmc * Cs / Amc
     dAme = dAme * Cs / Ame
     dAm  = dAm  * Cs / Am
  endif
  
end subroutine fAmdAm

!**********

END MODULE cable_canopy_module
