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
!        : Vanessa Haverd added new logical switch cable_user%litter. When 'true',
!          leaf litter suppresses soil evaporation
!        : Vanessa Haverd added new switch to enable SLI alternative to default soil
!          module. When cable_user%soil_struc=='sli', then SLI is used to compute
!          coupled transfers of heat and water in the soil and snow and at the surface
!          and an in-canopy stability correction is applied.
!        : See http://www.geosci-model-dev.net/9/3111/2016 for full documentation
!          of last 3 changes.
! ==============================================================================

MODULE cable_canopy_module

   use cable_data_module,    only: icanopy_type, point2constants
   use cable_common_module,  only: cable_user
   use cable_veg_hydraulics_module, only: optimisation, arrh, get_xylem_vulnerability, calc_plc
   USE cable_IO_vars_module, ONLY: logn
   USE cable_soil_hydraulics_module, ONLY : calc_soil_root_resistance, calc_swp, calc_psix
   implicit none

   public :: define_canopy, xvcmxt3, xejmxt3, ej3x, xrdt, xgmesT, &
      xvcmxt3_acclim, xejmxt3_acclim, light_inhibition

   private

   type(icanopy_type) :: C

CONTAINS

   SUBROUTINE define_canopy(ktau, ktau_tot, bal, rad, rough, air, met, dels, ssnow, soil, veg, canopy, climate, casapool, casabiome)

      USE cable_def_types_mod
      USE cable_radiation_module
      USE cable_air_module
      USE cable_common_module
      USE cable_roughness_module
      use cable_sli_main, only: sli_main
      use mo_utils,       only: eq
      USE casavariable
      implicit none
      INTEGER,                   INTENT(IN)    :: ktau
      INTEGER,                   INTENT(IN)    :: ktau_tot
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
      TYPE (casa_pool),  INTENT(IN)           :: casapool
      TYPE (casa_biome),  INTENT(IN)           :: casabiome
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
         ftemp, rlower_limit, &
         term1, term2, term3, term5

      REAL, DIMENSION(:), POINTER :: &
         cansat => null(),        & ! max canopy intercept. (mm)
         dsx => null(),           & ! leaf surface vpd
         dsxpsdo => null(),           & 
         dsy => null(),           & 
         dsypsdo => null(),           & 
         fwsoil => null(),        & ! soil water modifier of stom. cond
         fwsoilpsdo => null(),        &
         fwsoiltmp => null(),        & 
         fwsoiltmppsdo => null(),        & 
         tlfx => null(),          & ! leaf temp prev. iter (K)
         tlfxpsdo => null(),          &
         tlfy => null(),            & ! leaf temp (K)
         tlfypsdo => null()

      REAL(r_2), DIMENSION(mp) :: &
         gbvtop                   ! bnd layer cond. top leaf

      REAL(r_2), DIMENSION(:), POINTER :: &
         ecy => null(),           & ! lat heat fl dry big leaf
         ecypsdo => null(),           &
         hcy => null(),           & ! veg. sens heat
         hcypsdo => null(),           &
         rny => null(),           & ! net rad
         rnypsdo => null(),           &
         ghwet => null(), &           ! cond for heat for a wet canopy
         vpdpsdo => null()
      REAL(r_2) :: fwpsdo


      REAL(r_2), DIMENSION(:,:), POINTER :: &
         gbhu => null(),          & ! forcedConvectionBndryLayerCond
         gbhf => null(),          & ! freeConvectionBndryLayerCond
         csx => null(),           &   ! leaf surface CO2 concentration
         csxpsdo => null(),           &
         psily => null() ,        &
         psilypsdo => null() ,        &
         psilx => null(),         &
         psilxpsdo => null(),         &
         fwpsi => null(), &      ! psi leaf modifier of stom. cond zihan lu 02/11/2024
         fwpsipsdo => null(), &  
         wbpsdo => null()

      REAL :: rt_min
      REAL, DIMENSION(mp) :: zstar, rL, phist, csw, psihat,rt0bus

      INTEGER :: j

      INTEGER, SAVE :: call_number =0
      REAL, DIMENSION(ms) :: root_length

      ! END header

      call_number = call_number + 1

      ! assign local ptrs to constants defined in cable_data_module
      CALL point2constants(C)

      IF (.NOT. cable_runtime%um) canopy%cansto = canopy%oldcansto

      ALLOCATE(cansat(mp), gbhu(mp,mf))
      ALLOCATE(dsx(mp),dsxpsdo(mp), dsy(mp),dsypsdo(mp),fwsoil(mp), fwsoilpsdo(mp),fwsoiltmp(mp),fwsoiltmppsdo(mp), &
            tlfx(mp), tlfxpsdo(mp),tlfy(mp),tlfypsdo(mp))
      ALLOCATE(ecy(mp), ecypsdo(mp), hcy(mp), hcypsdo(mp), rny(mp), rnypsdo(mp))
      ALLOCATE(gbhf(mp,mf), csx(mp,mf), csxpsdo(mp,mf), psilx(mp,mf), psilxpsdo(mp,mf), psily(mp,mf), psilypsdo(mp,mf), &
       fwpsi(mp,mf), fwpsipsdo(mp,mf))
      ALLOCATE(ghwet(mp))
      ALLOCATE(vpdpsdo(mp))
      ALLOCATE(wbpsdo(mp,ms))

      canopy%N_neg = 0.0
      canopy%N_pos = 0.0
      canopy%N_neg_sw = 0.0
      canopy%N_pos_sw = 0.0
      ! BATS-type canopy saturation proportional to LAI:
      cansat = veg%canst1 * canopy%vlaiw

      !---compute surface wetness factor, update cansto, through
      ! print*, 'SS01.1 ', ssnow%wb(:,1)
      ! print*, 'SS01.2 ', soil%swilt
      ! print*, 'SS01.3 ', soil%sfc
      ! print*, 'SS01.4 ', ssnow%wbice(:,1)
      ! print*, 'SS01.5 ', ssnow%snowd(:)
      ! print*, 'SS01.6 ', ssnow%owetfac
      ! print*, 'SS01.7 ', ssnow%wetfac
      CALL surf_wetness_fact(cansat, canopy, ssnow, veg, met, soil, dels)
      ! print*, 'SS02.1 ', ssnow%wb(:,1)
      ! print*, 'SS02.2 ', soil%swilt
      ! print*, 'SS02.3 ', soil%sfc
      ! print*, 'SS02.4 ', ssnow%wbice(:,1)
      ! print*, 'SS02.5 ', ssnow%snowd(:)
      ! print*, 'SS02.6 ', ssnow%owetfac
      ! print*, 'SS02.7 ', ssnow%wetfac

      canopy%fevw_pot = 0.0
      canopy%gswx     = 1.0e-3     ! default stomatal conuctance
      gbhf            = 1.0e-3_r_2 ! default free convection boundary layer conductance
      gbhu            = 1.0e-3_r_2 ! default forced convection boundary layer conductance
      ssnow%evapfbl   = 0.0
      ssnow%rex       = 0.0_r_2
      ! Initialise in-canopy temperatures and humidity:
      csx = SPREAD(real(met%ca,r_2), 2, mf) ! initialise leaf surface CO2 concentration
      csxpsdo = csx
      if (ktau_tot==1) then
         !canopy%psi_can = SPREAD(real(ssnow%psi_rootzone,r_2), 2, mf)
         canopy%psi_can = SPREAD(real(canopy%psix,r_2), 2, mf)
      endif
      psilx = canopy%psi_can  ! SPREAD(real(ssnow%psi_rootzone,r_2), 2, mf)
      
      psily = canopy%psi_can  ! SPREAD(real(ssnow%psi_rootzone,r_2), 2, mf)
      
      if (INDEX(cable_user%FWSOIL_SWITCH,'LWP')>0) then
         fwpsi(:, 1) = (1.0_r_2 +exp(veg%slope_leaf(:) * veg%psi_50_leaf(:))) / &
         (1.0_r_2 + exp(veg%slope_leaf(:) * (veg%psi_50_leaf(:) - psilx(:, 1))))
         fwpsi(:, 2) = (1.0_r_2 +exp(veg%slope_leaf(:) * veg%psi_50_leaf(:))) / &
         (1.0_r_2 + exp(veg%slope_leaf(:) * (veg%psi_50_leaf(:) - psilx(:, 2))))
      else
         fwpsi = 1.0_r_2
      endif

      fwpsipsdo = fwpsi
      !print*, 'Entry psi ', psilx(:, 1)
      !print*, 'Entry fwpsi ', fwpsi(:, 1)

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
      dsy = dsx
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

      ! print*, 'SS03.1 ', niter, ssnow%rtsoil
      ! print*, 'SS03.2 ', ssnow%tss
      DO iter = 1, NITER
         ! AERODYNAMIC PROPERTIES: friction velocity us, thence turbulent
         ! resistances rt0, rt1 (elements of dispersion matrix):
         ! See CSIRO SCAM, Raupach et al 1997, eq. 3.46:
         CALL comp_friction_vel()

         ! print*, 'RR01 ', rough%zref_tq
         ! print*, 'RR02 ', rough%disp
         ! print*, 'RR03 ', rough%zruffs
         ! print*, 'RR04 ', rough%z0soilsn
         ! print*, 'RR05 ', rough%hruff
         ! E.Kowalczyk 2014
         IF (cable_user%l_new_roughness_soil) &
            CALL ruff_resist(veg, rough, ssnow, canopy)

         ! Turbulent aerodynamic resistance from roughness sublayer depth
         ! to reference height, x=1 if zref+disp>zruffs,
         ! 0 otherwise: thus rt1usc = 0 if zref+disp<zruffs
         ! See CSIRO SCAM, Raupach et al 1997, eq. 3.53:
         xx = 0.5 + sign(0.5, rough%zref_tq + rough%disp - rough%zruffs)
         ! print*, 'RR06 ', rough%zref_tq
         ! print*, 'RR07 ', rough%disp
         ! print*, 'RR08 ', rough%zruffs
         ! print*, 'RR09 ', rough%z0soilsn
         ! print*, 'RR10 ', rough%hruff

         ! correction by Ian Harman to the 2nd psis term
         rt1usc = xx * ( LOG( rough%zref_tq / MAX(rough%zruffs-rough%disp, &
            rough%z0soilsn) ) &
            - psis(canopy%zetar(:,iter)) &
            + psis(canopy%zetar(:,iter) * MAX(rough%zruffs-rough%disp, &
            rough%z0soilsn) / rough%zref_tq) ) / C%VONK
         rt_min = 5.0
         ! print*, 'RR11 ', canopy%zetar(:,iter)
         ! print*, 'RR12 ', psis(canopy%zetar(:,iter))

         !! vh_js !!
         IF (cable_user%soil_struc=='sli') THEN
            ! for stable conditions, update rough%rt0us & rough%rt1usa by
            ! replacing C%CSW by csw = cd/2* (U(hc)/ust)**2 according to Eqs 15 &
            ! 19 from notes by Ian Harman (9-9-2011)
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

               where (canopy%zetar(:,iter) > 1.0e-6)! stable conditions

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
            ! print*, 'RR13 ', rough%rt0us
            ! print*, 'RR14 ', canopy%us

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
         ! print*, 'RR15 ', rough%rt1usa
         ! print*, 'RR16 ', rough%rt1usb
         ! print*, 'RR17 ', rough%rt1

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
               gbhu(j,2) = 2.0_r_2 / real(rough%coexp(j),r_2) * gbvtop(j) * &
                  (1.0_r_2-EXP(-real( min(0.5*rough%coexp(j)*canopy%vlaiw(j), 20.0), r_2))) &
                  - gbhu(j,1)

               if (cable_user%amphistomatous) then
                  gbhu(j,1) = gbhu(j,1) * 2.0_r_2
                  gbhu(j,2) = gbhu(j,2) * 2.0_r_2
               endif

            ENDIF ! canopy%vlaiw(j) > C%LAI_THRESH

         ENDDO ! j=1, mp
         ! print*, 'RR19 ', air%cmolar
         ! print*, 'RR20 ', air%visc
         ! print*, 'RR21 ', rough%usuh
         ! print*, 'RR22 ', veg%dleaf
         ! print*, 'RR23 ', veg%shelrb
         ! print*, 'RR24 ', canopy%vlaiw
         ! print*, 'RR25 ', rough%coexp
         ! print*, 'RR26 ', rad%extkb

         rny = sum(real(rad%rniso,r_2), 2) ! init current estimate net rad
         rnypsdo = rny
         hcy = 0.0_r_2          ! init current estimate sensible heat
         hcypsdo = hcy
         ecy = rny - hcy        ! init current estimate latent heat
         ecypsdo = ecy

         sum_rad_rniso = sum(rad%rniso,2)
         ! !!!!!!!!!!!!!!! vpd = 0.6Kpa!!!!!!!!!!!!!!!!!!!!!
         ! if (iter==4) then
         !    vpdpsdo = SPREAD(real(600,r_2), 1, mp) ! 600Pa
         !    dsxpsdo = vpdpsdo
         !    dsypsdo = vpdpsdo
         !    psilxpsdo = psilx
         !    psilypsdo = psily
         !    fwsoilpsdo = fwsoil
         !    fwsoiltmppsdo = fwsoiltmp
         !    fwpsipsdo = fwpsi
         !    tlfxpsdo = tlfx
         !    tlfypsdo = tlfy
         !    ecypsdo = ecy
         !    hcypsdo = hcy
         !    rnypsdo = rny
         !    csxpsdo = csx
         !    CALL dryLeaf(ktau, ktau_tot,dels, rad, air, met,  &
         !    veg, canopy, soil, ssnow, casapool, dsxpsdo, dsypsdo, psilxpsdo, psilypsdo,&
         !    fwsoilpsdo, fwsoiltmppsdo, fwpsipsdo, tlfxpsdo, tlfypsdo, ecypsdo, hcypsdo,  &
         !    rnypsdo, gbhu, gbhf, csxpsdo, cansat,  &
         !    ghwet, iter, climate, vpdpsdo=vpdpsdo)
         !    dsxpsdo = dsx
         !    dsypsdo = dsy
         !    ! psilxpsdo = psilx
         !    ! psilypsdo = psily
         !    fwsoilpsdo = fwsoil
         !    ! fwsoiltmppsdo = fwsoiltmp
         !    ! fwpsipsdo = fwpsi
         !    tlfxpsdo = tlfx
         !    tlfypsdo = tlfy
         !    ecypsdo = ecy
         !    hcypsdo = hcy
         !    rnypsdo = rny
         !    csxpsdo = csx
         !    CALL dryLeaf_givengs(ktau, ktau_tot, dels, rad, air, met, &
         !       veg, canopy, soil, ssnow,casapool, dsxpsdo, dsypsdo, &
         !       tlfxpsdo, tlfypsdo, ecypsdo, hcypsdo,&
         !       rnypsdo, gbhu, gbhf, csxpsdo, &
         !       cansat, iter, climate,gspsdo = canopy%gsw_epotvpd)
         ! endif
         if (iter==4) then
            wbpsdo = SPREAD(real(soil%ssat,r_2), 2, ms) 
            DO j = 1, mp
               CALL calc_soil_root_resistance(ssnow, soil, veg, casapool, casabiome, root_length, j, wbpsdo)
               CALL calc_swp(ssnow, soil, j, wbpsdo)
            END DO
            dsxpsdo = dsx
            dsypsdo = dsy
            psilxpsdo = psilx
            psilypsdo = psily
            fwsoilpsdo = fwsoil
            fwsoiltmppsdo = fwsoiltmp
            fwpsipsdo = fwpsi
            tlfxpsdo = tlfx
            tlfypsdo = tlfy
            ecypsdo = ecy
            hcypsdo = hcy
            rnypsdo = rny
            csxpsdo = csx
            CALL dryLeaf(ktau, ktau_tot,dels, rad, air, met,  &
            veg, canopy, soil, ssnow, casapool, dsxpsdo, dsypsdo, psilxpsdo, psilypsdo,&
            fwsoilpsdo, fwsoiltmppsdo, fwpsipsdo, tlfxpsdo, tlfypsdo, ecypsdo, hcypsdo,  &
            rnypsdo, gbhu, gbhf, csxpsdo, cansat,  &
            ghwet, iter, climate, wbpsdo=wbpsdo)
            !!!!!!!!!!!!!!  wb = ssat, & vpd = 0.6Kpa !!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!! gs_ref
            wbpsdo = SPREAD(real(soil%ssat,r_2), 2, ms) 
            DO j = 1, mp
               CALL calc_soil_root_resistance(ssnow, soil, veg, casapool, casabiome, root_length, j, wbpsdo)
               CALL calc_swp(ssnow, soil, j, wbpsdo)
            END DO
            vpdpsdo = SPREAD(real(600,r_2), 1, mp) ! 600Pa
            dsxpsdo = vpdpsdo
            dsypsdo = vpdpsdo
            psilxpsdo = psilx
            psilypsdo = psily
            fwsoilpsdo = fwsoil
            fwsoiltmppsdo = fwsoiltmp
            fwpsipsdo = fwpsi
            tlfxpsdo = tlfx
            tlfypsdo = tlfy
            ecypsdo = ecy
            hcypsdo = hcy
            rnypsdo = rny
            csxpsdo = csx
            CALL dryLeaf(ktau, ktau_tot,dels, rad, air, met,  &
            veg, canopy, soil, ssnow, casapool, dsxpsdo, dsypsdo, psilxpsdo, psilypsdo,&
            fwsoilpsdo, fwsoiltmppsdo, fwpsipsdo, tlfxpsdo, tlfypsdo, ecypsdo, hcypsdo,  &
            rnypsdo, gbhu, gbhf, csxpsdo, cansat,  &
            ghwet, iter, climate, wbpsdo=wbpsdo,vpdpsdo=vpdpsdo)
            !!!!!!!!!!!!!! gs_coeff = 1 !!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!! gs_ref1
            wbpsdo = SPREAD(real(soil%ssat,r_2), 2, ms) 
            DO j = 1, mp
               CALL calc_soil_root_resistance(ssnow, soil, veg, casapool, casabiome, root_length, j, wbpsdo)
               CALL calc_swp(ssnow, soil, j, wbpsdo)
            END DO
            fwpsdo = 1.0_r_2
            dsxpsdo = dsx
            dsypsdo = dsy
            psilxpsdo = psilx
            psilypsdo = psily
            fwsoilpsdo = fwsoil
            fwsoiltmppsdo = fwsoiltmp
            fwpsipsdo = fwpsi
            tlfxpsdo = tlfx
            tlfypsdo = tlfy
            ecypsdo = ecy
            hcypsdo = hcy
            rnypsdo = rny
            csxpsdo = csx
            CALL dryLeaf(ktau, ktau_tot,dels, rad, air, met,  &
            veg, canopy, soil, ssnow, casapool, dsxpsdo, dsypsdo, psilxpsdo, psilypsdo,&
            fwsoilpsdo, fwsoiltmppsdo, fwpsipsdo, tlfxpsdo, tlfypsdo, ecypsdo, hcypsdo,  &
            rnypsdo, gbhu, gbhf, csxpsdo, cansat,  &
            ghwet, iter, climate, wbpsdo=wbpsdo,fwpsdo=fwpsdo)
         endif
         if (iter==4) then
            DO j = 1, mp
               ! reset psi_soil, soilR and rootR back to the realistic value
               CALL calc_soil_root_resistance(ssnow, soil, veg, casapool, casabiome, root_length, j)
               CALL calc_swp(ssnow, soil, j)
            END DO
         endif
         CALL dryLeaf(ktau, ktau_tot,dels, rad, air, met,  &
            veg, canopy, soil, ssnow, casapool, dsx, dsy, psilx, psily,&
            fwsoil, fwsoiltmp, fwpsi, tlfx, tlfy, ecy, hcy,  &
            rny, gbhu, gbhf, csx, cansat,  &
            ghwet, iter, climate)

         CALL wetLeaf( dels, rad, air, met, &
            canopy, cansat, tlfy,    &
            gbhu, gbhf, ghwet )


         ! Calculate latent heat from vegetation:
         ! Calculate sensible heat from vegetation:
         ! Calculate net rad absorbed by canopy:
         canopy%fev = real(canopy%fevc) + canopy%fevw
         canopy%fhv = (1.0-canopy%fwet) * real(hcy) + canopy%fhvw
         canopy%fnv = (1.0-canopy%fwet) * real(rny) + canopy%fevw + canopy%fhvw

         ! canopy rad. temperature calc from long-wave rad. balance
         sum_rad_gradis = SUM(rad%gradis,2)
         ! print*, 'RR27 ', canopy%fevc
         ! print*, 'RR28 ', canopy%fevw
         ! print*, 'RR29 ', canopy%fwet
         ! print*, 'RR30 ', canopy%fhvw
         ! print*, 'RR31 ', canopy%fevw
         ! print*, 'RR32 ', rad%gradis

         DO j=1,mp

            IF ( canopy%vlaiw(j) > C%LAI_THRESH .AND. &
               rough%hruff(j) > rough%z0soilsn(j) ) THEN

               rad%lwabv(j) = C%CAPP * C%rmair * ( tlfy(j) - met%tk(j) ) * &
                  sum_rad_gradis(j)

               !! vh_js !!
               if (  (rad%lwabv(j) / (2.0*(1.0-rad%transd(j)) &
                  * C%SBOLTZ*C%EMLEAF)+met%tk(j)**4) > 0.0) then
                  canopy%tv(j) = (rad%lwabv(j) / (2.0*(1.0-rad%transd(j)) &
                     * C%SBOLTZ*C%EMLEAF)+met%tk(j)**4)**0.25
               else
                  canopy%tv(j) = met%tk(j)
               endif

            ELSE ! sparse canopy

               canopy%tv(j) = met%tk(j)

            ENDIF

         ENDDO
         ! print*, 'RR33 ', rad%transd
         ! print*, 'RR34 ', canopy%tv

         ! Calculate net rad to soil:
         canopy%fns = rad%qssabs + rad%transd*met%fld + (1.0-rad%transd)*C%EMLEAF* &
            C%SBOLTZ*canopy%tv**4 - C%EMSOIL*C%SBOLTZ* tss4
         ! print*, 'RR35 ', rad%qssabs
         ! print*, 'RR36 ', rad%transd
         ! print*, 'RR37 ', canopy%tv
         ! print*, 'RR38 ', tss4
         ! print*, 'RR39 ', canopy%fns

         ! Saturation specific humidity at soil/snow surface temperature:
         call qsatfjh(ssnow%qstss, ssnow%tss-C%tfrz, met%pmb)
         ! print*, 'RR40 ', ssnow%tss
         ! print*, 'RR41 ', ssnow%qstss

         If (cable_user%soil_struc=='default') THEN

            IF (cable_user%ssnow_POTEV== "P-M") THEN

               !--- uses %ga from previous timestep
               ssnow%potev = Penman_Monteith(canopy%ga)

            ELSE !by default assumes Humidity Deficit Method

               ! Humidity deficit
               dq = ssnow%qstss - met%qv
               ssnow%potev =  Humidity_deficit_method(dq)

            ENDIF

            ! Soil latent heat:
            CALL latent_heat_flux()

            ! Calculate soil sensible heat:
            !canopy%fhs = air%rho*C%CAPP*(ssnow%tss - met%tk) /ssnow%rtsoil
            IF (cable_user%litter) THEN
               !! vh_js !! account for additional litter resistance to sensible heat transfer
               canopy%fhs = air%rho * C%CAPP * (ssnow%tss - met%tk) / &
                  ( ssnow%rtsoil + real(1-ssnow%isflag) * &
                  real(veg%clitt*0.003_r_2/canopy%kthLitt) / (air%rho*C%CAPP) )
            ELSE
               canopy%fhs = air%rho * C%CAPP * (ssnow%tss - met%tk) / ssnow%rtsoil
            ENDIF

         ELSEIF (cable_user%soil_struc=='sli') THEN
            ! SLI SEB to get canopy%fhs, canopy%fess, canopy%ga
            ! (Based on old Tsoil, new canopy%tv, new canopy%fns)
            CALL sli_main(1,dels,veg,soil,ssnow,met,canopy,air,rad,1)
         ENDIF

         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         ! ! SPECIAL for global offine met: set screen t and humidity as met inputs
         ! tstar = - (canopy%fhv+canopy%fhs) / ( air%rho*C%capp*canopy%us)
         ! qstar = - (canopy%fev+canopy%fes) / ( air%rho*air%rlam *canopy%us)
         ! zscrn = max(rough%z0m,1.8-rough%disp)
         ! denom = ( log(rough%zref_tq/zscrn)- psis(canopy%zetar(:,iter)) + &
         !      psis(canopy%zetar(:,iter) * zscrn / rough%zref_tq) ) /C%vonk
         !
         ! where (canopy%zetar(:,iter) > 0.7)
         !    zeta2=canopy%zetar(:,iter) * zscrn / rough%zref_tq
         !    denom =alpha1* ((canopy%zetar(:,iter)**beta1* &
         !         (1.0+gamma1*canopy%zetar(:,iter)**(1.0-beta1))) &
         !         - (zeta2**beta1*(1.0+gamma1*zeta2**(1.0-beta1)))) /C%vonk
         ! endwhere
         !
         ! where (abs( tstar*denom ).lt. 5.0)
         !    met%tk = canopy%tscrn + tstar*denom
         !    met%qv = max(canopy%qscrn + qstar*denom, 0.0)
         ! endwhere
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         CALL within_canopy( gbhu, gbhf )

         ! Saturation specific humidity at soil/snow surface temperature:
         call qsatfjh(ssnow%qstss, ssnow%tss-C%tfrz, met%pmb)

         IF (cable_user%soil_struc=='default') THEN

            IF(cable_user%ssnow_POTEV== "P-M") THEN

               !--- uses %ga from previous timestep
               ssnow%potev = Penman_Monteith(canopy%ga)

            ELSE !by default assumes Humidity Deficit Method

               ! Humidity deficit
               dq = ssnow%qstss - met%qvair
               ssnow%potev =  Humidity_deficit_method(dq)

            ENDIF

            ! Soil latent heat:
            CALL latent_heat_flux()

            ! Soil sensible heat:
            !canopy%fhs = air%rho*C%CAPP*(ssnow%tss - met%tvair) /ssnow%rtsoil
            IF (cable_user%litter) THEN
               !! vh_js !! account for additional litter resistance to sensible heat transfer
               canopy%fhs =  air%rho *C%CAPP *(ssnow%tss - met%tvair) / &
                  ( ssnow%rtsoil +  real(1-ssnow%isflag) * &
                  real(veg%clitt*0.003_r_2/canopy%kthLitt) / (air%rho*C%CAPP) )
            else
               canopy%fhs = air%rho*C%CAPP*(ssnow%tss - met%tvair) /ssnow%rtsoil
            ENDIF

            !! Ticket #90 ssnow%cls factor should be retained: required for energy balance
            canopy%ga = canopy%fns - canopy%fhs - real(canopy%fes) !*ssnow%cls

         ELSEIF (cable_user%soil_struc=='sli') THEN

            ! SLI SEB to get canopy%fhs, canopy%fess, canopy%ga
            ! (Based on old Tsoil, new canopy%tv, new canopy%fns)
            CALL sli_main(1,dels,veg,soil,ssnow,met,canopy,air,rad,1)

         ENDIF

         ! Set total latent heat:
         canopy%fe = canopy%fev + real(canopy%fes)

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
         where (eq(rlower_limit, 0.0)) rlower_limit = 1.e-7 !prevent from 0. by adding 1.e-7 (W/m2)

         canopy%wetfac_cs = max(0., min(1.0, canopy%fe / rlower_limit))

         DO j=1,mp
            IF ( canopy%wetfac_cs(j) .LE. 0. ) &
               canopy%wetfac_cs(j) = MAX( 0., MIN( 1., &
               MAX( canopy%fev(j)/canopy%fevw_pot(j), &
               real(canopy%fes(j))/ssnow%potev(j) ) ) )
         ENDDO

         CALL update_zetar()

         ! print*, 'SS04.1 ', iter, ssnow%rtsoil
         ! print*, 'SS04.2 ', ssnow%tss
         ! print*, 'SS04.3 ', ssnow%qstss
         ! print*, 'SS04.4 ', ssnow%potev
         ! print*, 'SS04.5 ', ssnow%cls
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

         IF (qtgnet(j) > 0. ) THEN
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
      canopy%dewmm = -(min(0.0,canopy%fevw) + min(0.0,real(canopy%fevc))) * &
         dels * 1.0e3 / (C%RHOW*air%rlam)

      ! Add dewfall to canopy water storage:
      canopy%cansto = canopy%cansto + canopy%dewmm

      ! Modify canopy water storage for evaporation:
      canopy%cansto = MAX(canopy%cansto-MAX(0.0,canopy%fevw)*dels &
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
            ssnow%dfh_dtg = air%rho * C%CAPP / ( ssnow%rtsoil + &
               real(1-ssnow%isflag) * real(veg%clitt*0.003_r_2/canopy%kthLitt) / &
               (air%rho*C%CAPP) )
            ssnow%dfe_ddq = ssnow%wetfac * air%rho * air%rlam * ssnow%cls/ &
               ( ssnow%rtsoil + real(1-ssnow%isflag) * &
               real(veg%clitt*0.003_r_2/canopy%DvLitt) )
         ELSE
            ssnow%dfh_dtg = air%rho * C%CAPP / ssnow%rtsoil
            ssnow%dfe_ddq = ssnow%wetfac * air%rho * air%rlam * ssnow%cls / ssnow%rtsoil
         ENDIF

         ssnow%ddq_dtg = (C%rmh2o/C%rmair) /met%pmb * C%TETENA*C%TETENB * C%TETENC &
            / ( ( C%TETENC + ssnow%tss-C%tfrz )**2 )*EXP( C%TETENB * &
            ( ssnow%tss-C%tfrz ) / ( C%TETENC + ssnow%tss-C%tfrz ) )
         canopy%dgdtg = real(ssnow%dfn_dtg - ssnow%dfh_dtg - ssnow%dfe_ddq * &
            ssnow%ddq_dtg, r_2)
      ENDIF

      bal%drybal = REAL(ecy+hcy) - SUM(rad%rniso,2) &
         + C%CAPP*C%rmair*(tlfy-met%tk)*SUM(rad%gradis,2)  ! YP nov2009

      bal%wetbal = canopy%fevw + canopy%fhvw - SUM(rad%rniso,2) * canopy%fwet &
         + C%CAPP*C%rmair * (tlfy-met%tk) * SUM(rad%gradis,2) * &
         canopy%fwet  ! YP nov2009
      CALL Penman_Monteith_canopy(gbhu,gbhf)
      DEALLOCATE(cansat,gbhu)
      DEALLOCATE(dsx,dsxpsdo,dsy, dsypsdo, fwsoil, fwsoilpsdo, &
       fwsoiltmp,fwsoiltmppsdo, fwpsi,fwpsipsdo, tlfx, tlfxpsdo, tlfy, tlfypsdo)
      DEALLOCATE(ecy,ecypsdo, hcy,hcypsdo, rny, rnypsdo)
      DEALLOCATE(gbhf, csx, csxpsdo)
      DEALLOCATE(ghwet)
      DEALLOCATE(wbpsdo, vpdpsdo)
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
            sss,              & ! var for Penman-Monteith soil evap
            cc1,              & ! var for Penman-Monteith soil evap
            cc2,              & ! var for Penman-Monteith soil evap
            qsatfvar            !

         ! Penman-Monteith formula
         sss = air%dsatdk
         cc1 = sss / (sss + air%psyc)
         cc2 = air%psyc / (sss + air%psyc)

         call qsatfjh(qsatfvar, met%tvair-c%tfrz, met%pmb)

         if (cable_user%litter) then
            !! vh_js !!
            ssnowpotev = cc1 * (canopy%fns - ground_H_flux) + &
               cc2 * air%rho * air%rlam * (qsatfvar - met%qvair) / &
               (ssnow%rtsoil + real(1-ssnow%isflag) * &
               real(veg%clitt*0.003_r_2/canopy%DvLitt))
         else
            ssnowpotev = cc1 * (canopy%fns - ground_H_flux) + &
               cc2 * air%rho * air%rlam * (qsatfvar - met%qvair) / ssnow%rtsoil
         endif
         ! print*, 'PP01 ', air%dsatdk
         ! print*, 'PP02 ', air%psyc
         ! print*, 'PP03 ', canopy%fns
         ! print*, 'PP04 ', air%rho
         ! print*, 'PP05 ', air%rlam
         ! print*, 'PP06 ', ssnow%rtsoil
         ! print*, 'PP07 ', ground_h_flux

      END FUNCTION Penman_Monteith

      SUBROUTINE Penman_Monteith_canopy(gbhu,gbhf) 
         !USE cable_def_types_mod, only : mp
         ! FAO 56
         implicit none
         REAL(r_2), INTENT(IN), DIMENSION(:,:) :: &
            gbhu, & ! forcedConvectionBndryLayerCond
            gbhf    ! freeConvectionBndryLayerCond
         !real, dimension(mp), intent(in)  :: ground_H_flux
         real(r_2), dimension(mp)              :: canopypot1,canopypot2 ! returned result of functio
         real(r_2), dimension(mp,mf)              :: gh,ghr,gs,gw,psycst,psycstr ! returned result of function
         
         real(r_2), dimension(mp)   :: RhoA, gsw_ratio
         INTEGER :: i
         RhoA = met%pmb * 100.0 / (8.314 * (met%tk)) ! air density [molA/m3]
         ! gsw_ratio = canopy%gswx(:,2) / canopy%gswx(:,1)
         gs(:,1) =  (0.5_r_2 * veg%vlai) / 100.0_r_2 * RhoA! ! m s-1 to mol m-2 s-1
         ! gs(:,2) = gs(:,1) * gsw_ratio
         gs(:,2) = gs(:,1)
         ! print*, "gs in penman_Monteith: ", gs
         DO i =1, mp
            
            gh(i,:) = real(2.0_r_2 * (gbhu(i,:) + gbhf(i,:)))
            ghr(i,:) = gh(i,:) + rad%gradis(i,:) ! mol m-2 s-1 
            
            gw(i,:) = 1.0 / ( 1.0 / gs(i,:) + &
            1.0 / ( 1.075 * real(gbhu(i,:) + gbhf(i,:)) ) )
            
            psycst(i,:) = air%psyc(i) * REAL( gh(i,:)  / gw(i,:) ) ! for gh
            psycstr(i,:) = air%psyc(i) * REAL( ghr(i,:)  / gw(i,:) ) ! for ghr

            canopypot1(i) = real( ( air%dsatdk(i) * ( rad%rniso(i,1) - C%capp * C%rmair &
            * ( met%tvair(i) - met%tk(i) ) * rad%gradis(i,1) ) &
            + C%capp * C%rmair * met%dva(i) * ghr(i,1) ) &
            / ( air%dsatdk(i) + psycstr(i,1) ),r_2) + real( ( air%dsatdk(i) * ( rad%rniso(i,2) - C%capp * C%rmair &
            * ( met%tvair(i) - met%tk(i) ) * rad%gradis(i,2) ) &
            + C%capp * C%rmair * met%dva(i) * ghr(i,2) ) &
            / ( air%dsatdk(i) + psycstr(i,2) ),r_2) 

            canopypot2(i) = real( ( air%dsatdk(i) * ( rad%rniso(i,1) )  &
            + C%capp * C%rmair * met%dva(i) * gh(i,1) ) &
            / ( air%dsatdk(i) + psycst(i,1) ),r_2) + real( ( air%dsatdk(i) * ( rad%rniso(i,2) - C%capp * C%rmair &
            * ( met%tvair(i) - met%tk(i) ) * rad%gradis(i,2) ) &
            + C%capp * C%rmair * met%dva(i) * gh(i,2) ) &
            / ( air%dsatdk(i) + psycst(i,2) ),r_2) 
            canopy%epotcan1(i) = canopypot1(i) *  (1.0_r_2-real(canopy%fwet(i), r_2))
            canopy%epotcan2(i) = canopypot2(i) * (1.0_r_2-real(canopy%fwet(i), r_2))
         END DO
         ! Penman-Monteith formula


      END SUBROUTINE Penman_Monteith_canopy
      ! ------------------------------------------------------------------------------

      ! method alternative to P-M formula above
      FUNCTION humidity_deficit_method(dq) RESULT(ssnowpotev)

         use cable_def_types_mod, only: mp
         use mo_utils,            only: eq

         implicit none

         real, dimension(mp), intent(in) :: dq
         real, dimension(mp)             :: ssnowpotev

         real, dimension(mp) :: idq
         integer :: j

         idq = dq
         do j=1, mp
            if ((ssnow%snowd(j)>1.0) .or. eq(ssnow%tgg(j,1), c%tfrz)) &
               idq(j) = max(-0.1e-3, idq(j))
         enddo

         if (cable_user%litter) then
            !! vh_js !!
            ssnowpotev = air%rho * air%rlam * idq / &
               (ssnow%rtsoil + real(1-ssnow%isflag) * real(veg%clitt*0.003_r_2/canopy%DvLitt))
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

            if ((ssnow%snowd(j) < 0.1) .and. (canopy%fess(j) > 0.0_r_2)) then
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

         USE cable_def_types_mod, only: mp, r_2
         use mo_utils,            only: eq

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
                  IF ( eq((met%fsd(j,1)+met%fsd(j,2)), 0.0) ) &
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

    ! not used
    ! FUNCTION qsatf(j,tair,pmb) RESULT(r)
    !   ! MRR, 1987
    !   ! AT TEMP tair (DEG C) AND PRESSURE pmb (MB), GET SATURATION SPECIFIC
    !   ! HUMIDITY (KG/KG) FROM TETEN FORMULA

    !   REAL, INTENT(IN) :: &
    !        tair,         & ! air temperature (C)
    !        pmb             ! pressure PMB (mb)
    !   INTEGER, INTENT(IN) :: j

    !   REAL           :: r    ! result; sat sp humidity

    !   r = (C%RMH2o/C%rmair) * (C%TETENA*EXP(C%TETENB*tair/(C%TETENC+tair))) / pmb

    ! END FUNCTION qsatf

    ! -----------------------------------------------------------------------------

    SUBROUTINE qsatfjh(var, tair, pmb)

      implicit none

      REAL, INTENT(OUT), DIMENSION(:) :: var ! result; sat sp humidity
      REAL, INTENT(IN),  DIMENSION(:) :: &
           tair, & ! air temperature (C)
           pmb     ! pressure (mbar)

      INTEGER :: j, mp

      mp = size(var, 1)
      
      DO j=1, mp
         var(j) = (C%RMH2O / C%rmair) * ( C%TETENA * &
              EXP( C%TETENB * tair(j) / (C%TETENC + tair(j)) ) ) &
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

      use mo_constants, only: pi => pi_sp

      implicit none

      real, dimension(:), intent(in) :: zeta !
      real, dimension(size(zeta, 1)) :: r    ! function result

      real, dimension(size(zeta, 1)) :: &
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

      z = 0.5 + sign(0.5, zeta) ! z=1 in stable, 0 in unstable

      ! Beljaars and Holtslag (1991) for stable
      stable = -a * zeta - b * (zeta - xc / d) * exp(-d * zeta) - b * xc / d
      x = (1.0 + gu * abs(zeta))**0.25
      unstable = log((1.0 + x * x) * (1.0 + x)**2 / 8.) - 2.0 * atan(x) + pi * 0.5
      r = z * stable + (1.0 - z) * unstable

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

      ! not used
      ! ELEMENTAL FUNCTION rplant(rpconst, rpcoef, tair) result(z)

      !   implicit none

      !   real, intent(in) :: rpconst
      !   real, intent(in) :: rpcoef
      !   real, intent(in) :: tair
      !   real             :: z

      !   z = rpconst * exp(rpcoef * tair)

      ! END FUNCTION rplant

      ! -----------------------------------------------------------------------------

      SUBROUTINE wetLeaf( dels, rad, air, met, canopy, cansat, tlfy, &
         gbhu, gbhf, ghwet )

         USE cable_def_types_mod

         implicit none

         REAL,                     INTENT(IN)    :: dels ! integration time step (s)
         TYPE(radiation_type),     INTENT(INOUT) :: rad
         TYPE(air_type),           INTENT(INOUT) :: air
         TYPE(met_type),           INTENT(INOUT) :: met
         TYPE(canopy_type),        INTENT(INOUT) :: canopy
         REAL,INTENT(IN), DIMENSION(:) :: &
            cansat, & ! max canopy intercept. (mm)
            tlfy      ! leaf temp (K) - assC%UMINg the temperature of
         ! wet leaf is equal that of dry leaf ="tlfy"
         REAL(r_2), INTENT(IN), DIMENSION(:,:) :: &
            gbhu,          & ! forcedConvectionBndryLayerCond
            gbhf             ! freeConvectionBndryLayerCond
         REAL(r_2), INTENT(OUT), DIMENSION(:) :: &
            ghwet            ! cond for heat for a wet canopy

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

         ghwet = 1.0e-3_r_2
         gwwet = 1.0e-3
         ghrwet= 1.0e-3
         canopy%fevw = 0.0
         canopy%fhvw = 0.0
         sum_gbh        = real(SUM((gbhu+gbhf),2))
         sum_rad_rniso  = SUM(rad%rniso,2)
         sum_rad_gradis = SUM(rad%gradis,2)

         ! print*, 'DD48 ', rad%rniso
         ! print*, 'DD49 ', rad%gradis
         ! print*, 'DD50 ', canopy%cansto
         ! print*, 'DD51 ', canopy%fwet
         ! print*, 'DD52 ', met%tvair
         ! print*, 'DD53 ', met%dva
         DO j=1,mp

            IF (canopy%vlaiw(j) > C%LAI_THRESH) THEN

               ! VEG SENSIBLE & LATENT HEAT FLUXES fevw, fhvw (W/m2) for a wet canopy
               ! calculate total thermal resistance, rthv in s/m
               ghwet(j)  = 2.0_r_2 * real(sum_gbh(j), r_2)
               gwwet(j)  = 1.075 * sum_gbh(j)
               ghrwet(j) = sum_rad_gradis(j) + real(ghwet(j))

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


   SUBROUTINE surf_wetness_fact(cansat, canopy, ssnow, veg, met, soil, dels)

      USE cable_common_module
      USE cable_def_types_mod

      ! max canopy intercept. (mm)
      REAL, DIMENSION(:),        INTENT(IN)    :: cansat
      TYPE(canopy_type),         INTENT(INOUT) :: canopy
      TYPE(soil_snow_type),      intent(inout) :: ssnow
      TYPE(veg_parameter_type),  INTENT(INOUT) :: veg
      TYPE(met_type),            INTENT(INOUT) :: met
      TYPE(soil_parameter_type), intent(inout) :: soil
      ! integration time setp (s)
      REAL,                      INTENT(IN)    :: dels

      !local variables
      REAL, DIMENSION(mp) :: lower_limit, upper_limit, ftemp
      INTEGER :: j

      ! Rainfall variable is limited so canopy interception is limited,
      ! used to stabilise latent fluxes.
      ! to avoid excessive direct canopy evaporation (EK nov2007, snow scheme)
      upper_limit = 4.0 * MIN(dels,1800.0) / (60.0 * 1440.0 )
      ftemp = MIN(met%precip-met%precip_sn, upper_limit )
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
      ssnow%wetfac = 0.5 * (ssnow%wetfac + ssnow%owetfac)

   END SUBROUTINE surf_wetness_fact


   ! -----------------------------------------------------------------------------


   SUBROUTINE dryLeaf(ktau, ktau_tot, dels, rad, air, met, &
      veg, canopy, soil, ssnow,casapool, dsx, dsy, psilx, psily, &
      fwsoil, fwsoiltmp, fwpsi, tlfx, tlfy, ecy, hcy, &
      rny, gbhu, gbhf, csx, &
      cansat, ghwet, iter, climate, wbpsdo, vpdpsdo,fwpsdo)

      use cable_def_types_mod
      use cable_common_module
      USE casavariable
#ifdef __MPI__
      use mpi, only: MPI_Abort
#endif

      implicit none
      INTEGER,                   INTENT(IN)    :: ktau
      INTEGER,                   INTENT(IN)    :: ktau_tot
      real,                      intent(in)    :: dels ! integration time step (s)
      type(radiation_type),      intent(inout) :: rad
      type(air_type),            intent(inout) :: air
      type(met_type),            intent(inout) :: met
      type(veg_parameter_type),  intent(inout) :: veg
      type(canopy_type),         intent(inout) :: canopy
      type(soil_parameter_type), intent(inout) :: soil
      type(soil_snow_type),      intent(inout) :: ssnow
      TYPE (casa_pool),  INTENT(IN)           :: casapool
      real,      dimension(:),   intent(inout) :: &
         dsx,        & ! leaf surface vpd
         dsy,        & ! leaf surface vpd
         fwsoil,     & ! soil water modifier of stom. cond
         fwsoiltmp,  &
         tlfx,       & ! leaf temp prev. iter (K)
         tlfy  ! leaf temp (K)
      real(r_2), dimension(:),   intent(inout) :: &
         ecy,        & ! lat heat fl dry big leaf
         hcy,        & ! veg. sens heat
         rny         !& !
    
      real(r_2), dimension(:,:), intent(inout) :: &
         gbhu,       & ! forcedConvectionBndryLayerCond
         gbhf,       & ! freeConvectionBndryLayerCond
         csx,        &  ! leaf surface CO2 concentration
         psilx,      &
         psily,      &
         fwpsi   
      real,      dimension(:),   intent(in)    :: cansat
      real(r_2), dimension(:),   intent(out)   :: ghwet  ! cond for heat for a wet canopy
      integer,                   intent(in)    :: iter
      type(climate_type),        intent(in)    :: climate
      real(r_2), dimension(:), intent(in), optional :: vpdpsdo
      real(r_2), dimension(:,:), intent(in), optional :: wbpsdo
      real(r_2), intent(in), optional :: fwpsdo
      !local variables (JK: move parameters to different routine at some point)
      real, parameter :: jtomol = 4.6e-6  ! Convert from J to Mol for light

      real, dimension(mp) :: &
         conkct,        & ! Michaelis Menton const.
         conkot,        & ! Michaelis Menton const.
         cx1,           & ! "d_{3}" in Wang and Leuning,
         cx2,           & !     1998, appendix E
         cx2y,         &
         tdiff,         & ! leaf air temp diff.
         tlfxx,         & ! leaf temp of current iteration (K)
         tlfxm,         & ! 
         dsxx,         & ! vpd at leaf surface of current iteration 
         dsxm,         &
         abs_deltlf,    & ! ABS(deltlf)
         abs_deltds,    &
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
         temp_c3,       & !
         temp_c4,       & !
         temp_sun_c3,   & !
         temp_shade_c3, & !
         temp_sun_c4,   & !
         temp_shade_c4, &    !
         fwsoil_nongs 

      real(r_2), dimension(mp)  :: &
         ecx,        & ! lat. hflux big leaf
         hcx,        & ! sens heat fl big leaf prev iteration
         rnx, &          ! net rad prev timestep
         ecxs, &  ! lat. hflux big leaf (sap flux)
         psixx, &
         kplantx, &
         total_est_evap
         
      real(r_2) :: psixxi, kplantxi
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
         vcmxt3y,     & ! vcmax big leaf C3
         vcmxt4y,     & ! vcmax big leaf C4
         vx3,        & ! carboxylation C3 plants
         vx4,        & ! carboxylation C4 plants
         co2cp,      & ! CO2 compensation point (needed for Leuing stomatal conductance)
      ! Ticket #56, xleuning is no longer used, we replace it with gs_coeff,
      ! which is computed differently based on the new GS_SWITCH. If GS_SWITCH
      ! is "leuning", it's the same, if "medlyn", then the new Medlyn model
      ! xleuning,   & ! leuning stomatal coeff
         gs_coeff,   & ! stom coeff, Ticket #56
         gs_coeffy,   & 
         psycst,     & ! modified pych. constant
         frac42,     & ! 2D frac4
         temp2,      &
         anrubiscox, & ! net photosynthesis (rubisco limited)
         anrubpx,    & ! net photosynthesis (rubp limited)
         ansinkx,    & ! net photosynthesis (sink limited)
         anrubiscoy, & ! net photosynthesis (rubisco limited)
         anrubpy,    & ! net photosynthesis (rubp limited)
         kc4,        & ! An-Ci slope in Collatz 1992 model
         gmes,       & ! mesophyll conductance of big leaf
         gswy,       &
         csy

      real(r_2), dimension(mp,mf)  ::  dAnrubiscox, & ! CO2 elasticity of net photosynthesis (rubisco limited)
         dAnrubpx, &   !  (rubp limited)
         dAnsinkx, &  ! (sink limited)
         dAnx, & !(actual rate)
         dAnrubiscoy, & ! CO2 elasticity of net photosynthesis (rubisco limited)
         dAnrubpy, &   !  (rubp limited)
         dAnsinky, &  ! (sink limited)
         dAn_y, & !(actual rate)
         eta_y, &
         eta_x, &
         ex,    &  !transpiration, kg m-2 s-1
         fwpsiy, &
         psilxx, &
         psilxm, &
         csxx, &
         abs_deltpsil , &
         abs_deltcs, &
         abs_deltfwpsi

      real ::  gam0,    &
         conkc0,  &
         conko0,  &
         egam,    &
         ekc,     &
         eko,     &
         qs,      &
         qm,      &
         qb

      ! real, dimension(:,:), pointer :: gswmin => null() ! min stomatal conductance
      ! real, dimension(:,:), allocatable :: gswmin ! min stomatal conductance
      real, dimension(mp,mf) :: gswmin ! min stomatal conductance
      real(r_2), dimension(mp,mf) :: fwpsixx, fwpsi1_tmp
      real(r_2) :: fw, inc
      real, dimension(mp,2) ::  gsw_term, lower_limit2  ! local temp var
      real(r_2), dimension(mp,ms) :: wbtmp
      real(r_2), dimension(mp):: vpdtmp
      real, dimension(mp,mf):: g0, g0xx
      integer :: i, j, k, kk, h, iter_ini ! iteration count
      integer :: NN,m,kmax,Mtag,Numtag
      integer, allocatable :: nktau(:), allktau(:), nktau_end(:)
      real :: vpd, g1, ktot, refill  ! Ticket #56
      REAL, PARAMETER :: & ! Ref. params from Bernacchi et al. (2001)
         co2cp325 = 42.75, & ! CO2 compensation pt C3 at 25 degrees, umol mol-1
         Eaco2cp325 = 37830. ! activation energy for the CO2 compensation pt
      REAL :: co2cp3 ! CO2 compensation pt C3, N.B.: this is set to 0 by default in CABLE

      REAL :: trans_mmol, press
      REAL, PARAMETER :: KG_2_G = 1000.0
      REAL, PARAMETER :: G_WATER_TO_MOL = 1.0 / 18.01528
      REAL, PARAMETER :: MOL_2_MMOL = 1000.0
      REAL, PARAMETER :: MMOL_2_MOL = 1.0 / MOL_2_MMOL
      REAL, PARAMETER :: MB_TO_PA = 100.
      INTEGER, PARAMETER :: resolution = 800 ! jumps in Ci ~ 0.5 umol mol-1
      REAL, DIMENSION(2) :: an_canopy
      REAL :: e_canopy
      REAL(r_2), DIMENSION(resolution) :: p

      REAL :: MOL_WATER_2_G_WATER, G_TO_KG, UMOL_TO_MOL, MB_TO_KPA, PA_TO_KPA
      REAL, DIMENSION(mp) :: kcmax, avg_kcan
      REAL :: new_plc_sat, new_plc_stem, new_plc_can
      REAL :: MOL_TO_UMOL, J_TO_MOL, dc
      CHARACTER(LEN=200) :: txtname,num_str
      logical :: NonStoLim = .True.
      logical :: modify_rule0, modify_rule
  
#ifdef __MPI__
      integer :: ierr
#endif

      ! END header
      ! allocate(gswmin(mp,mf))
    
      g0 = spread(veg%g0, 2, mf)
      g0xx = g0
      ! Soil water limitation on stomatal conductance:
      iter_ini = 1
      if (present(wbpsdo)) then
         iter_ini = 4
         wbtmp = wbpsdo
         total_est_evap = ssnow%total_est_evap_sat
      else
         wbtmp = ssnow%wb
         total_est_evap = ssnow%total_est_evap
      endif
      if (present(vpdpsdo)) then
         iter_ini = 4
         vpdtmp = vpdpsdo
      else
         vpdtmp = met%dva
      endif
      if (iter==iter_ini) then 
         if ((cable_user%soil_struc=='default') .and. (INDEX(cable_user%FWSOIL_SWITCH, 'Haverd2013') <= 0)) then
            if (cable_user%fwsoil_switch == 'standard') then
               call fwsoil_calc_std(fwsoil, soil, ssnow, veg, wbtmp)
            elseif (cable_user%fwsoil_switch == 'non-linear extrapolation') then
               !EAK, 09/10 - replace linear approx by polynomial fitting
               call fwsoil_calc_non_linear(fwsoil, soil, ssnow, veg, wbtmp)
            elseif (cable_user%fwsoil_switch == 'Lai and Katul 2000') then
               call fwsoil_calc_Lai_Katul(fwsoil, soil, ssnow, veg, wbtmp)
            elseif (cable_user%FWSOIL_SWITCH == 'profitmax') then
               fwsoil = 1.0
            elseif (cable_user%FWSOIL_SWITCH == 'constant1') then
               fwsoil = 0.98
            elseif (INDEX(cable_user%FWSOIL_SWITCH, 'LWP') > 0) then
               fwsoil = 1.0

            else
               write(*,*) 'fwsoil_switch failed.'
#ifdef __MPI__
               call MPI_Abort(0, 126, ierr) ! Do not know comm nor rank here
#else
               stop 126
#endif
            endif
            canopy%fwsoil = real(fwsoil, r_2)
            canopy%fwsoiltmp = real(fwsoil, r_2)
         elseif ((cable_user%soil_struc=='sli') .or. (INDEX(cable_user%FWSOIL_SWITCH, 'Haverd2013') > 0)) then
            fwsoil = real(canopy%fwsoil)
            canopy%fwsoiltmp = real(fwsoil, r_2)
         endif
      
      endif
      ! that could be changed later  zihanlu 
      if (INDEX(cable_user%FWSOIL_SWITCH, 'LWP') > 0) then
         !print*, 'fwsoil switch:  LWP'
         ! fwsoil = 1.0
         ! canopy%fwsoil = real(fwsoil, r_2)
         ! canopy%fwsoiltmp = real(fwsoil, r_2)
         canopy%fwpsi = real(fwpsi, r_2)
      end if
      if ((INDEX(cable_user%FWSOIL_SWITCH, 'LWP') > 0) .or. (.NOT. NonStoLim)) Then
         fwsoil_nongs = 1.0
      else 
         fwsoil_nongs = fwsoil
      endif

      psixx = real(canopy%psix,r_2)
      kplantx = real(canopy%kplant,r_2)
     ! write(logn,*) 'fwsoil of ', cable_user%fwsoil_switch, ': ', fwsoil(1)
      MOL_TO_UMOL = 1E6
      J_TO_MOL = 4.6E-6  ! Convert from J to Mol for light
      MOL_WATER_2_G_WATER = 18.02
      G_TO_KG = 1E-03
      UMOL_TO_MOL = 1E-06
      MB_TO_KPA = 0.1
      PA_TO_KPA = 1E-03
      ! weight min stomatal conductance by C3 an C4 plant fractions
      frac42       = spread(veg%frac4, 2, mf) ! frac C4 plants
      gsw_term     = spread(veg%gswmin, 2, mf)
      lower_limit2 = rad%scalex * gsw_term
      gswmin       = max(1.e-6, lower_limit2)
      ! print*, 'DD02 ', veg%frac4
      ! print*, 'DD03 ', veg%gswmin
      ! print*, 'DD04 ', rad%scalex

      gw          = 1.0e-3 ! default values of conductance
      gh          = 1.0e-3
      ghr         = 1.0e-3
      rdx         = 0.0
      anx         = 0.0
      anrubiscox  = 0.0
      anrubpx     = 0.0
      ansinkx     = 0.0
      dAnrubiscox = 0.0_r_2
      dAnrubpx    = 0.0_r_2
      dAnsinkx    = 0.0_r_2
      dAnx        = 0.0
      eta_x       = 0.0_r_2
      rnx         = SUM(real(rad%rniso,r_2),2)
      abs_deltlf  = 999.0
      abs_deltds  = 999.0
      abs_deltcs(:,:)   = SPREAD(SPREAD(999.0_r_2,1,mp),2,mf)
      vcmxt3      = 0.0
      vcmxt4      = 0.0
      ! print*, 'DD05 ', rad%rniso

      gras  = 1.0e-6
      an_y  = 0.0
      hcx   = 0.0_r_2              ! init sens heat iteration memory variable
      hcy   = 0.0_r_2
      rdy   = 0.0
      vcmxt3y      = 0.0
      vcmxt4y      = 0.0
      ecx   = SUM(real(rad%rniso,r_2),2) ! init lat heat iteration memory variable
      tlfxx = tlfx
      if (present(fwpsdo)) then
      fwpsixx = fwpsdo
      else
      fwpsixx = fwpsi
      endif
      dsxx = dsx
      csxx = csx
      psilxx = psilx
      psycst(:,:)   = SPREAD(air%psyc,2,mf)
      canopy%fevc   = 0.0_r_2
      ssnow%evapfbl = 0.0

      ecxs = 0.0_r_2
      canopy%fevcs = 0.0_r_2

      ghwet = 1.0e-3_r_2
      gwwet = 1.0e-3
      ghrwet= 1.0e-3
      canopy%fevw = 0.0
      canopy%fhvw = 0.0
      sum_gbh        = real(SUM((gbhu+gbhf),2))
      sum_rad_rniso  = SUM(rad%rniso,2)
      sum_rad_gradis = SUM(rad%gradis,2)
      IF (cable_user%GS_SWITCH == 'tuzet' .AND. &
                  INDEX(cable_user%FWSOIL_SWITCH,'LWP')>0) then
         abs_deltpsil(:,:)   = SPREAD(SPREAD(999.0_r_2,1,mp),2,mf)
         abs_deltfwpsi(:,:)   = SPREAD(SPREAD(999.0_r_2,1,mp),2,mf)
      else
         abs_deltpsil(:,:)   = SPREAD(SPREAD(0.0_r_2,1,mp),2,mf)
         abs_deltfwpsi(:,:)   = SPREAD(SPREAD(0.0_r_2,1,mp),2,mf)
      endif
      ! print*, 'DD06 ', rad%gradis

      ! default for variable d_3 of RuBP-limited photosynthesis of
      ! Wang and Leuning (1998), which is 2 * Gamma^star
      if (cable_user%explicit_gm) then
         if (trim(cable_user%Rubisco_parameters) == 'Bernacchi_2002') then
            gam0 = C%gam0cc
            egam = C%egamcc
         else if (trim(cable_user%Rubisco_parameters) == 'Walker_2013') then
            gam0 = C%gam0ccw
            egam = C%egamccw
         endif
      else
         gam0 = C%gam0
         egam = C%egam
      endif
      cx2(:) = 2.0 * gam0 * exp( egam / (C%rgas * C%trefk) &
         * (1.0 - C%trefk / tlfx(:)) )
      !a = k25 * exp((Ea * (Tk - 298.15)) / (298.15 * RGAS * Tk))

      DO kk=1,mp

         IF(canopy%vlaiw(kk) <= C%LAI_THRESH) THEN
            rnx(kk) = 0.0_r_2 ! intialise
            ecx(kk) = 0.0_r_2 ! intialise
            ecy(kk) = ecx(kk) ! store initial values
            ecxs(kk) = 0.0_r_2 ! initialise
            abs_deltlf(kk) = 0.0
            abs_deltds(kk) = 0.0
            abs_deltpsil(kk,:) = 0.0
            abs_deltfwpsi(kk,:) = 0.0
            abs_deltcs(kk,:) = 0.0
            rny(kk) = rnx(kk) ! store initial values
            ! calculate total thermal resistance, rthv in s/m
         END IF

      ENDDO

      deltlfy = abs_deltlf
      

      !kdcorbin, 08/10 - doing all points all the time'
      !allocate(nktau(4), nktau_end(4))
      nktau=[192212,238148,385692]
      NN=14
      nktau_end = nktau + NN - 1
      
      m = size(nktau) * NN
      allocate(allktau(m))
      do i = 1, size(nktau)
         do k = 0, NN-1
             j = (i - 1) * NN + k
             allktau(j+1) = nktau(i) + k
         end do
     end do
      !write(num_str, '(I0)') nktau
      ! txtname = trim(filename%path) // '/testIteration_cable_out_' // trim(num_str) &
      ! // '.txt'
      if (present(wbpsdo) .and. present(vpdpsdo)) then
         txtname = trim(filename%path) // 'testIteration_refww_cable_out.txt'
      elseif (present(wbpsdo) .and. present(fwpsdo)) then
      elseif (present (wbpsdo)) then
         txtname = trim(filename%path) // 'testIteration_wbpsdo_cable_out.txt'
      elseif (present (vpdpsdo)) then
         txtname = trim(filename%path) // 'testIteration_vpdpsdo_cable_out.txt'
      elseif (.not. present(wbpsdo) .and. .not. present(vpdpsdo) .and. .not. present(fwpsdo)) then
         txtname = trim(filename%path) // 'testIteration_cable_out.txt'
      endif
      ! if (any(nktau == ktau_tot) .and. iter==1) then
      !    if (ktau_tot == nktau(1)) then
      !        ! Open the file for overwrite if k is the first element
      !        open(unit=134, file=txtname, status="unknown", action="write")
      !    else
      !        ! Open the file for appending otherwise
      !        open(unit=134, file=txtname, status="unknown", position="append", action="write")
      !    end if
      ! end if
      if (present(wbpsdo) .and. present(vpdpsdo)) then
         if (ktau_tot == nktau(1) .and. iter==4) then
            ! Open the file for overwrite if k is the first element
            open(unit=137, file=txtname, status="unknown", action="write")
         end if
      elseif (present(wbpsdo) .and. present(fwpsdo)) then
      elseif (present (wbpsdo)) then
         if (ktau_tot == nktau(1) .and. iter==4) then
            ! Open the file for overwrite if k is the first element
            open(unit=135, file=txtname, status="unknown", action="write")
         end if
      elseif (present (vpdpsdo)) then
         if (ktau_tot == nktau(1) .and. iter==4) then
            ! Open the file for overwrite if k is the first element
            open(unit=136, file=txtname, status="unknown", action="write")
         end if
      elseif (.not. present(wbpsdo) .and. .not. present(vpdpsdo) .and. .not. present(fwpsdo)) then
         if (ktau_tot == nktau(1) .and. iter==1) then
            ! Open the file for overwrite if k is the first element
            open(unit=134, file=txtname, status="unknown", action="write")
         end if
      endif
      ! if (ktau_tot==nktau .and. iter==1) then
      ! open(unit=134, file=txtname)
      ! !print*, 'write iteration file '
      ! end if
      
      if (present(wbpsdo) .or. present(vpdpsdo)) then
         kmax = veg%Nmax(1)
      else
       !kmax = C%MAXITER
       kmax = veg%Nmax(1)
      endif
      

         ! print*, 'DD07 ', k, C%MAXITER, canopy%cansto
         ! print*, 'DD08 ', canopy%vlaiw
         ! print*, 'DD09 ', air%rlam
         ! print*, 'DD10 ', met%tvair
         ! print*, 'DD11 ', veg%dleaf
         ! print*, 'DD12 ', rad%fvlai
         ! print*, 'DD13 ', air%cmolar
         ! print*, 'DD14 ', veg%vcmax_sun
         ! print*, 'DD15 ', veg%vcmax_shade
         ! print*, 'DD16 ', veg%frac4
         ! print*, 'DD17 ', climate%mtemp
         ! print*, 'DD18 ', veg%ejmax_sun
         ! print*, 'DD19 ', veg%ejmax_shade
         ! print*, 'DD20 ', veg%alpha
         ! print*, 'DD21 ', veg%convex
         ! print*, 'DD22 ', rad%qcan
         ! print*, 'DD23 ', veg%cfrd
         ! print*, 'DD24 ', climate%qtemp_max_last_year
         ! print*, 'DD25 ', veg%a1gs
         ! print*, 'DD26 ', veg%d0gs
         ! print*, 'DD27 ', veg%g0
      DO i=1,mp 
         k = 0
         Numtag = 3
         DO WHILE (k < kmax)
            k = k + 1
            ! modify_rule0 = (canopy%vlaiw(i) > C%LAI_THRESH .AND. &
            !    (abs_deltlf(i) > 0.1 .OR. ANY(abs_deltpsil(i,:) > 0.1) .OR. ANY(abs_deltfwpsi(i,:) > 0.02)) &
            !    .OR. (Numtag > 0))
            modify_rule0 = (canopy%vlaiw(i) > C%LAI_THRESH .AND. &
               (abs_deltlf(i) > 0.1 .OR. ANY(abs_deltpsil(i,:) > 0.1) .OR. ANY(abs_deltfwpsi(i,:) > 0.02)))
            IF (modify_rule0) THEN

               ghwet(i)  = 2.0_r_2 * real(sum_gbh(i),r_2)
               gwwet(i)  = 1.075 * sum_gbh(i)
               ghrwet(i) = sum_rad_gradis(i) + real(ghwet(i))

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
               gbhf(i,1) = real( rad%fvlai(i,1) * air%cmolar(i) * 0.5*C%dheat &
                  * ( gras(i)**0.25 ) / veg%dleaf(i), r_2)
               gbhf(i,2) = real( rad%fvlai(i,2) * air%cmolar(i) * 0.5 * C%dheat &
                  * ( gras(i)**0.25 ) / veg%dleaf(i), r_2)
               gbhf(i,:) = max(1.e-6_r_2, gbhf(i,:))

               ! Conductance for heat:
               gh(i,:) = real(2.0_r_2 * (gbhu(i,:) + gbhf(i,:)))

               if (cable_user%amphistomatous) then
                  gbhf(i,:) = 2.0_r_2 * gbhf(i,:)
                  gh(i,:)   = real(gbhu(i,:) + gbhf(i,:))
               endif

               ! Conductance for heat and longwave radiation:
               ghr(i,:) = rad%gradis(i,:)+gh(i,:)

               ! Choose Ci-based (implicit gm) or Cc-based (explicit gm) parameters
               if (cable_user%explicit_gm) then
                  if (trim(cable_user%Rubisco_parameters) == 'Bernacchi_2002') then
                     gam0    = C%gam0cc
                     conkc0  = C%conkc0cc
                     conko0  = C%conko0cc
                     egam    = C%egamcc
                     ekc     = C%ekccc
                     eko     = C%ekocc
                  else if (trim(cable_user%Rubisco_parameters) == 'Walker_2013') then
                     gam0    = C%gam0ccw
                     conkc0  = C%conkc0ccw
                     conko0  = C%conko0ccw
                     egam    = C%egamccw
                     ekc     = C%ekcccw
                     eko     = C%ekoccw
                  endif
                  qs      = C%qs
                  qm      = C%qm
                  qb      = C%qb
               else
                  gam0   = C%gam0
                  conkc0 = C%conkc0
                  conko0 = C%conko0
                  egam   = C%egam
                  ekc    = C%ekc
                  eko    = C%eko
                  qs     = 0.5
                  qm     = 0.0     ! not used
                  qb     = 1.0
               endif
               ! JK: veg%vcmax_sun and veg%vcmax_shade are Cc-based if cable_user%explicit_gm = TRUE
               !     and Ci-based otherwise. If cable_user%coordinate_photosyn = FALSE,
               !     veg%vcmax_sun = veg%vcmax_shade = veg%vcmaxcc if cable_user%explicit_gm = TRUE
               !     and veg%vcmax_sun = veg%vcmax_shade = veg%vcmax otherwise.
               !     See also Subroutine 'casa_feedback' in casa_cable.F90. Same applies to ejmax.

               ! Leuning 2002 (PCE) equation for temperature response
               ! used for Vcmax for C3 plants:
               if (.not. cable_user%acclimate_photosyn) then
                  temp_sun_c3(i)   = xvcmxt3(tlfx(i)) * veg%vcmax_sun(i) * (1.0-veg%frac4(i))
                  temp_shade_c3(i) = xvcmxt3(tlfx(i)) * veg%vcmax_shade(i) * (1.0-veg%frac4(i))
                  temp_sun_c4(i)   = xvcmxt4(tlfx(i)) * veg%vcmax_sun(i) * veg%frac4(i)
                  temp_shade_c4(i) = xvcmxt4(tlfx(i)) * veg%vcmax_shade(i) * veg%frac4(i)
               else
                  call xvcmxt3_acclim(tlfx(i), climate%mtemp(i) , temp_c3(i))
                  call xvcmxt4_acclim(tlfx(i), climate%mtemp(i) , temp_c4(i))
                  temp_sun_c3(i)   = temp_c3(i) * veg%vcmax_sun(i) * (1.0-veg%frac4(i))
                  temp_shade_c3(i) = temp_c3(i) * veg%vcmax_shade(i) * (1.0-veg%frac4(i))
                  temp_sun_c4(i)   = temp_c4(i) * veg%vcmax_sun(i) * veg%frac4(i)
                  temp_shade_c4(i) = temp_c4(i) * veg%vcmax_shade(i) * veg%frac4(i)
               endif
               vcmxt3(i,1) = rad%scalex(i,1) * temp_sun_c3(i) * fwsoil_nongs(i)**qb
               vcmxt3(i,2) = rad%scalex(i,2) * temp_shade_c3(i) * fwsoil_nongs(i)**qb
               vcmxt4(i,1) = rad%scalex(i,1) * temp_sun_c4(i) * fwsoil_nongs(i)**qb
               vcmxt4(i,2) = rad%scalex(i,2) * temp_shade_c4(i) * fwsoil_nongs(i)**qb
               ! temperature response of the CO2 compensation point/gamma star used
               ! for C3 plants (Bernacchi et al. 2001), umol mol-1
               ! N.B.: this is always zero by default in CABLE, which is wrong.
               co2cp3 = arrh(co2cp325, Eaco2cp325, tlfx(i)) ! CABLE's gamma_star

               ! apply same scaling for k as for Vcmax in C4 plants
               if (.not. cable_user%explicit_gm) then
                  kc4(i,1) = veg%c4kci(i) * vcmxt4(i,1)/veg%vcmax_sun(i)
                  kc4(i,2) = veg%c4kci(i) * vcmxt4(i,2)/veg%vcmax_shade(i)
               else
                  kc4(i,1) = veg%c4kcc(i) * vcmxt4(i,1)/veg%vcmax_sun(i)
                  kc4(i,2) = veg%c4kcc(i) * vcmxt4(i,2)/veg%vcmax_shade(i)
               endif

               ! Leuning 2002 (PCE) equation for temperature response
               ! used for Jmax for C3 plants:
               if (.not.cable_user%acclimate_photosyn) then
                  temp_sun_c3(i) = xejmxt3(tlfx(i)) * &
                     veg%ejmax_sun(i) * (1.0-veg%frac4(i))
                  temp_shade_c3(i) = xejmxt3(tlfx(i)) * &
                     veg%ejmax_shade(i) * (1.0-veg%frac4(i))
               else
                  call xejmxt3_acclim(tlfx(i), climate%mtemp(i), climate%mtemp_max20(i), temp_c3(i))
                  temp_sun_c3(i)   = temp_c3(i) * veg%ejmax_sun(i) * (1.0-veg%frac4(i))
                  temp_shade_c3(i) = temp_c3(i) * veg%ejmax_shade(i) * (1.0-veg%frac4(i))
               endif
               ejmxt3(i,1) = rad%scalex(i,1) * temp_sun_c3(i) * fwsoil_nongs(i)**qb
               ejmxt3(i,2) = rad%scalex(i,2) * temp_shade_c3(i) * fwsoil_nongs(i)**qb

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
               ! Store dsx zihanlu
               dsxx(i) = dsx(i)
               ! "d_{3}" in Wang and Leuning, 1998, appendix E:
               cx1(i) = conkct(i) * (1.0+0.21/conkot(i))
               !cx2(i) = 2.0 * C%gam0 * ( 1.0 + C%gam1 * tdiff(i) &
               !                             + C%gam2 * tdiff(i) * tdiff(i) )
               cx2(i) = 2.0 * gam0 * EXP( ( egam / (C%rgas*C%trefk) ) &
                  * ( 1.0 - C%trefk/tlfx(i) ) )

               ! All equations below in appendix E in Wang and Leuning 1998 are
               ! for calculating anx, csx and gswx for Rubisco limited,
               ! RuBP limited, sink limited.
               ! JK: vx4 changed to correspond to formulation in Collatz et al. 1992
               temp2(i,:) = rad%qcan(i,:,1) * jtomol * (1.0-veg%frac4(i))
               vx3(i,:)   = ej3x(temp2(i,:), veg%alpha(i), veg%convex(i), ejmxt3(i,:))
               vx4(i,:)   = veg%alpha(i) * rad%qcan(i,:,1) * jtomol * veg%frac4(i) * fwsoil_nongs(i)**qb
               !vx4(i,:)   = ej4x(temp2(i,:), veg%alpha(i), veg%convex(i), vcmxt4(i,:))
               rdx(i,:)   = veg%cfrd(i)*Vcmxt3(i,:) + veg%cfrd(i)*vcmxt4(i,:)

               !Vanessa - the trunk does not contain xleauning as of Ticket#56 inclusion
               !as well as other inconsistencies here that need further investigation. In the
               !interests of getting this into the trunk ASAP just isolate this code for now
               !default side of this condition is to use trunk version

               !#ifdef VanessasCanopy


               if (cable_user%CALL_climate) then
                  ! if (veg%iveg(i).eq.1) then
                  !    temp2(i,1) = rad%qcan(i,1,1) * jtomol * (1.0-veg%frac4(i))
                  !    temp2(i,2) = rad%qcan(i,2,1) * jtomol * (1.0-veg%frac4(i))
                  !    vx3(i,1)  = ej3x(temp2(i,1),climate%frec(i)*veg%alpha(i), &
                  !         veg%convex(i),ejmxt3(i,1))
                  !    vx3(i,2)  = ej3x(temp2(i,2),climate%frec(i)*veg%alpha(i), &
                  !         veg%convex(i),ejmxt3(i,2))
                  !    temp(i) =  xvcmxt3(tlfx(i)) * veg%vcmax(i) * (1.0-veg%frac4(i))
                  !    vcmxt3(i,1) = rad%scalex(i,1) * temp(i) * climate%frec(i)
                  !    vcmxt3(i,2) = rad%scalex(i,2) * temp(i) * climate%frec(i)
                  ! endif

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

                  if (veg%iveg(i) == 2) then ! evergreen broadleaf forest

                     rdx(i,1) = 1.0*(1.2818e-6+0.0116*veg%vcmax(i)- &
                        0.0334*climate%qtemp_max_last_year(i)*1e-6)

                     rdx(i,2) = rdx(i,1)

                  elseif (veg%iveg(i) == 4) then ! decid broadleaf forest

                     rdx(i,1) = 1.0*(1.2818e-6+0.0116*veg%vcmax(i)- &
                        0.0334*climate%qtemp_max_last_year(i)*1e-6)

                     rdx(i,2) = rdx(i,1)

                  elseif (veg%iveg(i) == 1) then ! evergreen needleleaf forest

                     rdx(i,1) = 1.0*(1.2877e-6+0.0116*veg%vcmax(i)- &
                        0.0334*climate%qtemp_max_last_year(i)*1e-6)

                     rdx(i,2) = rdx(i,1)

                  elseif (veg%iveg(i) == 3) then ! decid needleleaf forest

                     rdx(i,1) = 1.0*(1.2877e-6+0.0116*veg%vcmax(i)- &
                        0.0334*climate%qtemp_max_last_year(i)*1e-6)

                     rdx(i,2) = rdx(i,1)

                  elseif ((veg%iveg(i) == 6) .or. (veg%iveg(i) == 8) .or. &
                     (veg%iveg(i) == 9)) then ! C3 grass, tundra, crop

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
                  rdx(i,1) = rdx(i,1) * xrdt(tlfx(i)) * rad%scalex(i,1) * fwsoil_nongs(i)**qb
                  rdx(i,2) = rdx(i,2) * xrdt(tlfx(i)) * rad%scalex(i,2) * fwsoil_nongs(i)**qb

                  ! reduction of daytime leaf dark-respiration to account for
                  !photo-inhibition
                  rdx(i,1) = rdx(i,1) * light_inhibition(rad%qcan(i,1,1)*jtomol*1.0e6)
                  rdx(i,2) = rdx(i,2) * light_inhibition(rad%qcan(i,2,1)*jtomol*1.0e6)

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

               ! calculate canopy-level mesophyll conductance
               if (cable_user%explicit_gm) then
                  gmes(i,1) = veg%gm(i) * rad%scalex(i,1) * xgmesT(tlfx(i)) * max(0.15,fwsoil(i)**qm)
                  gmes(i,2) = veg%gm(i) * rad%scalex(i,2) * xgmesT(tlfx(i)) * max(0.15,fwsoil(i)**qm)

                  !gmes(i,1) = gmes(i,1) * MAX(0.15_r_2,real(fwsoil(i)**qm,r_2))
                  !gmes(i,2) = gmes(i,2) * MAX(0.15_r_2,real(fwsoil(i)**qm,r_2))
               else ! gmes = 0.0 easier to debug than inf
                  gmes(i,1) = 0.0
                  gmes(i,2) = 0.0
               endif

               ! Ticket #56 added switch for Belinda Medlyn's model
               IF (cable_user%GS_SWITCH == 'leuning') THEN
                  ! vh: added calculation of CO2 compensation point
                  if ((vcmxt3(i,1)+vcmxt4(i,1)) > 1.0e-8) then
                     co2cp(i,1) = (cx2(i)/2.0 +  conkct(i) * (1.0 + 0.21/conkot(i))* &
                        rdx(i,1)/(vcmxt3(i,1)+vcmxt4(i,1)))/ &
                        (1.0 - rdx(i,1)/(vcmxt3(i,1)+vcmxt4(i,1)))
                  else
                     co2cp(i,1) = 0.0
                  endif

                  if ((vcmxt3(i,2)+vcmxt4(i,2)) > 1.0e-8) then
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
                  if (INDEX(cable_user%FWSOIL_SWITCH,'LWP')>0) then
                     fwpsi(i,1) = (1.0_r_2 +exp(veg%slope_leaf(i) * veg%psi_50_leaf(i))) / &
                      (1.0_r_2+exp(veg%slope_leaf(i) * (veg%psi_50_leaf(i)-psilx(i,1))))
                     fwpsi(i,2) = (1.0_r_2 +exp(veg%slope_leaf(i) * veg%psi_50_leaf(i))) / &
                      (1.0_r_2+exp(veg%slope_leaf(i) * (veg%psi_50_leaf(i)-psilx(i,2))))

                     
                     gs_coeff(i,1) = ( fwpsi(i,1)**qs / ( real(csx(i,1)) - co2cp(i,1) ) ) &
                        * ( 1.0 / ( 1.0/veg%a1gs(i) + (vpd/(veg%d0gs(i)*veg%a1gs(i)))))

                     gs_coeff(i,2) = ( fwpsi(i,2)**qs / ( real(csx(i,2)) - co2cp(i,2) ) ) &
                        * ( 1.0 / ( 1.0/veg%a1gs(i) + (vpd/(veg%d0gs(i)*veg%a1gs(i)))))
                     
                  else
                     gs_coeff(i,1) = ( fwsoil(i)**qs / ( real(csx(i,1)) - co2cp(i,1) ) ) &
                        * ( 1.0 / ( 1.0/veg%a1gs(i) + (vpd/(veg%d0gs(i)*veg%a1gs(i)))))

                     gs_coeff(i,2) = ( fwsoil(i)**qs / ( real(csx(i,2)) - co2cp(i,2) ) ) &
                        * ( 1.0 / ( 1.0/veg%a1gs(i) + (vpd/(veg%d0gs(i)*veg%a1gs(i)))))
                  endif
                  if (present(fwpsdo)) then
                     gs_coeff(i,1) = ( veg%a1gs(i) / ( real(csx(i,1)) - co2cp(i,1) ) ) &
                        * fwpsdo
                     gs_coeff(i,2) = ( veg%a1gs(i) / ( real(csx(i,2)) - co2cp(i,2) ) ) &
                        * fwpsdo
                  endif
                  ! vh re-write so that a1 and d0 are not correlated


                  ! gs_coeff(i,1) = ( fwsoil(i)**qs / ( csx(i,1) - co2cp3 ) ) &
                  !      * ( veg%a1gs(i) / ( 1.0 + dsx(i)/veg%d0gs(i)))

                  ! gs_coeff(i,2) = ( fwsoil(i)**qs / ( csx(i,2) - co2cp3 ) ) &
                  !      * ( veg%a1gs(i) / ( 1.0 + dsx(i)/veg%d0gs(i)))

                  ! Medlyn BE et al (2011) Global Change Biology 17: 2134-2144.
               ELSEIF(cable_user%GS_SWITCH == 'medlyn' .AND. &
                  cable_user%FWSOIL_SWITCH /= 'profitmax') THEN

                  gswmin(i,1) = veg%g0(i) * rad%scalex(i,1)
                  gswmin(i,2) = veg%g0(i) * rad%scalex(i,2)

                  IF (dsx(i) < 50.0) THEN
                     vpd  = 0.05 ! kPa
                  ELSE
                     vpd = dsx(i) * 1E-03 ! Pa -> kPa
                     !vpd = met%dva(i) * 1E-03 ! Pa -> kPa
                  END IF

                  g1 = veg%g1(i)

                  !gs_coeff(i,1) = (1.0* fwsoil(i)**qs + (g1 * fwsoil(i)**qs) / SQRT(vpd)) / csx(i,1)
                  !gs_coeff(i,2) = (1.0* fwsoil(i)**qs + (g1 * fwsoil(i)**qs) / SQRT(vpd)) / csx(i,2)

                  ! gs_coeff for CO2, hence no 1.6 factor as in the original model
                  if (INDEX(cable_user%FWSOIL_SWITCH,'LWP')>0) then
                     fwpsi(i,1) = (1.0_r_2 +exp(veg%slope_leaf(i) * veg%psi_50_leaf(i))) / &
                      (1.0_r_2+exp(veg%slope_leaf(i) * (veg%psi_50_leaf(i)-psilx(i,1))))
                     fwpsi(i,2) = (1.0_r_2 +exp(veg%slope_leaf(i) * veg%psi_50_leaf(i))) / &
                      (1.0_r_2+exp(veg%slope_leaf(i) * (veg%psi_50_leaf(i)-psilx(i,2))))

                     if (fwpsi(i,1) .LE. 0.05) then
                        gs_coeff(i,1) = (1.0* fwpsi(i,1)**qs + (g1 * fwsoil(i)**qs) / SQRT(vpd)) / real(csx(i,1))
                        gs_coeff(i,2) = (1.0* fwpsi(i,2)**qs + (g1 * fwsoil(i)**qs) / SQRT(vpd)) / real(csx(i,2))
                     else
                        gs_coeff(i,1) = (1.0 + (g1 * fwpsi(i,1)**qs) / SQRT(vpd)) / real(csx(i,1))
                        gs_coeff(i,2) = (1.0 + (g1 * fwpsi(i,2)**qs) / SQRT(vpd)) / real(csx(i,2))
                     endif
                  else
                        if (fwsoil(i) .LE. 0.05) then
                           gs_coeff(i,1) = (1.0* fwsoil(i)**qs + (g1 * fwsoil(i)**qs) / SQRT(vpd)) / real(csx(i,1))
                           gs_coeff(i,2) = (1.0* fwsoil(i)**qs + (g1 * fwsoil(i)**qs) / SQRT(vpd)) / real(csx(i,2))
                        else
                           gs_coeff(i,1) = (1.0 + (g1 * fwsoil(i)**qs) / SQRT(vpd)) / real(csx(i,1))
                           gs_coeff(i,2) = (1.0 + (g1 * fwsoil(i)**qs) / SQRT(vpd)) / real(csx(i,2))
                        endif
                  endif
                  if (present(fwpsdo)) then
                     gs_coeff(i,1) = ( veg%a1gs(i) / ( real(csx(i,1)) - co2cp(i,1) ) ) &
                        * fwpsdo
                     gs_coeff(i,2) = ( veg%a1gs(i) / ( real(csx(i,2)) - co2cp(i,2) ) ) &
                        * fwpsdo
                  endif
               ELSE IF (cable_user%GS_SWITCH == 'tuzet' .AND. &
                  INDEX(cable_user%FWSOIL_SWITCH,'LWP')>0) THEN

                  if (present(fwpsdo)) then
                     gs_coeff(i,1) = g1 / real(csx(i,1)) * fwpsdo
                     gs_coeff(i,2) = g1 / real(csx(i,2)) * fwpsdo
                  else
                     fwpsi(i,1) = (1.0_r_2 +exp(veg%slope_leaf(i) * veg%psi_50_leaf(i))) / &
                      (1.0_r_2+exp(veg%slope_leaf(i) * (veg%psi_50_leaf(i)-psilx(i,1))))
                     fwpsi(i,2) = (1.0_r_2 +exp(veg%slope_leaf(i) * veg%psi_50_leaf(i))) / &
                      (1.0_r_2+exp(veg%slope_leaf(i) * (veg%psi_50_leaf(i)-psilx(i,2))))
                     !print*, 1.0, real(1.0, r_2), 1.0 / 3.0_r_2, 1.0_r_2 / 3.0  
                     !print *, '!!!!!!!!!!!!!!! fwpsi:', fwpsi(i,1),psilx(i,1)
                     gs_coeff(i,1) =fwpsi(i,1) * g1 / real(csx(i,1))
                     gs_coeff(i,2) =fwpsi(i,2) * g1 / real(csx(i,2))
                     fwpsixx = fwpsi
                  endif
                     !gswmin(i,:) = veg%g0(i) * rad%scalex(i,:) * fwpsi(i,:)
                     gswmin(i,:) = g0xx(i,:) * rad%scalex(i,:)
                     !gswmin(i,1) = g0(i) * rad%scalex(i,1) * fwpsi(i,1)
                     !gswmin(i,2) = g0(i) * rad%scalex(i,2) * fwpsi(i,1)
                     g1 = veg%g1tuzet(i)
                     psilxx(i,:) = psilx(i,:)
                     

               ELSE IF (cable_user%GS_SWITCH == 'profitmax' .AND. &
                  cable_user%FWSOIL_SWITCH == 'profitmax') THEN
                  ! write(logn,*) '1:gs_coeff of 1', gs_coeff(1,1)
                  ! write(logn,*) '1:gs_coeff of 2', gs_coeff(1,2)
                  ! Profix-Max hydraulics model
                  vpd = dsx(i) * PA_TO_KPA ! this is leaf vpd
                  press = met%pmb(i) * MB_TO_KPA

                  IF (vpd < 0.05) THEN
                     ecx(i) = 0.0_r_2
                     ecxs(i) = 0.0_r_2
                     anx(i, 1) = 0.0 - rdx(i, 1)  !MC Why A=0 if VPD low?
                     anx(i, 2) = 0.0 - rdx(i, 2)
                  ELSE
                     !MC Photosynthesis is calculated in here: to change for "our"
                     !   photosynthesis model
                     !  CALL optimisation(canopy, rad, vpd, press, tlfx(i), &
                     !       csx, ssnow%psi_rootzone(i), &
                     !       kcmax, veg%kmax(i), veg%PLCcrit(i), &
                     !       veg%b_plant(i), veg%c_plant(i), resolution, vcmxt3, &
                     !       ejmxt3, rdx, vx3, cx1(i), an_canopy, e_canopy, &
                     !       avg_kcan, co2cp3, p, i)
                     ! zhanlu: difference between optimisation and photosynthesis_gm to
                     ! calculate Ac, Aj is co2cp3 Vs gammast. change co2cp3 to gammast = cx2/2

                     CALL optimisation(canopy, rad, vpd, press, tlfx(i), &
                        csx, ssnow%psi_rootzone(i), &
                        kcmax(i), veg%kmax(i), veg%PLCcrit(i), &
                        veg%b_plant(i), veg%c_plant(i), resolution, vcmxt3, &
                        ejmxt3, rdx, vx3, cx1(i), an_canopy, e_canopy, &
                        avg_kcan, cx2(i)/2, p, i, gswmin(i,:))

                     ! fix units for CABLE and pack into arrays
                     anx(i,1) = an_canopy(1) * UMOL_TO_MOL
                     anx(i,2) = an_canopy(2) * UMOL_TO_MOL

                     ! fix units for CABLE and pack sap flux E into arrays
                     !MC check for single/double precision
                     IF (e_canopy > 0.0) THEN
                        ecxs(i) = e_canopy * air%rlam(i) * MOL_WATER_2_G_WATER * G_TO_KG
                     ELSE
                        ecxs(i) = 0.0_r_2
                     END IF
                  END IF

               ELSE
                  PRINT *, cable_user%GS_SWITCH, cable_user%FWSOIL_SWITCH, veg%iveg(i)

                  write(*,*) 'gs_model_switch failed.'
#ifdef __MPI__
                  call MPI_Abort(0, 127, ierr) ! Do not know comm nor rank here
#else
                  stop 127
#endif
               ENDIF ! IF (cable_user%GS_SWITCH == 'leuning') THEN
               !#endif


            ENDIF !IF (canopy%vlaiw(i) > C%LAI_THRESH .AND. (abs_deltlf(i) > 0.1 or Any(abs_deltpsi(i,:) > 0.1))
               ! gmes is 0.0 if explicit_gm = FALSE (easier to debug)
            
            IF (cable_user%GS_SWITCH /= 'profitmax') THEN
            if (present(wbpsdo)) then
               CALL photosynthesis_gm( csx(:,:), &
                  spread(cx1(:),2,mf), &
                  spread(cx2(:),2,mf), &
                  gswmin(:,:), rdx(:,:), vcmxt3(:,:), &
                  vcmxt4(:,:), vx3(:,:), vx4(:,:), &
               ! Ticket #56, xleuning replaced with gs_coeff here
                  gs_coeff(:,:), rad%fvlai(:,:), modify_rule0, &
                  anx(:,:), fwsoil(:), qs, gmes(:,:), kc4(:,:), &
                  anrubiscox(:,:), anrubpx(:,:), ansinkx(:,:), eta_x(:,:), dAnx(:,:),i,ktau_tot,k )
            else
               CALL photosynthesis_gm( csx(:,:), &
                  spread(cx1(:),2,mf), &
                  spread(cx2(:),2,mf), &
                  gswmin(:,:), rdx(:,:), vcmxt3(:,:), &
                  vcmxt4(:,:), vx3(:,:), vx4(:,:), &
               ! Ticket #56, xleuning replaced with gs_coeff here
                  gs_coeff(:,:), rad%fvlai(:,:), modify_rule0, &
                  anx(:,:), fwsoil(:), qs, gmes(:,:), kc4(:,:), &
                  anrubiscox(:,:), anrubpx(:,:), ansinkx(:,:), eta_x(:,:), dAnx(:,:),i )
            endif
            ENDIF
            !DO i=1,mp


            IF (modify_rule0) Then
               DO kk=1, mf

                  IF (rad%fvlai(i,kk)>C%LAI_THRESH) THEN
                     !  write(logn,*) '2: gs_coeff of ',kk, ': ',gs_coeff(i,kk)


                     csx(i,kk) = real(met%ca(i),r_2) - real(C%RGBWC*anx(i,kk),r_2) / &
                        (gbhu(i,kk) + gbhf(i,kk))
                     csx(i,kk) = max(1.0e-4_r_2, csx(i,kk))
                     if (cable_user%FWSOIL_SWITCH /= 'profitmax') THEN

                        ! Ticket #56, xleuning replaced with gs_coeff here
                        if (cable_user%g0_switch == 'default') then
                           canopy%gswx(i,kk) = MAX( 1.e-3, gswmin(i,kk)*fwsoil(i)**qs + &
                              MAX( 0.0, C%RGSWC * gs_coeff(i,kk) * &
                              anx(i,kk) ) )
                        elseif (cable_user%g0_switch == 'maximum') then
                           ! set gsw to maximum of g0*fwsoil and humidity-dependent term,
                           ! according to third formulation suggested by Lombardozzi et al.,
                           ! GMD 10, 321-331, 2017, and applied in CLM
                           canopy%gswx(i,kk) = MAX( 1.e-3, max(gswmin(i,kk)*fwsoil(i)**qs, &
                              MAX( 0.0, C%RGSWC * gs_coeff(i,kk) * &
                              anx(i,kk) )) )

                        endif


                     endif
                     ! write(logn,*) '3: gsw of ',kk, ': ', canopy%gswx(i,kk)
                     !Recalculate conductance for water:
                     gw(i,kk) = 1.0 / ( 1.0 / canopy%gswx(i,kk) + &
                        1.0 / ( 1.075 * real(gbhu(i,kk) + gbhf(i,kk)) ) )
                     gw(i,kk) = max(gw(i,kk), 0.00001)

                     ! Modified psychrometric constant
                     ! (Monteith and Unsworth, 1990)
                     psycst(i,kk) = air%psyc(i) * REAL( ghr(i,kk) / gw(i,kk) )

                  ENDIF

               ENDDO
               ex(i,1) = real( ( air%dsatdk(i) * ( rad%rniso(i,1) - C%capp * C%rmair &
               * ( met%tvair(i) - met%tk(i) ) * rad%gradis(i,1) ) &
               + C%capp * C%rmair * vpdtmp(i) * ghr(i,1) ) &
               / ( air%dsatdk(i) + psycst(i,1) ),r_2)
               ex(i,2) = real( ( air%dsatdk(i) &
               * ( rad%rniso(i,2) - C%capp * C%rmair * ( met%tvair(i) - &
               met%tk(i) ) * rad%gradis(i,2) ) + C%capp * C%rmair * &
               vpdtmp(i) * ghr(i,2) ) / &
               ( air%dsatdk(i) + psycst(i,2) ), r_2)  
               ecx(i) = sum(ex(i,:))
               ex(i,:) = ex(i,:) * (1.0_r_2-real(canopy%fwet(i), r_2)) / real(air%rlam(i), r_2) 
               ! ecx(i) = real( ( air%dsatdk(i) * ( rad%rniso(i,1) - C%capp * C%rmair &
               !    * ( met%tvair(i) - met%tk(i) ) * rad%gradis(i,1) ) &
               !    + C%capp * C%rmair * met%dva(i) * ghr(i,1) ) &
               !    / ( air%dsatdk(i) + psycst(i,1) ) + ( air%dsatdk(i) &
               !    * ( rad%rniso(i,2) - C%capp * C%rmair * ( met%tvair(i) - &
               !    met%tk(i) ) * rad%gradis(i,2) ) + C%capp * C%rmair * &
               !    met%dva(i) * ghr(i,2) ) / &
               !    ( air%dsatdk(i) + psycst(i,2) ), r_2)
               
               IF (INDEX(cable_user%FWSOIL_SWITCH, 'LWP') > 0) then

                     
                  ! ex(i,:) = ex(i,:) * (1.0_r_2-real(canopy%fwet(i), r_2)) / real(air%rlam(i), r_2) 
                  ! convert from kg m-2 ground s-1 to mmol m-2 leaf s-1*
                  !ex(i,:)= ex(i,:) * 1.0e6_r_2/18.0_r_2  
                  ! psilx(i,1) = ssnow%psi_rootzone(i) - ex(i,1) / canopy%kplant(i)
                  ! psilx(i,2) = ssnow%psi_rootzone(i) - ex(i,2) / canopy%kplant(i)
                  
                  CALL calc_psix(ssnow, soil, canopy, veg, casapool,max(sum(real(ex(i,:),r_2)), 0.0_r_2),psixxi,kplantxi,i)
                  !if psixx(i)<-10
                  psixx(i) = psixxi
                  kplantx(i) = kplantxi
                  psilx(i,1) = psixx(i) - max(ex(i,1),0.0_r_2) / kplantx(i)
                  psilx(i,2) = psixx(i) - max(ex(i,2),0.0_r_2) / kplantx(i)
                  if (any(psilx(i,:) < -20.0_r_2) ) then
                     ! if (present(wbpsdo)) then
                     !    print*, 'saturation: psilx(i,:)=', psilx(i,:), iter, k
                     ! else
                     !  !  print*, 'psilx(i,:)=', psilx(i,:), iter, k
                     ! endif
                     where (psilx(i,:) < -10.0_r_2)
                     psilx(i,:) = -10.0_r_2
                     end where
                     !g0(i) = g0(i)*0.5
                     ! if (present(wbpsdo)) then
                     !    print*, 'saturation: ex and total_est_evap', sum(real(ex(i,:),r_2)), total_est_evap(i)
                     ! else
                     !   ! print*, 'ex and total_est_evap', sum(real(ex(i,:),r_2)), total_est_evap(i)
                     ! endif
                     if (sum(real(ex(i,:),r_2))>total_est_evap(i)) then
                     ex(i,1) = total_est_evap(i) * ex(i,1)/sum(real(ex(i,:),r_2))
                     ex(i,2) = total_est_evap(i) * ex(i,2)/sum(real(ex(i,:),r_2))
                     !ecx(i) = sum(ex(i,:))
                     ! if (present(wbpsdo)) then
                     !    print*, 'saturation: ex_modified: ',ex(i,:), sum(real(ex(i,:),r_2))
                     ! else
                     !   ! print*, 'ex and total_est_evap', sum(real(ex(i,:),r_2)), total_est_evap(i)
                     ! endif
                     CALL calc_psix(ssnow, soil, canopy, veg, casapool,max(sum(real(ex(i,:),r_2)), 0.0_r_2),psixxi,kplantxi,i)
                     !if psixx(i)<-10
                     psixx(i) = psixxi
                     kplantx(i) = kplantxi
                     psilx(i,1) = psixx(i) - max(ex(i,1),0.0_r_2) / kplantx(i)
                     psilx(i,2) = psixx(i) - max(ex(i,2),0.0_r_2) / kplantx(i)

                     g0xx(i,:) = ex(i,:) * 1000 /18.0 /(vpdtmp(i) /100.0 /met%pmb(i) /C%rmair * C%rmh2o) 
                     !met%dva = (qstvair - met%qvair) *  C%rmair/C%rmh2o * met%pmb * 100.0 
                     ! if (present(wbpsdo)) then
                     !    print*, 'saturation: psix_modified: ',psixx(i)
                     !    print*, 'saturation: kplant_modified: ',kplantx(i)
                     !    print*, 'saturation: psil_modified: ',psilx(i,:)
                     !    print*, 'saturation: g0_modified: ',g0xx(i,:)
                     ! else
                     !   ! print*, 'ex and total_est_evap', sum(real(ex(i,:),r_2)), total_est_evap(i)
                     ! endif
                     
                     endif
                     ! fwpsi(i,1) = (1.0_r_2 +exp(veg%slope_leaf(i) * veg%psi_50_leaf(i))) / &
                     !       (1.0_r_2+exp(veg%slope_leaf(i) * (veg%psi_50_leaf(i)-psilx(i,1))))
                     ! fwpsi(i,2) = (1.0_r_2 +exp(veg%slope_leaf(i) * veg%psi_50_leaf(i))) / &
                     !       (1.0_r_2+exp(veg%slope_leaf(i) * (veg%psi_50_leaf(i)-psilx(i,2))))
                     ! g0xx(i,:) = veg%g0(i) * fwpsi(i,:)
                  endif
         
                  ! if (any(allktau == ktau_tot) .and. (iter==4)) then
                  ! print*,'ktau_tot, i and k',ktau_tot, iter, k
                  ! if (present(wbpsdo)) then
                  !    print*, '    psdo ex, psix and kplant: ',max(sum(real(ex(i,:),r_2)), 0.0_r_2), psixx(i), kplantx(i)
                  !    print*, '    psdo psilx: ', psilx(i,1)
                  ! else

                  !    print*, '         ex, psix and kplant: ',max(sum(real(ex(i,:),r_2)), 0.0_r_2), psixx(i), kplantx(i)
                  !    print*, '         psilx: ', psilx(i,1)
                  ! endif
                  ! endif

                  ! !!!!!!!!!!!!!!!!!!! another option!!!!!!!!!!!!!!!!!!!!!!
                  ! psilx(i,1) = psixx(i) - max(ex(i,1),0.0_r_2) / kplantx(i)
                  ! psilx(i,2) = psixx(i) - max(ex(i,2),0.0_r_2) / kplantx(i)
                  ! if (any(allktau == ktau_tot)) then
                  !    print*,'ktau_tot, i and k',ktau_tot, iter, k
                  !    print*, '         ex, psilx, Lastkplant: ', max(sum(real(ex(i,:),r_2)), 0.0_r_2), psilx(i,1), kplantx(i)
                  ! endif
                  ! CALL calc_psix(ssnow, soil, canopy, veg, casapool,max(sum(real(ex(i,:),r_2)), 0.0_r_2),psixxi,kplantxi,i)
                  ! psixx(i) = psixxi
                  ! kplantx(i) = kplantxi
                  ! if (any(allktau == ktau_tot)) then
                  !    print*, '         psixx and kplant: ', psixx(i), kplantx(i)
                  ! endif

               ENDIF
               IF (cable_user%SOIL_SCHE == 'Haverd2013' .or. (INDEX(cable_user%FWSOIL_SWITCH, 'Haverd2013') > 0)) then
                  ! avoid root-water extraction when fwsoil is zero
                  if (fwsoil(i) < 1e-6) then
                     anx(i,:) = -rdx(i,:)
                     ecx(i)   = 0.0_r_2
                  endif

                  canopy%fevc(i) = ecx(i) * (1.0_r_2-real(canopy%fwet(i),r_2))

                  call getrex_1d(wbtmp(i,:)-real(ssnow%wbice(i,:),r_2), ssnow%rex(i,:), &
                     canopy%fwsoil(i), &
                     real(veg%froot(i,:),r_2), SPREAD(real(soil%ssat(i),r_2),1,ms), &
                     SPREAD(real(soil%swilt(i),r_2),1,ms), &
                     max(canopy%fevc(i)/real(air%rlam(i),r_2)/1000.0_r_2, 0.0_r_2), &
                     veg%gamma(i), &
                     real(soil%zse, r_2), real(dels,r_2), veg%zr(i))
                  canopy%fwsoiltmp(i) = real(canopy%fwsoil(i))
                  where (ssnow%rex(i,:) > tiny(1.0_r_2)) &
                  ssnow%evapfbl(i,:) = real(ssnow%rex(i,:))*dels*1000. ! mm water &
                  IF (INDEX(cable_user%FWSOIL_SWITCH, 'Haverd2013') > 0) then
                     fwsoil(i) = real(canopy%fwsoil(i))
                     !(root water extraction) per time step
   
                     if (cable_user%Cumberland_soil) then
                        canopy%fwsoil(i) = max(canopy%fwsoil(i), 0.6_r_2)
                        fwsoil(i) = real(canopy%fwsoil(i))
                     endif
                  ENDIF 


               ELSE IF (cable_user%FWSOIL_SWITCH == 'profitmax' .AND. cable_user%SOIL_SCHE == 'hydraulics') THEN


                  IF (ecx(i) > 0.0 .AND. canopy%fwet(i) < 1.0) THEN
                     evapfb(i) = ( 1.0 - canopy%fwet(i)) * REAL( ecx(i) ) *dels  &
                        / air%rlam(i)

                     DO kk = 1,ms

                        !ssnow%evapfbl(i,kk) = MIN(evapfb(i) * &
                        !                          ssnow%fraction_uptake(i,kk),  &
                        !                          MAX(0.0, &
                        !                              REAL(ssnow%wb(i,kk)) -    &
                        !                              soil%swilt(i)) * &
                        !                          soil%zse(kk) * 1000.0)

                        ! ms8355: no bounding by swilt in this version, or the
                        ! "beneits" from hydraulics are cancelled
                        ssnow%evapfbl(i,kk) = evapfb(i) * &
                           ssnow%fraction_uptake(i,kk)



                     ENDDO
                     canopy%fevc(i) = SUM(ssnow%evapfbl(i,:))*air%rlam(i)/dels

                     ecx(i) = canopy%fevc(i) / (1.0-canopy%fwet(i))
                  ENDIF

               ENDIF

               IF (cable_user%SOIL_SCHE .ne.'Haverd2013') THEN

                  if (ecx(i) > 0.0_r_2 .and. canopy%fwet(i) < 1.0) then
                     evapfb(i) = ( 1.0 - canopy%fwet(i)) * real(ecx(i)) *dels &
                        / air%rlam(i)

                     DO kk = 1,ms

                        ssnow%evapfbl(i,kk) = MIN( evapfb(i) * veg%froot(i,kk), &
                           MAX( 0.0, REAL( wbtmp(i,kk) ) - &
                           1.1 * soil%swilt(i) ) * &
                           soil%zse(kk) * 1000.0 )

                     ENDDO
                     IF (cable_user%soil_struc=='default') then
                        canopy%fevc(i) = real(sum(ssnow%evapfbl(i,:))*air%rlam(i)/dels, r_2)
                        ecx(i) = canopy%fevc(i) / (1.0_r_2-real(canopy%fwet(i),r_2))
                     ELSEIF (cable_user%soil_struc=='sli') then
                        canopy%fevc(i) = ecx(i) * (1.0_r_2-real(canopy%fwet(i),r_2))
                     ENDIF

                  ENDIF

               ENDIF
               ! Update canopy sensible heat flux:
               hcx(i) = ( SUM(real(rad%rniso(i,:),r_2)) - ecx(i) &
                  - real(C%capp*C%rmair*(met%tvair(i)-met%tk(i)), r_2) &
                  * SUM(real(rad%gradis(i,:),r_2)) ) &
                  * SUM(real(gh(i,:),r_2))/ SUM(real(ghr(i,:),r_2))

               ! Update leaf temperature:
               tlfx(i) = met%tvair(i)+REAL(hcx(i))/(C%capp*C%rmair*SUM(gh(i,:)))

               ! Update net radiation for canopy:
               rnx(i) = real( SUM( rad%rniso(i,:)) - &
                  C%CAPP * C%rmair *( tlfx(i)-met%tk(i) ) * &
                  SUM( rad%gradis(i,:) ), r_2)

               ! Update leaf surface vapour pressure deficit:
               dsx(i) = vpdtmp(i) + air%dsatdk(i) * (tlfx(i)-met%tvair(i))
               if (cable_user%perturb_dva_by_T) then
                  dsx(i) = vpdtmp(i) + air%dsatdk(i) * (tlfx(i) -met%tvair(i) + cable_user%dva_T_perturbation )
               endif

               dsx(i)=  max(dsx(i),0.0)

               ! Store change in leaf temperature between successive iterations:
               deltlf(i) = tlfxx(i)-tlfx(i)
               
               if (cable_user%GS_SWITCH == 'tuzet' .AND. &
                  INDEX(cable_user%FWSOIL_SWITCH,'LWP')>0) then
                  abs_deltpsil(i,:) = ABS(psilxx(i,:)-psilx(i,:))
               fwpsi(i,1) = (1.0_r_2 +exp(veg%slope_leaf(i) * veg%psi_50_leaf(i))) / &
                     (1.0_r_2+exp(veg%slope_leaf(i) * (veg%psi_50_leaf(i)-psilx(i,1))))
               fwpsi(i,2) = (1.0_r_2 +exp(veg%slope_leaf(i) * veg%psi_50_leaf(i))) / &
                     (1.0_r_2+exp(veg%slope_leaf(i) * (veg%psi_50_leaf(i)-psilx(i,2))))
                  abs_deltfwpsi(i,:) = ABS(fwpsixx(i,:)-fwpsi(i,:))
               endif
               abs_deltlf(i) = ABS(deltlf(i))
               abs_deltds(i) = ABS(dsxx(i)-dsx(i))
               abs_deltcs(i,:) = ABS(csxx(i,:)-csx(i,:))
            ENDIF !modify_rule0

            !ENDDO !i=1,mp
            psilxm(i,:) = 0.0_r_2
            tlfxm(i) = 0.0_r_2
            dsxm(i) = 0.0_r_2
            Mtag = 0
            modify_rule = (abs_deltlf(i) > 0.1 .or. Any(abs_deltpsil(i,:) > 0.1) .or. Any(abs_deltfwpsi(i,:) > 0.02)) &
            .AND. (k < kmax) 
            !if ( (abs_deltlf(i) > 0.1 .or. Any(abs_deltpsil(i,:) > 0.1)) .AND. k > 5 .AND. k < kmax ) then
            if ( modify_rule ) then
               Numtag = 3
               ! after 4 iterations, take mean of current & previous estimates
               ! as the next estimate of leaf temperature, to avoid oscillation
               dc = veg%dc(i)
               ! tlfx(i) = ( 0.5 * ( MAX( 0, k-5 ) / ( k - 4.9999 ) ) ) *tlfxx(i) + &
               !    ( 1.0 - ( 0.5 * ( MAX( 0, k-5 ) / ( k - 4.9999 ) ) ) ) &
               !    * tlfx(i)
               tlfxm(i) = dc *tlfxx(i) + ( 1.0 - dc ) * tlfx(i)
               dsxm(i) = dc *dsxx(i) + ( 1.0 - dc ) * dsx(i)
               csx(i,1) = dc *csxx(i,1) + ( 1.0 - dc ) * csx(i,1)
               csx(i,2) = dc *csxx(i,2) + ( 1.0 - dc ) * csx(i,2)
               psilxm(i,1) = dc *psilxx(i,1) + ( 1.0 - dc ) * psilx(i,1)
               psilxm(i,2) = dc *psilxx(i,2) + ( 1.0 - dc ) * psilx(i,2)
               fwpsi1_tmp(i,1) = 0.90 *fwpsixx(i,1) + ( 1.0 - 0.90 ) * fwpsi(i,1)
               fwpsi1_tmp(i,2) = 0.90 *fwpsixx(i,2) + ( 1.0 - 0.90 ) * fwpsi(i,2)
               ! if (fwpsixx(i,1)<0.45 .AND. fwpsixx(i,1)>0.3) then
               !    if (fwpsi(i,1)-fwpsixx(i,1)>0.0_r_2) then
               !       fwpsi1_tmp(i,1) = fwpsixx(i,1)+0.01
               !    else
               !       fwpsi1_tmp(i,1) = fwpsixx(i,1)-0.01
               !    endif
               !    ! fwpsi1_tmp(i,1) = 0.90 *fwpsixx(i,1) + ( 1.0 - 0.90 ) * fwpsi(i,1)
               !    ! fwpsi1_tmp(i,2) = 0.90 *fwpsixx(i,2) + ( 1.0 - 0.90 ) * fwpsi(i,2)
               ! else
               !    fwpsi1_tmp(i,1) = dc *fwpsixx(i,1) + ( 1.0 - dc ) * fwpsi(i,1)
               ! endif
               ! if (fwpsixx(i,2)<0.45 .AND. fwpsixx(i,2)>0.3) then
               !    if (fwpsi(i,2)-fwpsixx(i,2)>0.0_r_2) then
               !       fwpsi1_tmp(i,2) = fwpsixx(i,2)+0.01
               !    else
               !       fwpsi1_tmp(i,2) = fwpsixx(i,2)-0.01
               !    endif
               !    ! fwpsi1_tmp(i,1) = 0.90 *fwpsixx(i,1) + ( 1.0 - 0.90 ) * fwpsi(i,1)
               !    ! fwpsi1_tmp(i,2) = 0.90 *fwpsixx(i,2) + ( 1.0 - 0.90 ) * fwpsi(i,2)
               ! else
               !    fwpsi1_tmp(i,2) = dc *fwpsixx(i,2) + ( 1.0 - dc ) * fwpsi(i,2)
               ! endif
               ! psilxm(i,1) = veg%psi_50_leaf(i) - (1.0_r_2 / veg%slope_leaf(i)) * &
               ! log( (1.0_r_2 + exp(veg%slope_leaf(i)*veg%psi_50_leaf(i)) - fwpsi1_tmp(i,1)) / fwpsi1_tmp(i,1) )
               ! psilxm(i,2) = veg%psi_50_leaf(i) - (1.0_r_2 / veg%slope_leaf(i)) * &
               ! log( (1.0_r_2 + exp(veg%slope_leaf(i)*veg%psi_50_leaf(i)) - fwpsi1_tmp(i,2)) / fwpsi1_tmp(i,2) )
               Mtag=1
               ! inc = 0.1_r_2
               ! if (k>10) then
               !    if (psilx(i,1)-psilxx(i,1)<-0.5_r_2) then
               !       psilxm(i,1) = psilxx(i,1) - inc
               !    elseif (psilx(i,1)-psilxx(i,1)>0.5_r_2) then
               !       psilxm(i,1) = psilxx(i,1) + inc
               !    endif
               !    if (psilx(i,2)-psilxx(i,2)<-0.5_r_2) then
               !       psilxm(i,2) = psilxx(i,2) - inc
               !    elseif (psilx(i,2)-psilxx(i,2)>0.5_r_2) then
               !       psilxm(i,2) = psilxx(i,2) + inc
               !    endif
               ! !elseif (k>=16) then
                  
               ! endif
      
               ! ! calculate the new fwpsi1_tmp based on modified psilx
               ! fwpsi1_tmp(i,1) = (1.0_r_2 +exp(veg%slope_leaf(i) * veg%psi_50_leaf(i))) / &
               !        (1.0_r_2+exp(veg%slope_leaf(i) * (veg%psi_50_leaf(i)-psilxm(i,1))))
               ! fwpsi1_tmp(i,2) = (1.0_r_2 +exp(veg%slope_leaf(i) * veg%psi_50_leaf(i))) / &
               !        (1.0_r_2+exp(veg%slope_leaf(i) * (veg%psi_50_leaf(i)-psilxm(i,2))))
               ! if (fwpsi1_tmp(i,1)-fwpsixx(i,1)<-0.1_r_2) then
               !    Mtag=2
               !    if (present(wbpsdo)) then
               !       canopy%N_neg_sw = canopy%N_neg_sw + 1.0
               !    else
               !       canopy%N_neg = canopy%N_neg + 1.0
               !    endif
               !    inc = 0.1_r_2
               !    fwpsi1_tmp(i,1) = fwpsixx(i,1) - inc
               !    psilxm(i,1) = veg%psi_50_leaf(i) - (1.0_r_2 / veg%slope_leaf(i)) * &
               !    log( (1.0_r_2 + exp(veg%slope_leaf(i)*veg%psi_50_leaf(i)) - fwpsi1_tmp(i,1)) / fwpsi1_tmp(i,1) )
               !    fwpsi1_tmp(i,2) = fwpsixx(i,2) - inc
               !    psilxm(i,2) = veg%psi_50_leaf(i) - (1.0_r_2 / veg%slope_leaf(i)) * &
               !    log( (1.0_r_2 + exp(veg%slope_leaf(i)*veg%psi_50_leaf(i)) - fwpsi1_tmp(i,2)) / fwpsi1_tmp(i,2) )
               ! elseif (fwpsi1_tmp(i,1)-fwpsixx(i,1)>0.1_r_2) then
               !    Mtag=3
               !    if (present(wbpsdo)) then
               !       canopy%N_pos_sw = canopy%N_pos_sw + 1.0
               !    else
               !       canopy%N_pos = canopy%N_pos + 1.0
               !    endif
               !    fwpsi1_tmp(i,1) = fwpsixx(i,1) + inc
               !    psilxm(i,1) = veg%psi_50_leaf(i) - (1.0_r_2 / veg%slope_leaf(i)) * &
               !    log( (1.0_r_2 + exp(veg%slope_leaf(i)*veg%psi_50_leaf(i)) - fwpsi1_tmp(i,1)) / fwpsi1_tmp(i,1) )
               !    fwpsi1_tmp(i,2) = fwpsixx(i,2) + inc
               !    psilxm(i,2) = veg%psi_50_leaf(i) - (1.0_r_2 / veg%slope_leaf(i)) * &
               !    log( (1.0_r_2 + exp(veg%slope_leaf(i)*veg%psi_50_leaf(i)) - fwpsi1_tmp(i,2)) / fwpsi1_tmp(i,2) )   
               ! endif
            !elseif ( k > 5 .AND. k < kmax ) then
            elseif ( k < kmax ) then
               Numtag = Numtag - 1
            endif
            IF ( abs_deltlf(i) < ABS( deltlfy(i) ) ) THEN

               deltlfy(i)       = deltlf(i)
               tlfy(i)          = tlfx(i)
               dsy(i) = dsx(i)
               !print*,'when dT<last dT, tlfy: ', tlfy(i)
               psily(i,:)       = psilx(i,:)
               rny(i)           = rnx(i)
               hcy(i)           = hcx(i)
               ecy(i)           = ecx(i)
               rdy(i,:)         = rdx(i,:)
               an_y(i,:)        = anx(i,:)
               anrubiscoy(i,:)  = anrubiscox(i,:)
               anrubpy(i,:)     = anrubpx(i,:)
               dAn_y(i,:)       = dAnx(i,:)
               eta_y(i,:)       = eta_x(i,:)
               dAnrubiscoy(i,:) = dAnrubiscox(i,:)
               dAnrubpy(i,:)    = dAnrubpx(i,:)
               dAnsinky(i,:)    = dAnsinkx(i,:)

               ! save last values calculated for ssnow%evapfbl
               oldevapfbl(i,:) = ssnow%evapfbl(i,:)
               gswy(i,:)        = canopy%gswx(i,:)
               fwpsiy(i,:)     = fwpsi(i,:)
               csy(i,:) = csx(i,:)
               gs_coeffy(i,:) = gs_coeff(i,:)
               vcmxt3y(i,:) = vcmxt3(i,:)
               vcmxt4y(i,:) = vcmxt4(i,:)
               cx2y(i) = cx2(i)
            ENDIF
            ! if (ktau>=5184) then
            ! print*, 'check y=x ktau & k= ',ktau,k
            ! end if
            ! if (cable_user%GS_SWITCH == 'tuzet' .AND. &
            !       INDEX(cable_user%FWSOIL_SWITCH,'LWP') > 0 .AND. &
            !       Any(abs_deltpsil(i,:) > 0.5) ) then
            !    ! after 4 iterations, take mean of current & previous estimates
            !    ! as the next estimate of leaf temperature, to avoid oscillation
            !    psilx(i,1) = ( 0.5 * ( MAX( 0, k-5 ) / ( k - 4.9999 ) ) ) *psilxx(i,1) + &
            !       ( 1.0 - ( 0.5 * ( MAX( 0, k-5 ) / ( k - 4.9999 ) ) ) ) &
            !       * psilx(i,1)
            !    psilx(i,2) = ( 0.5 * ( MAX( 0, k-5 ) / ( k - 4.9999 ) ) ) *psilxx(i,2) + &
            !       ( 1.0 - ( 0.5 * ( MAX( 0, k-5 ) / ( k - 4.9999 ) ) ) ) &
            !       * psilx(i,2)
            ! endif
            ! if (ktau>=5184) then
            ! print*, 'check tlfx ktau & k= ',ktau,k
            ! end if
            IF (k==1) THEN
               ! take the first iterated estimates as the defaults
               tlfy(i) = tlfx(i)
               dsy(i) = dsx(i)
               !print*,'when k=1, tlfy: ', tlfy(i)
               psily(i,:) = psilx(i,:)
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
               gswy(i,:)       = canopy%gswx(i,:)
               fwpsiy(i,:)     = fwpsi(i,:)
               csy(i,:) = csx(i,:)
               gs_coeffy(i,:) = gs_coeff(i,:)
               vcmxt3y(i,:) = vcmxt3(i,:)
               vcmxt4y(i,:) = vcmxt4(i,:)
               cx2y(i) = cx2(i)
            END IF
            !print*, 'check after k==1 ',ktau,k
            !if (ktau_tot>=nktau .and. ktau_tot<=(nktau+NN-1)) then

            ! if (ktau>=5184) then
            ! print*, 'write 134 ',ktau,k
            ! end if
            if (any(allktau == ktau_tot)) then
               if (present(wbpsdo) .and. present(vpdpsdo)) then
                  write(137,*) ktau_tot, iter, i, k, tlfxx(i),tlfx(i),tlfxm(i), deltlf(i), &
                  dsxx(i),dsx(i),dsxm(i),abs_deltds(i),Numtag, Mtag, &
                  psilxx(i,1),psilxx(i,2),psilx(i,1), psilx(i,2),psilxm(i,1), psilxm(i,2), &
                  abs_deltpsil(i,1),abs_deltpsil(i,2),abs_deltfwpsi(i,1),abs_deltfwpsi(i,2), &
                  fwpsixx(i,1),fwpsixx(i,2),fwpsi(i,1),fwpsi(i,2),fwpsi1_tmp(i,1),fwpsi1_tmp(i,2), &
                  psixx(i), csx(i,1), csx(i,2),abs_deltcs(i,1), abs_deltcs(i,2),anx(i,1), anx(i,2), anrubiscox(i,1), &
                  anrubiscox(i,2), anrubpx(i,1),anrubpx(i,2),ansinkx(i,1),ansinkx(i,2), &
                  canopy%gswx(i,1), canopy%gswx(i,2),canopy%gswx(i,1), canopy%gswx(i,2), &
                  vcmxt3(i,1),vcmxt3(i,2),gs_coeff(i,1),gs_coeff(i,2),rdx(i,1),rdx(i,2),ex(i,1),ex(i,2), &
                  ssnow%rootR(i,1),ssnow%rootR(i,2),ssnow%rootR(i,3),ssnow%rootR(i,4),ssnow%rootR(i,5),ssnow%rootR(i,6), &
                  ssnow%soilR(i,1),ssnow%soilR(i,2),ssnow%soilR(i,3),ssnow%soilR(i,4),ssnow%soilR(i,5),ssnow%soilR(i,6), &
                  kplantx(i),g0xx(i,1),g0xx(i,2),fwsoil(i)
               elseif (present(wbpsdo) .and. present(fwpsdo)) then
               elseif (present (wbpsdo)) then 
                  write(135,*) ktau_tot, iter, i, k, tlfxx(i),tlfx(i),tlfxm(i), deltlf(i), &
                  dsxx(i),dsx(i),dsxm(i),abs_deltds(i),Numtag, Mtag, &
                  psilxx(i,1),psilxx(i,2),psilx(i,1), psilx(i,2),psilxm(i,1), psilxm(i,2), &
                  abs_deltpsil(i,1),abs_deltpsil(i,2),abs_deltfwpsi(i,1),abs_deltfwpsi(i,2), &
                  fwpsixx(i,1),fwpsixx(i,2),fwpsi(i,1),fwpsi(i,2),fwpsi1_tmp(i,1),fwpsi1_tmp(i,2), &
                  psixx(i), csx(i,1), csx(i,2),abs_deltcs(i,1), abs_deltcs(i,2),anx(i,1), anx(i,2), anrubiscox(i,1), &
                  anrubiscox(i,2), anrubpx(i,1),anrubpx(i,2),ansinkx(i,1),ansinkx(i,2), &
                  canopy%gswx(i,1), canopy%gswx(i,2),canopy%gswx(i,1), canopy%gswx(i,2), &
                  vcmxt3(i,1),vcmxt3(i,2),gs_coeff(i,1),gs_coeff(i,2),rdx(i,1),rdx(i,2),ex(i,1),ex(i,2), &
                  ssnow%rootR(i,1),ssnow%rootR(i,2),ssnow%rootR(i,3),ssnow%rootR(i,4),ssnow%rootR(i,5),ssnow%rootR(i,6), &
                  ssnow%soilR(i,1),ssnow%soilR(i,2),ssnow%soilR(i,3),ssnow%soilR(i,4),ssnow%soilR(i,5),ssnow%soilR(i,6), &
                  kplantx(i),g0xx(i,1),g0xx(i,2),fwsoil(i)
               elseif (present (vpdpsdo)) then 
                  write(136,*) ktau_tot, iter, i, k, tlfxx(i),tlfx(i),tlfxm(i), deltlf(i), &
                  dsxx(i),dsx(i),dsxm(i),abs_deltds(i),Numtag, Mtag, &
                  psilxx(i,1),psilxx(i,2),psilx(i,1), psilx(i,2),psilxm(i,1), psilxm(i,2), &
                  abs_deltpsil(i,1),abs_deltpsil(i,2), abs_deltfwpsi(i,1),abs_deltfwpsi(i,2), &
                  fwpsixx(i,1),fwpsixx(i,2),fwpsi(i,1),fwpsi(i,2),fwpsi1_tmp(i,1),fwpsi1_tmp(i,2), &
                  psixx(i), csx(i,1), csx(i,2),abs_deltcs(i,1), abs_deltcs(i,2),anx(i,1), anx(i,2), anrubiscox(i,1), &
                  anrubiscox(i,2), anrubpx(i,1),anrubpx(i,2),ansinkx(i,1),ansinkx(i,2), &
                  canopy%gswx(i,1), canopy%gswx(i,2),canopy%gswx(i,1), canopy%gswx(i,2), &
                  vcmxt3(i,1),vcmxt3(i,2),gs_coeff(i,1),gs_coeff(i,2),rdx(i,1),rdx(i,2),ex(i,1),ex(i,2), &
                  ssnow%rootR(i,1),ssnow%rootR(i,2),ssnow%rootR(i,3),ssnow%rootR(i,4),ssnow%rootR(i,5),ssnow%rootR(i,6), &
                  ssnow%soilR(i,1),ssnow%soilR(i,2),ssnow%soilR(i,3),ssnow%soilR(i,4),ssnow%soilR(i,5),ssnow%soilR(i,6), &
                  kplantx(i),g0xx(i,1),g0xx(i,2),fwsoil(i)
               else
                  write(134,*) ktau_tot, iter, i, k, tlfxx(i),tlfx(i),tlfxm(i), deltlf(i), &
                  dsxx(i),dsx(i),dsxm(i),abs_deltds(i), Numtag,Mtag, &
                  psilxx(i,1),psilxx(i,2),psilx(i,1), psilx(i,2),psilxm(i,1), psilxm(i,2), &
                  abs_deltpsil(i,1),abs_deltpsil(i,2), abs_deltfwpsi(i,1),abs_deltfwpsi(i,2), &
                  fwpsixx(i,1),fwpsixx(i,2),fwpsi(i,1),fwpsi(i,2),fwpsi1_tmp(i,1),fwpsi1_tmp(i,2), &
                  psixx(i), csx(i,1), csx(i,2),abs_deltcs(i,1), abs_deltcs(i,2),anx(i,1), anx(i,2), anrubiscox(i,1), &
                  anrubiscox(i,2), anrubpx(i,1),anrubpx(i,2),ansinkx(i,1),ansinkx(i,2), &
                  canopy%gswx(i,1), canopy%gswx(i,2),canopy%gswx(i,1), canopy%gswx(i,2), &
                  vcmxt3(i,1),vcmxt3(i,2),gs_coeff(i,1),gs_coeff(i,2),rdx(i,1),rdx(i,2),ex(i,1),ex(i,2), &
                  ssnow%rootR(i,1),ssnow%rootR(i,2),ssnow%rootR(i,3),ssnow%rootR(i,4),ssnow%rootR(i,5),ssnow%rootR(i,6), &
                  ssnow%soilR(i,1),ssnow%soilR(i,2),ssnow%soilR(i,3),ssnow%soilR(i,4),ssnow%soilR(i,5),ssnow%soilR(i,6), &
                  kplantx(i),g0xx(i,1),g0xx(i,2),fwsoil(i)

               endif

            END IF
            !if ( (abs_deltlf(i) > 0.1 .or. Any(abs_deltpsil(i,:) > 0.1)) .AND. k > 5 .AND. k < kmax ) then
            if ( modify_rule ) then
               psilx = psilxm
               tlfx = tlfxm
               dsx = dsxm
            endif
         END DO  ! DO WHILE k < C%MAXITER)
      ENDDO !i=1,mp


         ! print*, 'DD28 ', rad%fvlai
         ! print*, 'DD29 ', met%ca
         ! print*, 'DD30 ', canopy%gswx
         ! print*, 'DD31 ', air%dsatdk
         ! print*, 'DD32 ', ssnow%wbice
         ! print*, 'DD33 ', ssnow%rex
         ! print*, 'DD34 ', veg%froot
         ! print*, 'DD35 ', soil%ssat
         ! print*, 'DD36 ', soil%swilt
         ! print*, 'DD37 ', canopy%fevc
         ! print*, 'DD38 ', veg%gamma
         ! print*, 'DD39 ', air%rlam
         ! print*, 'DD40 ', veg%zr
         ! print*, 'DD41 ', soil%zse
         ! print*, 'DD42 ', ssnow%rex
         ! print*, 'DD43 ', ssnow%evapfbl
         ! if (ktau>=5184) then
         ! print*, 'check after getrex_1d: ktau & k= ',ktau,k
         ! end if
         ! Where leaf temp change b/w iterations is significant, and
         ! difference is smaller than the previous iteration, store results:
         



         
         !if (k == C%MAXITER .AND. ktau>=5760 .AND. ktau<=11664) then
      ! if (k == C%MAXITER .AND. ktau_tot>=nktau .and. ktau_tot<=(nktau+NN-1)) then
      ! !if (k == C%MAXITER .AND. ktau>=5760 .AND. ktau<=11664) then
      !        print*, 'value x:',anrubiscox(:,1),anx(:,1)
      ! endif
      
      !print*,'when k end, tlfy: ', tlfy(1)
      ! if (ktau_tot>=nktau .and. ktau_tot<=(nktau+NN-1)) then
      ! i = 1
      ! write(134,*) ktau, iter, i, 21, tlfy(i), deltlfy(i), &
      ! dsx(i), psily(i,1), psily(i,2),fwpsiy(:,1),fwpsiy(:,2),canopy%fwpsi(i,1),canopy%fwpsi(i,2),&
      ! csy(i,1), csy(i,2),csx(i,1), csx(i,2), an_y(i,1), an_y(i,2), &
      ! gswy(i,1), gswy(i,2),canopy%gswx(i,1), canopy%gswx(i,2), &
      ! vcmxt3(i,1),vcmxt3(i,2), gs_coeff(i,1),gs_coeff(i,2),rdy(i,1),rdy(i,2)
      ! end if 
      !if (ktau_tot==(nktau+NN-1)  .and. iter==4) THEN
      !if (any(nktau_end == ktau_tot) .and. iter==4) THEN
      if (ktau_tot == nktau_end(size(nktau_end)) .and. iter==4) THEN
         if (present(wbpsdo) .and. present(vpdpsdo)) then
            close(137)
         elseif (present(wbpsdo) .and. present(fwpsdo)) then
         elseif (present (wbpsdo)) then 
            close(135)
         elseif (present (vpdpsdo)) then 
            close(136)
         else
            close(134)
         endif
      END IF
      deallocate(nktau) 
      deallocate(nktau_end) 
      deallocate(allktau) 

      
      ! print*, 'End k loop: ktau & iter & k= ',ktau,iter,k
      ! dry canopy flux
      canopy%fevc = (1.0_r_2-real(canopy%fwet,r_2)) * ecy
      if (present(wbpsdo) .and. present(vpdpsdo)) then
         canopy%epotref = canopy%fevc
         canopy%gsw_ref = canopy%gswx
         canopy%abs_deltpsil_ref = abs_deltpsil
         canopy%abs_deltcs_ref = abs_deltcs * 1.0e6_r_2
         canopy%abs_deltlf_ref = abs_deltlf
         canopy%abs_deltds_ref = abs_deltds
      elseif (present(wbpsdo) .and. present(fwpsdo)) then
         canopy%epotref1 = canopy%fevc
         canopy%gsw_ref1 = canopy%gswx
         canopy%abs_deltpsil_ref1 = abs_deltpsil
         canopy%abs_deltcs_ref1 = abs_deltcs * 1.0e6_r_2
         canopy%abs_deltlf_ref1 = abs_deltlf
         canopy%abs_deltds_ref1 = abs_deltds
      elseif (present (wbpsdo)) then
         canopy%epotcan3 = canopy%fevc
         canopy%frday_epotcan3 = 12.0 * SUM(rdy, 2)
         !! vh !! inserted min to avoid -ve values of GPP
         canopy%fpn_epotcan3 = min(-12.0 * SUM(an_y, 2), canopy%frday)
         canopy%gsw_epotcan3 = canopy%gswx
         canopy%abs_deltpsil_sw = abs_deltpsil
         canopy%abs_deltcs_sw = abs_deltcs * 1.0e6_r_2
         canopy%abs_deltlf_sw = abs_deltlf
         canopy%abs_deltds_sw = abs_deltds
      elseif (present (vpdpsdo)) then
         
         canopy%gsw_epotvpd = canopy%gswx / C%RGSWC
         canopy%abs_deltpsil_vpd = abs_deltpsil
         canopy%abs_deltcs_vpd = abs_deltcs * 1.0e6_r_2
         canopy%abs_deltlf_vpd = abs_deltlf
         canopy%abs_deltds_vpd = abs_deltds
      else
         canopy%abs_deltpsil = abs_deltpsil
         canopy%abs_deltcs = abs_deltcs * 1.0e6_r_2
         canopy%abs_deltlf = abs_deltlf
         canopy%abs_deltds = abs_deltds


         !write(logn,*) 'fevc', canopy%fevc(1)
         ! print*, 'DD45 ', canopy%fwet
         ! print*, 'DD46 ', canopy%fevc
         IF (cable_user%fwsoil_switch.NE.'Haverd2013' .AND. cable_user%fwsoil_switch.NE.'profitmax') THEN

            ! Recalculate ssnow%evapfbl as ecy may not be updated with the ecx
            ! calculated in the last iteration.
            ! DO NOT use simple scaling as there are times that ssnow%evapfbl is zero.
            ! ** ssnow%evapfbl(i,:) = ssnow%evapfbl(i,:) * ecy(i) / ecx(i) **
            DO i = 1, mp

               IF( ecy(i) > 0.0_r_2 .AND. canopy%fwet(i) < 1.0 ) THEN

                  IF( ABS( ecy(i) - ecx(i) ) > 1.0e-6_r_2 ) THEN

                     IF( ABS( canopy%fevc(i) - real(SUM(oldevapfbl(i,:)) * air%rlam(i) &
                        /dels, r_2) ) > 1.0e-4_r_2 ) THEN

                        WRITE(*,*) 'ERROR! oldevapfbl not right.', ktau_gl, i
                        WRITE(*,*) 'ecx, ecy = ', ecx(i), ecy(i)
                        WRITE(*,*) 'or in mm = ', ecx(i) * ( 1.0 - canopy%fwet(i) ) &
                           / air%rlam(i) * dels, &
                           ecy(i) * ( 1.0 - canopy%fwet(i) ) / &
                           air%rlam(i) * dels

                        WRITE(*,*)'fevc = ', canopy%fevc(i), SUM( oldevapfbl(i,:) ) * &
                           air%rlam(i) / dels
                        WRITE(*,*) 'fwet = ', canopy%fwet(i)
                        WRITE(*,*) 'oldevapfbl = ', oldevapfbl(i,:)
                        WRITE(*,*) 'ssnow%evapfbl before rescaling: ', &
                           ssnow%evapfbl(i,:)
                        ! STOP

                     ELSE

                        ssnow%evapfbl(i,:) = oldevapfbl(i,:)

                     END IF

                  END IF

               END IF

            END DO

         ENDIF
         ! print*, 'Recalculate ssnow%evapfbl: ktau: ',ktau
         canopy%psi_can = real(psily, r_2)
         canopy%frday = 12.0 * SUM(rdy, 2)
         !! vh !! inserted min to avoid -ve values of GPP
         canopy%fpn = min(-12.0 * SUM(an_y, 2), canopy%frday)
         ! change canopy%gswx from gswx to gswy, zihanlu, 19/12/2024
         canopy%gswx = gswy
         canopy%fwpsi = fwpsiy
         ! additional diagnostic variables for assessing contributions of rubisco and rubp limited photosynthesis to
         ! net photosynthesis in sunlit and shaded leaves.
         canopy%A_sh = real(an_y(:,2), r_2)
         canopy%A_sl = real(an_y(:,1), r_2)

         canopy%GPP_sh = real(an_y(:,2) + rdy(:,2), r_2)
         canopy%GPP_sl = real(an_y(:,1) + rdy(:,1), r_2)

         canopy%A_shC = real(anrubiscoy(:,2), r_2)
         canopy%A_shJ = real(anrubpy(:,2), r_2)

         where (anrubiscoy(:,2) > an_y(:,2)) canopy%A_shC = 0.0_r_2
         where (anrubpy(:,2) > an_y(:,2))    canopy%A_shJ = 0.0_r_2
         
         canopy%A_slC = real(anrubiscoy(:,1), r_2)
         canopy%A_slJ = real(anrubpy(:,1), r_2)
         !if (ktau>=5760 .AND. ktau<=11664) then
         ! if (ktau_tot>=nktau .and. ktau_tot<=(nktau+NN-1)) then
         !    print*, 'value y:',canopy%A_slC, canopy%A_sl
         ! endif
         where (anrubiscoy(:,1) > an_y(:,1)) canopy%A_slC = 0.0_r_2
         where (anrubpy(:,1)    > an_y(:,1)) canopy%A_slJ = 0.0_r_2

         canopy%eta_A_cs = canopy%A_sh * min(eta_y(:,2),5.0_r_2) + canopy%A_sl * min(eta_y(:,1),5.0_r_2)
         where ((canopy%A_sl > 0.0_r_2) .and. (canopy%A_sh > 0.0_r_2))
            canopy%eta_GPP_cs = canopy%GPP_sh  * min(eta_y(:,2),5.0_r_2) * canopy%GPP_sh/canopy%A_sh + &
               canopy%GPP_sl * min(eta_y(:,1),5.0_r_2)* canopy%GPP_sl/canopy%A_sl
         elsewhere ((canopy%A_sl .le. 0.0_r_2) .and. (canopy%A_sh > 0.0_r_2))
            canopy%eta_GPP_cs = canopy%GPP_sh * min(eta_y(:,2),5.0_r_2) * canopy%GPP_sh/canopy%A_sh
         elsewhere
            canopy%eta_GPP_cs = canopy%eta_A_cs
         end where
         ! change from gs_coeff to gs_coeffy
         where (canopy%A_sl > 0.0_r_2)
            canopy%eta_A_cs_sl    =  min(eta_y(:,1),5.0_r_2)
            canopy%eta_fevc_cs_sl = (min(eta_y(:,1),5.0_r_2) - 1.0_r_2) * &
               real(max(0.0, gs_coeffy(:,1)*an_y(:,1)), r_2) / real(canopy%gswx(:,1), r_2)
         elsewhere
            canopy%eta_A_cs_sl    = 0.0_r_2
            canopy%eta_fevc_cs_sl = 0.0_r_2
         endwhere
         ! change from gs_coeff to gs_coeffy
         where (canopy%A_sh > 0.0_r_2)
            canopy%eta_A_cs_sh    =  min(eta_y(:,2),5.0_r_2)
            canopy%eta_fevc_cs_sh = (min(eta_y(:,2),5.0_r_2) - 1.0_r_2) * &
               real(max(0.0, gs_coeffy(:,2)*an_y(:,2)), r_2) / real(canopy%gswx(:,2), r_2)
         elsewhere
            canopy%eta_A_cs_sh    = 0.0_r_2
            canopy%eta_fevc_cs_sh = 0.0_r_2
         endwhere

         where ((rad%fvlai(:,1)+rad%fvlai(:,2)) > 0.01)
            canopy%eta_fevc_cs = ( canopy%eta_fevc_cs_sl *  real(rad%fvlai(:,1), r_2) + &
               canopy%eta_fevc_cs_sh * real(rad%fvlai(:,2), r_2) ) * canopy%fevc / &
               real(rad%fvlai(:,1)+rad%fvlai(:,2), r_2)
         elsewhere
            canopy%eta_fevc_cs = 0.0_r_2
         endwhere
         ! canopy%cs_sl = csx(:,1) * 1.0e6_r_2
         ! canopy%cs_sh = csx(:,2) * 1.0e6_r_2
         ! change canopy%cs from csx to csy, zihanlu, 19/12/2024
         canopy%cs_sl = csy(:,1) * 1.0e6_r_2
         canopy%cs_sh = csy(:,2) * 1.0e6_r_2

         canopy%dAdcs = canopy%A_sl * dAn_y(:,1) + canopy%A_sh * dAn_y(:,2)
         canopy%cs    = canopy%A_sl * canopy%cs_sl * 1.0e6_r_2 + canopy%A_sh * canopy%cs_sh * 1.0e6_r_2

         canopy%tlf   = real(tlfy, r_2)
         ! change canopy%dlf from dsx to dsy, zihanlu, 19/12/2024
         !canopy%dlf   = real(dsx, r_2)
         canopy%dlf   = real(dsy, r_2)


         ! print*, 'DD47 ', ssnow%evapfbl
         canopy%evapfbl = ssnow%evapfbl

         ! 13C
         canopy%An        = real(an_y, r_2)
         canopy%Rd        = real(rdy, r_2)
         canopy%isc3      = (1.0-veg%frac4) > epsilon(1.0)
         ! change from vcmxt3 to vcmxt3y, 
         canopy%vcmax     = real(merge(vcmxt3y, vcmxt4y, spread(canopy%isc3,2,mf)), r_2)
         ! change from cx2 to cx2y, 
         canopy%gammastar = real(spread(cx2y*0.5,2,mf), r_2)
         canopy%gsc       = real(canopy%gswx / C%rgswc, r_2)
         canopy%gbc       = (gbhu + gbhf) / real(C%rgbwc, r_2)
         canopy%gac       = huge(1.0_r_2)
         ! replace 1.0_r_2/canopy%gac by tiny(1.0_r_2) to avoid underflow
         canopy%ci        = real(spread(met%ca, 2, mf), r_2) - &
            canopy%An / &
            (1.0_r_2 / (1.0_r_2/canopy%gbc + 1.0_r_2/canopy%gsc + tiny(1.0_r_2)))
         IF (cable_user%FWSOIL_SWITCH == 'profitmax') THEN
            canopy%fevcs = (1.0-canopy%fwet) * ecxs  ! trans. from sapflux

            DO i = 1, mp

               new_plc_sat = calc_plc(kcmax(i), veg%kmax(i), veg%PLCcrit(i))

               IF (canopy%psi_can_opt(i) < ssnow%psi_rootzone(i) .AND. ecxs(i) > 1E-9) THEN
                  IF (avg_kcan(i) > 1E-9 .AND. avg_kcan(i) < veg%kmax(i)) THEN

                     ! Partitioning ratios from Drake et al. (2015), PC&E, 38: 1628-1636.
                     new_plc_stem = calc_plc(5. / (2. / avg_kcan(i) + 3. / kcmax(i)), &
                        veg%kmax(i), veg%PLCcrit(i)) ! 2:5 to 3:5

                     new_plc_can = calc_plc(avg_kcan(i), veg%kmax(i), veg%PLCcrit(i))

                  ELSE
                     new_plc_stem = 0.0
                     new_plc_can = 0.0
                  ENDIF

               ELSE
                  new_plc_stem = 0.0
                  new_plc_can = 0.0
               ENDIF

               IF (new_plc_sat > canopy%day_plc_sat(i)) THEN
                  canopy%day_plc_sat(:) = new_plc_sat
               ENDIF

               IF (new_plc_stem > canopy%day_plc_stem(i)) THEN
                  canopy%day_plc_stem(:) = new_plc_stem
               ENDIF

               IF (new_plc_can > canopy%day_plc_can(i)) THEN
                  canopy%day_plc_can(:) = new_plc_can
               ENDIF

               ! set the day's PLC and reset the PLC tracker here
               IF (met%hod(i) >= 24.0 - dels / (60.0 * 60.0) - 1E-6) THEN
                  canopy%day_plc_sat(:) = 0.0
                  canopy%day_plc_stem(:) = 0.0
                  canopy%day_plc_can(:) = 0.0

               ELSE IF (met%hod(i) >= 24.0 - 2.0 * dels / (60.0 * 60.0) - 1E-6) THEN
                  canopy%plc_sat(i) = canopy%day_plc_sat(i)
                  canopy%plc_stem(i) = canopy%day_plc_stem(i)
                  canopy%plc_can(i) = canopy%day_plc_can(i)
               ENDIF

            END DO
         ENDIF
         ! deallocate( gswmin )
      endif ! if not present wbpsdo

   END SUBROUTINE dryLeaf

SUBROUTINE dryLeaf_givengs(ktau, ktau_tot, dels, rad, air, met, &
      veg, canopy, soil, ssnow,casapool, dsx, dsy, &
      tlfx, tlfy, ecy, hcy, &
      rny, gbhu, gbhf, csx, &
      cansat, iter, climate,gspsdo)

      use cable_def_types_mod
      use cable_common_module
      USE casavariable
#ifdef __MPI__
      use mpi, only: MPI_Abort
#endif

      implicit none
      INTEGER,                   INTENT(IN)    :: ktau
      INTEGER,                   INTENT(IN)    :: ktau_tot
      real,                      intent(in)    :: dels ! integration time step (s)
      type(radiation_type),      intent(in) :: rad
      type(air_type),            intent(in) :: air
      type(met_type),            intent(in) :: met
      type(veg_parameter_type),  intent(in) :: veg
      type(canopy_type),         intent(in) :: canopy
      type(soil_parameter_type), intent(in) :: soil
      type(soil_snow_type),      intent(in) :: ssnow
      TYPE (casa_pool),  INTENT(IN)           :: casapool
      real,      dimension(:),   intent(inout) :: &
         dsx,        & ! leaf surface vpd
         dsy,        & ! leaf surface vpd

         ! fwsoiltmp,  &
         tlfx,       & ! leaf temp prev. iter (K)
         tlfy  ! leaf temp (K)
      real(r_2), dimension(:),   intent(inout) :: &
         ecy,        & ! lat heat fl dry big leaf
         hcy,        & ! veg. sens heat
         rny         !& !
    
      real(r_2), dimension(:,:), intent(inout) :: &
         gbhu,       & ! forcedConvectionBndryLayerCond
         gbhf,       & ! freeConvectionBndryLayerCond
         csx   ! leaf surface CO2 concentration
         ! psilx,      &
         ! psily,      &
         ! fwpsi   
      real,      dimension(:),   intent(in)    :: cansat
      integer,                   intent(in)    :: iter
      type(climate_type),        intent(in)    :: climate
      real(r_2), dimension(:,:), intent(in), optional :: gspsdo

      !local variables (JK: move parameters to different routine at some point)
      real, parameter :: jtomol = 4.6e-6  ! Convert from J to Mol for light

      real, dimension(mp) :: &
         conkct,        & ! Michaelis Menton const.
         conkot,        & ! Michaelis Menton const.
         cx1,           & ! "d_{3}" in Wang and Leuning,
         cx2,           & !     1998, appendix E
         cx2y,         &
         tdiff,         & ! leaf air temp diff.
         tlfxx,         & ! leaf temp of current iteration (K)
         dsxx,         & ! vpd at leaf surface of current iteration 
         abs_deltlf,    & ! ABS(deltlf)
         abs_deltds,    &
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
         temp_c3,       & !
         temp_c4,       & !
         temp_sun_c3,   & !
         temp_shade_c3, & !
         temp_sun_c4,   & !
         temp_shade_c4, &    !
         fwsoil

      real(r_2), dimension(mp)  :: &
         ecx,        & ! lat. hflux big leaf
         hcx,        & ! sens heat fl big leaf prev iteration
         rnx, &          ! net rad prev timestep
         ecxs, &  ! lat. hflux big leaf (sap flux)
         psixx, &
         kplantx, &
         ghwet
         
      real(r_2) :: psixxi, kplantxi
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
         vcmxt3y,     & ! vcmax big leaf C3
         vcmxt4y,     & ! vcmax big leaf C4
         vx3,        & ! carboxylation C3 plants
         vx4,        & ! carboxylation C4 plants
         co2cp,      & ! CO2 compensation point (needed for Leuing stomatal conductance)
      ! Ticket #56, xleuning is no longer used, we replace it with gs_coeff,
      ! which is computed differently based on the new GS_SWITCH. If GS_SWITCH
      ! is "leuning", it's the same, if "medlyn", then the new Medlyn model
      ! xleuning,   & ! leuning stomatal coeff
         gs_coeff,   & ! stom coeff, Ticket #56
         gs_coeffy,   & 
         psycst,     & ! modified pych. constant
         frac42,     & ! 2D frac4
         temp2,      &
         anrubiscox, & ! net photosynthesis (rubisco limited)
         anrubpx,    & ! net photosynthesis (rubp limited)
         ansinkx,    & ! net photosynthesis (sink limited)
         anrubiscoy, & ! net photosynthesis (rubisco limited)
         anrubpy,    & ! net photosynthesis (rubp limited)
         kc4,        & ! An-Ci slope in Collatz 1992 model
         gmes,       & ! mesophyll conductance of big leaf
         gswy,       &
         csy

      real(r_2), dimension(mp,mf)  ::  dAnrubiscox, & ! CO2 elasticity of net photosynthesis (rubisco limited)
         dAnrubpx, &   !  (rubp limited)
         dAnsinkx, &  ! (sink limited)
         dAnx, & !(actual rate)
         dAnrubiscoy, & ! CO2 elasticity of net photosynthesis (rubisco limited)
         dAnrubpy, &   !  (rubp limited)
         dAnsinky, &  ! (sink limited)
         dAn_y, & !(actual rate)
         eta_y, &
         eta_x, &
         ex,    &  !transpiration, kg m-2 s-1
         fwpsiy, &
         psilxx, &
         csxx, &
         abs_deltpsil , &
         abs_deltcs

      real ::  gam0,    &
         conkc0,  &
         conko0,  &
         egam,    &
         ekc,     &
         eko,     &
         qs,      &
         qm,      &
         qb

      ! real, dimension(:,:), pointer :: gswmin => null() ! min stomatal conductance
      ! real, dimension(:,:), allocatable :: gswmin ! min stomatal conductance
      real, dimension(mp,mf) :: gswmin ! min stomatal conductance

      real, dimension(mp,2) ::  gsw_term, lower_limit2  ! local temp var
      real(r_2), dimension(mp,ms) :: wbtmp
      real(r_2), dimension(mp):: vpdtmp
      integer :: i, j, k, kk, h,iter_ini ! iteration count
      integer :: NN,m
      integer, allocatable :: nktau(:), allktau(:), nktau_end(:)
      real :: vpd, g1, ktot, refill  ! Ticket #56
      REAL, PARAMETER :: & ! Ref. params from Bernacchi et al. (2001)
         co2cp325 = 42.75, & ! CO2 compensation pt C3 at 25 degrees, umol mol-1
         Eaco2cp325 = 37830. ! activation energy for the CO2 compensation pt
      REAL :: co2cp3 ! CO2 compensation pt C3, N.B.: this is set to 0 by default in CABLE

      REAL :: trans_mmol, press
      REAL, PARAMETER :: KG_2_G = 1000.0
      REAL, PARAMETER :: G_WATER_TO_MOL = 1.0 / 18.01528
      REAL, PARAMETER :: MOL_2_MMOL = 1000.0
      REAL, PARAMETER :: MMOL_2_MOL = 1.0 / MOL_2_MMOL
      REAL, PARAMETER :: MB_TO_PA = 100.
      INTEGER, PARAMETER :: resolution = 800 ! jumps in Ci ~ 0.5 umol mol-1
      REAL, DIMENSION(2) :: an_canopy
      REAL :: e_canopy
      REAL(r_2), DIMENSION(resolution) :: p

      REAL :: MOL_WATER_2_G_WATER, G_TO_KG, UMOL_TO_MOL, MB_TO_KPA, PA_TO_KPA
      REAL, DIMENSION(mp) :: kcmax, avg_kcan
      REAL :: new_plc_sat, new_plc_stem, new_plc_can
      REAL :: MOL_TO_UMOL, J_TO_MOL, dc
      CHARACTER(LEN=200) :: txtname,num_str
  
#ifdef __MPI__
      integer :: ierr
#endif
      iter_ini = 4
      if (iter==iter_ini) then 
         if ((cable_user%soil_struc=='default') .and. (INDEX(cable_user%FWSOIL_SWITCH, 'Haverd2013') <= 0)) then
            if (cable_user%fwsoil_switch == 'standard') then
               call fwsoil_calc_std(fwsoil, soil, ssnow, veg, wbtmp)
            elseif (cable_user%fwsoil_switch == 'non-linear extrapolation') then
               !EAK, 09/10 - replace linear approx by polynomial fitting
               call fwsoil_calc_non_linear(fwsoil, soil, ssnow, veg, wbtmp)
            elseif (cable_user%fwsoil_switch == 'Lai and Katul 2000') then
               call fwsoil_calc_Lai_Katul(fwsoil, soil, ssnow, veg, wbtmp)
            elseif (cable_user%FWSOIL_SWITCH == 'profitmax') then
               fwsoil = 1.0
            elseif (cable_user%FWSOIL_SWITCH == 'constant1') then
               fwsoil = 0.98
            elseif (INDEX(cable_user%FWSOIL_SWITCH, 'LWP') > 0) then
               fwsoil = 1.0

            else
               write(*,*) 'fwsoil_switch failed.'
#ifdef __MPI__
               call MPI_Abort(0, 126, ierr) ! Do not know comm nor rank here
#else
               stop 126
#endif
            endif
            ! canopy%fwsoil = real(fwsoil, r_2)
            ! canopy%fwsoiltmp = real(fwsoil, r_2)
         elseif ((cable_user%soil_struc=='sli') .or. (INDEX(cable_user%FWSOIL_SWITCH, 'Haverd2013') > 0)) then
            fwsoil = real(canopy%fwsoil)
            ! canopy%fwsoiltmp = real(fwsoil, r_2)
         endif
      
      endif
      ! END header
      ! allocate(gswmin(mp,mf))


      wbtmp = ssnow%wb

      vpdtmp = met%dva

     ! write(logn,*) 'fwsoil of ', cable_user%fwsoil_switch, ': ', fwsoil(1)
      MOL_TO_UMOL = 1E6
      J_TO_MOL = 4.6E-6  ! Convert from J to Mol for light
      MOL_WATER_2_G_WATER = 18.02
      G_TO_KG = 1E-03
      UMOL_TO_MOL = 1E-06
      MB_TO_KPA = 0.1
      PA_TO_KPA = 1E-03
      ! weight min stomatal conductance by C3 an C4 plant fractions
      frac42       = spread(veg%frac4, 2, mf) ! frac C4 plants
      gsw_term     = spread(veg%gswmin, 2, mf)
      lower_limit2 = rad%scalex * gsw_term
      gswmin       = max(1.e-6, lower_limit2)
      ! print*, 'DD02 ', veg%frac4
      ! print*, 'DD03 ', veg%gswmin
      ! print*, 'DD04 ', rad%scalex

      gw          = 1.0e-3 ! default values of conductance
      gh          = 1.0e-3
      ghr         = 1.0e-3
      rdx         = 0.0
      anx         = 0.0
      anrubiscox  = 0.0
      anrubpx     = 0.0
      ansinkx     = 0.0
      dAnrubiscox = 0.0_r_2
      dAnrubpx    = 0.0_r_2
      dAnsinkx    = 0.0_r_2
      dAnx        = 0.0
      eta_x       = 0.0_r_2
      rnx         = SUM(real(rad%rniso,r_2),2)
      abs_deltlf  = 999.0
      abs_deltds  = 999.0
      abs_deltcs(:,:)   = SPREAD(SPREAD(999.0_r_2,1,mp),2,mf)
      vcmxt3      = 0.0
      vcmxt4      = 0.0
      ! print*, 'DD05 ', rad%rniso

      gras  = 1.0e-6
      an_y  = 0.0
      hcx   = 0.0_r_2              ! init sens heat iteration memory variable
      hcy   = 0.0_r_2
      rdy   = 0.0
      vcmxt3y      = 0.0
      vcmxt4y      = 0.0
      ecx   = SUM(real(rad%rniso,r_2),2) ! init lat heat iteration memory variable
      tlfxx = tlfx
      dsxx = dsx
      csxx = csx
      ! psilxx = psilx
      psycst(:,:)   = SPREAD(air%psyc,2,mf)
      canopy%fevc   = 0.0_r_2
      ssnow%evapfbl = 0.0

      ecxs = 0.0_r_2
      canopy%fevcs = 0.0_r_2

      ghwet = 1.0e-3_r_2
      gwwet = 1.0e-3
      ghrwet= 1.0e-3
      canopy%fevw = 0.0
      canopy%fhvw = 0.0
      sum_gbh        = real(SUM((gbhu+gbhf),2))
      sum_rad_rniso  = SUM(rad%rniso,2)
      sum_rad_gradis = SUM(rad%gradis,2)

      ! print*, 'DD06 ', rad%gradis

      ! default for variable d_3 of RuBP-limited photosynthesis of
      ! Wang and Leuning (1998), which is 2 * Gamma^star
      if (cable_user%explicit_gm) then
         if (trim(cable_user%Rubisco_parameters) == 'Bernacchi_2002') then
            gam0 = C%gam0cc
            egam = C%egamcc
         else if (trim(cable_user%Rubisco_parameters) == 'Walker_2013') then
            gam0 = C%gam0ccw
            egam = C%egamccw
         endif
      else
         gam0 = C%gam0
         egam = C%egam
      endif
      cx2(:) = 2.0 * gam0 * exp( egam / (C%rgas * C%trefk) &
         * (1.0 - C%trefk / tlfx(:)) )
      !a = k25 * exp((Ea * (Tk - 298.15)) / (298.15 * RGAS * Tk))

      DO kk=1,mp

         IF(canopy%vlaiw(kk) <= C%LAI_THRESH) THEN
            rnx(kk) = 0.0_r_2 ! intialise
            ecx(kk) = 0.0_r_2 ! intialise
            ecy(kk) = ecx(kk) ! store initial values
            ecxs(kk) = 0.0_r_2 ! initialise
            abs_deltlf(kk) = 0.0
            abs_deltds(kk) = 0.0
            ! abs_deltpsil(kk,:) = 0.0
            abs_deltcs(kk,:) = 0.0
            rny(kk) = rnx(kk) ! store initial values
            ! calculate total thermal resistance, rthv in s/m
         END IF

      ENDDO

      deltlfy = abs_deltlf
      

      !kdcorbin, 08/10 - doing all points all the time'
      !allocate(nktau(4), nktau_end(4))
      nktau=[79029,80861,82966,96215]
      NN=2
      nktau_end = nktau + NN - 1
      
      m = size(nktau) * NN
      allocate(allktau(m))
      do i = 1, size(nktau)
         do k = 0, NN-1
             j = (i - 1) * NN + k
             allktau(j+1) = nktau(i) + k
         end do
     end do
      !write(num_str, '(I0)') nktau
      ! txtname = trim(filename%path) // '/testIteration_cable_out_' // trim(num_str) &

      txtname = trim(filename%path) // 'testIteration_givengs_cable_out.txt'
 
      ! if (any(nktau == ktau_tot) .and. iter==1) then
      !    if (ktau_tot == nktau(1)) then
      !        ! Open the file for overwrite if k is the first element
      !        open(unit=134, file=txtname, status="unknown", action="write")
      !    else
      !        ! Open the file for appending otherwise
      !        open(unit=134, file=txtname, status="unknown", position="append", action="write")
      !    end if
      ! end if

 
         if (ktau_tot == nktau(1) .and. iter==4) then
            ! Open the file for overwrite if k is the first element
            open(unit=138, file=txtname, status="unknown", action="write")
         end if
   
      ! if (ktau_tot==nktau .and. iter==1) then
      ! open(unit=134, file=txtname)
      ! !print*, 'write iteration file '
      ! end if
      k = 0
      DO WHILE (k < C%MAXITER)
         k = k + 1
         ! print*, 'DD07 ', k, C%MAXITER, canopy%cansto
         ! print*, 'DD08 ', canopy%vlaiw
         ! print*, 'DD09 ', air%rlam
         ! print*, 'DD10 ', met%tvair
         ! print*, 'DD11 ', veg%dleaf
         ! print*, 'DD12 ', rad%fvlai
         ! print*, 'DD13 ', air%cmolar
         ! print*, 'DD14 ', veg%vcmax_sun
         ! print*, 'DD15 ', veg%vcmax_shade
         ! print*, 'DD16 ', veg%frac4
         ! print*, 'DD17 ', climate%mtemp
         ! print*, 'DD18 ', veg%ejmax_sun
         ! print*, 'DD19 ', veg%ejmax_shade
         ! print*, 'DD20 ', veg%alpha
         ! print*, 'DD21 ', veg%convex
         ! print*, 'DD22 ', rad%qcan
         ! print*, 'DD23 ', veg%cfrd
         ! print*, 'DD24 ', climate%qtemp_max_last_year
         ! print*, 'DD25 ', veg%a1gs
         ! print*, 'DD26 ', veg%d0gs
         ! print*, 'DD27 ', veg%g0
         DO i=1,mp

            IF (canopy%vlaiw(i) > C%LAI_THRESH .AND. (abs_deltlf(i) > 0.1 )) THEN

               ghwet(i)  = 2.0_r_2 * real(sum_gbh(i),r_2)
               gwwet(i)  = 1.075 * sum_gbh(i)
               ghrwet(i) = sum_rad_gradis(i) + real(ghwet(i))

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
               gbhf(i,1) = real( rad%fvlai(i,1) * air%cmolar(i) * 0.5*C%dheat &
                  * ( gras(i)**0.25 ) / veg%dleaf(i), r_2)
               gbhf(i,2) = real( rad%fvlai(i,2) * air%cmolar(i) * 0.5 * C%dheat &
                  * ( gras(i)**0.25 ) / veg%dleaf(i), r_2)
               gbhf(i,:) = max(1.e-6_r_2, gbhf(i,:))

               ! Conductance for heat:
               gh(i,:) = real(2.0_r_2 * (gbhu(i,:) + gbhf(i,:)))

               if (cable_user%amphistomatous) then
                  gbhf(i,:) = 2.0_r_2 * gbhf(i,:)
                  gh(i,:)   = real(gbhu(i,:) + gbhf(i,:))
               endif

               ! Conductance for heat and longwave radiation:
               ghr(i,:) = rad%gradis(i,:)+gh(i,:)

               ! Choose Ci-based (implicit gm) or Cc-based (explicit gm) parameters
               if (cable_user%explicit_gm) then
                  if (trim(cable_user%Rubisco_parameters) == 'Bernacchi_2002') then
                     gam0    = C%gam0cc
                     conkc0  = C%conkc0cc
                     conko0  = C%conko0cc
                     egam    = C%egamcc
                     ekc     = C%ekccc
                     eko     = C%ekocc
                  else if (trim(cable_user%Rubisco_parameters) == 'Walker_2013') then
                     gam0    = C%gam0ccw
                     conkc0  = C%conkc0ccw
                     conko0  = C%conko0ccw
                     egam    = C%egamccw
                     ekc     = C%ekcccw
                     eko     = C%ekoccw
                  endif
                  qs      = C%qs
                  qm      = C%qm
                  qb      = C%qb
               else
                  gam0   = C%gam0
                  conkc0 = C%conkc0
                  conko0 = C%conko0
                  egam   = C%egam
                  ekc    = C%ekc
                  eko    = C%eko
                  qs     = 0.5
                  qm     = 0.0     ! not used
                  qb     = 1.0
               endif
               ! JK: veg%vcmax_sun and veg%vcmax_shade are Cc-based if cable_user%explicit_gm = TRUE
               !     and Ci-based otherwise. If cable_user%coordinate_photosyn = FALSE,
               !     veg%vcmax_sun = veg%vcmax_shade = veg%vcmaxcc if cable_user%explicit_gm = TRUE
               !     and veg%vcmax_sun = veg%vcmax_shade = veg%vcmax otherwise.
               !     See also Subroutine 'casa_feedback' in casa_cable.F90. Same applies to ejmax.

               ! Leuning 2002 (PCE) equation for temperature response
               ! used for Vcmax for C3 plants:
               if (.not. cable_user%acclimate_photosyn) then
                  temp_sun_c3(i)   = xvcmxt3(tlfx(i)) * veg%vcmax_sun(i) * (1.0-veg%frac4(i))
                  temp_shade_c3(i) = xvcmxt3(tlfx(i)) * veg%vcmax_shade(i) * (1.0-veg%frac4(i))
                  temp_sun_c4(i)   = xvcmxt4(tlfx(i)) * veg%vcmax_sun(i) * veg%frac4(i)
                  temp_shade_c4(i) = xvcmxt4(tlfx(i)) * veg%vcmax_shade(i) * veg%frac4(i)
               else
                  call xvcmxt3_acclim(tlfx(i), climate%mtemp(i) , temp_c3(i))
                  call xvcmxt4_acclim(tlfx(i), climate%mtemp(i) , temp_c4(i))
                  temp_sun_c3(i)   = temp_c3(i) * veg%vcmax_sun(i) * (1.0-veg%frac4(i))
                  temp_shade_c3(i) = temp_c3(i) * veg%vcmax_shade(i) * (1.0-veg%frac4(i))
                  temp_sun_c4(i)   = temp_c4(i) * veg%vcmax_sun(i) * veg%frac4(i)
                  temp_shade_c4(i) = temp_c4(i) * veg%vcmax_shade(i) * veg%frac4(i)
               endif
               vcmxt3(i,1) = rad%scalex(i,1) * temp_sun_c3(i) * fwsoil(i)**qb
               vcmxt3(i,2) = rad%scalex(i,2) * temp_shade_c3(i) * fwsoil(i)**qb
               vcmxt4(i,1) = rad%scalex(i,1) * temp_sun_c4(i) * fwsoil(i)**qb
               vcmxt4(i,2) = rad%scalex(i,2) * temp_shade_c4(i) * fwsoil(i)**qb
               ! temperature response of the CO2 compensation point/gamma star used
               ! for C3 plants (Bernacchi et al. 2001), umol mol-1
               ! N.B.: this is always zero by default in CABLE, which is wrong.
               co2cp3 = arrh(co2cp325, Eaco2cp325, tlfx(i)) ! CABLE's gamma_star

               ! apply same scaling for k as for Vcmax in C4 plants
               if (.not. cable_user%explicit_gm) then
                  kc4(i,1) = veg%c4kci(i) * vcmxt4(i,1)/veg%vcmax_sun(i)
                  kc4(i,2) = veg%c4kci(i) * vcmxt4(i,2)/veg%vcmax_shade(i)
               else
                  kc4(i,1) = veg%c4kcc(i) * vcmxt4(i,1)/veg%vcmax_sun(i)
                  kc4(i,2) = veg%c4kcc(i) * vcmxt4(i,2)/veg%vcmax_shade(i)
               endif

               ! Leuning 2002 (PCE) equation for temperature response
               ! used for Jmax for C3 plants:
               if (.not.cable_user%acclimate_photosyn) then
                  temp_sun_c3(i) = xejmxt3(tlfx(i)) * &
                     veg%ejmax_sun(i) * (1.0-veg%frac4(i))
                  temp_shade_c3(i) = xejmxt3(tlfx(i)) * &
                     veg%ejmax_shade(i) * (1.0-veg%frac4(i))
               else
                  call xejmxt3_acclim(tlfx(i), climate%mtemp(i), climate%mtemp_max20(i), temp_c3(i))
                  temp_sun_c3(i)   = temp_c3(i) * veg%ejmax_sun(i) * (1.0-veg%frac4(i))
                  temp_shade_c3(i) = temp_c3(i) * veg%ejmax_shade(i) * (1.0-veg%frac4(i))
               endif
               ejmxt3(i,1) = rad%scalex(i,1) * temp_sun_c3(i) * fwsoil(i)**qb
               ejmxt3(i,2) = rad%scalex(i,2) * temp_shade_c3(i) * fwsoil(i)**qb

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
               ! Store dsx zihanlu
               dsxx(i) = dsx(i)
               ! "d_{3}" in Wang and Leuning, 1998, appendix E:
               cx1(i) = conkct(i) * (1.0+0.21/conkot(i))
               !cx2(i) = 2.0 * C%gam0 * ( 1.0 + C%gam1 * tdiff(i) &
               !                             + C%gam2 * tdiff(i) * tdiff(i) )
               cx2(i) = 2.0 * gam0 * EXP( ( egam / (C%rgas*C%trefk) ) &
                  * ( 1.0 - C%trefk/tlfx(i) ) )

               ! All equations below in appendix E in Wang and Leuning 1998 are
               ! for calculating anx, csx and gswx for Rubisco limited,
               ! RuBP limited, sink limited.
               ! JK: vx4 changed to correspond to formulation in Collatz et al. 1992
               temp2(i,:) = rad%qcan(i,:,1) * jtomol * (1.0-veg%frac4(i))
               vx3(i,:)   = ej3x(temp2(i,:), veg%alpha(i), veg%convex(i), ejmxt3(i,:))
               vx4(i,:)   = veg%alpha(i) * rad%qcan(i,:,1) * jtomol * veg%frac4(i) * fwsoil(i)**qb
               !vx4(i,:)   = ej4x(temp2(i,:), veg%alpha(i), veg%convex(i), vcmxt4(i,:))
               rdx(i,:)   = veg%cfrd(i)*Vcmxt3(i,:) + veg%cfrd(i)*vcmxt4(i,:)

               !Vanessa - the trunk does not contain xleauning as of Ticket#56 inclusion
               !as well as other inconsistencies here that need further investigation. In the
               !interests of getting this into the trunk ASAP just isolate this code for now
               !default side of this condition is to use trunk version

               !#ifdef VanessasCanopy


               if (cable_user%CALL_climate) then
                  ! if (veg%iveg(i).eq.1) then
                  !    temp2(i,1) = rad%qcan(i,1,1) * jtomol * (1.0-veg%frac4(i))
                  !    temp2(i,2) = rad%qcan(i,2,1) * jtomol * (1.0-veg%frac4(i))
                  !    vx3(i,1)  = ej3x(temp2(i,1),climate%frec(i)*veg%alpha(i), &
                  !         veg%convex(i),ejmxt3(i,1))
                  !    vx3(i,2)  = ej3x(temp2(i,2),climate%frec(i)*veg%alpha(i), &
                  !         veg%convex(i),ejmxt3(i,2))
                  !    temp(i) =  xvcmxt3(tlfx(i)) * veg%vcmax(i) * (1.0-veg%frac4(i))
                  !    vcmxt3(i,1) = rad%scalex(i,1) * temp(i) * climate%frec(i)
                  !    vcmxt3(i,2) = rad%scalex(i,2) * temp(i) * climate%frec(i)
                  ! endif

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

                  if (veg%iveg(i) == 2) then ! evergreen broadleaf forest

                     rdx(i,1) = 1.0*(1.2818e-6+0.0116*veg%vcmax(i)- &
                        0.0334*climate%qtemp_max_last_year(i)*1e-6)

                     rdx(i,2) = rdx(i,1)

                  elseif (veg%iveg(i) == 4) then ! decid broadleaf forest

                     rdx(i,1) = 1.0*(1.2818e-6+0.0116*veg%vcmax(i)- &
                        0.0334*climate%qtemp_max_last_year(i)*1e-6)

                     rdx(i,2) = rdx(i,1)

                  elseif (veg%iveg(i) == 1) then ! evergreen needleleaf forest

                     rdx(i,1) = 1.0*(1.2877e-6+0.0116*veg%vcmax(i)- &
                        0.0334*climate%qtemp_max_last_year(i)*1e-6)

                     rdx(i,2) = rdx(i,1)

                  elseif (veg%iveg(i) == 3) then ! decid needleleaf forest

                     rdx(i,1) = 1.0*(1.2877e-6+0.0116*veg%vcmax(i)- &
                        0.0334*climate%qtemp_max_last_year(i)*1e-6)

                     rdx(i,2) = rdx(i,1)

                  elseif ((veg%iveg(i) == 6) .or. (veg%iveg(i) == 8) .or. &
                     (veg%iveg(i) == 9)) then ! C3 grass, tundra, crop

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
                  rdx(i,1) = rdx(i,1) * xrdt(tlfx(i)) * rad%scalex(i,1) * fwsoil(i)**qb
                  rdx(i,2) = rdx(i,2) * xrdt(tlfx(i)) * rad%scalex(i,2) * fwsoil(i)**qb

                  ! reduction of daytime leaf dark-respiration to account for
                  !photo-inhibition
                  rdx(i,1) = rdx(i,1) * light_inhibition(rad%qcan(i,1,1)*jtomol*1.0e6)
                  rdx(i,2) = rdx(i,2) * light_inhibition(rad%qcan(i,2,1)*jtomol*1.0e6)

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

               ! calculate canopy-level mesophyll conductance
               if (cable_user%explicit_gm) then
                  gmes(i,1) = veg%gm(i) * rad%scalex(i,1) * xgmesT(tlfx(i)) * max(0.15,fwsoil(i)**qm)
                  gmes(i,2) = veg%gm(i) * rad%scalex(i,2) * xgmesT(tlfx(i)) * max(0.15,fwsoil(i)**qm)

                  !gmes(i,1) = gmes(i,1) * MAX(0.15_r_2,real(fwsoil(i)**qm,r_2))
                  !gmes(i,2) = gmes(i,2) * MAX(0.15_r_2,real(fwsoil(i)**qm,r_2))
               else ! gmes = 0.0 easier to debug than inf
                  gmes(i,1) = 0.0
                  gmes(i,2) = 0.0
               endif



            ENDIF !IF (canopy%vlaiw(i) > C%LAI_THRESH .AND. abs_deltlf(i) > 0.1)

         ENDDO !i=1,mp
         ! gmes is 0.0 if explicit_gm = FALSE (easier to debug)
         IF (cable_user%GS_SWITCH /= 'profitmax') THEN
            ! CALL photosynthesis_gm( csx(:,:), &
            !    spread(cx1(:),2,mf), &
            !    spread(cx2(:),2,mf), &
            !    gswmin(:,:), rdx(:,:), vcmxt3(:,:), &
            !    vcmxt4(:,:), vx3(:,:), vx4(:,:), &
            ! ! Ticket #56, xleuning replaced with gs_coeff here
            !    gs_coeff(:,:), rad%fvlai(:,:), &
            !    spread(abs_deltlf,2,mf), abs_deltpsil, &
            !    anx(:,:), fwsoil(:), qs, gmes(:,:), kc4(:,:), &
            !    anrubiscox(:,:), anrubpx(:,:), ansinkx(:,:), eta_x(:,:), dAnx(:,:) )
         CALL photosynthesis_gm_givengs(csx(:,:),gspsdo(:,:), &
              spread(cx1(:),2,mf), &
               spread(cx2(:),2,mf), &
               rdx(:,:), vcmxt3(:,:), &
               vcmxt4(:,:), vx3(:,:), vx4(:,:), &
               rad%fvlai(:,:), spread(abs_deltlf,2,mf), &
               anx(:,:),gmes(:,:), kc4(:,:), &
               anrubiscox(:,:), anrubpx(:,:), ansinkx(:,:) )

         ENDIF

         ! print*, 'DD28 ', rad%fvlai
         ! print*, 'DD29 ', met%ca
         ! print*, 'DD30 ', canopy%gswx
         ! print*, 'DD31 ', air%dsatdk
         ! print*, 'DD32 ', ssnow%wbice
         ! print*, 'DD33 ', ssnow%rex
         ! print*, 'DD34 ', veg%froot
         ! print*, 'DD35 ', soil%ssat
         ! print*, 'DD36 ', soil%swilt
         ! print*, 'DD37 ', canopy%fevc
         ! print*, 'DD38 ', veg%gamma
         ! print*, 'DD39 ', air%rlam
         ! print*, 'DD40 ', veg%zr
         ! print*, 'DD41 ', soil%zse
         ! print*, 'DD42 ', ssnow%rex
         ! print*, 'DD43 ', ssnow%evapfbl
         DO i=1,mp


            IF (canopy%vlaiw(i) > C%LAI_THRESH .AND. (abs_deltlf(i) > 0.1)) Then
               DO kk=1, mf

                  IF (rad%fvlai(i,kk)>C%LAI_THRESH) THEN
                     !  write(logn,*) '2: gs_coeff of ',kk, ': ',gs_coeff(i,kk)


                     csx(i,kk) = real(met%ca(i),r_2) - real(C%RGBWC*anx(i,kk),r_2) / &
                        (gbhu(i,kk) + gbhf(i,kk))
                     csx(i,kk) = max(1.0e-4_r_2, csx(i,kk))

                     ! write(logn,*) '3: gsw of ',kk, ': ', canopy%gswx(i,kk)
                     !Recalculate conductance for water:
                     gw(i,kk) = 1.0 / ( 1.0 / (gspsdo(i,kk) * C%RGSWC) + &
                        1.0 / ( 1.075 * real(gbhu(i,kk) + gbhf(i,kk)) ) )
                     gw(i,kk) = max(gw(i,kk), 0.00001)

                     ! Modified psychrometric constant
                     ! (Monteith and Unsworth, 1990)
                     psycst(i,kk) = air%psyc(i) * REAL( ghr(i,kk) / gw(i,kk) )

                  ENDIF

               ENDDO
               ex(i,1) = real( ( air%dsatdk(i) * ( rad%rniso(i,1) - C%capp * C%rmair &
               * ( met%tvair(i) - met%tk(i) ) * rad%gradis(i,1) ) &
               + C%capp * C%rmair * vpdtmp(i) * ghr(i,1) ) &
               / ( air%dsatdk(i) + psycst(i,1) ),r_2)
               ex(i,2) = real( ( air%dsatdk(i) &
               * ( rad%rniso(i,2) - C%capp * C%rmair * ( met%tvair(i) - &
               met%tk(i) ) * rad%gradis(i,2) ) + C%capp * C%rmair * &
               vpdtmp(i) * ghr(i,2) ) / &
               ( air%dsatdk(i) + psycst(i,2) ), r_2)  
               ecx(i) = sum(ex(i,:))
               ex(i,:) = ex(i,:) * (1.0_r_2-real(canopy%fwet(i), r_2)) / real(air%rlam(i), r_2) 
               ! ecx(i) = real( ( air%dsatdk(i) * ( rad%rniso(i,1) - C%capp * C%rmair &
               !    * ( met%tvair(i) - met%tk(i) ) * rad%gradis(i,1) ) &
               !    + C%capp * C%rmair * met%dva(i) * ghr(i,1) ) &
               !    / ( air%dsatdk(i) + psycst(i,1) ) + ( air%dsatdk(i) &
               !    * ( rad%rniso(i,2) - C%capp * C%rmair * ( met%tvair(i) - &
               !    met%tk(i) ) * rad%gradis(i,2) ) + C%capp * C%rmair * &
               !    met%dva(i) * ghr(i,2) ) / &
               !    ( air%dsatdk(i) + psycst(i,2) ), r_2)
               !!!!!!!!!!!!!!! I don't think here we need use soil water to limit the fevc or ecx -zihanlu !!!!!!!!!!
   
               ! IF (cable_user%SOIL_SCHE == 'Haverd2013' .or. (INDEX(cable_user%FWSOIL_SWITCH, 'Haverd2013') > 0)) then
               !    ! avoid root-water extraction when fwsoil is zero
               !    if (fwsoil(i) < 1e-6) then
               !       anx(i,:) = -rdx(i,:)
               !       ecx(i)   = 0.0_r_2
               !    endif

               !    canopy%fevc(i) = ecx(i) * (1.0_r_2-real(canopy%fwet(i),r_2))

               !    call getrex_1d(wbtmp(i,:)-real(ssnow%wbice(i,:),r_2), ssnow%rex(i,:), &
               !       canopy%fwsoil(i), &
               !       real(veg%froot(i,:),r_2), SPREAD(real(soil%ssat(i),r_2),1,ms), &
               !       SPREAD(real(soil%swilt(i),r_2),1,ms), &
               !       max(canopy%fevc(i)/real(air%rlam(i),r_2)/1000.0_r_2, 0.0_r_2), &
               !       veg%gamma(i), &
               !       real(soil%zse, r_2), real(dels,r_2), veg%zr(i))
               !    canopy%fwsoiltmp(i) = real(canopy%fwsoil(i))
               !    where (ssnow%rex(i,:) > tiny(1.0_r_2)) &
               !    ssnow%evapfbl(i,:) = real(ssnow%rex(i,:))*dels*1000. ! mm water &
               !    IF (INDEX(cable_user%FWSOIL_SWITCH, 'Haverd2013') > 0) then
               !       fwsoil(i) = real(canopy%fwsoil(i))
               !       !(root water extraction) per time step
   
               !       if (cable_user%Cumberland_soil) then
               !          canopy%fwsoil(i) = max(canopy%fwsoil(i), 0.6_r_2)
               !          fwsoil(i) = real(canopy%fwsoil(i))
               !       endif
               !    ENDIF 


               ! ELSE IF (cable_user%FWSOIL_SWITCH == 'profitmax' .AND. cable_user%SOIL_SCHE == 'hydraulics') THEN


               !    IF (ecx(i) > 0.0 .AND. canopy%fwet(i) < 1.0) THEN
               !       evapfb(i) = ( 1.0 - canopy%fwet(i)) * REAL( ecx(i) ) *dels  &
               !          / air%rlam(i)

               !       DO kk = 1,ms

               !          !ssnow%evapfbl(i,kk) = MIN(evapfb(i) * &
               !          !                          ssnow%fraction_uptake(i,kk),  &
               !          !                          MAX(0.0, &
               !          !                              REAL(ssnow%wb(i,kk)) -    &
               !          !                              soil%swilt(i)) * &
               !          !                          soil%zse(kk) * 1000.0)

               !          ! ms8355: no bounding by swilt in this version, or the
               !          ! "beneits" from hydraulics are cancelled
               !          ssnow%evapfbl(i,kk) = evapfb(i) * &
               !             ssnow%fraction_uptake(i,kk)



               !       ENDDO
               !       canopy%fevc(i) = SUM(ssnow%evapfbl(i,:))*air%rlam(i)/dels

               !       ecx(i) = canopy%fevc(i) / (1.0-canopy%fwet(i))
               !    ENDIF

               ! ENDIF

               ! IF (cable_user%SOIL_SCHE .ne.'Haverd2013') THEN

               !    if (ecx(i) > 0.0_r_2 .and. canopy%fwet(i) < 1.0) then
               !       evapfb(i) = ( 1.0 - canopy%fwet(i)) * real(ecx(i)) *dels &
               !          / air%rlam(i)

               !       DO kk = 1,ms

               !          ssnow%evapfbl(i,kk) = MIN( evapfb(i) * veg%froot(i,kk), &
               !             MAX( 0.0, REAL( wbtmp(i,kk) ) - &
               !             1.1 * soil%swilt(i) ) * &
               !             soil%zse(kk) * 1000.0 )

               !       ENDDO
               !       IF (cable_user%soil_struc=='default') then
               !          canopy%fevc(i) = real(sum(ssnow%evapfbl(i,:))*air%rlam(i)/dels, r_2)
               !          ecx(i) = canopy%fevc(i) / (1.0_r_2-real(canopy%fwet(i),r_2))
               !       ELSEIF (cable_user%soil_struc=='sli') then
               !          canopy%fevc(i) = ecx(i) * (1.0_r_2-real(canopy%fwet(i),r_2))
               !       ENDIF

               !    ENDIF

               ! ENDIF
               ! Update canopy sensible heat flux:
               hcx(i) = ( SUM(real(rad%rniso(i,:),r_2)) - ecx(i) &
                  - real(C%capp*C%rmair*(met%tvair(i)-met%tk(i)), r_2) &
                  * SUM(real(rad%gradis(i,:),r_2)) ) &
                  * SUM(real(gh(i,:),r_2))/ SUM(real(ghr(i,:),r_2))

               ! Update leaf temperature:
               tlfx(i) = met%tvair(i)+REAL(hcx(i))/(C%capp*C%rmair*SUM(gh(i,:)))

               ! Update net radiation for canopy:
               rnx(i) = real( SUM( rad%rniso(i,:)) - &
                  C%CAPP * C%rmair *( tlfx(i)-met%tk(i) ) * &
                  SUM( rad%gradis(i,:) ), r_2)

               ! Update leaf surface vapour pressure deficit:
               dsx(i) = vpdtmp(i) + air%dsatdk(i) * (tlfx(i)-met%tvair(i))
               if (cable_user%perturb_dva_by_T) then
                  dsx(i) = vpdtmp(i) + air%dsatdk(i) * (tlfx(i) -met%tvair(i) + cable_user%dva_T_perturbation )
               endif

               dsx(i)=  max(dsx(i),0.0)

               ! Store change in leaf temperature between successive iterations:
               deltlf(i) = tlfxx(i)-tlfx(i)
               
  
               abs_deltlf(i) = ABS(deltlf(i))
               abs_deltds(i) = ABS(dsxx(i)-dsx(i))
               abs_deltcs(i,:) = ABS(csxx(i,:)-csx(i,:))
            ENDIF !lai/abs_deltlf

         ENDDO !i=1,mp
         ! if (ktau>=5184) then
         ! print*, 'check after getrex_1d: ktau & k= ',ktau,k
         ! end if
         ! Where leaf temp change b/w iterations is significant, and
         ! difference is smaller than the previous iteration, store results:
         DO i=1,mp

            IF ( abs_deltlf(i) < ABS( deltlfy(i) ) ) THEN

               deltlfy(i)       = deltlf(i)
               tlfy(i)          = tlfx(i)
               dsy(i) = dsx(i)
               !print*,'when dT<last dT, tlfy: ', tlfy(i)
               ! psily(i,:)       = psilx(i,:)
               rny(i)           = rnx(i)
               hcy(i)           = hcx(i)
               ecy(i)           = ecx(i)
               rdy(i,:)         = rdx(i,:)
               an_y(i,:)        = anx(i,:)
               anrubiscoy(i,:)  = anrubiscox(i,:)
               anrubpy(i,:)     = anrubpx(i,:)
               ! save last values calculated for ssnow%evapfbl
               ! oldevapfbl(i,:) = ssnow%evapfbl(i,:)
               gswy(i,:)        = canopy%gswx(i,:)

               csy(i,:) = csx(i,:)
              
               vcmxt3y(i,:) = vcmxt3(i,:)
               vcmxt4y(i,:) = vcmxt4(i,:)
               cx2y(i) = cx2(i)
            ENDIF
            ! if (ktau>=5184) then
            ! print*, 'check y=x ktau & k= ',ktau,k
            ! end if
            ! if (cable_user%GS_SWITCH == 'tuzet' .AND. &
            !       INDEX(cable_user%FWSOIL_SWITCH,'LWP') > 0 .AND. &
            !       Any(abs_deltpsil(i,:) > 0.5) ) then
            !    ! after 4 iterations, take mean of current & previous estimates
            !    ! as the next estimate of leaf temperature, to avoid oscillation
            !    psilx(i,1) = ( 0.5 * ( MAX( 0, k-5 ) / ( k - 4.9999 ) ) ) *psilxx(i,1) + &
            !       ( 1.0 - ( 0.5 * ( MAX( 0, k-5 ) / ( k - 4.9999 ) ) ) ) &
            !       * psilx(i,1)
            !    psilx(i,2) = ( 0.5 * ( MAX( 0, k-5 ) / ( k - 4.9999 ) ) ) *psilxx(i,2) + &
            !       ( 1.0 - ( 0.5 * ( MAX( 0, k-5 ) / ( k - 4.9999 ) ) ) ) &
            !       * psilx(i,2)
            ! endif
            ! if (ktau>=5184) then
            ! print*, 'check tlfx ktau & k= ',ktau,k
            ! end if
            IF (k==1) THEN
               ! take the first iterated estimates as the defaults
               tlfy(i) = tlfx(i)
               dsy(i) = dsx(i)
               !print*,'when k=1, tlfy: ', tlfy(i)
               ! psily(i,:) = psilx(i,:)
               rny(i) = rnx(i)
               hcy(i) = hcx(i)
               ecy(i) = ecx(i)
               rdy(i,:) = rdx(i,:)
               an_y(i,:) = anx(i,:)
               anrubiscoy(i,:) = anrubiscox(i,:)
               anrubpy(i,:) = anrubpx(i,:)

               ! save last values calculated for ssnow%evapfbl
               ! oldevapfbl(i,:) = ssnow%evapfbl(i,:)
               gswy(i,:)       = canopy%gswx(i,:)
         
               csy(i,:) = csx(i,:)
   
               vcmxt3y(i,:) = vcmxt3(i,:)
               vcmxt4y(i,:) = vcmxt4(i,:)
               cx2y(i) = cx2(i)
            END IF
            !print*, 'check after k==1 ',ktau,k
            !if (ktau_tot>=nktau .and. ktau_tot<=(nktau+NN-1)) then
            if (any(allktau == ktau_tot)) then
          
                  write(138,*) ktau_tot, iter, i, k, tlfx(i), deltlf(i), &
                  dsx(i),abs_deltds(i), 0, 0,0,0,0,0, &
                  0, csx(i,1), csx(i,2),abs_deltcs(i,1), abs_deltcs(i,2),anx(i,1), anx(i,2),anrubiscox(i,1), &
                  anrubiscox(i,2), anrubpx(i,1),anrubpx(i,2),ansinkx(i,1),ansinkx(i,2), &
                  canopy%gswx(i,1), canopy%gswx(i,2),canopy%gswx(i,1), canopy%gswx(i,2), &
                  vcmxt3(i,1),vcmxt3(i,2),gs_coeff(i,1),gs_coeff(i,2),rdx(i,1),rdx(i,2),ex(i,1),ex(i,2), &
                  ssnow%rootR(i,1),ssnow%rootR(i,2),ssnow%rootR(i,3),ssnow%rootR(i,4), &
                  ssnow%soilR(i,1),ssnow%soilR(i,2),ssnow%soilR(i,3),ssnow%soilR(i,4), &
                  kplantx(i)

         

            END IF
            ! if (ktau>=5184) then
            ! print*, 'write 134 ',ktau,k
            ! end if
            if ( abs_deltlf(i) > 0.1 .AND. k > 5 .AND. k < C%MAXITER ) then
               ! after 4 iterations, take mean of current & previous estimates
               ! as the next estimate of leaf temperature, to avoid oscillation
               dc = veg%dc(i)
               ! tlfx(i) = ( 0.5 * ( MAX( 0, k-5 ) / ( k - 4.9999 ) ) ) *tlfxx(i) + &
               !    ( 1.0 - ( 0.5 * ( MAX( 0, k-5 ) / ( k - 4.9999 ) ) ) ) &
               !    * tlfx(i)
               tlfx(i) = dc *tlfxx(i) + ( 1.0 - dc ) * tlfx(i)
               dsx(i) = dc *dsxx(i) + ( 1.0 - dc ) * dsx(i)
               csx(i,1) = dc *csxx(i,1) + ( 1.0 - dc ) * csx(i,1)
               csx(i,2) = dc *csxx(i,2) + ( 1.0 - dc ) * csx(i,2)

            endif

         END DO !over mp
         !if (k == C%MAXITER .AND. ktau>=5760 .AND. ktau<=11664) then
      ! if (k == C%MAXITER .AND. ktau_tot>=nktau .and. ktau_tot<=(nktau+NN-1)) then
      ! !if (k == C%MAXITER .AND. ktau>=5760 .AND. ktau<=11664) then
      !        print*, 'value x:',anrubiscox(:,1),anx(:,1)
      ! endif
      END DO  ! DO WHILE (ANY(abs_deltlf > 0.1) .AND.  k < C%MAXITER)
      !print*,'when k end, tlfy: ', tlfy(1)
      ! if (ktau_tot>=nktau .and. ktau_tot<=(nktau+NN-1)) then
      ! i = 1
      ! write(134,*) ktau, iter, i, 21, tlfy(i), deltlfy(i), &
      ! dsx(i), psily(i,1), psily(i,2),fwpsiy(:,1),fwpsiy(:,2),canopy%fwpsi(i,1),canopy%fwpsi(i,2),&
      ! csy(i,1), csy(i,2),csx(i,1), csx(i,2), an_y(i,1), an_y(i,2), &
      ! gswy(i,1), gswy(i,2),canopy%gswx(i,1), canopy%gswx(i,2), &
      ! vcmxt3(i,1),vcmxt3(i,2), gs_coeff(i,1),gs_coeff(i,2),rdy(i,1),rdy(i,2)
      ! end if 
      !if (ktau_tot==(nktau+NN-1)  .and. iter==4) THEN
      !if (any(nktau_end == ktau_tot) .and. iter==4) THEN
      if (ktau_tot == nktau_end(size(nktau_end)) .and. iter==4) THEN

         close(138)
      
      END IF
      deallocate(nktau) 
      deallocate(nktau_end) 
      deallocate(allktau) 

      
      ! print*, 'End k loop: ktau & iter & k= ',ktau,iter,k
      ! dry canopy flux
      canopy%epotvpd = (1.0_r_2-real(canopy%fwet,r_2)) * ecy
      canopy%abs_deltcs_vpd = abs_deltcs * 1.0e6_r_2
      canopy%abs_deltlf_vpd = abs_deltlf
      canopy%abs_deltds_vpd = abs_deltds
      IF (.NOT. PRESENT(gspsdo)) THEN
         canopy%abs_deltcs = abs_deltcs * 1.0e6_r_2
         canopy%abs_deltlf = abs_deltlf
         canopy%abs_deltds = abs_deltds


         !write(logn,*) 'fevc', canopy%fevc(1)
         ! print*, 'DD45 ', canopy%fwet
         ! print*, 'DD46 ', canopy%fevc
           !!!!!!!!!!!!!!!! since we don't have the evapfbl calculation before, we don't need that anymore -zihanLu !!!!!!!!!!!!!!!!!!
         ! IF (cable_user%fwsoil_switch.NE.'Haverd2013' .AND. cable_user%fwsoil_switch.NE.'profitmax') THEN

            ! Recalculate ssnow%evapfbl as ecy may not be updated with the ecx
            ! calculated in the last iteration.
            ! DO NOT use simple scaling as there are times that ssnow%evapfbl is zero.
            ! ** ssnow%evapfbl(i,:) = ssnow%evapfbl(i,:) * ecy(i) / ecx(i) **

         !    DO i = 1, mp

         !       IF( ecy(i) > 0.0_r_2 .AND. canopy%fwet(i) < 1.0 ) THEN

         !          IF( ABS( ecy(i) - ecx(i) ) > 1.0e-6_r_2 ) THEN

         !             IF( ABS( canopy%fevc(i) - real(SUM(oldevapfbl(i,:)) * air%rlam(i) &
         !                /dels, r_2) ) > 1.0e-4_r_2 ) THEN

         !                WRITE(*,*) 'ERROR! oldevapfbl not right.', ktau_gl, i
         !                WRITE(*,*) 'ecx, ecy = ', ecx(i), ecy(i)
         !                WRITE(*,*) 'or in mm = ', ecx(i) * ( 1.0 - canopy%fwet(i) ) &
         !                   / air%rlam(i) * dels, &
         !                   ecy(i) * ( 1.0 - canopy%fwet(i) ) / &
         !                   air%rlam(i) * dels

         !                WRITE(*,*)'fevc = ', canopy%fevc(i), SUM( oldevapfbl(i,:) ) * &
         !                   air%rlam(i) / dels
         !                WRITE(*,*) 'fwet = ', canopy%fwet(i)
         !                WRITE(*,*) 'oldevapfbl = ', oldevapfbl(i,:)
         !                WRITE(*,*) 'ssnow%evapfbl before rescaling: ', &
         !                   ssnow%evapfbl(i,:)
         !                ! STOP

         !             ELSE

         !                ssnow%evapfbl(i,:) = oldevapfbl(i,:)

         !             END IF

         !          END IF

         !       END IF

         !    END DO

         ! ENDIF
         ! print*, 'Recalculate ssnow%evapfbl: ktau: ',ktau
         ! canopy%psi_can = real(psily, r_2)
         canopy%frday = 12.0 * SUM(rdy, 2)
         !! vh !! inserted min to avoid -ve values of GPP
         canopy%fpn = min(-12.0 * SUM(an_y, 2), canopy%frday)
         ! change canopy%gswx from gswx to gswy, zihanlu, 19/12/2024
         canopy%gswx = gswy
         canopy%fwpsi = fwpsiy
         ! additional diagnostic variables for assessing contributions of rubisco and rubp limited photosynthesis to
         ! net photosynthesis in sunlit and shaded leaves.
         canopy%A_sh = real(an_y(:,2), r_2)
         canopy%A_sl = real(an_y(:,1), r_2)

         canopy%GPP_sh = real(an_y(:,2) + rdy(:,2), r_2)
         canopy%GPP_sl = real(an_y(:,1) + rdy(:,1), r_2)

         canopy%A_shC = real(anrubiscoy(:,2), r_2)
         canopy%A_shJ = real(anrubpy(:,2), r_2)

         where (anrubiscoy(:,2) > an_y(:,2)) canopy%A_shC = 0.0_r_2
         where (anrubpy(:,2) > an_y(:,2))    canopy%A_shJ = 0.0_r_2
         
         canopy%A_slC = real(anrubiscoy(:,1), r_2)
         canopy%A_slJ = real(anrubpy(:,1), r_2)
         !if (ktau>=5760 .AND. ktau<=11664) then
         ! if (ktau_tot>=nktau .and. ktau_tot<=(nktau+NN-1)) then
         !    print*, 'value y:',canopy%A_slC, canopy%A_sl
         ! endif
         where (anrubiscoy(:,1) > an_y(:,1)) canopy%A_slC = 0.0_r_2
         where (anrubpy(:,1)    > an_y(:,1)) canopy%A_slJ = 0.0_r_2

         canopy%eta_A_cs = canopy%A_sh * min(eta_y(:,2),5.0_r_2) + canopy%A_sl * min(eta_y(:,1),5.0_r_2)
         where ((canopy%A_sl > 0.0_r_2) .and. (canopy%A_sh > 0.0_r_2))
            canopy%eta_GPP_cs = canopy%GPP_sh  * min(eta_y(:,2),5.0_r_2) * canopy%GPP_sh/canopy%A_sh + &
               canopy%GPP_sl * min(eta_y(:,1),5.0_r_2)* canopy%GPP_sl/canopy%A_sl
         elsewhere ((canopy%A_sl .le. 0.0_r_2) .and. (canopy%A_sh > 0.0_r_2))
            canopy%eta_GPP_cs = canopy%GPP_sh * min(eta_y(:,2),5.0_r_2) * canopy%GPP_sh/canopy%A_sh
         elsewhere
            canopy%eta_GPP_cs = canopy%eta_A_cs
         end where
         ! change from gs_coeff to gs_coeffy
         where (canopy%A_sl > 0.0_r_2)
            canopy%eta_A_cs_sl    =  min(eta_y(:,1),5.0_r_2)
            canopy%eta_fevc_cs_sl = (min(eta_y(:,1),5.0_r_2) - 1.0_r_2) * &
               real(max(0.0, gs_coeffy(:,1)*an_y(:,1)), r_2) / real(canopy%gswx(:,1), r_2)
         elsewhere
            canopy%eta_A_cs_sl    = 0.0_r_2
            canopy%eta_fevc_cs_sl = 0.0_r_2
         endwhere
         ! change from gs_coeff to gs_coeffy
         where (canopy%A_sh > 0.0_r_2)
            canopy%eta_A_cs_sh    =  min(eta_y(:,2),5.0_r_2)
            canopy%eta_fevc_cs_sh = (min(eta_y(:,2),5.0_r_2) - 1.0_r_2) * &
               real(max(0.0, gs_coeffy(:,2)*an_y(:,2)), r_2) / real(canopy%gswx(:,2), r_2)
         elsewhere
            canopy%eta_A_cs_sh    = 0.0_r_2
            canopy%eta_fevc_cs_sh = 0.0_r_2
         endwhere

         where ((rad%fvlai(:,1)+rad%fvlai(:,2)) > 0.01)
            canopy%eta_fevc_cs = ( canopy%eta_fevc_cs_sl *  real(rad%fvlai(:,1), r_2) + &
               canopy%eta_fevc_cs_sh * real(rad%fvlai(:,2), r_2) ) * canopy%fevc / &
               real(rad%fvlai(:,1)+rad%fvlai(:,2), r_2)
         elsewhere
            canopy%eta_fevc_cs = 0.0_r_2
         endwhere
         ! canopy%cs_sl = csx(:,1) * 1.0e6_r_2
         ! canopy%cs_sh = csx(:,2) * 1.0e6_r_2
         ! change canopy%cs from csx to csy, zihanlu, 19/12/2024
         canopy%cs_sl = csy(:,1) * 1.0e6_r_2
         canopy%cs_sh = csy(:,2) * 1.0e6_r_2

         canopy%dAdcs = canopy%A_sl * dAn_y(:,1) + canopy%A_sh * dAn_y(:,2)
         canopy%cs    = canopy%A_sl * canopy%cs_sl * 1.0e6_r_2 + canopy%A_sh * canopy%cs_sh * 1.0e6_r_2

         canopy%tlf   = real(tlfy, r_2)
         ! change canopy%dlf from dsx to dsy, zihanlu, 19/12/2024
         !canopy%dlf   = real(dsx, r_2)
         canopy%dlf   = real(dsy, r_2)


         ! print*, 'DD47 ', ssnow%evapfbl
         canopy%evapfbl = ssnow%evapfbl

         ! 13C
         canopy%An        = real(an_y, r_2)
         canopy%Rd        = real(rdy, r_2)
         canopy%isc3      = (1.0-veg%frac4) > epsilon(1.0)
         ! change from vcmxt3 to vcmxt3y, 
         canopy%vcmax     = real(merge(vcmxt3y, vcmxt4y, spread(canopy%isc3,2,mf)), r_2)
         ! change from cx2 to cx2y, 
         canopy%gammastar = real(spread(cx2y*0.5,2,mf), r_2)
         canopy%gsc       = real(canopy%gswx / C%rgswc, r_2)
         canopy%gbc       = (gbhu + gbhf) / real(C%rgbwc, r_2)
         canopy%gac       = huge(1.0_r_2)
         ! replace 1.0_r_2/canopy%gac by tiny(1.0_r_2) to avoid underflow
         canopy%ci        = real(spread(met%ca, 2, mf), r_2) - &
            canopy%An / &
            (1.0_r_2 / (1.0_r_2/canopy%gbc + 1.0_r_2/canopy%gsc + tiny(1.0_r_2)))
         IF (cable_user%FWSOIL_SWITCH == 'profitmax') THEN
            canopy%fevcs = (1.0-canopy%fwet) * ecxs  ! trans. from sapflux

            DO i = 1, mp

               new_plc_sat = calc_plc(kcmax(i), veg%kmax(i), veg%PLCcrit(i))

               IF (canopy%psi_can_opt(i) < ssnow%psi_rootzone(i) .AND. ecxs(i) > 1E-9) THEN
                  IF (avg_kcan(i) > 1E-9 .AND. avg_kcan(i) < veg%kmax(i)) THEN

                     ! Partitioning ratios from Drake et al. (2015), PC&E, 38: 1628-1636.
                     new_plc_stem = calc_plc(5. / (2. / avg_kcan(i) + 3. / kcmax(i)), &
                        veg%kmax(i), veg%PLCcrit(i)) ! 2:5 to 3:5

                     new_plc_can = calc_plc(avg_kcan(i), veg%kmax(i), veg%PLCcrit(i))

                  ELSE
                     new_plc_stem = 0.0
                     new_plc_can = 0.0
                  ENDIF

               ELSE
                  new_plc_stem = 0.0
                  new_plc_can = 0.0
               ENDIF

               IF (new_plc_sat > canopy%day_plc_sat(i)) THEN
                  canopy%day_plc_sat(:) = new_plc_sat
               ENDIF

               IF (new_plc_stem > canopy%day_plc_stem(i)) THEN
                  canopy%day_plc_stem(:) = new_plc_stem
               ENDIF

               IF (new_plc_can > canopy%day_plc_can(i)) THEN
                  canopy%day_plc_can(:) = new_plc_can
               ENDIF

               ! set the day's PLC and reset the PLC tracker here
               IF (met%hod(i) >= 24.0 - dels / (60.0 * 60.0) - 1E-6) THEN
                  canopy%day_plc_sat(:) = 0.0
                  canopy%day_plc_stem(:) = 0.0
                  canopy%day_plc_can(:) = 0.0

               ELSE IF (met%hod(i) >= 24.0 - 2.0 * dels / (60.0 * 60.0) - 1E-6) THEN
                  canopy%plc_sat(i) = canopy%day_plc_sat(i)
                  canopy%plc_stem(i) = canopy%day_plc_stem(i)
                  canopy%plc_can(i) = canopy%day_plc_can(i)
               ENDIF

            END DO
         ENDIF
         ! deallocate( gswmin )
      endif ! if not present gspsdo

   END SUBROUTINE dryLeaf_givengs
   ! ------------------------------------------------------------------------------


   ! JK: subroutine photosynthesis_gm now used with and without explicit gm (cable_user%explicit_gm)
   SUBROUTINE photosynthesis_gm( csxz, cx1z, cx2z, gswminz, &
      rdxz, vcmxt3z, vcmxt4z, vx3z, &
      vx4z, gs_coeffz, vlaiz,photo_rule,anxz, fwsoilz, qs, &
      gmes, kc4, anrubiscoz, anrubpz, ansinkz, eta, dA,imp,ktau,k)

    use cable_def_types_mod, only: r_2
    use cable_common_module, only: cable_user

      implicit none

    real(r_2), dimension(:, :), intent(in) :: csxz
    real,      dimension(:, :), intent(in) :: gmes
    real,      dimension(:, :), intent(in) :: &
         cx1z,       & !
         cx2z,       & !
         rdxz,       & !
         vcmxt3z,    & !
         vcmxt4z,    & !
         vx4z,       & !
         vx3z,       & !
         gs_coeffz,  & ! Ticket #56, xleuningz repalced with gs_coeffz
         vlaiz,      & !
         !deltlfz, &
         kc4           !
    real,      dimension(:),    intent(in)    :: fwsoilz
    real,                       intent(in)    :: qs
    real,      dimension(:, :), intent(inout) :: gswminz
    real,      dimension(:, :), intent(inout) :: anxz, anrubiscoz, anrubpz, ansinkz
    real(r_2), dimension(:, :), intent(out)   :: eta, dA
    integer,   intent(in)    :: imp
    integer,   intent(in), optional ::ktau, k
    logical,   intent(in)    :: photo_rule

    ! local variables
    real(r_2), dimension(size(csxz, 1), size(csxz, 2)) :: dAmp, dAme, dAmc, eta_p, eta_e, eta_c
    !real, dimension(mp) :: fwsoilz  ! why was this local??
    real(r_2) :: gamma, beta, gammast, gm, g0, X, Rd, cs
    real(r_2) :: cc, gsm
    real(r_2) :: Am,a,b,ctmp
    integer :: i, j, mp, mf
    ! for minimum of 3 rates and corresponding elasticities
    real,      dimension(3) :: tmp3
    real(r_2), dimension(3) :: dtmp3
    integer,   dimension(1) :: ii

    mp = size(csxz, 1)
    mf = size(csxz, 2)
    do i=imp, imp

       do j=1, mf

          if (vlaiz(i,j) > C%lai_thresh) then

             !if ((deltlfz(i,j) > 0.1 .or. deltpsiz(i,j) > 0.1 .or. deltfwpsiz(i,j) > 0.02) .or. (Numtagz>0)) then
             if (photo_rule) then

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
                if ((vcmxt3z(i,j) > 1.0e-8) .and. (gs_coeffz(i,j) > 100.)) then
                   cs      = csxz(i,j)
                   g0      = real(gswminz(i,j) * fwsoilz(i)**qs / C%RGSWC, r_2)
                   X       = real(gs_coeffz(i,j), r_2)
                   gamma   = real(vcmxt3z(i,j), r_2)
                   beta    = real(cx1z(i,j), r_2)
                   gammast = real(cx2z(i,j) / 2.0, r_2)
                   Rd      = real(rdxz(i,j), r_2)
                   gm      = gmes(i,j)

                   if (trim(cable_user%g0_switch) == 'default') then
                      ! get partial derivative of A wrt cs

                     if (cable_user%explicit_gm) then
                         call fAmdAm_c3(cs, g0, X*cs, gamma, beta, gammast, Rd, &
                              gm, Am, dAmc(i,j))
                     else
                           call fAndAn_c3(cs, g0, X*cs, gamma, beta, gammast, Rd, &
                              Am, dAmc(i,j),a,b,ctmp)
                     endif

                      
                     elseif (trim(cable_user%g0_switch) == 'maximum') then
                        ! set g0 to zero initially
                        if (cable_user%explicit_gm) then
                           call fAmdAm_c3(cs, 0.0_r_2, X*cs, gamma, beta, gammast, Rd, &
                              gm, Am, dAmc(i,j))
                           if (g0 > Am*X) then ! repeat calculation if g0 > A*X
                              call fAmdAm_c3(cs, g0, 0.1e-4_r_2, gamma, beta, gammast, Rd, &
                                 gm, Am, dAmc(i,j))
                           endif
                        else
                           call fAndAn_c3(cs, 0.0_r_2, X*cs, gamma, beta, gammast, Rd, &
                              Am, dAmc(i,j),a,b,ctmp)
                           if (present(ktau)) then
                              if (ktau == 192213) then
                                 print *, '1.  Rubusco limited ----------------'
                                 print *, '  ktau,k,a,b,c: ', ktau, k, a, b, ctmp
                                 print *, '  An: ', Am
                           call fAn_c3(a, b, ctmp, Amtmp)
                           print*, a, b, ctmp, Amtmp
                           print*, b**2 - 4.0_r_2 * a * ctmp, (-b-sqrt(b**2 - 4.0_r_2 * a * ctmp)), 2.0_r_2*a, (-b-sqrt(b**2 - 4.0_r_2 * a * ctmp))/(2.0_r_2*a)

                                 print *, 'cs, X, gamma, beta, gammast, Rd, ', cs, X, gamma, beta, gammast, Rd
                                 print*, '1-x', 1-X*cs
                                 print*,'(1-x)*(Rd-gamma)',(1-X*cs)*(Rd-gamma)
                                 print*,'(1-x)*(Rd-gamma)*cs',(1-X*cs)*(Rd-gamma)*cs
                                 print*,'-(gamma*gammast+rd*beta)',-(gamma*gammast+Rd*beta)
                                 print*,'-x(gamma*gammast+rd*beta)',-X*cs*(gamma*gammast+Rd*beta)
                              end if
                           end if
                           if (g0 > Am*X) then ! repeat calculation if g0 > A*X
                              call fAndAn_c3(cs, g0, 0.1e-12_r_2, gamma, beta, gammast, Rd, &
                                 Am, dAmc(i,j),a,b,ctmp)
                           if (present(ktau)) then
                              if (ktau == 192213) then
                                 print *, '1.1  g0 > Am*X ----------------'
                                 print *, '  ktau,k,a,b,c: ', ktau, k, a, b, ctmp
                                 print *, '  An: ', Am
                              end if
                           end if
                           endif
                        endif
                     endif

                     anrubiscoz(i,j) = real(Am)

                     if (cable_user%explicit_gm) then
                        gsm = (gm * (g0+X*Am)) / (gm + (g0+X*Am))
                        cc  = cs - Am / max(gsm, 1.0e-4_r_2)
                     endif
                     if (Am > 0.0_r_2) eta_c(i,j) = dAmc(i,j) * cs / Am
                  endif  ! C3

                  ! C4, Rubisco limited, accounting for explicit mesophyll conductance
                  if (vcmxt4z(i,j) > 1.0e-8) then
                     anrubiscoz(i,j) = vcmxt4z(i,j)- rdxz(i,j)
                     dAmc(i,j)  = 0.0_r_2
                     eta_c(i,j) = 0.0_r_2
                  endif

                  ! C3, RuBP regeneration-limited, accounting for explicit mesophyll conductance
                  if ( (vcmxt3z(i,j) > 0.0) .and. (gs_coeffz(i,j) > 100.) .and. &
                     (vx3z(i,j) > 1.0e-8) ) then
                     cs      = csxz(i,j)
                     g0      = real(gswminz(i,j) * fwsoilz(i)**qs / C%RGSWC, r_2)
                     X       = real(gs_coeffz(i,j), r_2)
                     gamma   = real(vx3z(i,j), r_2)
                     beta    = real(cx2z(i,j), r_2)
                     gammast = real(cx2z(i,j) / 2.0, r_2)
                     Rd      = real(rdxz(i,j), r_2)
                     gm      = gmes(i,j)

                     if (trim(cable_user%g0_switch) == 'default') then
                     
                        if (cable_user%explicit_gm) then
                           call fAmdAm_c3(cs, g0, X*cs, gamma, beta, gammast, Rd, &
                              gm, Am, dAme(i,j))
                        else
                           call fAndAn_c3(cs, g0, X*cs, gamma, beta, gammast, Rd, &
                              Am, dAme(i,j),a,b,ctmp)
                        endif

                     elseif (trim(cable_user%g0_switch) == 'maximum') then
                        if (cable_user%explicit_gm) then
                           call fAmdAm_c3(cs, 0.0_r_2, X*cs, gamma, beta, gammast, Rd, &
                              gm, Am, dAme(i,j))
                           ! repeat calculation if g0 > A*X
                           if (g0 > Am*X) then
                              call fAmdAm_c3(cs, g0, 0.1e-4_r_2, gamma, beta, gammast, Rd, &
                                 gm, Am, dAme(i,j))
                           endif
                        else
                           call fAndAn_c3(cs, 0.0_r_2, X*cs, gamma, beta, gammast, Rd, &
                              Am, dAme(i,j),a,b,ctmp)
                           ! if (present(ktau)) then
                           !    if (ktau == 192213) then
                           !       print *, '2.  Rubp limited ----------------'
                           !       print *, '  ktau,k,a,b,c: ', ktau, k, a, b, ctmp
                           !       print *, '  An: ', Am
                           !    end if
                           ! end if
                           ! repeat calculation if g0 > A*X
                           if (g0 > Am*X) then
                              call fAndAn_c3(cs, g0, 0.1e-4_r_2, gamma, beta, gammast, Rd, &
                                 Am, dAme(i,j),a,b,ctmp)
                           ! if (present(ktau)) then
                           !    if (ktau == 192213) then
                           !       print *, '2.1  g0 > Am*X ----------------'
                           !       print *, '  ktau,k,a,b,c: ', ktau, k, a, b, ctmp
                           !       print *, '  An: ', Am
                           !    end if
                           ! end if
                           endif
                        endif
                     endif
                     anrubpz(i,j) = real(Am)

                     if (cable_user%explicit_gm) then
                        gsm = (gm * (g0+X*Am)) / (gm + (g0+X*Am))
                        cc  = cs - Am / max(gsm, 1.0e-4_r_2)
                     endif
                     if (Am > 0.0_r_2) eta_e(i,j) = dAme(i,j) * cs / Am
                  endif  ! C3

                  ! C4, RuBP limited, accounting for explicit mesophyll conductance
                  if (vx4z(i,j) > 1.0e-8) then
                     anrubpz(i,j) = vx4z(i,j) - rdxz(i,j)
                     dAme(i,j)  = 0.0_r_2
                     eta_e(i,j) = 0.0_r_2
                  endif

                  ! Sink limited, accounting for explicit mesophyll conductance
                  if (vcmxt3z(i,j) > 1.0e-10) then ! C3

                     ansinkz(i,j) = 0.5 * vcmxt3z(i,j) - rdxz(i,j)
                     dAmp(i,j)  = 0.0_r_2
                     eta_p(i,j) = 0.0_r_2

                  elseif ((vcmxt4z(i,j) > 1.0e-10) .and. (gs_coeffz(i,j) > 100.)) then ! C4
                     cs         = csxz(i,j)
                     g0         = real(gswminz(i,j) * fwsoilz(i)**qs / C%RGSWC, r_2)
                     X          = real(gs_coeffz(i,j), r_2)
                     gamma      = real(kc4(i,j),r_2) ! k in Collatz 1992
                     beta       = 0.0_r_2
                     gammast    = 0.0_r_2
                     Rd         = real(rdxz(i,j), r_2)
                     gm         = gmes(i,j)

                     if (trim(cable_user%g0_switch) == 'default') then
                        if (cable_user%explicit_gm) then
                           call fAmdAm_c4(cs, g0, X*cs, gamma, beta, gammast, Rd, gm, &
                              Am, dAmp(i,j))
                        else
                           call fAndAn_c4(cs, g0, X*cs, gamma, beta, gammast, Rd, &
                              Am, dAmp(i,j))
                        endif
                     elseif (trim(cable_user%g0_switch) == 'maximum') then
                        if (cable_user%explicit_gm) then
                           call fAmdAm_c4(cs, 0.0_r_2, X*cs, gamma, beta, gammast, Rd, gm, &
                              Am, dAmp(i,j))
                           if (g0 > Am*X) then
                              call fAmdAm_c4(cs, g0, 0.1e-4_r_2, gamma, beta, gammast, Rd, gm, &
                                 Am, dAmp(i,j))
                           endif
                        else
                           call fAndAn_c4(cs, 0.0_r_2, X*cs, gamma, beta, gammast, Rd, &
                              Am, dAmp(i,j))
                           if (g0 > Am*X) then
                              call fAndAn_c4(cs, g0, 0.1e-4_r_2, gamma, beta, gammast, Rd, &
                                 Am, dAmp(i,j))
                           endif
                        endif
                     endif
                     ansinkz(i,j) = real(Am)

                     dAmp(i,j)  = 0.0_r_2
                     eta_p(i,j) = 0.0_r_2
                     if (Am > 0.0_r_2) eta_p(i,j) = dAmp(i,j) * cs / Am
                  endif ! C3/C4

            endif ! deltlfz > 0.1

               ! minimal of three limited rates
               tmp3      = (/ anrubiscoz(i,j), anrubpz(i,j), ansinkz(i,j) /)
               ii        = minloc(tmp3)
               anxz(i,j) = tmp3(ii(1))
               dtmp3     = (/ dAmc(i,j), dAme(i,j), dAmp(i,j) /)
               dA(i,j)   = dtmp3(ii(1))
               dtmp3     = (/ eta_c(i,j), eta_e(i,j), eta_p(i,j) /)
               eta(i,j)  = dtmp3(ii(1))

         else ! vlaiz(i,j) > C%lai_thresh

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

         endif ! vlaiz(i,j) > C%lai_thresh

      enddo ! j=1, mf

   enddo ! i=1, mp

   END SUBROUTINE photosynthesis_gm

   SUBROUTINE photosynthesis_gm_givengs( csxz, gsz, cx1z, cx2z, &
      rdxz, vcmxt3z, vcmxt4z, vx3z, &
      vx4z, vlaiz, deltlfz, anxz, &
      gmes, kc4, anrubiscoz, anrubpz, ansinkz )

    use cable_def_types_mod, only: r_2
    use cable_common_module, only: cable_user

      implicit none

    real(r_2), dimension(:, :), intent(in) :: csxz, gsz
    real,      dimension(:, :), intent(in) :: gmes
    real,      dimension(:, :), intent(in) :: &
         cx1z,       & !
         cx2z,       & !
         rdxz,       & !
         vcmxt3z,    & !
         vcmxt4z,    & !
         vx4z,       & !
         vx3z,       & !
         vlaiz,      & !
         deltlfz, &
         kc4           !

    real,      dimension(:, :), intent(inout) :: anxz, anrubiscoz, anrubpz, ansinkz

    ! local variables
    real(r_2), dimension(size(csxz, 1), size(csxz, 2)) :: dAmp, dAme, dAmc, eta_p, eta_e, eta_c
    !real, dimension(mp) :: fwsoilz  ! why was this local??
    real(r_2) :: gamma, beta, gammast, gm, g0, X, Rd, cs,gss
    real(r_2) :: Am
    integer :: i, j, mp, mf
    ! for minimum of 3 rates and corresponding elasticities
    real,      dimension(3) :: tmp3
    real(r_2), dimension(3) :: dtmp3
    integer,   dimension(1) :: ii

    mp = size(csxz, 1)
    mf = size(csxz, 2)
    do i=1, mp

      do j=1, mf

         if (vlaiz(i,j) > C%lai_thresh) then

             if (deltlfz(i,j) > 0.1 ) then

                anxz(i,j)       = -rdxz(i,j)
                anrubiscoz(i,j) = -rdxz(i,j)
                anrubpz(i,j)    = -rdxz(i,j)
                ansinkz(i,j)    = -rdxz(i,j)
               !  dAmc(i,j)       = 0.0_r_2
               !  dAme(i,j)       = 0.0_r_2
               !  dAmp(i,j)       = 0.0_r_2
               !  dA(i,j)         = 0.0_r_2
               !  eta_c(i,j)      = 0.0_r_2
               !  eta_e(i,j)      = 0.0_r_2
               !  eta_p(i,j)      = 0.0_r_2
               !  eta(i,j)        = 0.0_r_2

                ! C3, Rubisco limited, accounting for explicit mesophyll conductance
                if ((vcmxt3z(i,j) > 1.0e-8)) then
                   cs      = csxz(i,j)
                   gss = gsz(i,j)
                   gamma   = real(vcmxt3z(i,j), r_2)
                   beta    = real(cx1z(i,j), r_2)
                   gammast = real(cx2z(i,j) / 2.0, r_2)
                   Rd      = real(rdxz(i,j), r_2)
                   gm      = gmes(i,j)

                  if (trim(cable_user%g0_switch) == 'default') then
                      ! get partial derivative of A wrt cs
         
                        call fAndAn_c3_givengs(cs, gss, gamma, beta, gammast, Rd, Am)
                     
                  elseif (trim(cable_user%g0_switch) == 'maximum') then
                        ! ! set g0 to zero initially
                        ! if (cable_user%explicit_gm) then
                        !    call fAmdAm_c3(cs, 0.0_r_2, X*cs, gamma, beta, gammast, Rd, &
                        !       gm, Am, dAmc(i,j))
                        !    if (g0 > Am*X) then ! repeat calculation if g0 > A*X
                        !       call fAmdAm_c3(cs, g0, 0.1e-4_r_2, gamma, beta, gammast, Rd, &
                        !          gm, Am, dAmc(i,j))
                        !    endif
                        ! else
                        !    call fAndAn_c3(cs, 0.0_r_2, X*cs, gamma, beta, gammast, Rd, &
                        !       Am, dAmc(i,j))
                        !    if (g0 > Am*X) then ! repeat calculation if g0 > A*X
                        !       call fAndAn_c3(cs, g0, 0.1e-4_r_2, gamma, beta, gammast, Rd, &
                        !          Am, dAmc(i,j))
                        !    endif
                        ! endif
                        Am = 0.0_r_2 ! temporaray
                  endif
                     anrubiscoz(i,j) = real(Am)
                        
                  endif  ! C3

                  ! C4, Rubisco limited, accounting for explicit mesophyll conductance
                  if (vcmxt4z(i,j) > 1.0e-8) then
                     anrubiscoz(i,j) = vcmxt4z(i,j)- rdxz(i,j)
                     dAmc(i,j)  = 0.0_r_2
                     eta_c(i,j) = 0.0_r_2
                  endif

                  ! C3, RuBP regeneration-limited, accounting for explicit mesophyll conductance
                  if ( (vcmxt3z(i,j) > 0.0) .and. &
                     (vx3z(i,j) > 1.0e-8) ) then
                     cs      = csxz(i,j)
                     gss = gsz(i,j)
                     gamma   = real(vx3z(i,j), r_2)
                     beta    = real(cx2z(i,j), r_2)
                     gammast = real(cx2z(i,j) / 2.0, r_2)
                     Rd      = real(rdxz(i,j), r_2)
                     gm      = gmes(i,j)

                     if (trim(cable_user%g0_switch) == 'default') then
                        call fAndAn_c3_givengs(cs, gss, gamma, beta, gammast, Rd, Am)
                     elseif (trim(cable_user%g0_switch) == 'maximum') then
                        ! if (cable_user%explicit_gm) then
                        !    call fAmdAm_c3(cs, 0.0_r_2, X*cs, gamma, beta, gammast, Rd, &
                        !       gm, Am, dAme(i,j))
                        !    ! repeat calculation if g0 > A*X
                        !    if (g0 > Am*X) then
                        !       call fAmdAm_c3(cs, g0, 0.1e-4_r_2, gamma, beta, gammast, Rd, &
                        !          gm, Am, dAme(i,j))
                        !    endif
                        ! else
                        !    call fAndAn_c3(cs, 0.0_r_2, X*cs, gamma, beta, gammast, Rd, &
                        !       Am, dAme(i,j))
                        !    ! repeat calculation if g0 > A*X
                        !    if (g0 > Am*X) then
                        !       call fAndAn_c3(cs, g0, 0.1e-4_r_2, gamma, beta, gammast, Rd, &
                        !          Am, dAme(i,j))
                        !    endif
                        ! endif
                     Am = 0.0_r_2 ! temporaray
                     endif
                     anrubpz(i,j) = real(Am)

                  endif  ! C3

                  ! C4, RuBP limited, accounting for explicit mesophyll conductance
                  if (vx4z(i,j) > 1.0e-8) then
                     anrubpz(i,j) = vx4z(i,j) - rdxz(i,j)
                     ! dAme(i,j)  = 0.0_r_2
                     ! eta_e(i,j) = 0.0_r_2
                  endif

                  ! Sink limited, accounting for explicit mesophyll conductance
                  if (vcmxt3z(i,j) > 1.0e-10) then ! C3

                     ansinkz(i,j) = 0.5 * vcmxt3z(i,j) - rdxz(i,j)
                     ! dAmp(i,j)  = 0.0_r_2
                     ! eta_p(i,j) = 0.0_r_2

                  elseif ((vcmxt4z(i,j) > 1.0e-10) ) then ! C4
                     cs         = csxz(i,j)
                     gamma      = real(kc4(i,j),r_2) ! k in Collatz 1992
                     beta       = 0.0_r_2
                     gammast    = 0.0_r_2
                     Rd         = real(rdxz(i,j), r_2)
                     gm         = gmes(i,j)

                     if (trim(cable_user%g0_switch) == 'default') then
                        if (cable_user%explicit_gm) then
                        !    call fAmdAm_c4(cs, g0, X*cs, gamma, beta, gammast, Rd, gm, &
                        !       Am, dAmp(i,j))
                        ! else
                        !    call fAndAn_c4(cs, g0, X*cs, gamma, beta, gammast, Rd, &
                        !       Am, dAmp(i,j))
                        Am = 0.0_r_2 ! temporaray
                        endif
                     elseif (trim(cable_user%g0_switch) == 'maximum') then
                        if (cable_user%explicit_gm) then
                        !    call fAmdAm_c4(cs, 0.0_r_2, X*cs, gamma, beta, gammast, Rd, gm, &
                        !       Am, dAmp(i,j))
                        !    if (g0 > Am*X) then
                        !       call fAmdAm_c4(cs, g0, 0.1e-4_r_2, gamma, beta, gammast, Rd, gm, &
                        !          Am, dAmp(i,j))
                        !    endif
                        ! else
                        !    call fAndAn_c4(cs, 0.0_r_2, X*cs, gamma, beta, gammast, Rd, &
                        !       Am, dAmp(i,j))
                        !    if (g0 > Am*X) then
                        !       call fAndAn_c4(cs, g0, 0.1e-4_r_2, gamma, beta, gammast, Rd, &
                        !          Am, dAmp(i,j))
                        !    endif
                        ! endif
                        Am = 0.0_r_2 ! temporaray
                        endif
                     endif
                     ansinkz(i,j) = real(Am)

                     ! dAmp(i,j)  = 0.0_r_2
                     ! eta_p(i,j) = 0.0_r_2
                     ! if (Am > 0.0_r_2) eta_p(i,j) = dAmp(i,j) * cs / Am
                  endif ! C3/C4

               endif ! deltlfz > 0.1

               ! minimal of three limited rates
               tmp3      = (/ anrubiscoz(i,j), anrubpz(i,j), ansinkz(i,j) /)
               ii        = minloc(tmp3)
               anxz(i,j) = tmp3(ii(1))


         else ! vlaiz(i,j) > C%lai_thresh

               anxz(i,:)       = 0.0
               anrubiscoz(i,:) = 0.0
               anrubpz(i,:)    = 0.0
               ansinkz(i,:)    = 0.0


         endif ! vlaiz(i,j) > C%lai_thresh

      enddo ! j=1, mf

   enddo ! i=1, mp

   END SUBROUTINE photosynthesis_gm_givengs
   !---------------------------------------------------------------------------------------


   FUNCTION light_inhibition(APAR) RESULT(xrd)
      !Mercado, L. M., Huntingford, C., Gash, J. H. C., Cox, P. M.,
      ! and Jogireddy, V.:
      ! Improving the representation of radiation
      !interception and photosynthesis for climate model applications,
      !Tellus B, 59, 553-565, 2007.
      ! Equation 3
      ! (Brooks and Farquhar, 1985, as implemented by Lloyd et al., 1995).
      ! Rc = Rd 0 < Io < 10 umol quanta m-2 s-1
      ! Rc = [0.5 * 0.05 ln(Io)] Rd Io > 10 umol quanta m-2 s-1
      ! JK: note that APAR is incoming radation in the original formulations
      !     However, APAR considered a good proxy here
      implicit none

      real, intent(in) :: APAR  ! absorbed PAR in umol m-2 s-1
      real             :: xrd   ! light inhibition of Rd (0-1)

      if (APAR > 10.0) then
         xrd = 0.5 - 0.05 * log(APAR)
      else
         xrd = 1.0
      endif

   END FUNCTION light_inhibition



   ELEMENTAL PURE FUNCTION ej3x(parx, alpha, convex, x) RESULT(z)

      implicit none

      real, intent(in) :: parx
      real, intent(in) :: alpha
      real, intent(in) :: convex
      real, intent(in) :: x
      real             :: z

      z = max(0.0, &
         0.25*( ( alpha*parx+x-sqrt( (alpha*parx+x)**2 - &
         4.0*convex*alpha*parx*x ) ) / (2.0*convex) ) )

   END FUNCTION ej3x


   ! ------------------------------------------------------------------------------


   ELEMENTAL PURE FUNCTION ej4x(parx,alpha,convex,x) RESULT(z)

      implicit none

      real, intent(in) :: parx ! Q mol photon m-2 s-1
      real, intent(in) :: alpha ! quantum efficiency in mol C (mol photon)-1
      real, intent(in) :: convex ! convexity parameter
      real, intent(in) :: x  ! vcmax
      real             :: z  ! rubsico-limited gross photosynthesis

      z = max(0.0, &
         ( alpha*parx+x-sqrt( (alpha*parx+x)**2 - &
         4.0*convex*alpha*parx*x ) ) / (2.0*convex) )

   END FUNCTION ej4x


   ! ------------------------------------------------------------------------------


   FUNCTION xgmesT(x) RESULT(z)
      ! Temperature response of mesophyll conductance
      ! parameters as in Knauer et al. 2019 (from Bernacchi et al. 2002)

      implicit none

      real, intent(in) :: x
      real             :: z

      real :: EHa, EHd, Entrop
      real :: xgmes

      if (trim(cable_user%Rubisco_parameters) == 'Bernacchi_2002') then
         EHa    = 49.6e3  ! J/mol
         EHd    = 437.4e3 ! J/mol
         Entrop = 1.4e3   ! J/mol/K
      else if (trim(cable_user%Rubisco_parameters) == 'Walker_2013') then
         EHa    = 7.4e3   ! J/mol
         EHd    = 434.0e3 ! J/mol
         Entrop = 1.4e3   ! J/mol/K
      endif

      CALL point2constants(C)

      xgmes = exp(EHa * (x - C%TrefK) / (C%TrefK * C%Rgas * x )) * &
         (1.0 + exp((C%TrefK * Entrop - EHd) / (C%TrefK * C%Rgas))) / &
         (1.0 + exp((x * Entrop - EHd) / (x * C%Rgas)))

      z = max(0.0, xgmes)

   END FUNCTION xgmesT


   ! ------------------------------------------------------------------------------


   FUNCTION xvcmxt3(Tk) RESULT(z)
      !  Vcmax temperature response (Arrhenius function)

      implicit none

      real, intent(in) :: Tk
      real             :: xv, z
      !real :: xvcnum, xvcden

      !real, parameter :: EHaVc  = 73637.0  ! J/mol (Leuning 2002)
      !real, parameter :: EHdVc  = 149252.0 ! J/mol (Leuning 2002)
      !real, parameter :: EntropVc = 486.0  ! J/mol/K (Leuning 2002)
      !real, parameter :: xVccoef = 1.17461 ! derived parameter

      ! Parameters calculated from Kumarathunge et al. 2019 with Tgrowth = 15 degC
      real, parameter :: EaV = 59700.0  ! J/mol
      real, parameter :: EdV = 200000.0 ! J/mol
      real, parameter :: dSV = 639.43   ! J/mol/K

      ! xVccoef=1.0+exp((EntropJx*C%TREFK-EHdJx)/(Rconst*C%TREFK))
      call point2constants(C)

      !xvcnum = xvccoef * exp(( EHaVc / (C%Rgas * C%TREFK ) )* ( 1.0 - C%TREFK/Tk ) )
      !xvcden = 1.0 + exp( ( EntropVc * Tk - EHdVc) / ( C%rgas * Tk ) )
      !z = max(0.0, xvcnum / xvcden)

      xv = exp(EaV * (Tk - C%TrefK) / (C%TrefK * C%Rgas * Tk )) * &
         (1.0 + exp((C%TrefK * dSV - EdV) / (C%TrefK * C%Rgas))) / &
         (1.0 + exp((Tk * dSV - EdV) / (Tk * C%Rgas)))
      z = max(0.0, xv)

   END FUNCTION xvcmxt3


   ! ------------------------------------------------------------------------------


   ! Explicit array dimensions as temporary work around for NEC inlining problem
   FUNCTION xvcmxt4(Tk) RESULT(z)

      implicit none

      real, intent(in) :: Tk  ! leaf temperature in Kelvin
      real             :: xc4, z

      !real, parameter :: q10c4 = 2.0
      ! parameter values from Massad et al. 2007, PCE 30, 1191-1204, slightly adjusted to
      ! get optimum at ~36degC as in Yamori et al. 2014
      real, parameter :: EHa    = 66000.0  ! 67.29e3_r_2   ! J/mol
      real, parameter :: EHd    = 145000.0 ! 144.57e3_r_2  ! J/mol
      real, parameter :: Entrop = 469.4217 ! 0.472e3_r_2   ! J/mol/K

      ! Q10 function replaced with Arrhenius function for the sake of parameter
      ! identifiability/interpretability
      !z = q10c4**(0.1*x - 2.5) / &
      !     ( (1.0 + exp(0.3 * (13.0 - x))) * (1.0 + exp(0.2 * (x - 38.0))) )

      call point2constants(C)

      xc4 = exp(EHa * (Tk - C%TrefK) / (C%TrefK * C%Rgas * Tk )) * &
         (1.0 + exp((C%TrefK * Entrop - EHd) / (C%TrefK * C%Rgas))) / &
         (1.0 + exp((Tk * Entrop - EHd) / (Tk * C%Rgas)))

      z = max(0.0, xc4)

   END FUNCTION xvcmxt4


   ! ------------------------------------------------------------------------------


   Subroutine xvcmxt3_acclim(Tk, Tgrowth, trf)
      ! acclimated temperature response of Vcmax. Kumarathunge et al., New. Phyt., 2019,
      ! Eq 7 and Table 2

      implicit none

      real, intent(in)  :: Tk, Tgrowth  ! instantaneous T in K, growth T in degC
      real, intent(out) :: trf

      real :: xVccoef, EHaVc, EHdVc, EntropVc, xvcnum, xvcden

      EHaVc = (42.6 + 1.14*Tgrowth) * 1000.0
      entropvc = (645.13 -0.38 * Tgrowth)
      EHdVc = 200000.0

      call point2constants(C)

      xVccoef  = 1.0 + exp((entropvc * C%TREFK - EHdVc)/  ( C%Rgas * C%TREFK ) )
      xvcnum   = xVccoef * exp((EHaVc / (C%Rgas * C%TREFK ) )* ( 1.0 - C%TREFK/Tk ) )
      xvcden   = 1.0 + exp((entropvc * Tk- EHdVc) / (C%Rgas * Tk ) )

      trf = max( 0.0, xvcnum / xvcden )

   end subroutine xvcmxt3_acclim


   ! ------------------------------------------------------------------------------


   Subroutine xvcmxt4_acclim(Tk, Tgrowth, trf)
      ! acclimated temperature response of Vcmax for C4 plants.
      ! Data fitted to obtain Topt-Tgrowth response as shown in Yamori et al. 2014, Photosny Research Fig. 5

      implicit none

      real, intent(in)  :: Tk, Tgrowth  ! instantaneous T in K, growth T in degC
      real, intent(out) :: trf

      real :: xVccoef, EHaVc, EHdVc, EntropVc, xvcnum, xvcden

      EHaVc = (45.0 + 1.05 * Tgrowth) * 1000.0
      entropvc = (472.0 -0.128913 * Tgrowth)
      EHdVc = 145000.0

      call point2constants(C)

      xVccoef  = 1.0 + exp((entropvc * C%TREFK - EHdVc)/  (C%Rgas * C%TREFK) )
      xvcnum   = xVccoef * exp((EHaVc / (C%Rgas * C%TREFK ) ) * (1.0 - C%TREFK/Tk) )
      xvcden   = 1.0 + exp((entropvc * Tk - EHdVc ) / (C%Rgas * Tk ) )

      trf = max(0.0, xvcnum / xvcden )

   end Subroutine xvcmxt4_acclim


   ! ------------------------------------------------------------------------------


   ELEMENTAL PURE FUNCTION xrdt(x)
      !  Atkin et al. (Eq 1, New Phytologist (2015) 206: 614-636)
      !variable Q10 temperature of dark respiration
      ! Originally from Tjoelker et al. 2001

      implicit none

      real, intent(in) :: x
      real             :: xrdt

      xrdt = (3.09 - 0.043*((x-273.15)+25.)/2.0)**((x-273.15 -25.0)/10.0)

   END FUNCTION xrdt


   ! ------------------------------------------------------------------------------


   FUNCTION xejmxt3(Tk) RESULT(z)
      ! Temperature response of Jmax (Arrhenius function)

      implicit none

      real, intent(in) :: Tk  ! Leaf temperature in Kelvin
      real             :: xj, z

      !real :: xjxnum, xjxden
      !real, parameter :: EHaJx  = 50300.0  ! J/mol (Leuning 2002)
      !real, parameter :: EHdJx  = 152044.0 ! J/mol (Leuning 2002)
      !real, parameter :: EntropJx = 495.0  ! J/mol/K (Leuning 2002)
      !real, parameter :: xjxcoef = 1.16715 ! derived parameter

      ! Parameters calculated from Kumarathunge et al. 2019 with Thome = 25degC and Tgrowth = 15degC
      real, parameter :: EaJ = 40710  ! J/mol (Leuning 2002)
      real, parameter :: EdJ = 200000 ! J/mol (Leuning 2002)
      real, parameter :: dSJ = 642.97 ! J/mol/K (Leuning 2002)

      call point2constants(C)

      !xjxnum = xjxcoef*exp( ( EHaJx / ( C%Rgas * C%TREFK)) * (1.0 - C%TREFK / Tk ) )
      !xjxden = 1.0 + exp((EntropJx * Tk - EHdJx) / ( C%Rgas * Tk ) )
      !z = max(0.0, xjxnum/xjxden)

      xj = exp(EaJ * (Tk - C%TrefK) / (C%TrefK * C%Rgas * Tk )) * &
         (1.0 + exp((C%TrefK * dSJ - EdJ) / (C%TrefK * C%Rgas))) / &
         (1.0 + exp((Tk * dSJ - EdJ) / (Tk * C%Rgas)))
      z = max(0.0, xj)

   END FUNCTION xejmxt3


   ! ------------------------------------------------------------------------------


   Subroutine xejmxt3_acclim(Tk, Tgrowth, Thome, trf)

      ! acclimated temperature response of Jmax. Kumarathunge et al., New. Phyt., 2019,
      ! Eq 7 and Table 2
      REAL, INTENT(IN) :: Tk, Tgrowth, Thome  ! instantaneous T in K, home and growth T in degC
      REAL, INTENT(OUT) :: trf
      REAL:: xVccoef, EHaVc, EHdVc, EntropVc, xvcnum, xvcden


      EHaVc = 40.71 * 1000.0
      EHdVc  = 200000.0
      entropvc = (658.77 -0.84 * Thome) -0.52 * (Tgrowth - Thome)

      call point2constants(C)

      xVccoef = 1.0 + exp((entropvc * C%TREFK - EHdVc)/  (C%Rgas * C%TREFK) )
      xvcnum  = xVccoef * exp((EHaVc / (C%Rgas * C%TREFK ) ) * ( 1.0 - C%TREFK/Tk) )
      xvcden  = 1.0 + exp((entropvc * Tk - EHdVc ) / (C%Rgas * Tk) )
      trf = max(real(0.0), xvcnum / xvcden )

   end subroutine xejmxt3_acclim


   ! ------------------------------------------------------------------------------


   subroutine fwsoil_calc_std(fwsoil, soil, ssnow, veg, wb)

      use cable_def_types_mod
      use cable_common_module, only : cable_user

      implicit none

      ! soil water modifier of stom. cond
      real, dimension(:),        intent(out) :: fwsoil
      type(soil_parameter_type), intent(in)  :: soil
      type(soil_snow_type),      intent(in)  :: ssnow
      type(veg_parameter_type),  intent(in)  :: veg
      real(r_2), dimension(:,:),  intent(in)  :: wb
      real, dimension(mp) :: rwater ! soil water availability

      rwater = max(1.0e-9, &
         sum(veg%froot * max(1.0e-9, min(1.0, real(wb) - &
         spread(soil%swilt, 2, ms))),2) / (soil%sfc-soil%swilt))

      ! Remove vbeta #56
      if (cable_user%gs_switch == 'medlyn') then
         fwsoil = max(1.0e-4, min(1.0, rwater))
      else
         fwsoil = max(1.0e-9, min(1.0, veg%vbeta * rwater))
      endif

   end subroutine fwsoil_calc_std


   ! ------------------------------------------------------------------------------


   subroutine fwsoil_calc_non_linear(fwsoil, soil, ssnow, veg, wb)

      use cable_def_types_mod

      implicit none

      ! soil water modifier of stom. cond
      real, dimension(:),        intent(out) :: fwsoil
      type(soil_parameter_type), intent(in)  :: soil
      type(soil_snow_type),      intent(in)  :: ssnow
      type(veg_parameter_type),  intent(in)  :: veg
      real(r_2), dimension(:,:),  intent(in)  :: wb

      real, dimension(mp)    :: rwater ! soil water availability
      real, dimension(mp, 3) :: xi, ti, si
      integer :: j

      rwater = max(1.0e-9, &
         sum(veg%froot * max(0.0, min(1.0, real(wb) - &
         spread(soil%swilt, 2, ms))),2) / (soil%sfc-soil%swilt))

      fwsoil = 1.

      rwater = soil%swilt + rwater * (soil%sfc-soil%swilt)

      xi(:, 1) = soil%swilt
      xi(:, 2) = soil%swilt + (soil%sfc - soil%swilt)/2.0
      xi(:, 3) = soil%sfc

      ti(:, 1) = 0.
      ti(:, 2) = 0.9
      ti(:, 3) = 1.0

      si(:, 1) = (rwater - xi(:,2)) / ( xi(:,1) - xi(:,2)) * &
         (rwater - xi(:,3)) / ( xi(:,1) - xi(:,3))

      si(:, 2) = (rwater - xi(:,1)) / ( xi(:,2) - xi(:,1)) * &
         (rwater - xi(:,3)) / ( xi(:,2) - xi(:,3))

      si(:, 3) = (rwater - xi(:,1)) / ( xi(:,3) - xi(:,1)) * &
         (rwater - xi(:,2)) / ( xi(:,3) - xi(:,2))

      do j=1, mp
         if (rwater(j) < soil%sfc(j) - 0.02) &
            fwsoil(j) = max(0., min(1., ti(j,1)*si(j,1) + &
            ti(j,2)*si(j,2) + ti(j,3)*si(j,3)))
      enddo

   end subroutine fwsoil_calc_non_linear


   ! ------------------------------------------------------------------------------


   ! ypw 19/may/2010 soil water uptake efficiency (see Lai and Katul 2000)
   subroutine fwsoil_calc_Lai_Katul(fwsoil, soil, ssnow, veg, wb)

      use cable_def_types_mod

      implicit none

      ! soil water modifier of stom. cond
      real, dimension(:),        intent(out) :: fwsoil
      type(soil_parameter_type), intent(in)  :: soil
      type(soil_snow_type),      intent(in)  :: ssnow
      type(veg_parameter_type),  intent(in)  :: veg
      real(r_2), dimension(:,:),  intent(in)  :: wb

      integer :: ns
      ! original is 0.01, 
      real, parameter :: rootgamma = 0.01    ! (19may2010)
      real, dimension(mp) :: dummy, normfac
      !--- local level dependent rwater
      real, dimension(mp,ms) :: frwater

      fwsoil(:)  = 0.0
      normFac(:) = 0.0

      do ns=1, ms
         dummy(:) = rootgamma / max(1.0e-3,real(wb(:,ns))-soil%swilt(:))

         frwater(:,ns) = max( 1.0e-4, &
            ((real(wb(:,ns))-soil%swilt(:))/soil%ssat(:))**dummy )

         fwsoil(:) = min(1.0, max(fwsoil(:), frwater(:,ns)))

         normFac(:) = normFac(:) + frwater(:,ns) * veg%froot(:,ns)
      enddo

   end subroutine fwsoil_calc_Lai_Katul


   ! ------------------------------------------------------------------------------

   ! not used
   ! SUBROUTINE fwsoil_calc_sli(fwsoil, soil, ssnow, veg)

   !   use cable_def_types_mod

   !   implicit none

   !   real, dimension(:),        intent(out)   :: fwsoil ! soil water modifier of stom. cond
   !   type(soil_parameter_type), intent(inout) :: soil
   !   type(soil_snow_type),      intent(inout) :: ssnow
   !   type(veg_parameter_type),  intent(inout) :: veg

   !   real(r_2), dimension(mp,ms) :: tmp1, tmp2, delta_root, alpha2a_root, alpha2_root

   !   ! Lai and Katul formulation for root efficiency function - vh 17/07/09
   !   alpha2a_root = max(ssnow%wb-soil%swilt_vec, 0.001_r_2) / (soil%ssat_vec)
   !   tmp1         = ssnow%wb - soil%swilt_vec
   !   tmp2         = spread(veg%gamma,2,ms) / tmp1 * log(alpha2a_root)
   !   where ((tmp1 > 0.001_r_2) .and. (tmp2 > -10.0_r_2))
   !      alpha2_root = exp(tmp2)
   !   elsewhere
   !      alpha2_root = 0.0_r_2
   !   endwhere

   !   WHERE (veg%froot > 0.0)
   !      delta_root = 1.0_r_2
   !   ELSEWHERE
   !      delta_root = 0.0_r_2
   !   ENDWHERE

   !   fwsoil  = maxval(real(alpha2_root*delta_root), dim=2)
   !   fwsoil  = max(0.0, fwsoil)

   ! END SUBROUTINE fwsoil_calc_sli


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

   ! ------------------------------------------------------------------------------
   ! not used
   !


   ! SUBROUTINE cubic_root_solver(a0,a1,a2,x1,x2,x3)

   !   USE cable_def_types_mod, only: r_2

   !   REAL(r_2), INTENT(IN) :: a0,a1,a2
   !   REAL(r_2), INTENT(OUT) :: x1,x2,x3

   !   REAL(r_2) :: Q, R, theta, a, b, c
   !   real :: pi_c = 3.1415927

   !   a = a2
   !   b = a1
   !   c = a0

   !   Q = (a**2 - 3*b)/9
   !   R = (2*a**3 - 9*a*b + 27*c)/54

   !   if (R**2 .lt. Q**2) then
   !      theta = acos(R/(Q**3)**0.5)
   !      x1 = -2 * (Q**0.5)*cos(theta/3) - a/3
   !      x2 = -2 * (Q**0.5)*cos((theta+2*pi_c)/3) - a/3
   !      x3 = -2 * (Q**0.5)*cos((theta-2*pi_c)/3) - a/3
   !   else
   !      x1 = 9999.0
   !      x2 = 9999.0
   !      x3 = 9999.0
   !   endif

   ! END SUBROUTINE cubic_root_solver


   ! ------------------------------------------------------------------------------


   ! functions for implicit mesophyll conductance
   elemental pure subroutine fabc(Cs, g0, x, gamma, beta, Gammastar, Rd, a, b, c)

      use cable_def_types_mod, only: r_2

      implicit none

      real(r_2), intent(in)  :: Cs, g0, x, gamma, beta, Gammastar, Rd
      real(r_2), intent(out) :: a, b, c

      a = (1.0_r_2-x)*Cs - x*beta
      b = -g0*Cs**2 + ((1.0_r_2-x)*(Rd-gamma)-g0*beta)*Cs - x*(gamma*Gammastar+Rd*beta)
      c = -g0*(Rd-gamma)*Cs**2 - g0*(gamma*Gammastar+Rd*beta)*Cs

   end subroutine fabc
   elemental pure subroutine fabc_givengs(Cs, gs, gamma, beta, Gammastar, Rd, a, b, c)

      use cable_def_types_mod, only: r_2

      implicit none

      real(r_2), intent(in)  :: Cs, gamma, beta, Gammastar, Rd, gs
      real(r_2), intent(out) :: a, b, c

      a = 1.0_r_2
      b = Rd - gamma - Cs*gs - beta*gs  
      c = gamma*gs*(Cs-Gammastar) - gs*Rd*(Cs+beta)

   end subroutine fabc_givengs

   elemental pure subroutine fabc_c4(Cs, g0, x, gamma, beta, Gammastar, Rd, a, b, c)

      use cable_def_types_mod, only: r_2

      implicit none

      real(r_2), intent(in)  :: Cs, g0, x, gamma, beta, Gammastar, Rd
      real(r_2), intent(out) :: a, b, c

      a = x
      b = (g0+gamma*(1.0_r_2-x))*Cs - x*(beta-Rd)
      c =  -gamma*g0*Cs**2 - g0*(beta-Rd)*Cs

   end subroutine fabc_c4


   elemental pure subroutine fAn_c3(a, b, c, A2)

      use cable_def_types_mod, only: r_2

      implicit none

      real(r_2), intent(in) :: a,b,c
      real(r_2), intent(out) :: A2

      real(r_2) :: s2

      s2 = b**2 - 4.0_r_2*a*c
      A2 = (-b - sqrt(s2))/(2.0_r_2*a)

   end subroutine fAn_c3


   elemental pure subroutine fAn_c4(a, b, c, A2)

      use cable_def_types_mod, only: r_2

      implicit none

      real(r_2), intent(in)  :: a, b, c
      real(r_2), intent(out) :: A2

      real(r_2) :: s2

      s2 = b**2 - 4.0_r_2*a*c
      A2 = (-b + sqrt(s2))/(2.0_r_2*a)

   end subroutine fAn_c4


   elemental pure subroutine fdabc(Cs, g0, x, gamma, beta, Gammastar, Rd, da, db, dc)

      use cable_def_types_mod, only: r_2

      implicit none

      real(r_2), intent(in)  :: Cs, g0, x, gamma, beta, Gammastar, Rd
      real(r_2), intent(out) :: da, db, dc

      da = 1.0_r_2-x
      db = -2.0_r_2*g0*Cs + (1.0_r_2-x)*(Rd-gamma) - g0*beta
      dc = -2.0_r_2*g0*(Rd-gamma)*Cs - g0*(gamma*Gammastar+Rd*beta)

   end subroutine fdabc


   elemental pure subroutine fdabc_c4(Cs, g0, x, gamma, beta, Gammastar, Rd, da, db, dc)

      use cable_def_types_mod, only: r_2

      implicit none

      real(r_2), intent(in)  :: Cs, g0, x, gamma, beta, Gammastar, Rd
      real(r_2), intent(out) :: da, db, dc

      da = 0.
      db = g0 + gamma*(1.0_r_2-x)
      dc = -2.0_r_2*gamma*g0*Cs - g0*(beta-Rd)

   end subroutine fdabc_c4


   elemental pure subroutine fdAn_c3(a, b, c, da, db, dc, dA2)

      use cable_def_types_mod, only: r_2

      implicit none

      real(r_2), intent(in) :: a, b, c, da, db, dc
      real(r_2), intent(out) :: dA2

      real(r_2) :: s, p

      s = sqrt(b**2 - 4.0_r_2*a*c)
      p = (2.0_r_2*b*db - 4.0_r_2*c*da - 4.0_r_2*a*dc)/(2.0_r_2*s)
      dA2 = (-db - p)/(2.0_r_2*a) - (-b - s)/(2.0_r_2*a**2)*da

   end subroutine fdAn_c3


   elemental pure subroutine fdAn_c4(a, b, c, da, db, dc, dA2)

      use cable_def_types_mod, only: r_2

      implicit none

      real(r_2), intent(in)  :: a, b, c, da, db, dc
      real(r_2), intent(out) :: dA2

      real(r_2) :: s, p

      s = sqrt(b**2 - 4.0_r_2*a*c)
      p = (2.0_r_2*b*db - 4.0_r_2*c*da - 4.0_r_2*a*dc)/(2.0_r_2*s)
      dA2 = (-db + p)/(2.0_r_2*a) - (-b + s)/(2.0_r_2*a**2)*da

   end subroutine fdAn_c4


   elemental pure subroutine fAndAn_c3(Cs, g0, x, gamma, beta, Gammastar, Rd, An, dAn,a,b,c)

      use cable_def_types_mod, only: r_2

      implicit none

      real(r_2), intent(in)  :: Cs, g0, x, gamma, beta, Gammastar, Rd
      real(r_2), intent(out) :: An, dAn,a,b,c

      real(r_2) ::  da, db, dc

      call fabc(Cs, g0, x, gamma, beta, Gammastar, Rd,a,b,c)

      call fAn_c3(a, b, c, An)

      call fdabc(Cs, g0, x, gamma, beta, Gammastar, Rd, da, db, dc)
      call fdAn_c3(a, b, c, da, db, dc, dAn)

   end subroutine fAndAn_c3

   elemental pure subroutine fAndAn_c3_givengs(Cs, gs, gamma, beta, Gammastar, Rd, An)

      use cable_def_types_mod, only: r_2

      implicit none

      real(r_2), intent(in)  :: Cs, gs, gamma, beta, Gammastar, Rd
      real(r_2), intent(out) :: An

      real(r_2) :: a, b, c

      call fabc_givengs(Cs, gs, gamma, beta, Gammastar, Rd,a,b,c)
      call fAn_c3(a, b, c, An)

   end subroutine fAndAn_c3_givengs
   elemental pure subroutine fAndAn_c4(Cs, g0, x, gamma, beta, Gammastar, Rd, An, dAn)

      use cable_def_types_mod, only: r_2

      implicit none

      real(r_2), intent(in)  :: Cs, g0, x, gamma, beta, Gammastar, Rd
      real(r_2), intent(out) :: An, dAn

      real(r_2) :: a, b, c, da, db, dc

      call fabc_c4(Cs, g0, x, gamma, beta, Gammastar, Rd,a,b,c)
      call fAn_c4(a, b, c, An)
      call fdabc_c4(Cs, g0, x, gamma, beta, Gammastar, Rd, da, db, dc)
      call fdAn_c4(a, b, c, da, db, dc, dAn)

   end subroutine fAndAn_c4


   ! not used
   ! elemental pure subroutine fAndAn(Cs, g0, x, Gammastar, Rd, gammac, betac, gammae, betae, &
   !      flag_eps, &
   !      Anc, Ane, An, dAnc, dAne, dAn)

   !   use cable_def_types_mod, only: r_2

   !   implicit none

   !   real(r_2), intent(in)  :: Cs, g0, x, Gammastar, Rd, gammac, betac, gammae, betae
   !   logical,   intent(in)  :: flag_eps
   !   real(r_2), intent(out) :: Anc, Ane, An, dAnc, dAne, dAn

   !   call fAndAn1(Cs, g0, x, gammac, betac, Gammastar, Rd, Anc, dAnc)
   !   call fAndAn1(Cs, g0, x, gammae, betae, Gammastar, Rd, Ane, dAne)

   !   ! An = min(Anc,Ane)
   !   ! if (An==Anc) then
   !   !    dAn = dAnc
   !   ! else
   !   !    dAn = dAne
   !   ! endif
   !   if (Anc < Ane) then
   !      An  = Anc
   !      dAn = dAnc
   !   else
   !      An  = Ane
   !      dAn = dAne
   !   endif

   !   if (flag_eps) then !  return elasticity instead of derivative
   !      dAnc = dAnc * Cs / Anc
   !      dAne = dAne * Cs / Ane
   !      dAn  = dAn  * Cs / An
   !   endif

   ! end subroutine fAndAn


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


   elemental pure subroutine fAm_c3(a, b, c1, d, p, q, Am)

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

   end subroutine fAm_c3


   elemental pure subroutine fAm_c4(a, b, c1, Am)

      use cable_def_types_mod, only: r_2

      implicit none

      real(r_2), intent(in) :: a, b, c1
      real(r_2), intent(out) :: Am

      real(r_2) :: s2

      s2 = b**2 - 4.0_r_2*a*c1
      Am = (-b + sqrt(s2))/(2.0_r_2*a)

   end subroutine fAm_c4


   elemental pure subroutine fdAm_c3(a, b, c1, d, p, q, da, db, dc, dd, dp, dq, dAm)

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

   end subroutine fdAm_c3


   elemental pure subroutine fdAm_c4(a, b, c1, da, db, dc, dAm)

      use cable_def_types_mod, only: r_2

      implicit none

      real(r_2), intent(in) :: a, b, c1, da, db, dc
      real(r_2), intent(out) :: dAm

      real(r_2) :: s, p

      s = sqrt(b**2 - 4.0_r_2*a*c1)
      p = (2.0_r_2*b*db - 4.0_r_2*c1*da - 4.0_r_2*a*dc)/(2.0_r_2*s)
      dAm = (-db + p)/(2.0_r_2*a) - (-b + s)/(2.0_r_2*a**2)*da

   end subroutine fdAm_c4


   elemental pure subroutine fAmdAm_c3(Cs, g0, x, gamma, beta, Gammastar, Rd, gm, &
      Am, dAm)

      use cable_def_types_mod, only: r_2

      implicit none

      real(r_2), intent(in)  :: Cs, g0, x, gamma, beta, Gammastar, Rd, gm
      real(r_2), intent(out) :: Am, dAm

      real(r_2) :: a, b, c, d, p, q, da, db, dc, dd, dp, dq

      call fabcd(Cs, g0, x, gamma, beta, Gammastar, Rd, gm, a, b, c, d)
      call fpq(a, b, c, d, p, q)
      call fAm_c3(a, b, c, d, p, q, Am)
      call fdabcd(Cs, g0, x, gamma, beta, Gammastar, Rd, gm, da, db, dc, dd)
      call fdpq(a, b, c, d, da, db, dc, dd, dp, dq)
      call fdAm_c3(a, b, c, d, p, q, da, db, dc, dd, dp, dq, dAm)

   end subroutine fAmdAm_c3


   elemental pure subroutine fAmdAm_c4(Cs, g0, x, gamma, beta, Gammastar, Rd, gm, &
      Am, dAm)

      use cable_def_types_mod, only: r_2

      implicit none

      real(r_2), intent(in)  :: Cs, g0, x, gamma, beta, Gammastar, Rd, gm
      real(r_2), intent(out) :: Am, dAm

      real(r_2) :: a, b, c1, da, db, dc

      call fabcm(Cs, g0, x, gamma, beta, Gammastar, Rd, gm,a,b,c1)
      call fAm_c4(a, b, c1, Am)
      call fdabcm(Cs, g0, x, gamma, beta, Gammastar, Rd, gm, da, db, dc)
      call fdAm_c4(a, b, c1, da, db, dc, dAm)

   end subroutine fAmdAm_c4


   ! not used
   ! elemental pure subroutine  fAmdAm(Cs, g0, x, Gammastar, Rd, gm, gammac, betac, gammae, betae, &
   !      flag_eps, &
   !      Amc, Ame, Am, dAmc, dAme, dAm)

   !   use cable_def_types_mod, only: r_2

   !   implicit none

   !   real(r_2), intent(in)  :: Cs, g0, x, Gammastar, &
   !        Rd, gm, gammac, betac,gammae, betae
   !   logical,   intent(in)  :: flag_eps
   !   real(r_2), intent(out) :: Amc, Ame, Am, dAmc, dAme, dAm

   !   call fAmdAm1(Cs, g0, x, gammac, betac, Gammastar, Rd, gm, Amc, dAmc)
   !   call fAmdAm1(Cs, g0, x, gammae, betae, Gammastar, Rd, gm, Ame, dAme)

   !   ! Am = min(Amc,Ame)
   !   ! if (Am==Amc) then
   !   !    dAm = dAmc
   !   ! else
   !   !    dAm = dAme
   !   ! endif
   !   if (Amc < Ame) then
   !      Am  = Amc
   !      dAm = dAmc
   !   else
   !      Am  = Ame
   !      dAm = dAme
   !   endif

   !   if (flag_eps)  then !  return elasticity instead of derivative
   !      dAmc = dAmc * Cs / Amc
   !      dAme = dAme * Cs / Ame
   !      dAm  = dAm  * Cs / Am
   !   endif

   ! end subroutine fAmdAm

END MODULE cable_canopy_module
