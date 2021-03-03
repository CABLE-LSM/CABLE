
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
! Purpose: acclimation of ratio of Jmax to Vcmax, with corresponding optimisation
! of contributions of electron-tranport and carboxylation rates to rate of net photosynthesis
!
! Called from: SUBROUTINE bgcdriver in casa_cable.F90
!
! History: Vanessa Haverd Aug 2017

! ==============================================================================
MODULE cable_optimise_JV_module

  Use cable_def_types_mod, ONLY: met_type, climate_type, canopy_type, veg_parameter_type, &
       mp, r_2

  USE cable_data_module,   ONLY: icanopy_type, point2constants
  USE TypeDef,             ONLY: i4b, dp
  USE cable_common_module, ONLY: cable_user

  TYPE( icanopy_type ) :: C

  ! variables local to module
  REAL, ALLOCATABLE :: APAR(:), Dleaf(:), Tleaf(:), cs(:), scalex(:), fwsoil(:), g0(:)
  REAL :: Anet, vcmax00, bjv, g1, kc0, ko0, ekc, eko, gam0, egam, alpha, gm0, a1, D0, qs, qm, qb
  REAL :: convex, Neff, relcost_J, Rd0, Tgrowth, Thome
  INTEGER :: nt,kk
  !REAL, PARAMETER :: relcost_J = 1.6 ! Chen et al. Oecologia, 1993, 93: 63-69
  ! (use this value for forced co-ordination)
  ! now moved to icanopy_type and (if Cc-based) recalculated in the
  ! adjust_JV_gm Subroutine
  !REAL, PARAMETER :: relcost_J = 2.3  ! use this value for optimisation algorithm
  LOGICAL, PARAMETER :: coord = .FALSE.  ! adjust ratioJV to force co-oridnation.
  ! otherwise maximise photosynthesis

CONTAINS

  ! ==============================================================================

  SUBROUTINE optimise_JV(veg, climate, ktauday, bjv, relcostJ)

    IMPLICIT NONE

    TYPE(veg_parameter_type), INTENT(INOUT) :: veg      ! vegetation parameters
    TYPE(climate_type),       INTENT(IN)    :: climate  ! climate variables
    INTEGER,                  INTENT(IN)    :: ktauday
    REAL, DIMENSION(mp),      INTENT(IN)    :: bjv
    REAL, DIMENSION(mp),      INTENT(IN)    :: relcostJ

    INTEGER :: k
    REAL, DIMENSION(mp) :: bjv_new
    REAL :: Anet_cost
    REAL, PARAMETER :: l_bound = 0.5
    REAL, PARAMETER :: u_bound = 5.0

    ! assign local ptrs to constants defined in cable_data_module
    CALL point2constants(C)

    nt = ktauday * C%ndays_optim

    ALLOCATE(APAR(nt))
    ALLOCATE(Dleaf(nt))
    ALLOCATE(Tleaf(nt))
    ALLOCATE(cs(nt))
    ALLOCATE(scalex(nt))
    ALLOCATE(g0(nt))
    ALLOCATE(fwsoil(nt))

    DO k=1,mp
       if (veg%frac4(k) < 0.001) then ! not C4
          if (cable_user%explicit_gm) then
             if (trim(cable_user%Rubisco_parameters) == 'Bernacchi_2002') then
                gam0 = C%gam0cc
                Kc0  = C%conkc0cc
                Ko0  = C%conko0cc
                egam = C%egamcc
                ekc  = C%ekccc
                eko  = C%ekocc
             else if (trim(cable_user%Rubisco_parameters) == 'Walker_2013') then
                gam0 = C%gam0ccw
                Kc0  = C%conkc0ccw
                Ko0  = C%conko0ccw
                egam = C%egamccw
                ekc  = C%ekcccw
                eko  = C%ekoccw
             endif
             gm0  = veg%gm(k)
             vcmax00 = veg%vcmaxcc(k) ! vcmax at standard temperature (25degC)
             qs   = C%qs
             qm   = C%qm
             qb   = C%qb
          else
             Kc0  = C%conkc0
             Ko0  = C%conko0
             ekc  = C%ekc
             eko  = C%eko
             gam0 = C%gam0
             egam = C%egam
             gm0  = veg%gm(k)    ! used in code but not in solution if explicit_gm = false
             vcmax00 = veg%vcmax(k) ! vcmax at standard temperature (25degC)
             qs   = 0.5
             qm   = 0.0
             qb   = 1.0
          endif
          g1  = veg%g1(k)
          Rd0 = veg%cfrd(k) * veg%vcmax(k)
          ! soil-moisture modifier to stomatal conductance
          fwsoil = climate%fwsoil(k,:)
          alpha  = veg%alpha(k) ! quantum efficiency for electron transport
          convex = veg%convex(k)
          a1 = veg%a1gs(k)
          D0 = veg%d0gs(k)
          relcost_J = relcostJ(k)
          Neff = vcmax00 + relcost_J*bjv(k)*vcmax00/4. ! effective nitrogen amount
          !for distribution between e-limited and c-limited processes

          ! optimisation for shade leaves
          APAR = climate%APAR_leaf_shade(k,:)*1e-6;
          Dleaf = max(climate%Dleaf_shade(k,:), 50.0)*1e-3 ! Pa -> kPa
          Tleaf = climate%Tleaf_shade(k,:)
          cs = climate%cs_shade(k,:)*1e-6
          scalex = climate%scalex_shade(k,:)
          g0 = veg%g0(k) * scalex / C%RGSWC

          if (cable_user%acclimate_photosyn) then
             Tgrowth = climate%mtemp(k)
             Thome = climate%mtemp_max20(k)
          endif

          if (coord) then
             if(diff_Ac_Aj(l_bound)*diff_Ac_Aj(u_bound) < 0.0) then
                bjv_new(k) = rtbis(diff_Ac_Aj,l_bound, u_bound, 0.001)
             else
                bjv_new(k) = bjv(k)
             endif
             veg%vcmax_shade(k) = Neff / (1. + relcost_J * bjv_new(k) / 4.0)
             veg%ejmax_shade(k) = veg%vcmax_shade(k) * bjv_new(k)
          else
             if(total_photosynthesis_cost(bjv(k)).lt.total_photosynthesis_cost(l_bound) .and. &
                  total_photosynthesis_cost(bjv(k)).lt.total_photosynthesis_cost(u_bound)) then

                Anet_cost = golden(l_bound,bjv(k),u_bound,total_photosynthesis_cost,0.01,bjv_new(k))
                veg%vcmax_shade(k) = Neff/(1.+relcost_J*bjv_new(k)/4.0)
                veg%ejmax_shade(k) = veg%vcmax_shade(k)*bjv_new(k)
             else
                bjv_new(k) = bjv(k)
                veg%vcmax_shade(k) = veg%vcmax(k)
                veg%ejmax_shade(k) = veg%ejmax(k)
             endif
          endif

          ! optimisation for sun leaves
          APAR = climate%APAR_leaf_sun(k,:)*1e-6;
          Dleaf = max(climate%Dleaf_sun(k,:), 50.0)*1e-3 ! Pa -> kPa
          Tleaf = climate%Tleaf_sun(k,:)
          cs = climate%cs_sun(k,:)*1e-6
          scalex = climate%scalex_sun(k,:)

          if (coord) then
             if (diff_Ac_Aj(l_bound)*diff_Ac_Aj(u_bound)<0) then
                bjv_new(k) = rtbis(diff_Ac_Aj,l_bound,u_bound,0.001)
                !call total_An_Ac_Aj(bjv_new(k),An,Ac,Aj)

             else
                bjv_new(k) = bjv(k)
             endif
             veg%vcmax_sun(k) = Neff/(1.+relcost_J*bjv_new(k)/4.0)
             veg%ejmax_sun(k) = veg%vcmax_sun(k)*bjv_new(k)
             !call total_An_Ac_Aj(bjv_new(k),An,Ac,Aj)
          else
             if(total_photosynthesis_cost(bjv(k)).lt.total_photosynthesis_cost(l_bound).and. &
                  total_photosynthesis_cost(bjv(k)).lt.total_photosynthesis_cost(u_bound)) then

                Anet_cost = golden(l_bound,bjv(k),u_bound,total_photosynthesis_cost,0.01,bjv_new(k))
                veg%vcmax_sun(k) = Neff/(1.+relcost_J*bjv_new(k)/4.0)
                veg%ejmax_sun(k) = veg%vcmax_sun(k)*bjv_new(k)
             else
                bjv_new(k) = bjv(k)
                veg%vcmax_sun(k) = veg%vcmax(k)
                veg%ejmax_sun(k) = veg%ejmax(k)
             endif
             !call total_An_Ac_Aj(bjv_new(k),An,Ac,Aj)
          endif
       else !C4
          bjv_new(k) = bjv(k)
          veg%vcmax_shade(k) = veg%vcmax(k)
          veg%ejmax_shade(k) = veg%ejmax(k)

          veg%vcmax_sun(k) = veg%vcmax(k)
          veg%ejmax_sun(k) = veg%ejmax(k)
       endif
    ENDDO

    DEALLOCATE(APAR)
    DEALLOCATE(Dleaf)
    DEALLOCATE(Tleaf)
    DEALLOCATE(cs)
    DEALLOCATE(g0)
    DEALLOCATE(scalex)
    DEALLOCATE(fwsoil)

  END SUBROUTINE optimise_JV

  ! Copied from cable_canopy.f90 (is there a better way??)
  ! ------------------------------------------------------------------------------

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

  ! ------------------------------------------------------------------------------

  subroutine fAn_c3(a, b, c, A2)

    use cable_def_types_mod, only: r_2

    implicit none

    real(r_2), intent(in)  :: a, b, c
    real(r_2), intent(OUT) :: A2

    real(r_2) :: s2

    s2 = b**2 - 4.0_r_2 * a * c
    A2 = (-b - sqrt(s2))/(2.0_r_2 * a)

  end subroutine fAn_c3

  ! ------------------------------------------------------------------------------

  subroutine fabc(Cs, g0, x, gamm, beta, Gammastar, Rd, a, b, c)

    use cable_def_types_mod, only: r_2

    implicit none

    real,      intent(in)  :: Cs, g0, x, gamm, beta, Gammastar, Rd
    real(r_2), intent(out) :: a, b, c

    a = real((1.0-x)*Cs - x*beta, r_2)
    b = real(-g0*Cs**2 + ((1.0-x)*(Rd-gamm)-g0*beta)*Cs - x*(gamm*Gammastar+Rd*beta), r_2)
    c = real(-g0*(Rd-gamm)*Cs**2 - g0*(gamm*Gammastar+Rd*beta)*Cs, r_2)

  end subroutine fabc

  ! ------------------------------------------------------------------------------

  subroutine fabcd(Cs, g0, x, gamm, beta, Gammastar, Rd, gm, a, b, c1, d)

    use cable_def_types_mod, only: r_2

    implicit none

    real,      intent(in)  :: Cs, g0, x, gamm, beta, Gammastar, Rd, gm
    real(r_2), intent(out) :: a, b, c1, d

    a  = real(x, r_2)
    b  = real((gm+g0-gm*x)*Cs + x*(Rd-gamm) - gm*x*beta, r_2)
    c1 = real(-gm*g0*Cs**2 + ((gm+g0-gm*x)*(Rd-gamm)-gm*g0*beta)*Cs - &
         gm*x*(gamm*Gammastar+Rd*beta), r_2)
    d  = real(-gm*g0*(Rd-gamm)*Cs**2 - gm*g0*(gamm*Gammastar+Rd*beta)*Cs, r_2)

  end subroutine fabcd

  ! ------------------------------------------------------------------------------

  subroutine fpq(a, b, c, d, p, q)

    use cable_def_types_mod, only: r_2

    implicit none

    real(r_2), intent(in)  :: a, b, c, d
    real(r_2), intent(out) :: p, q

    p = (3.0_r_2*a*c - b**2)/(3.0_r_2*a**2)
    q = (2.0_r_2*b**3 - 9.0_r_2*a*b*c + 27.0_r_2*a**2*d)/(27.0_r_2*a**3)

  end subroutine fpq

  ! ------------------------------------------------------------------------------

  subroutine fAm_c3(a, b, c1, d, p, q, Am)

    use cable_def_types_mod, only: r_2
    use mo_constants, only: pi => pi_dp

    implicit none

    real(r_2), intent(in)  :: a, b, c1, d, p, q
    real(r_2), intent(out) :: Am

    real(r_2) :: p3, pq, k

    p3 = -p/3.0_r_2
    pq = min(3.0_r_2*q/(2.0_r_2*p)*sqrt(1.0_r_2/p3), 0.999999999999_r_2)
    k  = 1.0_r_2
    Am = 2.0_r_2*sqrt(p3)*cos(acos(pq)/3.0_r_2 - 2.0_r_2*pi*k/3.0_r_2) - b/(3.0_r_2*a)

  end subroutine fAm_c3

  ! ------------------------------------------------------------------------------

  REAL FUNCTION total_photosynthesis_cost(bjv)

    use cable_canopy_module, only: xvcmxt3, xejmxt3, ej3x, xrdt, &
         xgmesT, xvcmxt3_acclim, xejmxt3_acclim
    use cable_def_types_mod, only: r_2

    implicit none

    real, intent(in) :: bjv

    INTEGER   :: k, j
    REAL      :: kct, kot ! , tdiff
    REAL      :: x, gamm, beta, gammastar, Rd, gm, jmaxt, trf, vcmax0
    REAL(r_2) :: a, b, c1, d
    REAL(r_2) :: p, q  ! if cable_user%explicit_gm
    REAL(r_2) :: Anc, Ane
    REAL      :: An(nt), Ac(nt), Aj(nt)

    CALL point2constants(C)

    An = 0.0
    Ac = 0.0
    Aj = 0.0
    j = 1
    vcmax0 = Neff/(1.+ relcost_J*bjv/4.0)

    DO k=1, nt
       if (APAR(k) .gt. C%tAPAR_optim) then
          if (cable_user%GS_SWITCH == 'medlyn') THEN
             x = 1.0 + (g1 * fwsoil(k)**qs) / SQRT(Dleaf(k))
          elseif (cable_user%GS_SWITCH == 'leuning') THEN
             x = 1.0 / ( 1.0/a1 + (Dleaf(k)*1.0e-3/D0))
          endif
          if (cable_user%acclimate_photosyn) then
             CALL xvcmxt3_acclim(Tleaf(k), Tgrowth, trf)
             gamm = Vcmax0*scalex(k)*trf*fwsoil(k)**qb
          else
             gamm = Vcmax0*scalex(k)*xvcmxt3(Tleaf(k))*fwsoil(k)**qb
          endif
          gm = gm0 * scalex(k) * max(0.15,fwsoil(k)**qm) * xgmesT(Tleaf(k))
          ! tdiff = Tleaf(k) - C%Trefk
          ! gammastar = C%gam0 * ( 1.0 + C%gam1 * tdiff &
          !     + C%gam2 * tdiff * tdiff )
          gammastar = gam0 * EXP( ( egam / (C%rgas*C%trefk) ) &
               * (1.0 - C%trefk/Tleaf(k) ) )

          Rd  = Rd0*scalex(k)*xrdt(Tleaf(k))*fwsoil(k)**qb * light_inhibition(APAR(k)*1.0e6)
          kct = kc0 * EXP( ( ekc / (C%rgas*C%trefk) ) &
               * ( 1.0 - C%trefk/Tleaf(k) ) )
          kot = ko0 * EXP( ( eko / (C%rgas*C%trefk) ) &
               * ( 1.0 - C%trefk/Tleaf(k) ) )

          beta = kct * (1.0+0.21/kot)
          ! Rubisco-limited
          if (cable_user%explicit_gm) then
             if (TRIM(cable_user%g0_switch) == 'default') then
                CALL fabcd(cs(k), g0(k)*fwsoil(k)**qs, x, gamm, beta, gammastar, Rd, gm, a, b, c1, d)
                CALL fpq(a, b, c1, d, p, q)
                CALL fAm_c3(a, b, c1, d, p, q, Anc)
             elseif (TRIM(cable_user%g0_switch) == 'maximum') then
                CALL fabcd(cs(k), 0.0, x, gamm, beta, gammastar, Rd, gm, a, b, c1, d)
                CALL fpq(a, b, c1, d, p, q)
                CALL fAm_c3(a, b, c1, d, p, q, Anc)
                if (g0(k)*fwsoil(k)**qs .gt. Anc*x/cs(k)) then
                   CALL fabcd(cs(k), g0(k)*fwsoil(k)**qs, 0.1e-4, gamm, beta, gammastar, Rd, gm, a, b, c1, d)
                   CALL fpq(a, b, c1, d, p, q)
                   CALL fAm_c3(a, b, c1, d, p, q, Anc)
                endif
             endif
          else  ! infinite (implicit) gm
             if (TRIM(cable_user%g0_switch) == 'default') then
                CALL fabc(cs(k), g0(k)*fwsoil(k)**qs, x, gamm, beta, gammastar, Rd, a, b, c1)
                CALL fAn_c3(a, b, c1, Anc) ! rubisco-limited
             elseif (TRIM(cable_user%g0_switch) == 'maximum') then
                CALL fabc(cs(k), 0.0, x, gamm, beta, gammastar, Rd, a, b, c1)
                CALL fAn_c3(a, b, c1, Anc) ! rubisco-limited
                if (g0(k)*fwsoil(k)**qs .gt. Anc*x/cs(k)) then
                   CALL fabc(cs(k), g0(k)*fwsoil(k)**qs, 0.1e-4, gamm, beta, gammastar, Rd, a, b, c1)
                   CALL fAn_c3(a, b, c1, Anc) ! rubisco-limited
                endif
             endif
          endif

          if (cable_user%acclimate_photosyn) then
             call xejmxt3_acclim(Tleaf(k), Tgrowth, Thome, trf)
             jmaxt = bjv*Vcmax0*scalex(k)*trf*fwsoil(k)**qb
          else
             jmaxt = bjv*Vcmax0*scalex(k)*xejmxt3(Tleaf(k))*fwsoil(k)**qb
          endif
          gamm = ej3x(APAR(k), alpha, convex, jmaxt)
          beta = 2.0 * gammastar
          ! RuBP regeneration-limited
          if (cable_user%explicit_gm) then
             if (TRIM(cable_user%g0_switch) == 'default') then
                CALL fabcd(cs(k), g0(k)*fwsoil(k)**qs, x, gamm, beta, gammastar, Rd, gm, a, b, c1, d)
                CALL fpq(a, b, c1, d, p, q)
                CALL fAm_c3(a, b, c1, d, p, q, Ane)
             elseif (TRIM(cable_user%g0_switch) == 'maximum') then
                CALL fabcd(cs(k), 0.0, x, gamm, beta, gammastar, Rd, gm, a, b, c1, d)
                CALL fpq(a, b, c1, d, p, q)
                CALL fAm_c3(a, b, c1, d, p, q, Ane)
                if (g0(k)*fwsoil(k)**qs .gt. Ane*x/cs(k)) then
                   CALL fabcd(cs(k), g0(k)*fwsoil(k)**qs, 0.1e-4, gamm, beta, gammastar, Rd, gm, a, b, c1, d)
                   CALL fpq(a, b, c1, d, p, q)
                   CALL fAm_c3(a, b, c1, d, p, q, Ane)
                endif
             endif
          else ! infinite (implicit) gm
             if (TRIM(cable_user%g0_switch) == 'default') then
                CALL fabc(cs(k), g0(k)*fwsoil(k)**qs, x, gamm, beta, gammastar, Rd, a, b, c1)
                CALL fAn_c3(a,b,c1,Ane) ! e-transport limited
             elseif (TRIM(cable_user%g0_switch) == 'maximum') then
                CALL fabc(cs(k), 0.0, x, gamm, beta, gammastar, Rd, a, b, c1)
                CALL fAn_c3(a,b,c1,Ane) ! e-transport limited
                if (g0(k)*fwsoil(k)**qs .gt. Ane*x/cs(k)) then
                   CALL fabc(cs(k), g0(k)*fwsoil(k)**qs, 0.1e-4, gamm, beta, gammastar, Rd, a, b, c1)
                   CALL fAn_c3(a, b, c1, Ane) ! e-transport limited
                endif
             endif
          endif

          An(j) = real(min(Anc, Ane))
       else
          An(j) = 0.0
       endif
       j = j+1
    ENDDO


    if (sum(An) > 0.0) then
       total_photosynthesis_cost = Neff / sum(An)
       !total_photosynthesis_cost = (relcost_J*Vcmax0*bjv/4.0 + Vcmax0)/sum(An)
    else
       total_photosynthesis_cost = 0.0
    endif

  END FUNCTION total_photosynthesis_cost

  ! ------------------------------------------------------------------------------

  REAL FUNCTION total_photosynthesis(bjv)

    use cable_canopy_module, only: xvcmxt3, xejmxt3, ej3x, xrdt, &
         xgmesT, xvcmxt3_acclim, xejmxt3_acclim
    use cable_def_types_mod, only: r_2

    implicit none

    !TYPE( icanopy_type ) :: C
    real, intent(in) :: bjv

    INTEGER   :: k, j
    REAL      :: kct, kot, tdiff
    REAL      :: x, gamm,  beta, gammastar, Rd, jmaxt, gm, trf, vcmax0
    REAL(r_2) :: a, b, c1, d, p, q
    REAL(r_2) :: Anc, Ane
    REAL      :: An(nt)

    CALL point2constants(C)
    An = 0.0
    j  = 1
    vcmax0 = Neff/(1.+ relcost_J*bjv/4.0)

    DO k=1,nt
       if (APAR(k) .gt. C%tAPAR_optim) then
          if (cable_user%GS_SWITCH == 'medlyn') THEN
             x = 1.0 + (g1 * fwsoil(k)**qs) / SQRT(Dleaf(k))
          elseif (cable_user%GS_SWITCH == 'leuning') THEN
             x = 1.0 / ( 1.0/a1 + (Dleaf(k)*1.0e-3/D0))
          endif

          if (cable_user%acclimate_photosyn) then
             cALL xvcmxt3_acclim(Tleaf(k), Tgrowth, trf)
             gamm = Vcmax0*scalex(k)*trf*fwsoil(k)**qb
          else
             gamm = Vcmax0*scalex(k)*xvcmxt3(Tleaf(k))*fwsoil(k)**qb
          endif
          gm = gm0 * scalex(k) * max(0.15,fwsoil(k)**qm) * xgmesT(Tleaf(k))
          tdiff = Tleaf(k) - C%Trefk
          !gammastar = C%gam0 * ( 1.0 + C%gam1 * tdiff                  &
          !     + C%gam2 * tdiff * tdiff )
          gammastar = gam0 * EXP( ( egam / (C%rgas*C%trefk) ) &
               * (1.0 - C%trefk/Tleaf(k) ) )
          Rd  = Rd0*scalex(k)*xrdt(Tleaf(k))*fwsoil(k)**qb * light_inhibition(APAR(k)*1.0e6)
          kct = kc0 * EXP( ( ekc / (C%rgas*C%trefk) ) &
               * ( 1.0 - C%trefk/Tleaf(k) ) )
          kot = ko0 * EXP( ( eko / (C%rgas*C%trefk) ) &
               * ( 1.0 - C%trefk/Tleaf(k) ) )

          beta = kct * (1.0+0.21/kot)
          ! Rubisco-limited
          if (cable_user%explicit_gm) then
             if (TRIM(cable_user%g0_switch) == 'default') then
                CALL fabcd(cs(k), g0(k)*fwsoil(k)**qs, x, gamm, beta, gammastar, Rd, gm, a, b, c1, d)
                CALL fpq(a,b,c1,d,p,q)
                CALL fAm_c3(a,b,c1,d,p,q,Anc)
             elseif (TRIM(cable_user%g0_switch) == 'maximum') then
                CALL fabcd(cs(k),0.0, x, gamm, beta, gammastar, Rd, gm, a, b, c1, d)
                CALL fpq(a,b,c1,d,p,q)
                CALL fAm_c3(a,b,c1,d,p,q,Anc)
                if (g0(k)*fwsoil(k)**qs .gt. Anc*x/cs(k)) then
                   CALL fabcd(cs(k),g0(k)*fwsoil(k)**qs, 0.1e-4, gamm, beta, gammastar, Rd, gm, a, b, c1, d)
                   CALL fpq(a,b,c1,d,p,q)
                   CALL fAm_c3(a,b,c1,d,p,q,Anc)
                endif
             endif
          else ! infinite (implicit) gm
             if (TRIM(cable_user%g0_switch) == 'default') then
                CALL fabc(cs(k), g0(k)*fwsoil(k)**qs, x, gamm, beta, gammastar, Rd, a, b, c1)
                CALL fAn_c3(a,b,c1,Anc) ! rubisco-limited
             elseif (TRIM(cable_user%g0_switch) == 'maximum') then
                CALL fabc(cs(k),0.0, x, gamm, beta, gammastar, Rd, a, b, c1)
                CALL fAn_c3(a,b,c1,Anc) ! rubisco-limited
                if (g0(k)*fwsoil(k)**qs .gt. Anc*x/cs(k)) then
                   CALL fabc(cs(k), g0(k)*fwsoil(k)**qs, 0.1e-4, gamm, beta, gammastar, Rd, a, b, c1)
                   CALL fAn_c3(a,b,c1,Anc) ! rubisco-limited
                endif
             endif
          endif

          if (cable_user%acclimate_photosyn) then
             call xejmxt3_acclim(Tleaf(k), Tgrowth, Thome, trf)
             jmaxt = bjv*Vcmax0*scalex(k)*trf*fwsoil(k)**qb
          else
             jmaxt = bjv*Vcmax0*scalex(k)*xejmxt3(Tleaf(k))*fwsoil(k)**qb
          endif
          gamm = ej3x(APAR(k), alpha,convex, jmaxt)
          beta = 2.0 * gammastar
          ! RuBP regeneration-limited
          if (cable_user%explicit_gm) then
             if (TRIM(cable_user%g0_switch) == 'default') then
                CALL fabcd(cs(k), g0(k)*fwsoil(k)**qs, x, gamm, beta, gammastar, Rd, gm, a, b, c1, d)
                CALL fpq(a,b,c1,d,p,q)
                CALL fAm_c3(a,b,c1,d,p,q,Ane)
             elseif (TRIM(cable_user%g0_switch) == 'maximum') then
                CALL fabcd(cs(k),0.0, x, gamm, beta, gammastar, Rd, gm, a, b, c1, d)
                CALL fpq(a,b,c1,d,p,q)
                CALL fAm_c3(a,b,c1,d,p,q,Ane)
                if (g0(k)*fwsoil(k)**qs .gt. Ane*x/cs(k)) then
                   CALL fabcd(cs(k),g0(k)*fwsoil(k)**qs, 0.1e-4, gamm, beta, gammastar, Rd, gm, a, b, c1, d)
                   CALL fpq(a,b,c1,d,p,q)
                   CALL fAm_c3(a,b,c1,d,p,q,Ane)
                endif
             endif
          else ! infinite (implicit) gm
             if (TRIM(cable_user%g0_switch) == 'default') then
                CALL fabc(cs(k), g0(k)*fwsoil(k)**qs, x, gamm, beta, gammastar, Rd, a, b, c1)
                CALL fAn_c3(a,b,c1,Anc) ! rubisco-limited
             elseif (TRIM(cable_user%g0_switch) == 'maximum') then
                CALL fabc(cs(k), 0.0, x, gamm, beta, gammastar, Rd, a, b, c1)
                CALL fAn_c3(a,b,c1,Ane) ! e-transport limited
                if (g0(k)*fwsoil(k)**qs .gt. Ane*x/cs(k)) then
                   CALL fabc(cs(k), g0(k)*fwsoil(k)**qs, 0.1e-4, gamm, beta, gammastar, Rd, a, b, c1)
                   CALL fAn_c3(a,b,c1,Ane) ! e-transport limited
                endif
             endif
          endif
          An(j) = real(min(Anc, Ane))
       else
          An(j) = 0.0
       endif
       j = j+1
    ENDDO

    total_photosynthesis = sum(An)

  END FUNCTION total_photosynthesis

  ! ------------------------------------------------------------------------------

  REAL FUNCTION diff_Ac_Aj(bjv)

    use cable_canopy_module, only: xvcmxt3, xejmxt3, ej3x, xrdt, &
         xgmesT, xvcmxt3_acclim, xejmxt3_acclim
    use cable_def_types_mod, only: r_2

    implicit none

    !TYPE( icanopy_type ) :: C
    real, intent(in) :: bjv

    INTEGER   :: k, j  ! k is timestep!
    REAL      :: kct, kot, tdiff
    REAL      :: x, gamm,  beta, gammastar, Rd, jmaxt, gm, trf, vcmax0 ! c1 because C already taken
    REAL(r_2) :: a, b, c1, d, p, q  ! if cable_user%explicit_gm
    REAL(r_2) :: Anc, Ane
    REAL      :: An(nt), Ac(nt), Aj(nt)
    REAL      :: total_An, total_Ac, total_Aj

    CALL point2constants(C)
    An = 0.0
    Ac = 0.0
    Aj = 0.0
    j = 1
    vcmax0 = Neff/(1. + relcost_J*bjv/4.0)
    !bjv = (vcmax0/Neff -1.)*4.0/relcost_J

    DO k=1,nt
       if (APAR(k) .gt. C%tAPAR_optim) then
          Ac(j) = 0.0
          Aj(j) = 0.0
          if (cable_user%GS_SWITCH == 'medlyn') THEN
             x = 1.0 + (g1 * fwsoil(k)**qs) / SQRT(Dleaf(k))
          elseif (cable_user%GS_SWITCH == 'leuning') THEN
             x = 1.0 / ( 1.0/a1 + (Dleaf(k)*1.0e-3/D0))
          endif
          if (cable_user%acclimate_photosyn) then
             CALL xvcmxt3_acclim(Tleaf(k), Tgrowth, trf)
             gamm = Vcmax0*scalex(k)*trf*fwsoil(k)**qb
          else
             gamm = Vcmax0*scalex(k)*xvcmxt3(Tleaf(k))*fwsoil(k)**qb
          endif
          gm = gm0 * scalex(k) * max(0.15,fwsoil(k)**qm) * xgmesT(Tleaf(k))
          tdiff = Tleaf(k) - C%Trefk
          !gammastar = C%gam0 * ( 1.0 + C%gam1 * tdiff                  &
          !     + C%gam2 * tdiff * tdiff )
          gammastar = gam0 * EXP( ( egam / (C%rgas*C%trefk) ) &
               * (1.0 - C%trefk/Tleaf(k) ) )
          Rd  = Rd0 * scalex(k) * xrdt(Tleaf(k))*fwsoil(k)**qb * light_inhibition(APAR(k)*1.0e6)
          kct = kc0 * EXP( ( ekc / (C%rgas*C%trefk) ) &
               * ( 1.0 - C%trefk/Tleaf(k) ) )
          kot = ko0 * EXP( ( eko / (C%rgas*C%trefk) ) &
               * ( 1.0 - C%trefk/Tleaf(k) ) )
          beta = kct * (1.0+0.21/kot)

          ! Rubisco-limited
          if (cable_user%explicit_gm) then
             if (TRIM(cable_user%g0_switch) == 'default') then
                CALL fabcd(cs(k), g0(k)*fwsoil(k)**qs, x, gamm, beta, gammastar, Rd, gm, a, b, c1, d)
                CALL fpq(a,b,c1,d,p,q)
                CALL fAm_c3(a,b,c1,d,p,q,Anc)
             elseif (TRIM(cable_user%g0_switch) == 'maximum') then
                CALL fabcd(cs(k), 0.0, x, gamm, beta, gammastar, Rd, gm, a, b, c1, d)
                CALL fpq(a,b,c1,d,p,q)
                CALL fAm_c3(a,b,c1,d,p,q,Anc)
                if (g0(k)*fwsoil(k)**qs .gt. Anc*x/cs(k)) then
                   CALL fabcd(cs(k), g0(k)*fwsoil(k)**qs, 0.1e-4, gamm, beta, gammastar, Rd, gm, a, b, c1, d)
                   CALL fpq(a,b,c1,d,p,q)
                   CALL fAm_c3(a,b,c1,d,p,q,Anc)
                endif
             endif
          else ! infinite (implicit) gm
             if (TRIM(cable_user%g0_switch) == 'default') then
                CALL fabc(cs(k), g0(k)*fwsoil(k)**qs, x, gamm, beta, gammastar, Rd, a, b, c1)
                CALL fAn_c3(a,b,c1,Anc) ! rubisco-limited
             elseif (TRIM(cable_user%g0_switch) == 'maximum') then
                CALL fabc(cs(k), 0.0, x, gamm, beta, gammastar, Rd, a, b, c1)
                CALL fAn_c3(a,b,c1,Anc) ! rubisco-limited
                if (g0(k)*fwsoil(k)**qs .gt. Anc*x/cs(k)) then
                   CALL fabc(cs(k), g0(k)*fwsoil(k)**qs, 0.1e-4, gamm, beta, gammastar, Rd, a, b, c1)
                   CALL fAn_c3(a,b,c1,Anc) ! rubisco-limited
                endif
             endif
          endif

          if (cable_user%acclimate_photosyn) then
             call xejmxt3_acclim(Tleaf(k), Tgrowth, Thome, trf)
             jmaxt = bjv*Vcmax0*scalex(k)*trf*fwsoil(k)**qb
          else
             jmaxt = bjv*Vcmax0*scalex(k)*xejmxt3(Tleaf(k))*fwsoil(k)**qb
          endif
          gamm = ej3x(APAR(k), alpha,convex, jmaxt)
          beta  = 2.0 * gammastar

          ! RuBP regeneration-limited
          if (cable_user%explicit_gm) then
             if (TRIM(cable_user%g0_switch) == 'default') then
                CALL fabcd(cs(k), g0(k)*fwsoil(k)**qs, x, gamm, beta, gammastar, Rd, gm, a, b, c1, d)
                CALL fpq(a,b,c1,d,p,q)
                CALL fAm_c3(a,b,c1,d,p,q,Ane)
             elseif (TRIM(cable_user%g0_switch) == 'maximum') then
                CALL fabcd(cs(k), 0.0, x, gamm, beta, gammastar, Rd, gm, a, b, c1, d)
                CALL fpq(a,b,c1,d,p,q)
                CALL fAm_c3(a,b,c1,d,p,q,Ane)
                if (g0(k)*fwsoil(k)**qs .gt. Ane*x/cs(k)) then
                   CALL fabcd(cs(k), g0(k)*fwsoil(k)**qs, 0.1e-4, gamm, beta, gammastar, Rd, gm, a, b, c1, d)
                   CALL fpq(a,b,c1,d,p,q)
                   CALL fAm_c3(a,b,c1,d,p,q,Ane)
                endif
             endif
          else ! infinite (implicit) gm
             if (TRIM(cable_user%g0_switch) == 'default') then
                CALL fabc(cs(k), g0(k)*fwsoil(k)**qs, x, gamm, beta, gammastar, Rd, a, b, c1)
                CALL fAn_c3(a,b,c1,Ane) ! e-transport limited
             elseif (TRIM(cable_user%g0_switch) == 'maximum') then
                CALL fabc(cs(k), 0.0, x, gamm, beta, gammastar, Rd, a, b, c1)
                CALL fAn_c3(a,b,c1,Ane) ! e-transport limited
                if (g0(k)*fwsoil(k)**qs .gt. Ane*x/cs(k)) then
                   CALL fabc(cs(k), g0(k)*fwsoil(k)**qs, 0.1e-4, gamm, beta, gammastar, Rd, a, b, c1)
                   CALL fAn_c3(a,b,c1,Ane) ! e-transport limited
                endif
             endif
          endif
          An(j) = real(min(Anc, Ane))
          if (Anc < Ane) Ac(j) = real(Anc)
          if (Ane < Anc) Aj(j) = real(Ane)
       else
          An(j) = 0.0
       endif

       j = j+1

    ENDDO

    total_An = sum(An)
    total_Ac = sum(Ac)
    total_Aj = sum(Aj)
    diff_Ac_Aj = total_Ac-total_Aj

  END FUNCTION diff_Ac_Aj

  ! ------------------------------------------------------------------------------

  FUNCTION golden(ax,bx,cx,func,tol,xmin)

    IMPLICIT NONE

    REAL, INTENT(IN)  :: ax, bx, cx
    INTERFACE
       FUNCTION func(x)
         IMPLICIT NONE
         REAL, INTENT(IN) :: x
         REAL :: func
       END FUNCTION func
    END INTERFACE
    REAL, INTENT(IN)  :: tol
    REAL, INTENT(OUT) :: xmin
    REAL              :: golden

    REAL, PARAMETER :: R=0.61803399
    REAL, PARAMETER :: C=1.0-R
    REAL :: f1,f2,x0,x1,x2,x3

    x0 = ax
    x3 = cx
    if (abs(cx-bx) > abs(bx-ax)) then
       x1 = bx
       x2 = bx+C*(cx-bx)
    else
       x2 = bx
       x1 = bx-C*(bx-ax)
    end if
    f1 = func(x1)
    f2 = func(x2)
    do
       if (abs(x3-x0) <= tol*(abs(x1)+abs(x2))) exit
       if (f2 < f1) then
          call shft3(x0,x1,x2,R*x2+C*x3)
          call shft2(f1,f2,func(x2))
       else
          call shft3(x3,x2,x1,R*x1+C*x0)
          call shft2(f2,f1,func(x1))
       end if
    end do
    if (f1 < f2) then
       golden = f1
       xmin = x1
    else
       golden = f2
       xmin = x2
    end if

  CONTAINS

    SUBROUTINE shft2(a,b,c)
      REAL, INTENT(OUT)   :: a
      REAL, INTENT(INOUT) :: b
      REAL, INTENT(IN)    :: c
      a=b
      b=c
    END SUBROUTINE shft2

    SUBROUTINE shft3(a,b,c,d)
      REAL, INTENT(OUT) :: a
      REAL, INTENT(INOUT) :: b,c
      REAL, INTENT(IN) :: d
      a=b
      b=c
      c=d
    END SUBROUTINE shft3

  END FUNCTION golden

  ! ------------------------------------------------------------------------------

  SUBROUTINE total_An_Ac_Aj(bjv, total_An, total_Ac, total_Aj)

    use cable_canopy_module, only: xvcmxt3,xejmxt3, ej3x, xrdt, &
         xgmesT, xvcmxt3_acclim, xejmxt3_acclim
    use cable_def_types_mod, only: r_2

    implicit none

    !TYPE( icanopy_type ) :: C
    real, intent(in)  :: bjv
    real, intent(out) :: total_An, total_Ac, total_Aj

    INTEGER :: k, j
    REAL :: kct, kot, tdiff
    REAL :: x, gamm,  beta, gammastar, Rd, jmaxt, gm, trf, vcmax0
    REAL(r_2) :: a, b, c1, d, p, q  ! if cable_user%explicit_gm
    REAL(r_2) :: Anc, Ane
    REAL :: An(nt), Ac(nt), Aj(nt)

    CALL point2constants(C)

    An = 0.0
    Ac = 0
    Aj = 0
    j = 1
    vcmax0 = Neff/(1. + relcost_J*bjv/4.0)

    DO k=1,nt
       if (APAR(k) .gt. C%tAPAR_optim) then
          Ac(j) = 0.0
          Aj(j) = 0.0
          if (cable_user%GS_SWITCH == 'medlyn') THEN
             x = 1.0 + (g1 * fwsoil(k)**qs) / SQRT(Dleaf(k))
          elseif (cable_user%GS_SWITCH == 'leuning') THEN
             x = 1.0 / ( 1.0/a1 + (Dleaf(k)*1.0e-3/D0))
          endif

          if (cable_user%acclimate_photosyn) then
             CALL xvcmxt3_acclim(Tleaf(k), Tgrowth, trf)
             gamm = Vcmax0*scalex(k)*trf*fwsoil(k)**qb
          else
             gamm = Vcmax0*scalex(k)*xvcmxt3(Tleaf(k))*fwsoil(k)**qb
          endif
          gm = gm0 * scalex(k) * max(0.15,fwsoil(k)**qm) * xgmesT(Tleaf(k))
          tdiff = Tleaf(k) - C%Trefk
          !gammastar = C%gam0 * ( 1.0 + C%gam1 * tdiff                  &
          !                                 + C%gam2 * tdiff * tdiff )
          gammastar = gam0 * EXP( ( egam / (C%rgas*C%trefk) ) &
               * (1.0 - C%trefk/Tleaf(k) ) )
          Rd = Rd0*scalex(k)*xrdt(Tleaf(k))*fwsoil(k)**qb * light_inhibition(APAR(k)*1.0e6)
          kct = kc0 * EXP( ( ekc / (C%rgas*C%trefk) ) &
               * ( 1.0 - C%trefk/Tleaf(k) ) )
          kot = ko0 *EXP ( ( eko / (C%rgas*C%trefk) ) &
               * ( 1.0 - C%trefk/Tleaf(k) ) )

          beta = kct * (1.0+0.21/kot)

          ! Rubisco-limited
          if (cable_user%explicit_gm) then
             CALL fabcd(cs(k), 0.0, x, gamm, beta, gammastar, Rd, gm, a, b, c1, d)
             CALL fpq(a,b,c1,d,p,q)
             CALL fAm_c3(a,b,c1,d,p,q,Anc)
             if (g0(k)*fwsoil(k)**qs .gt. Anc*x/cs(k)) then
                CALL fabcd(cs(k), g0(k)*fwsoil(k)**qs, 0.1e-4, gamm, beta, gammastar, Rd, gm, a, b, c1, d)
                CALL fpq(a,b,c1,d,p,q)
                CALL fAm_c3(a,b,c1,d,p,q,Anc)
             endif
          else ! infinite (implicit) gm
             CALL fabc(cs(k), 0.0, x, gamm, beta, gammastar, Rd, a, b, c1)
             CALL fAn_c3(a,b,c1,Anc)
             if (g0(k)*fwsoil(k)**qs .gt. Anc*x/cs(k)) then
                CALL fabc(cs(k), g0(k)*fwsoil(k)**qs, 0.1e-4, gamm, beta, gammastar, Rd, a, b, c1)
                CALL fAn_c3(a,b,c1,Anc) ! rubisco-limited
             endif
          endif

          if (cable_user%acclimate_photosyn) then
             call xejmxt3_acclim(Tleaf(k), Tgrowth, Thome, trf)
             jmaxt = bjv*Vcmax0*scalex(k)*trf*fwsoil(k)**qb
          else
             jmaxt = bjv*Vcmax0*scalex(k)*xejmxt3(Tleaf(k))*fwsoil(k)**qb
          endif

          gamm = ej3x(APAR(k), alpha,convex, jmaxt)
          beta = 2.0 * gammastar
          ! RuBP regeneration-limited
          if (cable_user%explicit_gm) then
             CALL fabcd(cs(k), 0.0, x, gamm, beta, gammastar, Rd, gm, a, b, c1, d)
             CALL fpq(a,b,c1,d,p,q)
             CALL fAm_c3(a,b,c1,d,p,q,Ane)
             if (g0(k)*fwsoil(k)**qs .gt. Ane*x/cs(k)) then
                CALL fabcd(cs(k), g0(k)*fwsoil(k)**qs, 0.1e-4, gamm, beta, gammastar, Rd, gm, a, b, c1, d)
                CALL fpq(a,b,c1,d,p,q)
                CALL fAm_c3(a,b,c1,d,p,q,Ane)
             endif
          else ! infinite (implicit) gm
             CALL fabc(cs(k), 0.0, x, gamm, beta, gammastar, Rd, a, b, c1)
             CALL fAn_c3(a,b,c1,Ane)
             if (g0(k)*fwsoil(k)**qs .gt. Ane*x/cs(k)) then
                CALL fabc(cs(k), g0(k)*fwsoil(k)**qs, 0.1e-4, gamm, beta, gammastar, Rd, a, b, c1)
                CALL fAn_c3(a,b,c1,Ane) ! e-transport limited
             endif
          endif

          An(j) = real(min(Anc, Ane))
          if (Anc < Ane) Ac(j) = real(Anc)
          if (Ane < Anc) Aj(j) = real(Ane)
       else
          An(j) = 0.0
       endif
       j = j+1
    ENDDO

    total_An = sum(An)
    total_Ac = sum(Ac)
    total_Aj = sum(Aj)

  END SUBROUTINE total_An_Ac_Aj

  ! ------------------------------------------------------------------------------

  FUNCTION rtbis(func,x1,x2,xacc)

    !USE nrtype; USE nrutil, ONLY : nrerror
    use mo_utils, only: eq, le, ge
#ifdef __MPI__
    use mpi,      only: MPI_Abort
#endif

    IMPLICIT NONE

    INTERFACE
       FUNCTION func(x)
         IMPLICIT NONE
         REAL, INTENT(IN) :: x
         REAL :: func
       END FUNCTION func
    END INTERFACE
    REAL, INTENT(IN) :: x1, x2, xacc
    REAL :: rtbis

    INTEGER, PARAMETER :: MAXIT=40
    INTEGER :: j
    REAL :: dx,f,fmid,xmid
#ifdef __MPI__
    integer :: ierr
#endif

    fmid = func(x2)
    f    = func(x1)
    if (ge(f*fmid, 0.0)) then
       write(*,*) 'rtbis: root must be bracketed'
#ifdef __MPI__
       call MPI_Abort(0, 89, ierr) ! Do not know comm nor rank here
#else
       stop 89
#endif
    endif

    if (f < 0.0) then
       rtbis=x1
       dx=x2-x1
    else
       rtbis=x2
       dx=x1-x2
    end if

    do j=1,MAXIT
       dx = dx*0.5
       xmid = rtbis+dx
       fmid = func(xmid)
       if (le(fmid, 0.0)) rtbis=xmid
       if (abs(dx) < xacc .or. eq(fmid, 0.0)) RETURN
    end do

    write(*,*) 'rtbis: too many bisections'
#ifdef __MPI__
    call MPI_Abort(0, 90, ierr) ! Do not know comm nor rank here
#else
    stop 90
#endif
  END FUNCTION rtbis

  ! ==============================================================================

END MODULE cable_optimise_JV_module
