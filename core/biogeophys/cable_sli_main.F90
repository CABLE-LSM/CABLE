module cable_sli_main

  implicit none

contains

  SUBROUTINE sli_main(ktau, dt, veg, soil, ssnow, met, canopy, air, rad, SEB_only)

    ! Main subroutine for Soil-litter-iso soil model
    ! Vanessa Haverd, CSIRO Marine and Atmospheric Research
    ! and Matthias Cuntz, UFZ - Helmholtz Centre for Environmental Research, 2010
    ! Modified to operate for multiple veg tiles but a single soil column, March 2011
    ! Rewritten for same number of soil columns as veg tiles May 2012
    USE cable_def_types_mod,  ONLY: veg_parameter_type, soil_parameter_type, soil_snow_type, met_type, &
         canopy_type, air_type, radiation_type, ms, mp, r_2, i_d
    USE cable_common_module , ONLY: cable_user
    USE sli_numbers,          ONLY:  zero, half, one, two, four, thousand, & ! numbers
         Tzero, experiment, &                                       ! variables
         vars_met, vars, vars_snow, &                                  ! types
         MW, Rgas, Lambdas, lambdaf, csice, cswat, rhow, nsnow_max, e5, &
         freezefac, topmodel, alpha, botbc
    USE sli_utils,            ONLY: x, dx, par, setpar, setpar_Loetsch, printparams, &
         setx, plit, dxL, setlitterpar, esat, &
         esat_ice, slope_esat_ice, thetalmax, Tfrz,  hyofS, SEB, &
         zerovars, zerovars_met, zerovars_snow, &
         printvars, printvars_met, printvars_snow
    USE sli_roots,            ONLY: setroots, getrex
    USE sli_solve,            ONLY: solve
    USE cable_IO_vars_module, ONLY: wlogn

    IMPLICIT NONE

    REAL,                      INTENT(IN)    :: dt
    TYPE(veg_parameter_type),  INTENT(INOUT) :: veg     ! all r_1
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil    ! all r_1
    TYPE(soil_snow_type),      INTENT(INOUT) :: ssnow   ! r_1, r_2 desaster
    TYPE(met_type),            INTENT(INOUT) :: met     ! all r_1
    TYPE(canopy_type),         INTENT(INOUT) :: canopy  ! all r_1
    TYPE(air_type),            INTENT(INOUT) :: air     ! all r_1
    TYPE(radiation_type),      INTENT(IN)    :: rad
    INTEGER,                   INTENT(IN)    :: ktau ! integration step number
    INTEGER,                   INTENT(IN)    :: SEB_only ! integration step number

    REAL(r_2), PARAMETER :: emsoil = 0.97_r_2
    REAL(r_2), PARAMETER :: rhocp = 1.1822e3_r_2
    REAL(r_2), PARAMETER :: Dva   = 2.17e-5_r_2
    INTEGER(i_d) :: i, k, kk, setroot
    REAL(r_2)    :: ti, tf
    TYPE(vars_met),  DIMENSION(1:mp)      :: vmet ! Meteorology above soil
    TYPE(vars),      DIMENSION(1:mp)      :: vlit
    TYPE(vars),      DIMENSION(1:mp,1:ms) :: var
    TYPE(vars_snow), DIMENSION(1:mp)      :: vsnow
    INTEGER(i_d),    DIMENSION(1:mp)      :: nsteps
    REAL(r_2),       DIMENSION(1:mp,1:ms) :: Tsoil, S, thetai, Jsensible
    REAL(r_2),       DIMENSION(1:mp)      :: SL, TL, T0
    REAL(r_2),       DIMENSION(1:mp)      :: drn, evap, infil, qprec, qprec_snow,qprec_snow_tmp, runoff, runoff_sat
    REAL(r_2),       DIMENSION(1:mp)      :: win, wp, wpi, h0, deltah0, h0old,hsnowold, discharge
    REAL(r_2),       DIMENSION(1:mp)      :: ip, ipi ! volumetric ice content of profile (final and initial)
    REAL(r_2),       DIMENSION(1:mp,1:ms) :: wex, csoil, qex, kth, phi, thetal_max, Sliq, Ksat
    REAL(r_2),       DIMENSION(1:mp,1:ms) :: FS
    ! surface temperature (top of top soil layer or top of litter layer)
    REAL(r_2),       DIMENSION(1:mp)      :: Tsurface
    REAL(r_2),       DIMENSION(1:mp)      :: gr, grc
    REAL(r_2),       DIMENSION(1:mp)      :: Etrans
    REAL(r_2),       DIMENSION(1:mp)      :: gamm
    REAL(r_2),       DIMENSION(1:mp)      :: G0, H, lE, Epot
    REAL(r_2),       DIMENSION(1:mp,-nsnow_max:ms) :: qh, qvsig, qlsig, qvTsig, qvh
    REAL(r_2),       DIMENSION(1:mp)      :: deltaTa, lE_old !, SA, SB, wpAi, wpBi, wpA, wpB
    REAL(r_2),       DIMENSION(1:mp)      :: evap_pot, deltaice_cum_T, deltaice_cum_S, zdelta
    REAL(r_2),       DIMENSION(1:mp)      :: fws
    REAL(r_2),       DIMENSION(1:mp)      :: Qadvcum, Jcol_sensible, Jcol_latent_S, Jcol_latent_T
    REAL(r_2),       DIMENSION(1:mp)      :: tmp1d1, deltaEsnow
    REAL(r_2),       DIMENSION(1:mp)      :: hice, phie
    REAL(r_2)                             :: tmp1d1a, tmp1d2, tmp1d3, tmp1d4, &
         tmp1d5, tmp1d6, tmp1d7, tmp1d8, tmp1d9,tmp1d10, tmp1d11, &
         tmp1d12,tmp1d13, tmp1d14, tmp1d15, tmp1d16
    REAL(r_2) :: rbw, rbh, rrc ! resistances for output
    INTEGER(i_d), DIMENSION(1:mp) :: index

    ! Topmodel
    REAL(r_2), DIMENSION(1:mp) :: fsat ! topmodel saturated area
    REAL(r_2), DIMENSION(1:mp) :: qb   ! topmodel baseflow

    ! Model switches
    INTEGER(i_d), PARAMETER :: litter       = 2 ! which litter model
    ! 0: no litter
    ! 1: full litter
    ! 2: litter resistance
    INTEGER(i_d), PARAMETER :: advection    = 1 ! heat advection by water
    INTEGER(i_d), PARAMETER :: isotopologue = 0 ! which isotope
    ! 0: no isotope calculations
    ! 1: HDO
    ! 2: H218O
    ! 3: HDO & H218O
    ! 0: normal run
    INTEGER(i_d), PARAMETER :: condition    = 3 ! condition matrix before solving
    ! 0: no conditioning
    ! 1: condition columns
    ! 2: condition lines
    ! 3: condition first lines then columns
    LOGICAL, SAVE :: first = .true.
    INTEGER(i_d), SAVE  :: counter

    ! Error flag if nstep of SLI > nsteps_max: err=0 -> no error; err/=0 -> error
    INTEGER(i_d), DIMENSION(1:mp) :: err

    ! initialise cumulative variables
    ! Jcol_sensible = zero
    ! Jcol_latent_S = zero
    ! Jcol_latent_T = zero
    ! deltaice_cum_T = zero
    ! deltaice_cum_S = zero
    ! Jsensible = zero
    drn    = zero
    ! discharge = zero
    infil  = zero
    evap   = zero
    ! evap_pot  = zero
    runoff = zero
    ! qh = zerocanopy%fevc(i) = ecx(i)*(1.0-canopy%fwet(i))
    ! H      = zero
    ! G0      = zero
    ! lE     = zero
    ! csoil = zero
    ! kth = zero
    ! Qadvcum  = zero
    wex    = zero
    ! qlsig        = zero
    ! qvsig        = zero
    ! qvh        = zero
    ! qvtsig        = zero
    thetal_max = zero
    Sliq       = zero
    Ksat       = zero
    phie       = zero
    phi        = zero
    hice       = zero
    err        = 0

    vmet  = zerovars_met()
    vlit  = zerovars()
    var   = zerovars()
    vsnow = zerovars_snow()

    ! output files for testing purposes
    if (first) then
       ! open(unit=332, file="vh08.out", status="replace", position="rewind")
       ! open(unit=334, file="S.out", status="replace", position="rewind")
       ! open(unit=336, file="Tsoil.out", status="replace", position="rewind")
       ! open(unit=335, file="SEB.out", status="replace", position="rewind")
       ! open(unit=337, file="soil_log.out", status="replace", position="rewind")
       ! open(unit=338, file="thetai.out", status="replace", position="rewind")
       ! open(unit=340, file="snow.out", status="replace", position="rewind")
       ! open(unit=346, file="diags.out", status="replace", position="rewind")
       ! open(unit=369, file="vmet.out", status="replace", position="rewind", recl=20*20)
       ! open(unit=370, file="qex.out", status="replace", position="rewind")
       ! open(unit=371, file="q.out", status="replace", position="rewind")

       ! open(unit=339, file="latlong.out", status="replace", position="rewind")
       ! write(339,"(20000f8.2)") rad%latitude
       ! write(339,"(20000f8.2)") rad%longitude
       counter = 0
    endif

    counter = counter + 1

    ! Save soil / snow surface temperature from last time step:
    ssnow%otss = ssnow%tss

    ! set layer thicknesses
    if (.not. allocated(x)) call setx(mp, ms, soil)

    ! Set root density distribution (leave in for sli offline)
    setroot = 0  ! reset rooting depths
    if (setroot == 1) then
       call setroots(x*100.0_r_2, real(veg%F10,r_2), real(veg%ZR,r_2)*100.0_r_2, FS)
       veg%froot = real(FS)
    else
       FS = real(veg%froot,r_2)
    endif
    ! set required soil hydraulic params
    if (.not. allocated(par)) then
       index=(/(i,i=1,mp,1)/)
       if (experiment == 16) then
          call setpar_Loetsch(mp, ms, x-half*dx)
       else
          call setpar(mp, ms, soil, index)
       endif
    endif

    ! If we want solutes:
    ! if (.not. allocated(bd)) allocate(bd(soil%nhorizons(k)))

    ! Litter parameters:
    if (.not. allocated(plit)) then
       index=(/(i,i=1,mp,1)/)
       call setlitterpar(mp, veg, index)
    endif

    ! Met data above soil:
    if (first) then
       vmet%rha   = zero
       vmet%rrc   = zero
       vmet%Rn    = zero
       vmet%Rnsw  = zero
       vmet%cva   = zero
       vmet%civa  = zero
       vmet%phiva = zero
    endif

    vmet%Ta  = real(met%Tvair,r_2) - Tzero
    vmet%Da  = real(met%dva,r_2)
    vmet%rbh = ssnow%rtsoil
    vmet%rbw = vmet%rbh

    gr       = four * emsoil * (vmet%Ta+Tzero)**3 *5.67e-8_r_2 ! radiation conductance Wm-2K-1
    grc      = one/vmet%rbh   + gr/rhocp
    vmet%rrc = one/grc                            ! resistance to radiative and convective heat transfer

    vmet%rha   = max(min((esat(vmet%Ta)-vmet%Da)/esat(vmet%Ta),one),0.1_r_2)
    vmet%cva   = vmet%rha * esat(vmet%Ta)*0.018_r_2/thousand/8.314_r_2/(vmet%Ta+Tzero) ! m3 H2O (liq) m-3 (air)
    vmet%phiva = Dva * vmet%cva
    vmet%Rn    = canopy%fns
    ! vmet%Rnsw  = rad%qssabs  ! shortwave radiation absorbed
    vmet%Rnsw  = zero ! all radiation absorbed at snow surface
    ! vmet%Rnsw = vmet%Rn ! all radiation absorbed beneath snow surface
    Etrans     = max(canopy%fevc/air%rlam/thousand, zero) ! m s-1
    where (canopy%fevc .lt. zero)
       canopy%fevw = canopy%fevw + real(canopy%fevc)
       canopy%fevc = zero
    end where
    h0         = ssnow%h0

    ! zero runoff here, in case error is returned to avoid excessive runoff from previous time-step.
    ! (Runoff is multipled by dt in cable_driver.F90)
    ssnow%rnof1            = 0.0
    ssnow%rnof2            = 0.0
    ssnow%runoff           = 0.0
    ssnow%E_fusion_sn      = zero
    ssnow%E_sublimation_sn = zero
    ssnow%evap_liq_sn      = zero
    ssnow%surface_melt     = zero
    ! Set isotopes to zero
    vmet%civa = zero

    ! Initialisations:
    if (first) then
       do kk=1, mp
          vsnow(kk)%hsnow      = zero
          vsnow(kk)%depth      = zero
          vsnow(kk)%totdepth   = zero
          vsnow(kk)%wcol       = zero
          vsnow(kk)%dens       = 120._r_2 ! snow density kg m-3
          vsnow(kk)%tsn        = zero
          vsnow(kk)%kH         = 0.16_r_2 ! snow thermal cond (W m-2 K-1)
          vsnow(kk)%Dv         = Dva      ! m2 s-1
          vsnow(kk)%sl         = zero
          vsnow(kk)%kE         = zero
          vsnow(kk)%kth        = vsnow(kk)%kH
          vsnow(kk)%cv         = zero
          vsnow(kk)%hliq       = zero
          vsnow(kk)%melt       = zero
          vsnow(kk)%nsnow      = 0
          vsnow(kk)%nsnow_last = 0
          ssnow%cls(kk)        = one
          vsnow(kk)%J                      = zero
          vsnow(kk)%fsnowliq_max           = 0.03_r_2
          vsnow(kk)%deltaJlatent           = zero
          vsnow(kk)%deltaJsensible         = zero
          vsnow(kk)%Qadv_snow              = zero
          vsnow(kk)%Qadv_rain              = zero
          vsnow(kk)%Qadv_melt              = zero
          vsnow(kk)%Qadv_vap               = zero
          vsnow(kk)%Qcond_net              = zero
          vsnow(kk)%Qadv_transfer          = zero
          vsnow(kk)%Qmelt                  = zero
          vsnow(kk)%Qtransfer              = zero
          vsnow(kk)%FluxDivergence         = zero
          vsnow(kk)%deltaJ                 = zero
          vsnow(kk)%Qvap                   = zero
          vsnow(kk)%MoistureFluxDivergence = zero
          vsnow(kk)%Qprec                  = zero
          vsnow(kk)%Qevap                  = zero
          vsnow(kk)%deltawcol              = zero
          vsnow(kk)%Jsensible              = zero
          vsnow(kk)%Jlatent                = zero
          ssnow%latent_heat_sn(kk)         = zero
       enddo
       first = .false.
    endif

    Tsurface = ssnow%Tsurface
    T0       = ssnow%Tsurface
    deltaTa  = zero
    lE_old   = ssnow%lE
    zdelta   = ssnow%zdelta

    SL    = 0.5_r_2   ! degree of litter saturation
    Tsoil = ssnow%Tsoil
    TL(:) = Tsoil(:,1) ! litter T

    thetai      = ssnow%thetai
    var%thetai = thetai
    ssnow%smelt = zero
    ssnow%cls   = one
    S           = ssnow%S                ! degree of soil saturation

    ! ----------------------------------------------------------------
    ! Iinitialise phi where it is (frozen and saturated) and where (pond >zero)

    ! Test steep freezing curve
    where (Tsoil<Tfrz(S, par%he, one/(par%lambc*freezefac)))
       par%lam = par%lambc * freezefac
    elsewhere
       par%lam = par%lambc
    end where
    par%eta = two/par%lam + two + one

    ! where (Tsoil<Tfrz(S,par%he,one/par%lam))
    where (Tsoil<Tfrz(S, par%he, one/(par%lambc*freezefac)))
       thetal_max = thetalmax(Tsoil,S,par%he,one/par%lam,par%thre,par%the)
       Sliq       = (thetal_max - (par%the-par%thre))/par%thre
       Ksat       = par%Ke*exp(par%eta*log(Sliq))
    elsewhere
       Sliq       = S
       Ksat       = par%Ke
    endwhere

    ! where ((Tsoil<Tfrz(S,par%he,one/par%lam)) .and. (S>=1))
    where ((Tsoil<Tfrz(S, par%he, one/(par%lambc*freezefac))) .and. (S>=1))
       phi = par%phie*exp(-log(Sliq)/par%lam)*exp(par%eta*log(Sliq))
    endwhere

    where (h0>zero)
       hice(:)  = h0(:)*(S(:,1)-Sliq(:,1))
       phie(:)  = par(:,1)%phie*exp(-log(Sliq(:,1))/par(:,1)%lam)*exp(par(:,1)%eta*log(Sliq(:,1)))
       !phi(:,1) = max((phie -par(:,1)%he*Ksat(:,1)), (one+e5)*phie)+(h0(:)-hice(:))*Ksat(:,1)
       phi(:,1) = (one+e5)*phie + (h0(:)-hice(:))*Ksat(:,1)
    endwhere
    var%phi = phi

    ! ----------------------------------------------------------------

    do kk=1, mp
       vsnow(kk)%nsnow = ssnow%nsnow(kk)
       vsnow(kk)%wcol  = ssnow%snowd(kk)/thousand ! diagnostic SWE (with or without dedicated snow pack)
       if ( vsnow(kk)%nsnow > 0)  then ! dedicated snow pack
          ! define variables associated with snow
          ssnow%isflag(kk)   = 1
          vsnow(kk)%hsnow(:) = ssnow%smass(kk,1:nsnow_max)/thousand ! SWE (state variable)
          vsnow(kk)%depth(:) = ssnow%sdepth(kk,1:nsnow_max)  ! depth of snow pack (m)
          vsnow(kk)%dens(:)  = ssnow%ssdn(kk,1:nsnow_max)
          where (vsnow(kk)%dens <= 200._r_2)
             vsnow(kk)%fsnowliq_max = 0.03_r_2
          elsewhere
             vsnow(kk)%fsnowliq_max = 0.03_r_2 + (0.1_r_2 - 0.03_r_2)*(vsnow(kk)%dens-200._r_2)/vsnow(kk)%dens
          endwhere
          vsnow(kk)%fsnowliq_max = 0.1 !MC! ???
          vsnow(kk)%tsn(:) = real(ssnow%tggsn(kk,1:nsnow_max),r_2) - Tzero
          vsnow(kk)%kH(:)  = ssnow%sconds(kk,1:nsnow_max)
          vsnow(kk)%Dv(:)  = Dva*(max(real(ssnow%tggsn(kk,1:nsnow_max),r_2)/Tzero,0.0_r_2))**1.88_r_2 ! m2 s-1
          vsnow(kk)%sl     = slope_esat_ice(vsnow(kk)%tsn) * Mw/thousand/Rgas/(vsnow(kk)%tsn+Tzero)
          vsnow(kk)%kE     = vsnow(kk)%Dv*vsnow(kk)%sl*thousand*lambdaf
          vsnow(kk)%kth    = vsnow(kk)%kE + vsnow(kk)%kH
          vsnow(kk)%cv     = esat_ice(vsnow(kk)%tsn)*Mw/thousand/Rgas/(vsnow(kk)%tsn+Tzero) ! m3 m-3
          vsnow(kk)%hliq   = ssnow%snowliq(kk,1:nsnow_max)/thousand ! amount of liq snow water
          vsnow(kk)%melt   = zero ! amount of melted snow leaving each snowlayer (mm/dt)
       else
          ssnow%isflag(kk)       = 0
          vsnow(kk)%hsnow        = zero
          vsnow(kk)%depth        = zero
          vsnow(kk)%dens         = 120_r_2 ! snow density kg m-3
          vsnow(kk)%tsn          = zero
          vsnow(kk)%kH           = 0.16_r_2    ! snow thermal cond (W m-2 K-1)
          vsnow(kk)%Dv           = Dva ! *(Tzero/Tzero)**1.88_r_2 ! m2 s-1
          vsnow(kk)%sl           = zero
          vsnow(kk)%kE           = zero
          vsnow(kk)%kth          = vsnow(kk)%kH
          vsnow(kk)%cv           = zero
          vsnow(kk)%hliq         = zero
          vsnow(kk)%melt         = zero
          vsnow(kk)%Jlatent      = zero
          vsnow(kk)%Jsensible    = zero
          vsnow(kk)%J            = zero
          vsnow(kk)%nsnow        = 0
          vsnow(kk)%fsnowliq_max = 0.03_r_2
       endif
       ! heat stored in snowpack
       where (vsnow(kk)%hsnow > zero)
          vsnow(kk)%Jsensible = (vsnow(kk)%hsnow-vsnow(kk)%hliq)*rhow*csice*(vsnow(kk)%Tsn) + &
               vsnow(kk)%hliq*rhow*cswat*(vsnow(kk)%Tsn)
          vsnow(kk)%Jlatent   = (vsnow(kk)%hsnow-vsnow(kk)%hliq)*rhow*(-lambdaf)
       elsewhere
          vsnow(kk)%Jsensible = zero
          vsnow(kk)%Jlatent   = zero
       end where
       vsnow(kk)%totdepth               = sum(vsnow(kk)%depth)
       vsnow(kk)%J                      = sum(vsnow(kk)%Jsensible+vsnow(kk)%Jlatent)
       vsnow(kk)%deltaJ                 = zero
       vsnow(kk)%deltawcol              = zero
       vsnow(kk)%deltaJlatent           = zero
       vsnow(kk)%deltaJsensible         = zero
       vsnow(kk)%Qadv_snow              = zero
       vsnow(kk)%Qadv_rain              = zero
       vsnow(kk)%Qadv_melt              = zero
       vsnow(kk)%Qadv_vap               = zero
       vsnow(kk)%Qcond_net              = zero
       vsnow(kk)%Qadv_transfer          = zero
       vsnow(kk)%Qmelt                  = zero
       vsnow(kk)%Qtransfer              = zero
       vsnow(kk)%FluxDivergence         = zero
       vsnow(kk)%Qvap                   = zero
       vsnow(kk)%MoistureFluxDivergence = zero
       vsnow(kk)%Qprec                  = zero
       vsnow(kk)%Qevap                  = zero
       vsnow(kk)%nsnow_last             = vsnow(kk)%nsnow
       ssnow%osnowd(kk) = ssnow%snowd(kk)
    enddo

    deltaTa = zero
    lE_old  = ssnow%lE
    gamm    = real(veg%gamma,r_2)
    where (canopy%through>=met%precip_sn)
       qprec      = max((canopy%through-met%precip_sn)/thousand/dt , zero)             ! liq precip rate (m s-1)
       qprec_snow = (met%precip_sn)/thousand/dt
    elsewhere
       qprec = max(real(canopy%through,r_2), zero)
       qprec_snow = zero
    endwhere
    !if ( wlogn == 1011) then
    !     write(wlogn,*) 'prec', qprec(79), canopy%through(79), met%precip_sn(79)
    !  write(*,*) 'prec', qprec(79), canopy%through(79), met%precip_sn(79)
    !endif
    ! re-calculate qprec_snow and qprec based on total precip and air T (ref Jin et al. Table II, Hyd Proc, 1999
    ! qprec_tot = qprec + qprec_snow
    ! where (vmet%Ta > 2.5_r_2)
    !    qprec_snow = zero
    !    qprec = qprec_tot
    ! elsewhere ((vmet%Ta <= 2.5_r_2) .and. (vmet%Ta > 2.0_r_2))
    !    qprec_snow = 0.6_r_2 * qprec_tot
    !    qprec = qprec_tot - qprec_snow
    ! elsewhere ((vmet%Ta <= 2.0_r_2) .and. (vmet%Ta > zero))
    !    qprec_snow = (1._r_2 - (54.62_r_2 - 0.2_r_2 *(vmet%Ta + Tzero)))*qprec_tot
    !    qprec = qprec_tot - qprec_snow
    ! elsewhere (vmet%Ta <= zero)
    !    qprec = zero
    !    qprec_snow = qprec_tot
    ! endwhere
    qprec_snow_tmp = qprec_snow
    h0old = ssnow%h0 ! pond height

    do kk=1, mp
       hsnowold(kk) = sum(vsnow(kk)%hsnow(1:nsnow_max))
    enddo

    ! Heat balance variables

    ! Water balance variables:
    ipi = sum(ssnow%thetai*dx,2)  + h0*ssnow%thetai(:,1)/par(:,1)%thre       ! ice in profile initially
    ! water in profile initially
    wpi = sum((par%thr + (par%the-par%thr)*S)*dx,2)  + plit%thre*SL*dxL

    nsteps = 0
    ti     = zero
    tf     = dt ! initial and final times

    win    = zero ! water input (total precip)
    evap   = zero
    runoff = zero
    infil  = zero
    drn    = zero

    if (SEB_only == 0) then

       if (cable_user%fwsoil_switch.ne.'Haverd2013') then

          do kk=1, mp
             call getrex(ssnow%S(kk,:), ssnow%rex(kk,:), fws(kk), FS(kk,:), par(kk,:)%the, &
                  par(kk,:)%thw, Etrans(kk), gamm(kk), dx(kk,:), real(dt,r_2))
          enddo

       else
          fws = canopy%fwsoil

       endif

       do k=1,ms
          qex(:,k)= ssnow%rex(:,k)
       enddo

    endif

    ! topmodel - 0:no; 1: only sat; 2: only base; 3: sat and base
    zdelta     = zero
    fsat       = zero
    runoff_sat = zero
    qb         = zero
    if (topmodel > 0) then
       zdelta = max( sum((par%thr + (par%the-par%thr))*dx,2) - wpi, zero )
       if (mod(topmodel, 2) == 1) then
          ! saturated fraction
          fsat       = min( par(:,ms)%fsatmax * exp(-zdelta), one )
          runoff_sat = fsat * qprec * thousand * dt
          qprec      = qprec * (one-fsat)  ! precip available for infiltration
       endif
       if (topmodel > 1) then
          botbc = "zero flux"
          ! calculate topmodel base flow
          qb        = alpha * par(:,ms)%Ke*exp(-(par(:,ms)%zeta*zdelta))
          qex(:,ms) =  qex(:,ms) + qb
       endif
    endif

    if (SEB_only == 1) then
       do kk=1, mp
          ! call hyofS(S(kk,:), Tsoil(kk,:), par(kk,:), var(kk,:))
          call hyofS(S(kk,1), Tsoil(kk,1), par(kk,1), var(kk,1))
          CALL SEB(ms, par(kk,:), vmet(kk), vsnow(kk), var(kk,:), qprec(kk), qprec_snow(kk), dx(kk,:), &
               h0(kk), Tsoil(kk,:), &
               Tsurface(kk), G0(kk), lE(kk), Epot(kk), &
               tmp1d1a, tmp1d2, tmp1d3, tmp1d4, &
               tmp1d5, tmp1d6, tmp1d7, tmp1d8, tmp1d9,tmp1d10, tmp1d11, &
               tmp1d12,tmp1d13, tmp1d14, tmp1d15, tmp1d16)
       enddo
       print*, 'SEB01 ', mp, ms
       do kk=1, mp
          call printvars_met(vmet(kk))
          call printvars_snow(vsnow(kk))
          do i=1, ms
             print*, 'i/ms ', i, ms
             call printparams(par(kk,i))
             call printvars(var(kk,i))
          enddo
          print*, 'qprec ', qprec(kk)
          print*, 'qprec_snow ', qprec_snow(kk)
          print*, 'dx ', dx(kk,:)
          print*, 'h0 ', h0(kk)
          print*, 'Tsoil ', Tsoil(kk,:)
          print*, 'Tsurface ', Tsurface(kk)
          print*, 'G0 ', G0(kk)
          print*, 'lE ', lE(kk)
          print*, 'Epot ', Epot(kk)
       enddo
       canopy%ga   = real(G0)
       canopy%fes  = lE
       canopy%fhs  = canopy%fns - canopy%ga - real(canopy%fes)
       ssnow%tss   = real(Tsurface + Tzero)
       ssnow%potev = real(Epot)
    else ! full SLI
       ! save for output, because they get changed with litter in solve
       rbw = vmet(1)%rbw
       rbh = vmet(1)%rbh
       rrc = vmet(1)%rrc
       !write(*,*), 'b4 solve', ktau
       call hyofS(S, Tsoil, par, var)
       call solve( ti, tf, ktau, mp, qprec, qprec_snow, ms, dx, &
            h0, S, thetai, Jsensible, Tsoil, evap, &
            evap_pot, runoff, infil, drn, discharge, qh, &
            nsteps, vmet, vlit, vsnow, var, csoil, kth, phi, T0, Tsurface, &
            H, lE, G0, Qadvcum, Jcol_sensible, &
            Jcol_latent_S, Jcol_latent_T, deltaice_cum_T, &
            deltaice_cum_S, dxL, zdelta, SL, TL, &
            plit, par, qex=qex, &
            wex=wex, qvsig=qvsig, qlsig=qlsig, qvTsig=qvTsig, qvh=qvh, &
            deltaTa=deltaTa, lE_old=lE_old, &
            dolitter=litter, doisotopologue=isotopologue, docondition=condition, &
            doadvection=advection, err=err)

       qprec_snow = qprec_snow_tmp
       H             = (H/(tf-ti))
       lE            = lE/(tf-ti)
       G0            = G0/(tf-ti)
       Jcol_latent_S = Jcol_latent_S/(tf-ti)
       Jcol_latent_T = Jcol_latent_T/(tf-ti)
       Jcol_sensible = Jcol_sensible/(tf-ti)
       Qadvcum       = Qadvcum/(tf-ti)

       do kk=1, mp
          tmp1d1(kk) = (sum(vsnow(kk)%Jsensible) + sum(vsnow(kk)%Jlatent))
          ! heat stored in snowpack
          where (vsnow(kk)%hsnow > zero)
             vsnow(kk)%Jsensible = (vsnow(kk)%hsnow-vsnow(kk)%hliq)*rhow*csice*(vsnow(kk)%Tsn) + &
                  vsnow(kk)%hliq*rhow*cswat*(vsnow(kk)%Tsn)
             vsnow(kk)%Jlatent = (vsnow(kk)%hsnow-vsnow(kk)%hliq)*rhow*(-lambdaf)
          elsewhere
             vsnow(kk)%Jsensible = zero
             vsnow(kk)%Jlatent = zero
          endwhere
          deltaEsnow(kk) = sum(vsnow(kk)%Jsensible) + sum(vsnow(kk)%Jlatent) - tmp1d1(kk)
          deltah0(kk) = h0(kk)-h0old(kk)+sum(vsnow(kk)%hsnow(1:vsnow(kk)%nsnow))-hsnowold(kk)
       enddo
       do k=1, ms
          where (err(:)==0) ssnow%thetai(:,k) = thetai(:,k)
       end do
       ip  = sum(ssnow%thetai*dx,2) + h0*ssnow%thetai(:,1)/par(:,1)%thre   ! ice in profile at tf
       ! water at tf
       wp  = sum((par%thr + (par%the-par%thr)*S)*dx,2) + plit%thre*SL*dxL
       win = win + (qprec+qprec_snow)*(tf-ti)

       if (1 == 0 .and. wlogn == 1011) then
          k=79
          ! if(k==1) then
          !     write(*,"(100e16.6)") &
          !    (win(k)-(wp(k)-wpi(k)+deltah0(k)+runoff(k)+evap(k)+drn(k))-Etrans(k)*dt)*1000, &
          !    win(k)*1000, (wp(k)-wpi(k)+deltah0(k))*1000, runoff(k)*1000+drn(k)*1000, evap(k)*1000,  Etrans(k)*dt*1000
          ! endif
          write(332,"(i8,i8,18e16.6)") ktau, nsteps(k), wp(k)-wpi(k), infil(k)-drn(k), runoff(k), &
               win(k)-(wp(k)-wpi(k)+deltah0(k)+runoff(k)+evap(k)+drn(k))-Etrans(k)*dt, wp(k), &
               evap(k), evap_pot(k), infil(k), &
               drn(k), h0(k), Etrans(k)*dt, discharge(k), fws(k), (ip(k)-ipi(k)), fsat(k), runoff_sat(k), qb(k)
          write(334,"(100f15.6)") S(k,:), S(k,:)*par(k,:)%thre+par(k,:)%thr
          write(336,"(100f15.6)") Tsoil(k,:)
          write(335,"(100e20.12)") vmet(k)%Ta, Tsurface(k), zero, H(k), lE(k), &
               G0(k),Jcol_sensible(k),Jcol_latent_S(k), Jcol_latent_T(k), &
               vmet(k)%Rn, TL(k), SL(k), deltaice_cum_T(k), &
               deltaice_cum_S(k), zero, Tsurface(k), vmet(k)%rha, &
               Qadvcum(k), sum((Jsensible(k,:)-ssnow%gammzz(k,:)),1)
          write(338,"(100f18.6)") thetai(k,:)
          write(369,"(20e20.12)") vmet(k)%Ta, vmet(k)%rha, rbw, &
               rbh, rrc, vmet(k)%Rn, &
               vmet(k)%Da, vmet(k)%cva, vmet(k)%civa, &
               vmet(k)%phiva, Etrans(k), qprec(k), qprec_snow(k), rad%latitude(k),  rad%longitude(k)
          write(370,"(20e20.12)")  qex
          write(371,"(20e20.12)")  qvsig+qlsig
       endif
       if (1 == 0) then
          write(wlogn+100,"(100i8)")  nsteps
       endif

       ! Update variables for output:
       where (err(1:mp) == 0)
          ssnow%tss      = real(Tsurface + Tzero)
          ssnow%wbtot    = real(wp*thousand)
          canopy%ga      = real(G0)
          canopy%fes     = lE
          canopy%fhs     = canopy%fns - canopy%ga - real(canopy%fes)
          ssnow%rnof1    = real(runoff*thousand/dt )
          ssnow%rnof2    = real(drn*thousand/dt )
          ssnow%runoff   = ssnow%rnof1 + ssnow%rnof2
          ssnow%zdelta   = zdelta
          ssnow%SL       = SL
          ssnow%TL       = TL
          ssnow%delwcol  = (wp-wpi+deltah0)*thousand  ! includes change in snow pack via deltah0
          ssnow%Tsurface = Tsurface
          ssnow%lE       = lE
          ssnow%evap     = evap*thousand
          ssnow%nsteps   = real(nsteps)
          ssnow%h0       = h0
       endwhere
       do k=1, ms
          where (err(:) == 0)
             ssnow%tgg(:,k)    = real(Tsoil(:,k) + Tzero)
             ssnow%wb(:,k)     = real(S(:,k)*(par(:,k)%thr+(par(:,k)%the-par(:,k)%thr)))
             ssnow%wbice(:,k)  = thetai(:,k)
             ssnow%S(:,k)      = S(:,k)
             ssnow%Tsoil(:,k)  = Tsoil(:,k)
             ssnow%rex(:,k)    = wex(:,k)*thousand
             ssnow%kth(:,k)    = kth(:,k)
             ssnow%gammzz(:,k) = Jsensible(:,k)
          end where
       end do

       if (cable_user%fwsoil_switch.ne.'Haverd2013') then
          where (err(1:mp) == 0) canopy%fwsoil = fws
       endif

       if (litter==0) then
          ssnow%rlitt = zero
       else
          where (err(1:mp) == 0) ssnow%rlitt = dxL/vlit%Dv
       endif

       ! update CABLE snow variables
       do kk=1, mp
          if (err(kk) == 0) then
             ssnow%snowd(kk)               = real(vsnow(kk)%wcol*thousand)     ! amount of snow  (mm liq water eq)
             ! amount of snow in dedicated snow pack (mm liq water eq)
             ssnow%smass(kk,1:nsnow_max)   = real(vsnow(kk)%hsnow(:)*thousand)
             ssnow%sdepth(kk,1:nsnow_max)  = real(vsnow(kk)%depth(:))          ! depth of snow pack (m)
             ssnow%ssdn(kk,1:nsnow_max)    = real(vsnow(kk)%dens(:))           ! density of snow (kg m-3)
             ssnow%tggsn(kk,1:nsnow_max)   = real(vsnow(kk)%tsn(:) + Tzero)    ! abs T of snowpack
             ssnow%sconds(kk,1:nsnow_max)  = real(vsnow(kk)%kH(:))             ! thermal conductivty of snowpack
             ! amount of melted snow leaving bottom of snow pack (mm/dt)
             ssnow%smelt(kk)               = real(vsnow(kk)%Qmelt*thousand/dt)
             ssnow%nsnow(kk)               = vsnow(kk)%nsnow
             if (sum(ssnow%sdepth(kk,1:nsnow_max)) > zero) then
                ssnow%ssdnn(kk) = ssnow%snowd(kk)/sum(ssnow%sdepth(kk,1:nsnow_max))
             endif
             ! change in latent heat of snow-pack (excl snow-fall contribution)
             ssnow%E_fusion_sn(kk) = (sum(vsnow(kk)%Jlatent) - ssnow%latent_heat_sn(kk) - vsnow(kk)%Qadv_snow) / dt
             ssnow%E_sublimation_sn(kk) = zero
             if ((vsnow(kk)%tsn(1) < zero) .and. ((ssnow%tggsn(kk,1)-Tzero) < zero)) then
                ssnow%E_fusion_sn(kk) = zero
                ! assume change in latent heat occurs as sublimation
                ! when Tsnow at beginning and end of timestep is < 0
                ssnow%E_sublimation_sn(kk) = (sum(vsnow(kk)%Jlatent) - ssnow%latent_heat_sn(kk) - vsnow(kk)%Qadv_snow)/dt ! W m-2
                ssnow%E_sublimation_sn(kk) = ssnow%E_sublimation_sn(kk)/lambdas ! Wm-2 -> kgm-2s-1
             endif

             if (sum(vsnow(kk)%hliq(1:vsnow(kk)%nsnow)) > zero) then
                ssnow%evap_liq_sn(kk)  = evap(kk)*thousand
                ! net melt rate mm s-1
                ssnow%surface_melt(kk) = (sum(vsnow(kk)%hliq(1:vsnow(kk)%nsnow)*thousand) &
                     - sum(ssnow%snowliq(kk,1:nsnow_max)))/dt &
                     + ssnow%smelt(kk) - qprec(kk)*thousand
             endif
             ! amount of liq snow water
             ssnow%snowliq(kk,1:nsnow_max) = real(vsnow(kk)%hliq(:)*thousand)
             ssnow%latent_heat_sn(kk) = sum(vsnow(kk)%Jlatent)
             ssnow%Qadv_rain_sn(kk)   = vsnow(kk)%Qadv_rain/dt
          endif
       enddo

       where (err(1:mp) == 0) ssnow%isflag = 0

       ! snow output
       if (1 == 0) then
          k = 1
          write(340,"(100e16.6)") sum(vsnow(k)%hsnow(1:vsnow(k)%nsnow)), vsnow(k)%tsn(1),sum(vsnow(k)%hliq(1:vsnow(k)%nsnow)), &
               qprec_snow(k)*dt, vsnow(k)%Qmelt, qprec(k)*dt, &
               vsnow(k)%Qevap,vsnow(k)%Qvap,ssnow%albsoilsn(k,1), ssnow%albsoilsn(k,2), ssnow%sconds(k,1), &
               vsnow(k)%dens(1),sum(vsnow(k)%depth(1:vsnow(k)%nsnow)), vsnow(k)%J, &
               vsnow(k)%MoistureFluxDivergence, vsnow(k)%FluxDivergence, vsnow(k)%dens(nsnow_max), vsnow(k)%tsn(nsnow_max), &
               qh(k,0), vmet(k)%rbh
       endif

       where (err(1:mp) == 0)
          canopy%ofes = canopy%fes
          ! Update total latent heat to reflect updated soil component:
          canopy%fe = canopy%fev + real(canopy%fes)
          ! Update total sensible heat to reflect updated soil component:
          canopy%fh = real(canopy%fhv) + canopy%fhs
       end where

    endif ! SEB only

  END SUBROUTINE sli_main

end module cable_sli_main
