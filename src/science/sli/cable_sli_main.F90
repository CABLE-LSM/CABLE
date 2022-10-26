MODULE sli_main_mod

CONTAINS

  SUBROUTINE sli_main(ktau, dt, veg, soil, ssnow, met, canopy, air, rad, SEB_only)

    ! Main subroutine for Soil-litter-iso soil model
    ! Vanessa Haverd, CSIRO Marine and Atmospheric Research
    ! and Matthias Cuntz, UFZ - Helmholtz Centre for Environmental Research, 2010
    ! Modified to operate for multiple veg tiles but a single soil column, March 2011
    ! Rewritten for same number of soil columns as veg tiles May 2012
    USE cable_def_types_mod,       ONLY: veg_parameter_type, soil_parameter_type, soil_snow_type, met_type, &
         canopy_type, air_type, radiation_type, ms, mp, r_2, i_d
    USE cable_common_module , ONLY: cable_user
    USE sli_numbers,        ONLY:  zero, half, one, two, four, thousand, & ! numbers
         Tzero, experiment, &                                       ! variables
         vars_met, vars, params, vars_snow, &                                  ! types
         MW, snmin, Rgas, Lambdas, lambdaf, csice, cswat, rhow, nsnow_max, e5, &
         freezefac, topmodel, alpha, fsat_max, botbc
    USE sli_utils,          ONLY: x, dx, par, setpar, setpar_Loetsch, setx, plit, dxL, setlitterpar, esat, &
         esat_ice, slope_esat_ice, thetalmax, Tfrz,  hyofS, SEB
    USE sli_roots,          ONLY: setroots, getrex
    USE sli_solve,          ONLY: solve

    USE  cable_IO_vars_module, ONLY: wlogn, verbose

    IMPLICIT NONE
    !INTEGER, INTENT(IN)            :: wlogn
    !INTEGER :: wlogn = 10001   !use correct value from io_vars module
    REAL,                      INTENT(IN)    :: dt
    TYPE(veg_parameter_type),  INTENT(INOUT) :: veg     ! all r_1
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil    ! all r_1
    TYPE(soil_snow_type),      INTENT(INOUT) :: ssnow   ! r_1, r_2 desaster
    TYPE(met_type),            INTENT(INOUT) :: met     ! all r_1
    TYPE(canopy_type),         INTENT(INOUT) :: canopy  ! all r_1
    TYPE(air_type),            INTENT(INOUT) :: air     ! all r_1
    TYPE (radiation_type),     INTENT(IN)    :: rad
    INTEGER,                   INTENT(IN)    :: ktau ! integration step number
    INTEGER,                   INTENT(IN)    :: SEB_only ! integration step number

    REAL(r_2), PARAMETER :: emsoil=0.97
    REAL(r_2), PARAMETER :: rhocp=1.1822e3
    REAL(r_2), PARAMETER :: Dva = 2.17e-5
    INTEGER(i_d) :: i, k, kk, setroot
    REAL(r_2)    :: ti, tf
    TYPE(vars_met),  DIMENSION(1:mp)      :: vmet ! Meteorology above soil
    TYPE(vars),      DIMENSION(1:mp)      :: vlit
    TYPE(vars),      DIMENSION(1:mp,1:ms) :: var
    TYPE(vars_snow), DIMENSION(1:mp)      :: vsnow
    INTEGER(i_d),    DIMENSION(1:mp)      :: nsteps
    REAL(r_2),       DIMENSION(1:mp,1:ms) :: Tsoil, S, thetai, Jsensible
    REAL(r_2),       DIMENSION(1:mp)      :: SL, TL, T0
    REAL(r_2),       DIMENSION(1:mp)      :: drn, evap, infil, qprec, qprec_snow,qprec_snow_tmp, qprec_tot, runoff, runoff_sat
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
    INTEGER(i_d), PARAMETER :: septs        = 0 ! coupled or uncoupled energy and water calculation
    ! 0: coupled calc
    ! 1: uncoupled energy (T) and moisture (S)
    INTEGER(i_d), PARAMETER :: condition    = 3 ! condition matrix before solving
    ! 0: no conditioning
    ! 1: condition columns
    ! 2: condition lines
    ! 3: condition first lines then columns
    LOGICAL, SAVE :: first = .TRUE.
    INTEGER(i_d), SAVE  :: counter

    ! Error flag if nstep of SLI > nsteps_max: err=0 -> no error; err/=0 -> error
    INTEGER(i_d), DIMENSION(1:mp) :: err



    REAL(r_2), DIMENSION(ms) :: zmm,dzmm
    REAL(r_2), DIMENSION(mp) :: zaq

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

    ! output files for testing purposes
    IF (first .AND. verbose) THEN
       OPEN (unit=332,file="vh08.out",status="replace",position="rewind")
       OPEN (unit=334,file="S.out",status="replace",position="rewind")
       OPEN (unit=336,file="Tsoil.out",status="replace",position="rewind")
       OPEN (unit=335,file="SEB.out",status="replace",position="rewind")
       OPEN (unit=337,file="soil_log.out",status="replace",position="rewind")
       OPEN(unit=338, file="thetai.out", status="replace", position="rewind")
       OPEN(unit=340, file="snow.out", status="replace", position="rewind")
       OPEN(unit=345, file="diags.out",status="replace", position="rewind")
       OPEN(unit=369, file="vmet.out", status="replace", position="rewind", recl=20*20)
       OPEN(unit=370, file="qex.out",status="replace", position="rewind")
       OPEN(unit=371, file="q.out",status="replace", position="rewind")

       !open(unit=339, file="latlong.out",status="replace", position="rewind")
       ! write(339,"(20000f8.2)") rad%latitude
       !write(339,"(20000f8.2)") rad%longitude
       counter = 0
    ENDIF

    counter = counter + 1

    ! Save soil / snow surface temperature from last time step:
    ssnow%otss = ssnow%tss

    ! set layer thicknesses
    IF (.NOT. ALLOCATED(x)) CALL setx(mp, ms, soil)

    ! Set root density distribution (leave in for sli offline)
    setroot = 0  ! reset rooting depths
    IF (setroot == 1) THEN
       CALL setroots(x*100.0_r_2, REAL(veg%F10,r_2), REAL(veg%ZR,r_2)*100.0_r_2, FS)
       veg%froot = FS
    ELSE
       FS = REAL(veg%froot,r_2)
    ENDIF
    ! set required soil hydraulic params
    IF (.NOT. ALLOCATED(par)) THEN
       index=(/(i,i=1,mp,1)/)
       IF (experiment == 16) THEN
          CALL setpar_Loetsch(mp, ms, x-half*dx)
       ELSE
          CALL setpar(mp, ms, soil, index)
       ENDIF
    ENDIF

    ! If we want solutes:
!!$     if (.not. allocated(bd)) allocate(bd(soil%nhorizons(k)))

    ! Litter parameters:
    IF (.NOT. ALLOCATED(plit)) THEN
       index=(/(i,i=1,mp,1)/)
       CALL setlitterpar(mp, veg, index)
    ENDIF

    ! Met data above soil:
    IF (first) THEN
       vmet%rha   = zero
       vmet%rrc   = zero
       vmet%Rn    = zero
       vmet%Rnsw = zero
       vmet%cva   = zero
       vmet%civa  = zero
       vmet%phiva = zero
    ENDIF

    vmet%Ta  = REAL(met%Tvair,r_2) - Tzero
    vmet%Da  = REAL(met%dva,r_2)

    IF (cable_user%or_evap .AND. SEB_only==1) THEN
       vmet%rbh = ssnow%rtsoil +  (one - &
            ssnow%satfrac(:))*ssnow%rtevap_unsat(:) + ssnow%satfrac(:)*ssnow%rtevap_sat(:)
       vmet%rbw = vmet%rbh + canopy%sublayer_dz(:)/(0.27_r_2/1189.8)
    ELSE
       vmet%rbh = ssnow%rtsoil
       vmet%rbw = vmet%rbh
    END IF

    gr       = four * emsoil * (vmet%Ta+Tzero)**3 *5.67e-8_r_2 ! radiation conductance Wm-2K-1
    grc      = one/vmet%rbh   + gr/rhocp
    vmet%rrc = one/grc                            ! resistance to radiative and convective heat transfer

    vmet%rha   = MAX(MIN((esat(vmet%Ta)-vmet%Da)/esat(vmet%Ta),one),0.1_r_2)
    vmet%cva   = vmet%rha * esat(vmet%Ta)*0.018_r_2/thousand/8.314_r_2/(vmet%Ta+Tzero) ! m3 H2O (liq) m-3 (air)
    vmet%phiva = Dva * vmet%cva
    vmet%Rn    = canopy%fns
    ! vmet%Rnsw  = rad%qssabs  ! shortwave radiation absorbed
    vmet%Rnsw = zero ! all radiation absorbed at snow surface
    ! vmet%Rnsw = vmet%Rn ! all radiation absorbed beneath snow surface
    Etrans     = MAX(canopy%fevc/air%rlam/thousand, zero) ! m s-1
    WHERE (canopy%fevc .LT. zero)
       canopy%fevw = canopy%fevw+canopy%fevc
       canopy%fevc = zero
    END WHERE
    h0         = ssnow%h0

    ! zero runoff here, in case error is returned to avoid excessive runoff from previous time-step. (Runoff is multipled by dt in cable_driver.F90)
    ssnow%rnof1 = 0.0
    ssnow%rnof2 = 0.0
    ssnow%runoff = 0.0
    ! Set isotopes to zero
    vmet%civa = zero

    ! Initialisations:
    IF (first) THEN
       DO kk=1, mp
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
       ENDDO
       first = .FALSE.
    ENDIF


    Tsurface = ssnow%Tsurface
    T0       = ssnow%Tsurface
    deltaTa  = zero
    lE_old   = ssnow%lE
    zdelta   = ssnow%zdelta



    SL    = 0.5_r_2   ! degree of litter saturation
    Tsoil = ssnow%Tsoil
    TL(:) = Tsoil(:,1) ! litter T

    !thetai      = ssnow%thetai
    !var%thetai = thetai
    ssnow%smelt = zero
    ssnow%cls   = one

    !mrd561 UM changes
    thetai      =ssnow%wbice! ssnow%thetai
    var%thetai = thetai
    DO k=1,ms
       DO i=1,mp
          ssnow%S(i,k) = ssnow%wb(i,k)/(par(i,k)%thr+(par(i,k)%the-par(i,k)%thr))
       END DO
    END DO
    S           = ssnow%S                ! degree of soil saturation


    ! ----------------------------------------------------------------
    ! Iinitialise phi where it is (frozen and saturated) and where (pond >zero)

    ! Test steep freezing curve
    WHERE (Tsoil<Tfrz(S, par%he, one/(par%lambc*freezefac)))
       par%lam = par%lambc * freezefac
    ELSEWHERE
       par%lam = par%lambc
    END WHERE
    par%eta = two/par%lam + two + one

    ! where (Tsoil<Tfrz(S,par%he,one/par%lam))
    WHERE (Tsoil<Tfrz(S, par%he, one/(par%lambc*freezefac)))
       thetal_max = thetalmax(Tsoil,S,par%he,one/par%lam,par%thre,par%the)
       Sliq       = (thetal_max - (par%the-par%thre))/par%thre
       Ksat       = par%Ke*EXP(par%eta*LOG(Sliq))
    ELSEWHERE
       Sliq       = S
       Ksat       = par%Ke
    endwhere

    ! where ((Tsoil<Tfrz(S,par%he,one/par%lam)) .and. (S>=1))
    WHERE ((Tsoil<Tfrz(S, par%he, one/(par%lambc*freezefac))) .AND. (S>=1))
       phi = par%phie*EXP(-LOG(Sliq)/par%lam)*EXP(par%eta*LOG(Sliq))
    endwhere

    WHERE (h0>zero)
       hice(:)  = h0(:)*(S(:,1)-Sliq(:,1))
       phie(:)  = par(:,1)%phie*EXP(-LOG(Sliq(:,1))/par(:,1)%lam)*EXP(par(:,1)%eta*LOG(Sliq(:,1)))
       !phi(:,1) = max((phie -par(:,1)%he*Ksat(:,1)), (one+e5)*phie)+(h0(:)-hice(:))*Ksat(:,1)
       phi(:,1) = (one+e5)*phie + (h0(:)-hice(:))*Ksat(:,1)
    endwhere
    var%phi = phi

    ! ----------------------------------------------------------------

    DO kk=1, mp
       vsnow(kk)%nsnow = ssnow%nsnow(kk)
       vsnow(kk)%wcol  = ssnow%snowd(kk)/thousand ! diagnostic SWE (with or without dedicated snow pack)
       IF ( vsnow(kk)%nsnow > 0)  THEN ! dedicated snow pack
          ! define variables associated with snow
          ssnow%isflag(kk)   = 1
          vsnow(kk)%hsnow(:) = ssnow%smass(kk,1:nsnow_max)/thousand ! SWE (state variable)
          vsnow(kk)%depth(:) = ssnow%sdepth(kk,1:nsnow_max)  ! depth of snow pack (m)
          vsnow(kk)%dens(:)  = ssnow%ssdn(kk,1:nsnow_max)
          WHERE (vsnow(kk)%dens <= 200._r_2)
             vsnow(kk)%fsnowliq_max = 0.03_r_2
          ELSEWHERE
             vsnow(kk)%fsnowliq_max = 0.03_r_2 + (0.1_r_2 - 0.03_r_2)*(vsnow(kk)%dens-200._r_2)/vsnow(kk)%dens
          endwhere
          vsnow(kk)%fsnowliq_max = 0.1 !MC! ???
          vsnow(kk)%tsn(:) = REAL(ssnow%tggsn(kk,1:nsnow_max),r_2) - Tzero
          vsnow(kk)%kH(:)  = ssnow%sconds(kk,1:nsnow_max)
          vsnow(kk)%Dv(:)  = Dva*(MAX(REAL(ssnow%tggsn(kk,1:nsnow_max),r_2)/Tzero,0.0_r_2))**1.88_r_2 ! m2 s-1
          vsnow(kk)%sl     = slope_esat_ice(vsnow(kk)%tsn) * Mw/thousand/Rgas/(vsnow(kk)%tsn+Tzero)
          vsnow(kk)%kE     = vsnow(kk)%Dv*vsnow(kk)%sl*thousand*lambdaf
          vsnow(kk)%kth    = vsnow(kk)%kE + vsnow(kk)%kH
          vsnow(kk)%cv     = esat_ice(vsnow(kk)%tsn)*Mw/thousand/Rgas/(vsnow(kk)%tsn+Tzero) ! m3 m-3
          vsnow(kk)%hliq   = ssnow%snowliq(kk,1:nsnow_max)/thousand ! amount of liq snow water
          vsnow(kk)%melt   = zero ! amount of melted snow leaving each snowlayer (mm/dt)
       ELSE
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
       ENDIF
       ! heat stored in snowpack
       WHERE (vsnow(kk)%hsnow > zero)
          vsnow(kk)%Jsensible = (vsnow(kk)%hsnow-vsnow(kk)%hliq)*rhow*csice*(vsnow(kk)%Tsn) + &
               vsnow(kk)%hliq*rhow*cswat*(vsnow(kk)%Tsn)
          vsnow(kk)%Jlatent   = (vsnow(kk)%hsnow-vsnow(kk)%hliq)*rhow*(-lambdaf)
       ELSEWHERE
          vsnow(kk)%Jsensible = zero
          vsnow(kk)%Jlatent   = zero
       END WHERE
       vsnow(kk)%totdepth               = SUM(vsnow(kk)%depth)
       vsnow(kk)%J                      = SUM(vsnow(kk)%Jsensible+vsnow(kk)%Jlatent)
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
    ENDDO

    deltaTa = zero
    lE_old  = ssnow%lE
    gamm    = REAL(veg%gamma,r_2)

    IF (cable_user%test_new_gw) THEN
       WHERE (canopy%through>=met%precip_sn)
          qprec      = ssnow%fwtop/thousand            ! liq precip rate (m s-1)
          qprec_snow = (met%precip_sn)/thousand/dt
       ELSEWHERE
          qprec = MAX(ssnow%fwtop/thousand, zero)
          qprec_snow = zero
       endwhere
    ELSE
       WHERE (canopy%through>=met%precip_sn)
          qprec      = MAX((canopy%through-met%precip_sn)/thousand/dt , zero) ! liq precip rate (m s-1)
          qprec_snow = (met%precip_sn)/thousand/dt
       ELSEWHERE
          qprec = MAX(canopy%through, zero)
          qprec_snow = zero
       endwhere
    END IF

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

    DO kk=1, mp
       hsnowold(kk) = SUM(vsnow(kk)%hsnow(1:nsnow_max))
    ENDDO

    ! Heat balance variables

    ! Water balance variables:
    ipi = SUM(ssnow%thetai*dx,2)  + h0*ssnow%thetai(:,1)/par(:,1)%thre       ! ice in profile initially
    ! water in profile initially
    wpi = SUM((par%thr + (par%the-par%thr)*S)*dx,2)  + plit%thre*SL*dxL

    nsteps = 0
    ti     = zero
    tf     = dt ! initial and final times

    win    = zero ! water input (total precip)
    evap   = zero
    runoff = zero
    infil  = zero
    drn    = zero

    IF (SEB_only == 0) THEN

       IF (cable_user%fwsoil_switch.NE.'Haverd2013') THEN

          DO kk=1, mp
             CALL getrex(ssnow%S(kk,:), ssnow%rex(kk,:), fws(kk), FS(kk,:), par(kk,:)%the, &
                  par(kk,:)%thw, Etrans(kk), gamm(kk), dx(kk,:), REAL(dt,r_2))
          ENDDO

       ELSE
          fws = canopy%fwsoil

       ENDIF

       DO k=1,ms
          qex(:,k)= ssnow%rex(:,k)
       ENDDO

       IF (cable_user%test_new_gw) THEN

          qex(:,:) = qex(:,:) + ssnow%qhlev(:,1:ms)/thousand
          qex(:,ms) = qex(:,ms) + ssnow%Qrecharge(:)/thousand
          botbc = "zero flux"

       END IF


    ENDIF

    ! topmodel - 0:no; 1: only sat; 2: only base; 3: sat and base
    zdelta     = zero
    fsat       = zero
    runoff_sat = zero
    qb         = zero
    IF (topmodel > 0) THEN
       zdelta = MAX( SUM((par%thr + (par%the-par%thr))*dx,2) - wpi, zero )
       IF (MOD(topmodel, 2) == 1) THEN
          ! saturated fraction
          fsat       = MIN( par(:,ms)%fsatmax * EXP(-zdelta), one )
          runoff_sat = fsat * qprec * thousand * dt
          qprec      = qprec * (one-fsat)  ! precip available for infiltration
       ENDIF
       IF (topmodel > 1) THEN
          botbc = "zero flux"
          ! calculate topmodel base flow
          qb        = alpha * par(:,ms)%Ke*EXP(-(par(:,ms)%zeta*zdelta))
          qex(:,ms) =  qex(:,ms) + qb
       ENDIF
    ENDIF



    IF (SEB_only == 1) THEN

       DO kk=1, mp

          ! call hyofS(S(kk,:), Tsoil(kk,:), par(kk,:), var(kk,:))
          CALL hyofS(S(kk,1), Tsoil(kk,1), par(kk,1), var(kk,1))
          CALL SEB(ms, par(kk,:), vmet(kk), vsnow(kk), var(kk,:), qprec(kk), qprec_snow(kk), dx(kk,:), &
               h0(kk), Tsoil(kk,:), &
               Tsurface(kk), G0(kk), lE(kk),Epot(kk), &
               tmp1d1a, tmp1d2, tmp1d3, tmp1d4, &
               tmp1d5, tmp1d6, tmp1d7, tmp1d8, tmp1d9,tmp1d10, tmp1d11, &
               tmp1d12,tmp1d13, tmp1d14, tmp1d15, tmp1d16, ktau)

       ENDDO
       canopy%ga  = REAL(G0)
       canopy%fes = REAL(lE)
       canopy%fhs = canopy%fns - canopy%ga - REAL(canopy%fes)
       ssnow%tss  = REAL(Tsurface + Tzero)
       ssnow%potev  = REAL(Epot)

    ELSE ! full SLI
       ! save for output, because they get changed with litter in solve
       rbw = vmet(1)%rbw
       rbh = vmet(1)%rbh
       rrc = vmet(1)%rrc
       !write(*,*), 'b4 solve', ktau
       CALL solve( ti, tf, ktau, mp, qprec, qprec_snow, ms, dx, &
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

       DO kk=1, mp
          tmp1d1(kk) = (SUM(vsnow(kk)%Jsensible) + SUM(vsnow(kk)%Jlatent))
          ! heat stored in snowpack
          WHERE (vsnow(kk)%hsnow > zero)
             vsnow(kk)%Jsensible = (vsnow(kk)%hsnow-vsnow(kk)%hliq)*rhow*csice*(vsnow(kk)%Tsn) + &
                  vsnow(kk)%hliq*rhow*cswat*(vsnow(kk)%Tsn)
             vsnow(kk)%Jlatent = (vsnow(kk)%hsnow-vsnow(kk)%hliq)*rhow*(-lambdaf)
          ELSEWHERE
             vsnow(kk)%Jsensible = zero
             vsnow(kk)%Jlatent = zero
          endwhere
          deltaEsnow(kk) = SUM(vsnow(kk)%Jsensible) + SUM(vsnow(kk)%Jlatent) - tmp1d1(kk)
          deltah0(kk) = h0(kk)-h0old(kk)+SUM(vsnow(kk)%hsnow(1:vsnow(kk)%nsnow))-hsnowold(kk)
       ENDDO
       DO k=1, ms
          WHERE (err(:)==0) ssnow%thetai(:,k) = thetai(:,k)
       END DO
       ip  = SUM(ssnow%thetai*dx,2) + h0*ssnow%thetai(:,1)/par(:,1)%thre   ! ice in profile at tf
       ! water at tf
       wp  = SUM((par%thr + (par%the-par%thr)*S)*dx,2) + plit%thre*SL*dxL
       win = win + (qprec+qprec_snow)*(tf-ti)

       IF (verbose) THEN
          k=1
          WRITE(332,"(i8,i8,18e16.6)") ktau, nsteps(k), wp(k)-wpi(k), infil(k)-drn(k), runoff(k), &
               win(k)-(wp(k)-wpi(k)+deltah0(k)+runoff(k)+evap(k)+drn(k))-Etrans(k)*dt, wp(k), &
               evap(k), evap_pot(k), infil(k), &
               drn(k), h0(k), Etrans(k)*dt, discharge(k), fws(k), (ip(k)-ipi(k)), fsat(k), runoff_sat(k), qb(k)

!!$if (ktau==5) then
!!$
!!$ write(*,*) win(k), (wp(k)-wpi(k)), deltah0(k),runoff(k), evap(k), drn(k), Etrans(k)*dt, canopy%fwsoil(k), canopy%fevc(k)
!!$stop
!!$endif


          WRITE(334,"(100f15.6)") S(k,:), S(k,:)*par(k,:)%thre+par(k,:)%thr
          WRITE(336,"(100f15.6)") Tsoil(k,:)
          WRITE(335,"(100e20.12)") vmet(k)%Ta, Tsurface(k), zero, H(k), lE(k), &
               G0(k),Jcol_sensible(k),Jcol_latent_S(k), Jcol_latent_T(k), &
               vmet(k)%Rn, TL(k), SL(k), deltaice_cum_T(k), &
               deltaice_cum_S(k), zero, Tsurface(k), vmet(k)%rha, &
               Qadvcum(k), SUM((Jsensible(k,:)-ssnow%gammzz(k,:)),1)
          WRITE(338,"(100f18.6)") thetai(k,:)
          WRITE(369,"(20e20.12)") vmet(k)%Ta, vmet(k)%rha, rbw, &
               rbh, rrc, vmet(k)%Rn, &
               vmet(k)%Da, vmet(k)%cva, vmet(k)%civa, &
               vmet(k)%phiva, Etrans(k), qprec(k), qprec_snow(k), rad%latitude(k),  rad%longitude(k)
          WRITE(370,"(20e20.12)")  qex
          WRITE(371,"(20e20.12)")  qvsig+qlsig

       ENDIF

       IF (cable_user%test_new_gw) THEN

          DO i=1,mp
             ssnow%GWwb(i) = ssnow%GWwb(i)  + (ssnow%Qrecharge(i)-ssnow%qhlev(i,ms+1))*dt/soil%GWdz(i)/thousand

             IF (ssnow%GWwb(i) .GT. soil%GWssat_vec(i)) THEN
                ssnow%qhlev(i,ms+1) = ssnow%qhlev(i,ms+1)  + (ssnow%GWwb(i) - soil%GWssat_vec(i))*soil%GWdz(i)*thousand/dt
                ssnow%GWwb(i) = soil%GWssat_vec(i)
                ssnow%qhz(i) = SUM(ssnow%qhlev(i,1:ms+1),dim=1)
             END IF

             IF (ssnow%GWwb(i) .LT. soil%GWwatr(i)) THEN
                ssnow%qhlev(i,ms+1) = ssnow%qhlev(i,ms+1)  - (ssnow%GWwb(i) - soil%GWwatr(i))*soil%GWdz(i)*thousand/dt
                ssnow%GWwb(i) = soil%GWwatr(i)
                ssnow%qhz(i) = SUM(ssnow%qhlev(i,1:ms+1),dim=1)
             END IF
          END DO
       END IF


       ! Update variables for output:
       WHERE (err(1:mp) == 0)
          ssnow%tss      = REAL(Tsurface + Tzero)
          ssnow%wbtot    = REAL(wp*thousand)
          canopy%ga      = REAL(G0)
          canopy%fes     = lE
          canopy%fhs     = canopy%fns - canopy%ga - canopy%fes
          ssnow%rnof1    = REAL(runoff*thousand/dt )
          ssnow%rnof2    = REAL(drn*thousand/dt )
          ssnow%runoff = ssnow%rnof1 + ssnow%rnof2
          ssnow%zdelta   = zdelta
          ssnow%SL       = SL
          ssnow%TL       = TL
          ssnow%delwcol  = (wp-wpi+deltah0)*thousand  ! includes change in snow pack via deltah0
          ssnow%Tsurface = Tsurface
          ssnow%lE       = lE
          ssnow%evap     = evap*thousand
          ssnow%nsteps   = REAL(nsteps)
          ssnow%h0       = h0
       endwhere
       DO k=1, ms
          WHERE (err(:) == 0)
             ssnow%tgg(:,k)    = REAL(Tsoil(:,k) + Tzero)
             ssnow%wb(:,k)     = REAL(S(:,k)*(par(:,k)%thr+(par(:,k)%the-par(:,k)%thr)))
             ssnow%wbice(:,k)  = thetai(:,k)
             ssnow%S(:,k)      = S(:,k)
             ssnow%Tsoil(:,k)  = Tsoil(:,k)
             ssnow%rex(:,k)    = wex(:,k)*thousand
             ssnow%kth(:,k)    = kth(:,k)
             ssnow%gammzz(:,k) = Jsensible(:,k)
          END WHERE
       END DO

       IF (cable_user%fwsoil_switch.NE.'Haverd2013') THEN
          WHERE (err(1:mp) == 0) canopy%fwsoil = REAL(fws)
       ENDIF

       IF (litter==0) THEN
          ssnow%rlitt = zero
       ELSE
          WHERE (err(1:mp) == 0) ssnow%rlitt = dxL/vlit%Dv
       ENDIF

       ! update CABLE snow variables
       DO kk=1, mp
          IF (err(kk) == 0) THEN
             ssnow%snowd(kk)               = REAL(vsnow(kk)%wcol*thousand)     ! amount of snow  (mm liq water eq)
             ! amount of snow in dedicated snow pack (mm liq water eq)
             ssnow%smass(kk,1:nsnow_max)   = REAL(vsnow(kk)%hsnow(:)*thousand)
             ssnow%sdepth(kk,1:nsnow_max)  = REAL(vsnow(kk)%depth(:))          ! depth of snow pack (m)
             ssnow%ssdn(kk,1:nsnow_max)    = REAL(vsnow(kk)%dens(:))           ! density of snow (kg m-3)
             ssnow%tggsn(kk,1:nsnow_max)   = REAL(vsnow(kk)%tsn(:) + Tzero)    ! abs T of snowpack
             ssnow%sconds(kk,1:nsnow_max)  = REAL(vsnow(kk)%kH(:))             ! thermal conductivty of snowpack
             ssnow%snowliq(kk,1:nsnow_max) = REAL(vsnow(kk)%hliq(:)*thousand)  ! amount of liq snow water
             ! amount of melted snow leaving bottom of snow pack (mm/dt)
             ssnow%smelt(kk)               = REAL(vsnow(kk)%Qmelt*thousand/dt)
             ssnow%nsnow(kk)               = vsnow(kk)%nsnow
             IF (SUM(ssnow%sdepth(kk,1:nsnow_max)) > zero) THEN
                ssnow%ssdnn(kk) = ssnow%snowd(kk)/SUM(ssnow%sdepth(kk,1:nsnow_max))
             ENDIF
          ENDIF
       ENDDO

       WHERE (err(1:mp) == 0) ssnow%isflag = 0

       ! snow output
       IF (1 == 0) THEN
          k = 1
          WRITE(340,"(100e16.6)") SUM(vsnow(k)%hsnow(1:vsnow(k)%nsnow)), vsnow(k)%tsn(1),SUM(vsnow(k)%hliq(1:vsnow(k)%nsnow)), &
               qprec_snow(k)*dt, vsnow(k)%Qmelt, qprec(k)*dt, &
               vsnow(k)%Qevap,vsnow(k)%Qvap,ssnow%albsoilsn(k,1), ssnow%albsoilsn(k,2), ssnow%sconds(k,1), &
               vsnow(k)%dens(1),SUM(vsnow(k)%depth(1:vsnow(k)%nsnow)), vsnow(k)%J, &
               vsnow(k)%MoistureFluxDivergence, vsnow(k)%FluxDivergence, vsnow(k)%dens(nsnow_max), vsnow(k)%tsn(nsnow_max), &
               qh(k,0), vmet(k)%rbh
       ENDIF

       WHERE (err(1:mp) == 0)
          canopy%ofes = canopy%fes
          ! Update total latent heat to reflect updated soil component:
          canopy%fe = canopy%fev + REAL(canopy%fes)
          ! Update total sensible heat to reflect updated soil component:
          canopy%fh = REAL(canopy%fhv) + canopy%fhs
       END WHERE

    ENDIF ! SEB only

  END SUBROUTINE sli_main
END MODULE sli_main_mod
