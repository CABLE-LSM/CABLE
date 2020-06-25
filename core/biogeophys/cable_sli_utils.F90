MODULE sli_utils

  USE cable_def_types_mod, ONLY: r_2, i_d
  USE cable_def_types_mod, ONLY: soil_parameter_type, veg_parameter_type
  USE sli_numbers,         ONLY: &
       experiment, &
       zero, half, one, two, e3, pi, &
       Tzero, gravity, Rgas, thousand, Mw, rlambda, Dva, &
       params, vars_aquifer, vars, rapointer, &
       rhow, lambdaf, lambdas, csice, cswat, cpa,&
       dpmaxr, solve_type, &
       gf, hmin, rhmin, dsmmax, rhocp, vars_snow, vars_met, &
       freezefac, ithermalcond, rmair, Mw
  USE cable_common_module, ONLY: cable_user
  use mo_utils,            only: eq
#ifdef __MPI__
  use mpi,                 only: MPI_Abort
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: dx, dxL, par, plit, sol, x ! soil water parameters
  PUBLIC :: bd, dis, isopar, isotype   ! soil solute parameters
  PUBLIC :: aquifer_props, flux, generic_thomas, getfluxes_vp, getheatfluxes, hyofh, hyofS, isosub ! subroutines
  PUBLIC :: litter_props, massman_sparse, potential_evap, setlitterpar, setpar, setpar_Loetsch, setsol, setx, tri
  PUBLIC :: csat, csoil, dthetalmaxdT, dthetalmaxdTh, esat, esat_ice, gammln, igamma, phi, rh0_sol, rtbis_rh0 ! functions
  PUBLIC :: slope_csat, slope_esat,slope_esat_ice, Sofh, Tfrz, thetalmax, weight, zerovars, Tthetalmax, Tfrozen
  PUBLIC :: rtbis_Tfrozen, GTfrozen, JSoilLayer, forcerestore, SEB
  PUBLIC :: spline_b, mean, nse

  REAL(r_2),        DIMENSION(:,:),   ALLOCATABLE :: dx
  REAL(r_2),        DIMENSION(:),     ALLOCATABLE :: dxL
  TYPE(params),     DIMENSION(:,:),   ALLOCATABLE :: par
  TYPE(params),     DIMENSION(:),     ALLOCATABLE :: plit
  TYPE(solve_type), DIMENSION(:),     ALLOCATABLE :: sol
  REAL(r_2),        DIMENSION(:,:),   ALLOCATABLE :: x

  ! MC solute not done yet
  !REAL(r_2),        DIMENSION(:,:),   ALLOCATABLE         :: bd
  !REAL(r_2),        DIMENSION(:,:),   ALLOCATABLE         :: dis
  !TYPE(rapointer),  DIMENSION(:,:,:), ALLOCATABLE         :: isopar
  !CHARACTER(LEN=2), DIMENSION(:,:,:), ALLOCATABLE         :: isotype
  REAL(r_2),        DIMENSION(:),   ALLOCATABLE         :: bd ! function solute not done yet
  REAL(r_2),        DIMENSION(:),   ALLOCATABLE         :: dis
  TYPE(rapointer),  DIMENSION(:,:), ALLOCATABLE         :: isopar
  CHARACTER(LEN=2), DIMENSION(:,:), ALLOCATABLE         :: isotype
#ifdef __MPI__
  integer :: ierr
#endif

  ! Subroutine interfaces

  INTERFACE generic_thomas
     MODULE PROCEDURE generic_thomas_1d
     MODULE PROCEDURE generic_thomas_2d
  END INTERFACE generic_thomas

  INTERFACE getfluxes_vp
     MODULE PROCEDURE getfluxes_vp_1d
     MODULE PROCEDURE getfluxes_vp_2d
  END INTERFACE getfluxes_vp

  INTERFACE getheatfluxes
     MODULE PROCEDURE getheatfluxes_1d
     MODULE PROCEDURE getheatfluxes_2d
  END INTERFACE getheatfluxes

  INTERFACE massman_sparse
     MODULE PROCEDURE massman_sparse_1d
     MODULE PROCEDURE massman_sparse_2d
  END INTERFACE massman_sparse

  INTERFACE tri
     MODULE PROCEDURE tri_1d
     MODULE PROCEDURE tri_2d
  END INTERFACE tri

  ! Function interfaces

  INTERFACE mean
     MODULE PROCEDURE mean_1d, mean_2d
  END INTERFACE mean

  INTERFACE nse
     MODULE PROCEDURE nse_1d, nse_2d
  END INTERFACE nse

  !P.J. Ross 2005-2007:
  ! This module implements Brooks-Corey (BC) soil water retention and conductivity
  ! functions. It illustrates the structure required without the complicating
  ! detail of more flexible hydraulic property functions.
  ! Definitions of public entities (see above for default values):
  ! params  - type for water parameters. Params the, thre (=the-thr), he, lam,
  !           Ke and eta are for the BC functions. Params he, Ke, KSe, phie
  !           and phiSe are values of variables h, K, KS, phi and phiS at
  !           saturation (denoted by the "e" for "air entry"), needed by
  !           module flow (MF).
  ! vars    - type for water variables used by MF and returned by subroutine
  !           hyofS (except for isat, which is 0 for unsaturated layers and 1
  !           for saturated layers).
  ! gf      - gravity factor for flow direction (usually 1 for straight down).
  ! hmin    - minimum matric head h (used by MF).
  ! par(:)  - hydraulic property params for soil types (used by MF).
  ! allo    - subroutine to allocate parameter storage.
  ! hyofS   - subroutine to get water variable from saturation S (where S<1).
  ! hyofh   - subroutine to get some water variables from h.
  ! Sofh    - subroutine to get S from h.
  ! hypar   - subroutine to set soil hydraulic params.
  ! weight  - subroutine to get gravity flow conductivity weight w.

  !V. Haverd 2008:
  ! modified by vh to include isothermal vapour conductivity in phi
  ! change units from cm, h to m, s
  ! includes thermal properties in soil parameters
  ! litter_props - subroutine to define litter properties

  !**********************************************************************************************************************

CONTAINS

  !**********************************************************************************************************************
  ! SUBROUTINES - alphabetical
  !**********************************************************************************************************************

  ELEMENTAL PURE SUBROUTINE aquifer_props(v_aquifer)

    IMPLICIT NONE

    TYPE(vars_aquifer), INTENT(INOUT) :: v_aquifer

    v_aquifer%zzero     = 53.43_r_2  ! water table depth corresponding to Wa = zero
    v_aquifer%Sy        = 0.2_r_2  ! specific yield of aquifer
    ! initialise water content of aquifer
    v_aquifer%Wa        = v_aquifer%Sy*(v_aquifer%zzero-max(v_aquifer%zdelta,v_aquifer%zsoil))
    v_aquifer%isat      = 0
    if (v_aquifer%zdelta <= v_aquifer%zsoil) v_aquifer%isat = 1
    v_aquifer%f         = 1.25_r_2 ! multiplier in exponent of Rs (m-1)
    v_aquifer%Rsmax     = 4.5e-7_r_2  ! maximum discharge rate from aquifer (ms-1)
    v_aquifer%discharge = v_aquifer%Rsmax*exp(-v_aquifer%f*v_aquifer%zdelta)

  END SUBROUTINE aquifer_props

  !**********************************************************************************************************************

  ELEMENTAL PURE SUBROUTINE flux(parin, v1, v2, dz, q, qya, qyb, qTa, qTb)
    ! VH modified 25/05/10 to include t-dep component of liquid flux (in frozen soil)

    IMPLICIT NONE

    TYPE(params), INTENT(IN)  :: parin
    REAL(r_2),    INTENT(IN)  :: dz
    REAL(r_2),    INTENT(OUT) :: q, qya, qyb, qTa, qTb
    TYPE(vars),   INTENT(IN)  :: v1, v2
    ! Gets flux and partial derivs for specified flow path.
    ! Definitions of arguments:
    ! j   - soil type no.
    ! v1  - water vars at upper end of path.
    ! v2  - ditto at lower end.
    ! dz  - length of path.
    ! q   - flux.
    ! qya - partial deriv of flux wrt S (if unsat) or phi (if sat) at upper end.
    ! qyb - ditto at lower end.
    REAL(r_2) :: w, rdz

    ! gf is gravity factor (0 to 1) assumed available in module
    if (gf < zero) then
       !if ((v1%isat /= 0 .and. v2%isat /= 0) .or. v1%h-gf*(-dz) >= parin%he) then
       if ((v1%isat /= 0 .and. v2%isat /= 0) .or. v1%h-gf*(-dz) >= v1%he) then
          w = one
       else
          w = weight(parin, v1%h, v1%K*v1%macropore_factor, v1%phi, -dz)
          w = one-w
       end if
    else
       !if ((v1%isat /= 0 .and. v2%isat /= 0) .or. v2%h-gf*dz >= parin%he) then
       if ((v1%isat /= 0 .and. v2%isat /= 0) .or. v2%h-gf*dz >= v2%he) then
          w = zero
       else
          w = weight(parin, v2%h, v2%K*v2%macropore_factor, v2%phi, dz)
       end if
    end if

    rdz = one/dz
    q   = (v1%phi-v2%phi)*rdz + gf*(w*v1%K*v1%macropore_factor+(one-w)*v2%K*v2%macropore_factor)
    if (v1%isat==0) then
       qya = v1%phiS*rdz + gf*w*v1%KS*v1%macropore_factor
       qTa = v1%phiT*rdz + gf*w*v1%KT*v1%macropore_factor
    else
       qya = rdz
       qTa = zero
    end if

    if (v2%isat==0) then
       qyb = -v2%phiS*rdz + gf*(one-w)*v2%KS*v2%macropore_factor
       qTb = -v2%phiT*rdz + gf*(one-w)*v2%KT*v2%macropore_factor
    else
       qyb = -rdz
       qTb = zero
    end if

  END SUBROUTINE flux

  !*****************************************************************************************

  SUBROUTINE forcerestore(Tg0, Rnet0, lE0, dlEdTg, Ta, Tbar, d1, rrc, lambda, &
       dt, iice, Tg, G, H, lE)

    ! method applicable to multilayer soil or soil/snow column
    ! derived using Eq's 3-12 in Hirota et al. JGR 2002
    IMPLICIT NONE
    REAL(r_2), INTENT(IN)   :: Tg0 ! ground surface temp of previous time-step [deg C]
    REAL(r_2), INTENT(IN)   :: Rnet0 ! Rnet at current time step, assuming Tg of previous time-step [W m-2]
    REAL(r_2), INTENT(IN)   :: lE0 ! latent heat flux at current time step, assuming Tg of previous time-step  [W m-2]
    REAL(r_2), INTENT(IN)   :: dlEdTg ! derivative of latent heat flux wrt Tg   [W m-2 K-1]
    REAL(r_2), INTENT(IN)   :: Ta ! air temperature [deg C]
    REAL(r_2), INTENT(IN)   :: Tbar ! temperature at diurnal damping depth
    REAL(r_2), INTENT(IN)   :: d1 ! diurnal damping depth (m) = sqrt(2*lambda/c/omega) (Hirota et al. eq 40)
    REAL(r_2), INTENT(IN)   :: rrc ! resistance to sensible heat and radiation transfer at ground/air interface [m-1 s]
    REAL(r_2), INTENT(IN)   :: lambda ! thermal conductivity of soil or snow at surface [W m-1 K-1]
    REAL(r_2), INTENT(IN)   :: dt ! time step [s]
    INTEGER, INTENT(IN):: iice ! top layer frozen (1) or not (0)
    REAL(r_2), INTENT(OUT) :: Tg, G, H, lE

    ! local variables
    REAL(r_2), PARAMETER :: tau1 = 86400._r_2 ! period of diurnal forcing (seconds)
    REAL(r_2) :: c1, c2, omega, a, b

    a = Rnet0 - rhocp/rrc*(Tg0-Ta) -lE0 ! G0 at current time step, assuming surface T of previous time step

    omega = 2._r_2*pi/tau1  ! diurnal forcing frequency (s-1)
    c1 = omega * d1 / lambda
    c2 = omega
    b = -rhocp/rrc  - dlEdTg

    !Tg = (Tg0 + dt*(c1*a+c2*Tbar)) / (one + c2 *dt)

    Tg = (Tg0 + dt*(c1*(a-b*Tg0)+c2*Tbar)) / (one + c2 *dt - b*c1*dt)

    if (iice.eq.1) then
       Tg = min(zero, Tg)
    endif

    G = a + b*(Tg-Tg0)
    H = rhocp/rrc*(Tg-Ta)
    lE = lE0 + dlEdTg*(Tg-Tg0)

  END SUBROUTINE forcerestore

  !*****************************************************************************************

  ! not used
  ! SUBROUTINE forcerestore_Deardorff(Tg0, Rnet0, lE0, dlEdTg, Ta, Tbar, d1, rrc, rhos, &
  !      cs, dt, iice, Tg, G, H, lE)

  !   IMPLICIT NONE
  !   REAL(r_2), INTENT(IN)   :: Tg0 ! ground surface temp of previous time-step [deg C]
  !   REAL(r_2), INTENT(IN)   :: Rnet0 ! Rnet at current time step, assuming Tg of previous time-step [W m-2]
  !   REAL(r_2), INTENT(IN)   :: lE0 ! latent heat flux at current time step, assuming Tg of previous time-step  [W m-2]
  !   REAL(r_2), INTENT(IN)   :: dlEdTg ! derivative of latent heat flux wrt Tg   [W m-2 K-1]
  !   REAL(r_2), INTENT(IN)   :: Ta ! air temperature [deg C]
  !   REAL(r_2), INTENT(IN)   :: Tbar ! temperature at diurnal damping depth
  !   REAL(r_2), INTENT(IN)   :: d1 ! diurnal damping depth (m)
  !   REAL(r_2), INTENT(IN)   :: rrc ! resistance to sensible heat and radiation transfer at ground/air interface [m-1 s]
  !   REAL(r_2), INTENT(IN)   :: rhos ! density of soil or snow at surface [kg m-3]
  !   REAL(r_2), INTENT(IN)   :: cs ! heat capacity of soil or snow at surface [J kg-1 K-1]
  !   REAL(r_2), INTENT(IN)   :: dt ! time step [s]
  !   INTEGER, INTENT(IN):: iice ! top layer frozen (1) or not (0)
  !   REAL(r_2), INTENT(OUT) :: Tg, G, H, lE

  !   ! local variables
  !   REAL(r_2),    PARAMETER :: c1 = 3.72_r_2 ! Deardorff JGR (1978)
  !   REAL(r_2),    PARAMETER :: c2 = 7.4_r_2 ! Deardorff JGR (1978)
  !   REAL(r_2),    PARAMETER :: tau1 = 86400._r_2 ! period of diurnal forcing (seconds)
  !   REAL(r_2) :: a, b

  !   a = Rnet0 - rhocp/rrc*(Tg0-Ta) -lE0
  !   b = -rhocp/rrc  - dlEdTg

  !   Tg = Tg0 + (c1*a/(rhos*cs*d1) - c2/tau1*(Tg0-Tbar))/ &
  !        (1._r_2/dt - b*c1/(rhos*cs*d1) + c2/tau1)

  !   if (iice.eq.1) then
  !      Tg = min(zero, Tg)
  !   endif

  !   G = a + b*(Tg-Tg0)
  !   H = rhocp/rrc*(Tg-Ta)
  !   lE = lE0 + dlEdTg*(Tg-Tg0)

  ! END SUBROUTINE forcerestore_Deardorff

  !**********************************************************************************************************************

  ! Surface Energy Balance
  SUBROUTINE SEB(n, par, vmet, vsnow, var, qprec, qprec_snow, dx, h0, Tsoil,  &
       Tsurface, G0, lE0, Epot, qsurface, qevap, qliq, qv, &
       qyb, qTb, qlyb, qvyb, qlTb, qvTb, qh, qadv, qhyb, qhTb, qadvyb, qadvTb)

    IMPLICIT NONE

    INTEGER(i_d),                    INTENT(IN) :: n
    TYPE(params),    DIMENSION(1:n), INTENT(IN) :: par
    TYPE(vars_met),                  INTENT(IN) :: vmet
    TYPE(vars_snow),                 INTENT(IN) :: vsnow
    TYPE(vars),      DIMENSION(1:n), INTENT(IN) :: var
    REAL(r_2),                       INTENT(IN) :: qprec
    REAL(r_2),                       INTENT(IN) :: qprec_snow
    REAL(r_2),       DIMENSION(1:n), INTENT(IN) :: dx
    REAL(r_2),                       INTENT(IN) :: h0
    REAL(r_2),       DIMENSION(1:n), INTENT(IN) :: Tsoil

    REAL(r_2),                       INTENT(OUT)           :: Tsurface, G0, lE0, Epot ! SEB
    REAL(r_2),                       INTENT(OUT)           :: qsurface          ! water flux into surface
    REAL(r_2),                       INTENT(OUT)           :: qevap             ! evaporative water flux
    ! liquid and vapour components of water flux from surface into soil
    REAL(r_2),                       INTENT(OUT)           :: qliq, qv
    ! derivatives of water fluxes wrt moisture and T
    REAL(r_2),                       INTENT(OUT)           :: qyb, qTb, qlyb, qvyb, qlTb, qvTb
    ! total and advective components of heat flux into surface
    REAL(r_2),                       INTENT(OUT)           :: qh, qadv
    ! derivatives of heat fluxes wrt moiture and T
    REAL(r_2),                       INTENT(OUT)           :: qhyb, qhTb, qadvyb, qadvTb

    ! local variables
    INTEGER(i_d) :: surface_case
    REAL(r_2) :: Tsurface_pot,  Hpot, Gpot, dEdrha, dEdTs, dEdTsoil, dGdTa, dGdTsoil
    REAL(r_2) :: E_vap, dE_vapdT1, E_liq
    REAL(r_2) :: Kmin, Khmin, phimin
    REAL(r_2) :: Tqw, dtqwdtb, rhocp1, cs
    LOGICAL   :: isEpot

    if (vsnow%nsnow.eq.0) surface_case = 1
    if (vsnow%nsnow>0) surface_case = 2

    select case (surface_case)
    case (1) ! no snow
       call potential_evap(vmet%Rn, vmet%rbh, vmet%rbw, vmet%Ta, vmet%rha, &
            Tsoil(1), var(1)%kth, half*dx(1)+h0, var(1)%lambdav, Tsurface_pot, Epot, Hpot, &
            Gpot, dEdrha, dEdTs, dEdTsoil, dGdTa, dGdTsoil)
       if (var(1)%iice.eq.1.and.Tsurface_pot> zero) then
          Tsurface_pot = 0.0_r_2
          Tsurface = 0.0_r_2

          Epot = (esat(Tsurface)*0.018_r_2/thousand/8.314_r_2/(vmet%Ta+Tzero)  - & ! m3 H2O (liq) m-3 (air)
               vmet%cva)*rhow*var(1)%lambdav/vmet%rbw
          dEdTsoil = zero
          dGdTsoil = zero
          Hpot = rhocp*(Tsurface - vmet%Ta)/vmet%rbh
          Gpot = vmet%Rn - Hpot - Epot
          dEdTs = zero
       endif

       if (var(1)%isat.eq.1) then  ! saturated surface =>. potential evporation
          isEpot = .true.
          Tsurface = Tsurface_pot
          lE0 = Epot
          G0 = Gpot
          E_vap = zero
          dE_vapdT1 = zero
          E_liq = lE0
       else ! unsaturated surface: finite vapour transfer; surface flux may be supply limited
          call hyofh(hmin, par(1)%lam, par(1)%eta, par(1)%Ke, par(1)%he, &
               Kmin, Khmin, phimin) ! get phi at hmin
          E_liq = ((var(1)%phi-phimin)/(half*dx(1))-var(1)%K)*thousand*var(1)%lambdav
          if (var(1)%Dv > 1.e-12_r_2) then
             ! E = (cs-ca)/rbw = E_liq + E_vap = E_liq + (c1-cs)/Dv/dx/2
             cs = ( E_liq/var(1)%lambdav + var(1)%rh*csat(Tsoil(1))*var(1)%Dv/(half*dx(1)) + &
                  vmet%cva*thousand/vmet%rbw ) &
                  / (one/vmet%rbw + var(1)%Dv/(half*dx(1)))
             E_vap = (var(1)%rh*csat(Tsoil(1)) - cs) * var(1)%Dv/(half*dx(1)) * var(1)%lambdav
             ! dE_vapdT1 = var(1)%Dv/(half*dx(1))*var(1)%rh*slope_csat(Tsoil(1)) &
             !      * (one - vmet%rbw*var(1)%Dv/(half*dx(1))/(one+vmet%rbw*var(1)%Dv/(half*dx(1)))) &
             !      * var(1)%lambdav
             dE_vapdT1 = var(1)%Dv/(half*dx(1)) * var(1)%rh*slope_csat(Tsoil(1)) * var(1)%lambdav
          else
             E_vap     = zero
             dE_vapdT1 = zero
          endif
          if (Epot <= (E_vap+E_liq)) then
             isEpot = .true.
             lE0 = Epot
             cs  = csat(Tsurface_pot)
             E_vap = (var(1)%rh*csat(Tsoil(1)) - cs) * var(1)%Dv/(half*dx(1)) * var(1)%lambdav
             E_liq = lE0 - E_vap
             dE_vapdT1 = var(1)%Dv/(half*dx(1)) * var(1)%rh*slope_csat(Tsoil(1)) * var(1)%lambdav
          else
             isEpot = .false.
             lE0 = E_vap+E_liq
             dEdTs = zero
          endif
          ! lE0 = min(Epot, E_vap+E_liq) ! analytic approximation (See Haverd et al. 2013, Appxx)
          ! if (Epot .gt. (E_vap+E_liq)) dEdTs = zero
          Tsurface = (-half*dx(1)*lE0 + half*dx(1)*vmet%Rn + &
               var(1)%kth*Tsoil(1) + half*dx(1)*(one/vmet%rbh*rhocp)*vmet%Ta) &
               /(var(1)%kth + half*dx(1)*(one/vmet%rbh*rhocp))

          G0       = var(1)%kth/(half*dx(1))*(Tsurface-Tsoil(1))
          dGdTsoil  =  -var(1)%kth/(half*dx(1))

          if ((var(1)%iice.eq.1) .and. (Tsurface>zero)) then
             Tsurface = 0.0_r_2
             rhocp1 = rmair*101325._r_2/rgas/(vmet%Ta+Tzero)*cpa
             G0 = vmet%Rn - rhocp1*(Tsurface - vmet%Ta)/vmet%rbh - lE0
             dGdTsoil = 0.0_r_2
          endif

       endif
       ! write(*,*) var(1)%phi, phimin

       qevap = lE0/(thousand*var(1)%lambdav)
       qsurface  = qprec + qprec_snow - qevap

       ! derivatives
       ! q refers to moisture in numerator
       ! qh refers to heat in numerator
       ! y refers to moisture in denominator
       ! T refers to temperature in denominator
       ! a refers to the layer above
       ! b refers to the layer below

       ! initialise derivatives to zero
       qyb = zero
       qTb = zero
       qhyb = zero
       qhTb = zero
       ! liquid and vapour fluxes
       qlyb = zero
       qvyb = zero

       ! potential evap independent of S(1), dependent on T1
       if ((var(1)%isat.eq.1 .or. vsnow%nsnow.gt.0)) then
          qyb  = zero
          qTb  = -dEdTsoil/(thousand*var(1)%lambdav)
          qliq = -Epot/(thousand*var(1)%lambdav)
          qv   = zero
          qlyb = zero
          qvyb = zero
          qlTb = qTb
          qvTb = zero
       elseif (isEpot) then
          qTb  = -dE_vapdT1/(thousand*var(1)%lambdav)
          qyb  = zero
          qliq = -E_liq/(thousand*var(1)%lambdav)
          qv   = -E_vap/(thousand*var(1)%lambdav)
          qlyb = zero
          qvyb = zero
          qlTb = qTb
          qvTb = zero
       else ! supply limited
          qTb = -dE_vapdT1/(thousand*var(1)%lambdav)
          qyb = -(var(1)%phiS/(half*dx(1)) - var(1)%KS)  !!vh!! include vapour component??
          qliq = -E_liq/(thousand*var(1)%lambdav)
          qv   = -E_vap/(thousand*var(1)%lambdav)
          qlyb = -(var(1)%phiS/(half*dx(1)) - var(1)%KS)
          qvyb = zero
          qlTb = zero
          qvTb = qTb
       endif
       ! end of partial derivative evaluation

       ! advective component of heat flux
       qadv = rhow*cswat*qprec*(vmet%Ta) + rhow*csice*qprec_snow*(min(vmet%Ta,zero)) &
            - rhow*qprec_snow*lambdaf

       Tqw  = merge(vmet%Ta, Tsoil(1), (-qevap)>zero)
       dTqwdTb = merge(zero, one, (-qevap)>zero)
       qadv = qadv + rhow*cswat*Tqw*(-qevap)

       qadvTb = dTqwdTb + rhow*cswat*Tqw*qTb
       qadvyb =  rhow*cswat*qyb*Tqw
       qadvyb = 0


       qadvTb = 0 ! test vh!

       qh = qadv + G0
       qhyb = qadvyb
       qhTb = dGdTsoil + qadvTb

    case (2) !dedicated snow layer
       ! NB Only longwave component of net radiation directly affects SEB: sw component is absorbed internally

       ! SEB at snow/air interface
       if (vsnow%hliq(1)>zero) then
          Tsurface = 0.0_r_2
          !Epot = (esat(Tsurface)*0.018_r_2/thousand/8.314_r_2/(vmet%Ta+Tzero)  - & ! m3 H2O (liq) m-3 (air)
          !            vmet%cva)*rhow*rlambda/vmet%rbw !!vh check this !!
          Epot = (csat(Tsurface)/thousand - vmet%cva)/vmet%rbw *rlambda*rhow  ! m3 H2O (liq) m-3 (air) -> W/m2
          ! write(*,*) "Epot", vmet%rha, vmet%Ta,
          !            esat(Tsurface)*0.018_r_2/thousand/8.314_r_2/(vmet%Ta+Tzero)*rhow*rlambda/vmet%rbw, &
          !vmet%cva*rhow*rlambda/vmet%rbw, Epot
          dEdTsoil = zero
          dGdTsoil = zero
          rhocp1 = rmair*101325._r_2/rgas/(vmet%Ta+Tzero)*cpa
          Hpot = rhocp1*(Tsurface - vmet%Ta)/vmet%rbh
          Gpot = vmet%Rn-vmet%Rnsw - Hpot - Epot
          dEdTs = zero
          qevap = Epot/(rhow*rlambda)
          qTb = zero
          !     write(*,*) "Epot2", Tsurface, vmet%Ta, Epot, Hpot, vmet%rrc, rhocp1
       else
          !! vh !! use max snow depth of 20 cm in this calculation to avoid huge resistances
          !! leading to large negative surface temperatures when snow-pack is thick and
          !! Rn is large and negative (~-100 Wm-2)
          call potential_evap(vmet%Rn-vmet%Rnsw, vmet%rbh, vmet%rbw, vmet%Ta, vmet%rha, &
               vsnow%tsn(1), max(vsnow%kth(1),0.1_r_2), half*min(vsnow%depth(1),0.05_r_2), &
               lambdas, Tsurface, Epot, Hpot, &
               Gpot, dEdrha, dEdTs, dEdTsoil, dGdTa, dGdTsoil,iice=.TRUE.)
          if (Tsurface > zero) then ! temperature of frozen surface must be <= zero
             Tsurface = 0.0_r_2
             Epot = (esat(Tsurface)*0.018_r_2/thousand/8.314_r_2/(vmet%Ta+Tzero)  - & ! m3 H2O (liq) m-3 (air)
                  vmet%cva)*rhow*lambdas/vmet%rbw
             dEdTsoil = zero
             dGdTsoil = zero
             Hpot = rhocp*(Tsurface - vmet%Ta)/vmet%rbh
             Gpot = vmet%Rn-vmet%Rnsw - Hpot - Epot
             dEdTs= zero
             !   write(*,*) "Epot3", Tsurface, vmet%Ta, Epot, Hpot, vmet%rbh
          endif
!!$          elseif (abs(Tsurface - vmet%Ta).gt. 20) then
!!$             Tsurface = min(vmet%Ta, 0.0)
!!$             Epot = (esat(Tsurface)*0.018_r_2/thousand/8.314_r_2/(vmet%Ta+Tzero)  - & ! m3 H2O (liq) m-3 (air)
!!$                  vmet%cva)*rhow*lambdas/vmet%rbw
!!$             dEdTsoil = zero
!!$             dGdTsoil = zero
!!$             Hpot = rhocp*(Tsurface - vmet%Ta)/vmet%rbh
!!$             Gpot = vmet%Rn-vmet%Rnsw - Hpot - Epot
!!$             dEdTs= zero
!!$           endif
          qevap = Epot/(rhow*lambdas)
          qTb = -dEdTsoil/(thousand*lambdas)
       endif
       lE0 = Epot
       ! moisture flux at air/snow interface
       qsurface = qprec_snow+qprec-qevap
       qyb = zero

       ! conductive heat flux at air/snow interface
       G0 = Gpot
       qhTb = dGdTsoil
       if (vsnow%hliq(1)>zero) then
          qhTb = zero
          qTb = zero
       endif

       ! advective heat flux at air/snow interface
       qadv = rhow*(qprec_snow)*(csice*(min(vmet%Ta,zero))-lambdaf) + &
            rhow*(qprec)*cswat*(max(vmet%Ta,zero))
       Tqw  = merge(vmet%Ta, vsnow%tsn(1), -qevap>zero)
       dTqwdTb = merge(zero,one, -qevap>zero)
       if (vsnow%hliq(vsnow%nsnow)>zero) then
          qadv = qadv + rhow*(-qevap)*cswat*Tqw
          qadvTb = zero
       else
          qadv = qadv + rhow*(-qevap)*cswat*Tqw
          qadvTb = rhow*cswat*(-qevap)*dTqwdTb  +  rhow*cswat*Tqw*qTb
       endif

       ! add sw energy absorption to G0
       G0 = G0 + vmet%Rnsw
       qh     = qadv + G0
       qhyb   = zero
       qadvyb = zero
       qhyb   = qhyb + qadvyb
       qhTb   = qhTb + qadvTb


       qv   = -qevap
       qliq = zero
       qvyb = qyb
       qvTb = qTb
       qlyb = zero
       qlTb = zero

    end select ! surface_case

    ! finished all the surfaces

  END SUBROUTINE SEB

  !**********************************************************************************************************************

  ! not used
  ! ! Surface Energy Balance
  ! SUBROUTINE SEB_FR(n, par, vmet, vsnow, var, qprec, qprec_snow, &
  !      nsteps, dx, h0, Tsoil, dt, Tsurface0, &
  !      Tsurface, G0, lE0, TsurfaceFR, G0FR, lEFR, HFR, qsurface, qevap, qliq, qv, &
  !      qyb, qTb, qlyb, qvyb, qlTb, qvTb, qh, qadv, qhyb, qhTb, qadvyb, qadvTb, irec)

  !   IMPLICIT NONE

  !   INTEGER(i_d), INTENT(IN)              :: n
  !   TYPE(params),      DIMENSION(1:n), INTENT(IN)  :: par
  !   TYPE(vars_met),        INTENT(IN)           :: vmet
  !   TYPE(vars_snow),   INTENT(IN)           :: vsnow
  !   TYPE(vars),      DIMENSION(1:n),   INTENT(IN)           :: var
  !   REAL(r_2),   INTENT(IN)              :: qprec
  !   REAL(r_2),   INTENT(IN)              :: qprec_snow
  !   INTEGER(i_d), INTENT(IN)              ::  nsteps, irec
  !   REAL(r_2),      DIMENSION(1:n),   INTENT(IN)              :: dx
  !   REAL(r_2),          INTENT(IN)           :: h0
  !   REAL(r_2),      DIMENSION(1:n),   INTENT(IN)           :: Tsoil
  !   REAL(r_2),   INTENT(IN)              :: dt
  !   REAL(r_2),  INTENT(IN)           :: Tsurface0

  !   REAL(r_2),  INTENT(OUT)           :: Tsurface, G0, lE0  ! SEB (subdiurnal, uses T in top layer)
  !   REAL(r_2),  INTENT(OUT)           :: TsurfaceFR, G0FR, lEFR, HFR  ! SEB (Force-Restore)
  !   REAL(r_2),  INTENT(OUT)           :: qsurface ! water flux into surface
  !   REAL(r_2),  INTENT(OUT)           :: qevap ! evaporative water flux
  !   REAL(r_2),  INTENT(OUT)           :: qliq, qv ! liquid and vapour components of water flux from surface into soil
  !   REAL(r_2),  INTENT(OUT)           :: qyb, qTb, qlyb, qvyb, qlTb, qvTb ! derivatives of water fluxes wrt moisture and T
  !   REAL(r_2),  INTENT(OUT)           :: qh, qadv ! total and advective components of heat flux into surface
  !   REAL(r_2),  INTENT(OUT)           :: qhyb, qhTb, qadvyb, qadvTb ! derivatives of heat fluxes wrt moiture and T

  !   ! local variables
  !   INTEGER(i_d) :: surface_case, j
  !   REAL(r_2) :: Tsurface_pot, Epot, Hpot, Gpot, dEdrha, dEdTs, dEdTsoil, dGdTa, dGdTsoil
  !   REAL(r_2) :: E_vap, dE_vapdT1, E_liq
  !   REAL(r_2) :: Kmin, Khmin, phimin
  !   REAL(r_2) :: Tqw, dtqwdtb, d1, tmp1d2, Tbar, f, csnow

  !   if (vsnow%nsnow.eq.0) surface_case = 1
  !   if (vsnow%nsnow>0) surface_case = 2
  !   select case (surface_case)
  !   case (1)
  !      call potential_evap(vmet%Rn, vmet%rrc, vmet%rbw, vmet%Ta, vmet%rha, &
  !           Tsoil(1), var(1)%kth, half*dx(1)+h0, var(1)%lambdav, Tsurface_pot, Epot, Hpot, &
  !           Gpot, dEdrha, dEdTs, dEdTsoil, dGdTa, dGdTsoil)

  !      if (var(1)%iice.eq.1.and.Tsurface_pot> zero) then
  !         Tsurface_pot = 0.0_r_2
  !         Epot = (esat(Tsurface)*0.018_r_2/thousand/8.314_r_2/(vmet%Ta+Tzero)  - & ! m3 H2O (liq) m-3 (air)
  !              vmet%cva)*rhow*var(1)%lambdav/vmet%rbw
  !         dEdTsoil = zero
  !         dGdTsoil = zero
  !         Hpot = rhocp*(Tsurface - vmet%Ta)/vmet%rrc
  !         Gpot = vmet%Rn - Hpot - Epot
  !         dEdTs = zero
  !      endif

  !      if (var(1)%isat.eq.1) then  ! saturated surface =>. potential evporation
  !         Tsurface = Tsurface_pot
  !         lE0 = Epot
  !         G0 = Gpot
  !         E_vap = zero
  !         dE_vapdT1 = zero
  !         E_liq = lE0
  !      else ! unsaturated surface: finite vapour transfer; surface flux may be supply limited
  !         if (var(1)%Dv > 1.e-12_r_2) then
  !            E_vap = (var(1)%rh*csat(Tsoil(1))-vmet%cva*thousand)/(vmet%rbw + half*dx(1)/var(1)%Dv)*var(1)%lambdav
  !            dE_vapdT1 = (var(1)%rh*slope_csat(Tsoil(1)))/(vmet%rbw + half*dx(1)/var(1)%Dv)*var(1)%lambdav
  !         else
  !            E_vap = zero
  !            dE_vapdT1 = zero
  !         endif
  !         call hyofh(hmin, par(1)%lam, par(1)%eta, par(1)%Ke, par(1)%he, &
  !              Kmin, Khmin, phimin) ! get phi at hmin
  !         E_liq = ((var(1)%phi-phimin)/(half*dx(1))-var(1)%K)*thousand*var(1)%lambdav
  !         lE0 = min(Epot,E_vap+ E_liq) ! analytic approximation (See Haverd et al. 2013, Appxx)
  !         if (Epot.gt.(E_vap+ E_liq)) dEdTs = zero
  !         Tsurface = (-half*dx(1)*lE0 + half*dx(1)*vmet%Rn + &
  !              var(1)%kth*Tsoil(1) + half*dx(1)*(one/vmet%rrc*rhocp)*vmet%Ta) &
  !              /(var(1)%kth + half*dx(1)*(one/vmet%rrc*rhocp))
  !         if (var(1)%iice.eq.1.and.Tsurface> zero) Tsurface = 0.0_r_2
  !         G0       = var(1)%kth/(half*dx(1))*(Tsurface-Tsoil(1))
  !         dGdTsoil  =  -var(1)%kth/(half*dx(1))
  !      endif
  !      ! write(*,*) var(1)%phi, phimin

  !      qevap = lE0/(thousand*var(1)%lambdav)
  !      qsurface  = qprec + qprec_snow - qevap

  !      ! derivatives
  !      ! q refers to moisture in numerator
  !      ! qh refers to heat in numerator
  !      ! y refers to moisture in denominator
  !      ! T refers to temperature in denominator
  !      ! a refers to the layer above
  !      ! b refers to the layer below

  !      ! initialise derivatives to zero
  !      qyb = zero
  !      qTb = zero
  !      qhyb = zero
  !      qhTb = zero
  !      ! liquid and vapour fluxes
  !      qlyb = zero
  !      qvyb = zero

  !      ! potential evap independent of S(1), dependent on T1
  !      if ((Epot<=(E_vap+ E_liq)).or.(var(1)%isat.eq.1.or.vsnow%nsnow.gt.0)) then
  !         qyb = zero
  !         qTb = -dEdTsoil/(thousand*var(1)%lambdav)
  !         qliq = -Epot/(thousand*var(1)%lambdav)
  !         qv = zero
  !         qlyb = zero
  !         qvyb = zero
  !         qlTb = qTb
  !         qvTb = zero
  !      else ! supply limited
  !         qTb = -dE_vapdT1/(thousand*var(1)%lambdav)
  !         qyb = -(var(1)%phiS/(half*dx(1)) - var(1)%KS)  !!vh!! include vapour component??
  !         qliq = -E_liq/(thousand*var(1)%lambdav)
  !         qv = -E_vap/(thousand*var(1)%lambdav)
  !         qlyb = -(var(1)%phiS/(half*dx(1)) - var(1)%KS)
  !         qvyb = zero
  !         qlTb = zero
  !         qvTb = qTb
  !      endif
  !      ! end of partial derivative evaluation

  !      ! advective component of heat flux
  !      qadv = rhow*cswat*qprec*(vmet%Ta) + rhow*csice*qprec_snow*(min(vmet%Ta,zero)) &
  !           - rhow*qprec_snow*lambdaf

  !      Tqw  = merge(vmet%Ta, Tsoil(1), (-qevap)>zero)
  !      dTqwdTb = merge(zero, one, (-qevap)>zero)
  !      qadv = qadv + rhow*cswat*Tqw*(-qevap)

  !      qadvTb = dTqwdTb + rhow*cswat*Tqw*qTb
  !      qadvyb =  rhow*cswat*qyb*Tqw
  !      qh = qadv + G0
  !      qhyb = qadvyb
  !      qhTb = dGdTsoil + qadvTb

  !      if (nsteps.eq.-1) then
  !         j = 1
  !         tmp1d2 = zero
  !         do while (tmp1d2.lt.one)
  !            ! d1 = (var(j)%kth/(var(j)%csoileff)*86400./pi)**0.5
  !            d1 = sqrt(var(j)%kth/(var(j)%csoileff)*86400./pi)
  !            tmp1d2 = tmp1d2 + dx(j)/d1
  !            Tbar = Tsoil(j)
  !            j = j+1
  !         enddo
  !         j = j-1 ! integer corresponding to lowest layer contributing to soil col above damping depth
  !         ! damping depth = :sum(dx(k)) (k=1,j-1) +f dx(j)
  !         if (j.gt.1) then
  !            f = (one - (tmp1d2 -dx(j)/d1))*d1/dx(j) ! fraction of lowest layer contributing to soil col above damping depth
  !            d1 = sum(dx(1:j-1)) + f*dx(j)  ! multilayer damping depth
  !         else
  !            f = d1/dx(1)
  !            d1 = f*dx(1)
  !         endif
  !         write(21,"(2i8,3e16.6, i8, 8e16.6)") irec, j, f, d1, Tbar, vsnow%nsnow,Tsurface0, vmet%Rn, lE0, dEdTs, &
  !              vsnow%depth(1)!, vsnow%depth(2)

  !         call forcerestore(Tsurface0, vmet%Rn, lE0, dEdTs, vmet%Ta, &
  !              Tbar, d1, vmet%rrc, var(1)%kth, &
  !              var(1)%csoileff, dt, var(1)%iice, TsurfaceFR, G0FR, HFR, lEFR)
  !         ! write(37,"(16e16.6)") Tsurface0, vmet%Rn, lE0, dEdTs, vmet%Ta, &
  !         !                       Tbar, d1, vmet%rrc, var(1)%csoileff, dt, real(var(1)%iice), TsurfaceFR, G0FR, HFR, lEFR
  !      endif

  !   case (2) !dedicated snow layer

  !      ! SEB at snow/air interface
  !      if (vsnow%hliq(1)>zero) then
  !         Tsurface = 0.0_r_2
  !         Epot = (esat(Tsurface)*0.018_r_2/thousand/8.314_r_2/(vmet%Ta+Tzero)  - & ! m3 H2O (liq) m-3 (air)
  !              vmet%cva)*rhow*lambdaf/vmet%rbw
  !         dEdTsoil = zero
  !         dGdTsoil = zero
  !         Hpot = rhocp*(Tsurface - vmet%Ta)/vmet%rrc
  !         Gpot = vmet%Rn - Hpot - Epot
  !         dEdTs = zero
  !         qevap = Epot/(thousand*lambdaf)
  !         qTb = zero
  !      else

  !         call potential_evap(vmet%Rn, vmet%rrc, vmet%rbw, vmet%Ta, vmet%rha, &
  !              vsnow%tsn(1), vsnow%kth(1), half*vsnow%depth(1), &
  !              lambdas, Tsurface, Epot, Hpot, &
  !              Gpot, dEdrha, dEdTs, dEdTsoil, dGdTa, dGdTsoil)

  !         if (Tsurface > zero) then ! temperature of frozen surface must be <= zero
  !            Tsurface = 0.0_r_2
  !            Epot = (esat(Tsurface)*0.018_r_2/thousand/8.314_r_2/(vmet%Ta+Tzero)  - & ! m3 H2O (liq) m-3 (air)
  !                 vmet%cva)*rhow*lambdas/vmet%rbw
  !            dEdTsoil = zero
  !            dGdTsoil = zero
  !            Hpot = rhocp*(Tsurface - vmet%Ta)/vmet%rrc
  !            Gpot = vmet%Rn - Hpot - Epot
  !            dEdTs= zero
  !         endif
  !         qevap = Epot/(thousand*lambdas)
  !         qTb = -dEdTsoil/(thousand*lambdas)
  !      endif

  !      ! moisture flux at air/snow interface
  !      qsurface = qprec_snow+qprec-qevap
  !      qyb = zero

  !      ! conductive heat flux at air/snow interface
  !      G0 = Gpot
  !      qhTb = dGdTsoil
  !      if (vsnow%hliq(1)>zero) then
  !         qhTb = zero
  !         qTb = zero
  !      endif

  !      ! advective heat flux at air/snow interface
  !      qadv = rhow*(qprec_snow)*(csice*(min(vmet%Ta,zero))-lambdaf) + &
  !           rhow*(qprec)*cswat*(vmet%Ta)
  !      Tqw  = merge(vmet%Ta, vsnow%tsn(1), -qevap>zero)
  !      dTqwdTb = merge(zero,one, -qevap>zero)
  !      if (vsnow%hliq(vsnow%nsnow)>zero) then
  !         qadv = qadv + rhow*(-qevap)*cswat*Tqw
  !         qadvTb = zero
  !      else
  !         qadv = qadv + rhow*(-qevap)*cswat*Tqw
  !         qadvTb = rhow*cswat*(-qevap)*dTqwdTb  +  rhow*cswat*Tqw*qTb
  !      endif
  !      qh = qadv + G0
  !      qhyb = qhyb +  qadvyb
  !      qhTb = qhTb + qadvTb

  !      if (nsteps.eq.-1) then
  !         j = -vsnow%nsnow + 1
  !         tmp1d2 = zero
  !         do while (tmp1d2.lt.one)
  !            if (j.lt.1) then
  !               !write(*,*) 'chk1', -j+vsnow%nsnow, j, vsnow%nsnow
  !               csnow = (csice*(vsnow%hsnow(j+vsnow%nsnow)-vsnow%hliq(j+vsnow%nsnow))+cswat*vsnow%hliq(j+vsnow%nsnow))/ &
  !                    vsnow%depth(j+vsnow%nsnow)*rhow
  !               ! d1 = (vsnow%kth(j+vsnow%nsnow)/(csnow)*86400./pi)**0.5
  !               d1 = sqrt(vsnow%kth(j+vsnow%nsnow)/(csnow)*86400._r_2/pi)
  !               tmp1d2 = tmp1d2 + vsnow%depth(j+vsnow%nsnow)/d1  ! check snow index here!
  !               Tbar = vsnow%tsn(j+vsnow%nsnow)
  !            else
  !               ! d1 = (var(j)%kth/(var(j)%csoileff)*86400./pi)**0.5
  !               d1 = sqrt(var(j)%kth/(var(j)%csoileff)*86400._r_2/pi)
  !               tmp1d2 = tmp1d2 + dx(j)/d1
  !               Tbar = Tsoil(j)
  !            endif

  !            j = j+1
  !         enddo
  !         j = j-1 ! integer corresponding to lowest layer contributing to soil col above damping depth
  !         ! damping depth = sum(dx(k)) (k=1,j-1) +f dx(j)
  !         if (j.gt.(-vsnow%nsnow + 1).and.j.le.0) then ! damping depth within snowpack
  !            ! fraction of lowest layer contributing to soil col above damping depth
  !            f = (one - (tmp1d2 -vsnow%depth(j+vsnow%nsnow)/d1))*d1/vsnow%depth(j+vsnow%nsnow)
  !            d1 = sum(vsnow%depth(1:j+vsnow%nsnow-1)) + f*vsnow%depth(j+vsnow%nsnow) ! multilayer damping depth
  !         elseif (j.gt.0) then ! damping depth within soil column
  !            f = (one - (tmp1d2 -dx(j)/d1))*d1/dx(j) ! fraction of lowest layer contributing to soil col above damping depth
  !            if (j.gt.1) then
  !               d1 = sum(vsnow%depth(1:vsnow%nsnow))+sum(dx(1:j-1)) + f*dx(j)  ! multilayer damping depth
  !            else
  !               d1 = sum(vsnow%depth(1:vsnow%nsnow)) + f*dx(j)  ! multilayer damping depth
  !            endif
  !         elseif (j.eq.-vsnow%nsnow + 1) then ! damping depth within top layer of snowpack
  !            f = d1/vsnow%depth(j+vsnow%nsnow)
  !            d1 = f*vsnow%depth(j+vsnow%nsnow)
  !         endif

  !         csnow = (csice*(vsnow%hsnow(1)-vsnow%hliq(1))+cswat*vsnow%hliq(1))/vsnow%depth(1)*rhow
  !         write(21,"(2i8,3e16.6,i8,8e16.6)") irec,j, f, d1, Tbar, vsnow%nsnow, Tsurface0, vmet%Rn, lE0, dEdTs, &
  !              vsnow%depth(1)!, vsnow%depth(2)

  !         call forcerestore(Tsurface0, vmet%Rn, lE0, dEdTs, vmet%Ta, &
  !              Tbar, d1, vmet%rrc, vsnow%kth(1), &
  !              csnow, dt, 1, TsurfaceFR, G0FR, HFR, lEFR)
  !         ! write(37,"(16e16.6)") Tsurface0, vmet%Rn, lE0, dEdTs, vmet%Ta, &
  !         !                       Tbar, d1, vmet%rrc, var(1)%csoileff, dt, real(var(1)%iice), TsurfaceFR, G0FR, HFR, lEFR
  !      endif

  !   end select ! surface_case
  !   ! finished all the surfaces

  ! END SUBROUTINE SEB_FR

  !**********************************************************************************************************************

  SUBROUTINE generic_thomas_1d(n,A,B,C,r,u,err)

    USE sli_numbers,       ONLY: one

    IMPLICIT NONE

    ! in/out
    INTEGER(i_d),                      INTENT(IN)  :: n
    REAL(r_2), DIMENSION(1:n,1:2,1:2), INTENT(IN)  :: A, B, C
    REAL(r_2), DIMENSION(1:n,1:2),     INTENT(IN)  :: r
    REAL(r_2), DIMENSION(1:n,1:2),     INTENT(OUT) :: u
    INTEGER(i_d), OPTIONAL,            INTENT(OUT) :: err ! 0: no error; >0: error
    ! local
    REAL(r_2), DIMENSION(1:n,1:2,1:2) :: G
    REAL(r_2), DIMENSION(1:2,1:2)     :: bet
    REAL(r_2), DIMENSION(1:2)         :: d
    REAL(r_2)                         :: detbet, detbet1
    INTEGER(i_d)                      :: j

    if (present(err)) err = 0
    ! j=1
    bet(1:2,1:2) = B(1,1:2,1:2)
    detbet       = bet(1,1)*bet(2,2) - bet(1,2)*bet(2,1)
    if (abs(detbet) < epsilon(detbet)) then
       write(*,*) 'generic_thomas_1d error1: det = 0'
       if (present(err)) then
          err = 1
          return
       else
#ifdef __MPI__
          call MPI_Abort(0, 144, ierr) ! Do not know comm nor rank here
#else
          stop 144
#endif
       endif
    endif
    detbet1 = one/detbet
    d(1:2)  = r(1,1:2)
    u(1,1)  = (bet(2,2)*d(1) - bet(1,2)*d(2))*detbet1
    u(1,2)  = (bet(1,1)*d(2) - bet(2,1)*d(1))*detbet1
    ! j=2, n
    do j=2, n
       d(1:2)       = C(j-1,1:2,1)
       G(j,1,1)     = (bet(2,2)*d(1) - bet(1,2)*d(2))*detbet1
       G(j,2,1)     = (bet(1,1)*d(2) - bet(2,1)*d(1))*detbet1
       d(1:2)       = C(j-1,1:2,2)
       G(j,1,2)     = (bet(2,2)*d(1) - bet(1,2)*d(2))*detbet1
       G(j,2,2)     = (bet(1,1)*d(2) - bet(2,1)*d(1))*detbet1
       bet(1:2,1:2) = B(j,1:2,1:2) - matmul(A(j,1:2,1:2),G(j,1:2,1:2))
       detbet       = bet(1,1)*bet(2,2) - bet(1,2)*bet(2,1)
       if (abs(detbet) < epsilon(detbet)) then
          write(*,*) 'generic_thomas_1d error2: det = 0 at j=', j
          if (present(err)) then
             err = 1
             return
          else
#ifdef __MPI__
             call MPI_Abort(0, 145, ierr) ! Do not know comm nor rank here
#else
             stop 145
#endif
          endif
       endif
       detbet1      = one/detbet
       d(1:2)       = r(j,1:2) - matmul(A(j,1:2,1:2),u(j-1,1:2))
       u(j,1)       = (bet(2,2)*d(1) - bet(1,2)*d(2))*detbet1
       u(j,2)       = (bet(1,1)*d(2) - bet(2,1)*d(1))*detbet1
    end do
    ! back substitution
    do j=n-1, 1, -1
       u(j,1:2) = u(j,1:2) - matmul(G(j+1,1:2,1:2),u(j+1,1:2))
    end do
    !
  END SUBROUTINE generic_thomas_1d

  SUBROUTINE generic_thomas_2d(mp, n,A,B,C,r,u,err)

    USE sli_numbers,       ONLY: one

    IMPLICIT NONE

    ! in/out
    INTEGER(i_d),                           INTENT(IN)  :: mp
    INTEGER(i_d),                           INTENT(IN)  :: n
    REAL(r_2), DIMENSION(1:mp,1:n,1:2,1:2), INTENT(IN)  :: A, B, C
    REAL(r_2), DIMENSION(1:mp,1:n,1:2),     INTENT(IN)  :: r
    REAL(r_2), DIMENSION(1:mp,1:n,1:2),     INTENT(OUT) :: u
    INTEGER(i_d), OPTIONAL,                 INTENT(OUT) :: err ! 0: no error; >0: error
    ! local
    REAL(r_2), DIMENSION(1:mp,1:n,1:2,1:2) :: G
    REAL(r_2), DIMENSION(1:mp,1:2,1:2)     :: bet
    REAL(r_2), DIMENSION(1:mp,1:2)         :: d
    REAL(r_2), DIMENSION(1:mp)             :: detbet, detbet1
    INTEGER(i_d)                           :: j
    REAL(r_2), DIMENSION(1:mp,1:2)         :: tmp1d
    REAL(r_2), DIMENSION(1:mp,1:2,1:2)     :: tmp2d

    if (present(err)) err = 0
    ! j=1
    bet(1:mp,1:2,1:2) = B(1:mp,1,1:2,1:2)
    detbet(1:mp)      = bet(1:mp,1,1)*bet(1:mp,2,2) - bet(1:mp,1,2)*bet(1:mp,2,1)
    if (any(abs(detbet(1:mp)) < epsilon(detbet))) then
       write(*,*) 'generic_thomas_2d error1: det = 0'
       if (present(err)) then
          err = 1
          return
       else
#ifdef __MPI__
          call MPI_Abort(0, 146, ierr) ! Do not know comm nor rank here
#else
          stop 146
#endif
       endif
    endif
    detbet1(1:mp) = one/detbet(1:mp)
    d(1:mp,1:2)   = r(1:mp,1,1:2)
    u(1:mp,1,1)   = (bet(1:mp,2,2)*d(1:mp,1) - bet(1:mp,1,2)*d(1:mp,2))*detbet1(1:mp)
    u(1:mp,1,2)   = (bet(1:mp,1,1)*d(1:mp,2) - bet(1:mp,2,1)*d(1:mp,1))*detbet1(1:mp)
    ! j=2, n
    do j=2, n
       d(1:mp,1:2)       = C(1:mp,j-1,1:2,1)
       G(1:mp,j,1,1)     = (bet(1:mp,2,2)*d(1:mp,1) - bet(1:mp,1,2)*d(1:mp,2))*detbet1(1:mp)
       G(1:mp,j,2,1)     = (bet(1:mp,1,1)*d(1:mp,2) - bet(1:mp,2,1)*d(1:mp,1))*detbet1(1:mp)
       d(1:mp,1:2)       = C(1:mp,j-1,1:2,2)
       G(1:mp,j,1,2)     = (bet(1:mp,2,2)*d(1:mp,1) - bet(1:mp,1,2)*d(1:mp,2))*detbet1(1:mp)
       G(1:mp,j,2,2)     = (bet(1:mp,1,1)*d(1:mp,2) - bet(1:mp,2,1)*d(1:mp,1))*detbet1(1:mp)
       tmp2d(1:mp,1,1)   = A(1:mp,j,1,1)*G(1:mp,j,1,1) + A(1:mp,j,1,2)*G(1:mp,j,2,1)
       tmp2d(1:mp,1,2)   = A(1:mp,j,1,1)*G(1:mp,j,1,2) + A(1:mp,j,1,2)*G(1:mp,j,2,2)
       tmp2d(1:mp,2,1)   = A(1:mp,j,2,1)*G(1:mp,j,1,1) + A(1:mp,j,2,2)*G(1:mp,j,2,1)
       tmp2d(1:mp,2,2)   = A(1:mp,j,2,1)*G(1:mp,j,1,2) + A(1:mp,j,2,2)*G(1:mp,j,2,2)
       bet(1:mp,1:2,1:2) = B(1:mp,j,1:2,1:2) - tmp2d(1:mp,1:2,1:2)
       detbet(1:mp)      = bet(1:mp,1,1)*bet(1:mp,2,2) - bet(1:mp,1,2)*bet(1:mp,2,1)
       if (any(abs(detbet(1:mp)) < epsilon(detbet))) then
          write(*,*) 'generic_thomas_2d error2: det = 0 at j=', j
          if (present(err)) then
             err = 1
             return
          else
#ifdef __MPI__
             call MPI_Abort(0, 147, ierr) ! Do not know comm nor rank here
#else
             stop 147
#endif
          endif
       endif
       detbet1(1:mp) = one/detbet(1:mp)
       tmp1d(1:mp,1) = A(1:mp,j,1,1)*u(1:mp,j-1,1) + A(1:mp,j,1,2)*u(1:mp,j-1,2)
       tmp1d(1:mp,2) = A(1:mp,j,2,1)*u(1:mp,j-1,1) + A(1:mp,j,2,2)*u(1:mp,j-1,2)
       d(1:mp,1:2)   = r(1:mp,j,1:2) - tmp1d(1:mp,1:2)
       u(1:mp,j,1)   = (bet(1:mp,2,2)*d(1:mp,1) - bet(1:mp,1,2)*d(1:mp,2))*detbet1(1:mp)
       u(1:mp,j,2)   = (bet(1:mp,1,1)*d(1:mp,2) - bet(1:mp,2,1)*d(1:mp,1))*detbet1(1:mp)
    end do
    ! back substitution
    do j=n-1, 1, -1
       tmp1d(1:mp,1) = G(1:mp,j+1,1,1)*u(1:mp,j+1,1) + G(1:mp,j+1,1,2)*u(1:mp,j+1,2)
       tmp1d(1:mp,2) = G(1:mp,j+1,2,1)*u(1:mp,j+1,1) + G(1:mp,j+1,2,2)*u(1:mp,j+1,2)
       u(1:mp,j,1:2) = u(1:mp,j,1:2) - tmp1d(1:mp,1:2)
    end do
    !
  END SUBROUTINE generic_thomas_2d

  !**********************************************************************************************************************

  SUBROUTINE getfluxes_vp_1d(n, dx, vtop, vbot, parin, var, hint, phimin, q, qya, qyb, qTa, qTb, &
       ql, qlya, qlyb, qv, qvT, qvh, qvya, qvyb, &
       iflux, init, getq0, getqn, Tsoil, T0, nsat, nsatlast)

    IMPLICIT NONE

    INTEGER(i_d),                 INTENT(IN)    :: n
    REAL(r_2),    DIMENSION(1:n), INTENT(IN)    :: dx
    TYPE(vars),                   INTENT(IN)    :: vtop
    TYPE(vars),                   INTENT(IN)    :: vbot
    TYPE(params), DIMENSION(1:n), INTENT(IN)    :: parin
    TYPE(vars),   DIMENSION(1:n), INTENT(INOUT)    :: var
    REAL(r_2),    DIMENSION(1:n), INTENT(INOUT) :: hint
    REAL(r_2),    DIMENSION(1:n), INTENT(INOUT) :: phimin
    REAL(r_2),    DIMENSION(0:n), INTENT(INOUT) :: q
    REAL(r_2),    DIMENSION(0:n), INTENT(INOUT) :: qya
    REAL(r_2),    DIMENSION(0:n), INTENT(INOUT) :: qyb
    REAL(r_2),    DIMENSION(0:n), INTENT(INOUT) :: qTa
    REAL(r_2),    DIMENSION(0:n), INTENT(INOUT) :: qTb
    REAL(r_2),    DIMENSION(0:n), INTENT(INOUT)   :: ql
    REAL(r_2),    DIMENSION(0:n), INTENT(INOUT)   :: qlya
    REAL(r_2),    DIMENSION(0:n), INTENT(INOUT)   :: qlyb
    REAL(r_2),    DIMENSION(0:n), INTENT(INOUT)   :: qv
    REAL(r_2),    DIMENSION(0:n), INTENT(INOUT)   :: qvT
    REAL(r_2),    DIMENSION(0:n), INTENT(INOUT)   :: qvh
    REAL(r_2),    DIMENSION(0:n), INTENT(INOUT)   :: qvya
    REAL(r_2),    DIMENSION(0:n), INTENT(INOUT)   :: qvyb
    INTEGER(i_d),                 INTENT(IN)    :: iflux
    LOGICAL,                      INTENT(IN)    :: init
    LOGICAL,                      INTENT(IN)    :: getq0
    LOGICAL,                      INTENT(IN)    :: getqn
    REAL(r_2),    DIMENSION(1:n), INTENT(IN)    :: Tsoil
    REAL(r_2),                    INTENT(IN)    :: T0
    INTEGER(i_d),                 INTENT(IN)    :: nsat
    INTEGER(i_d),                 INTENT(IN)    :: nsatlast
    ! Gets fluxes q and partial derivs qya, qyb wrt S (if unsat) or phi (if sat).
    ! Fluxes at top and bottom of profile, and fluxes due to plant extraction of
    ! water are included.
    ! Definitions of arguments:
    ! k     - land point
    ! n     - no. of soil layers.
    ! jt(1:n)   - layer soil type nos.
    ! dx(1:n)   - layer thicknesses.
    ! dz(1:n-1)   - distances between layer centres.
    ! vtop    - water vars at soil surface.
    ! vbot    - water vars at bottom of profile.
    ! var(1:n)   - water vars at layer centres.
    ! hint(1:n)   - values of h at interfaces are stored sequentially in hint.
    ! phimin(1:n) - similarly for phi at hmin in layers above interfaces.
    ! q(0:n)   - fluxes; q(i), i=1,...,n-1 is flux from layer i to layer i+1.
    !    q(0) is surface flux and q(n) is flux at bottom of profile.
    ! qya(0:n)   - partial deriv of q(i), i=0,...,n, wrt the variable to be solved
    !    for (S, phi or h) at upper end of flow path.
    ! qyb(0:n)   - ditto for var at lower end.
    ! iflux    - if iflux/=1, get only fluxes involving sat layers.
    ! init    - true if hint and phimin to be initialised.
    ! getq0    - true if q(0) required.
    ! getqn    - true if q(n) required.
    LOGICAL               :: flag, limit
    INTEGER(i_d)          :: i, itmp, l
    REAL(r_2)             :: dphii1, dhi, h1, h2, hi, Khi1, Khi2, phii1, q2, qya2, qyb2, y, y1, y2
    REAL(r_2)             :: qTa2, qTb2
    TYPE(vars)            :: vi1, vi2
    REAL(r_2), DIMENSION(1:n-1) :: dz

    dz(:) = half*(dx(1:n-1)+dx(2:n))
    vi1 = zerovars()
    vi2 = zerovars()

    if ((iflux==1) .or. (var(1)%isat /= 0)) then ! get top flux if required
       if (getq0) then
          call flux(parin(1), vtop, var(1), half*dx(1), q(0), qya(0), qyb(0), qTa(0), qTb(0))

          q(0)  = q(0)+(T0-Tsoil(1))*(var(1)%kE)/thousand/var(1)%lambdav/dx(1)*two

          qTa(0) = zero
          qTb(0) = -(var(1)%kE)/thousand/var(1)%lambdav/dx(1)*two
          qv(0)  = (vtop%phiv-var(1)%phiv)/dx(1)*two +(T0-Tsoil(1))*(var(1)%kE)/thousand/var(1)%lambdav/dx(1)*two
          ql(0)  = q(0) - qv(0)

          qvya(0) = zero

          if (vtop%isat==0) then
             qvyb(0) = vtop%phivS/dx(1)*two
          else
             qvyb(0) = zero
          end if

          qlya(0) = qya(0) - qvya(0)
          qlyb(0) = qyb(0) - qvyb(0)
       end if
    end if
    ! otherwise undefined
    qvh(0) = zero
    qvT(0) = zero

    ! get other fluxes
    l = 0
    do i=1, n-1
       if (iflux==1 .or. var(i)%isat/=0 .or. var(i+1)%isat/=0 .or. nsat/=nsatlast) then ! get flux
          if (parin(i)%ishorizon == parin(i+1)%ishorizon) then ! same soil type, no interface
             call flux(parin(i), var(i), var(i+1), dz(i), q(i), qya(i), qyb(i), qTa(i), qTb(i))
          else ! interface
             l  = l+1
             if (init) then ! initialise
                call hyofh(hmin, parin(i)%lam, parin(i)%eta, parin(i)%Ke, parin(i)%he, vi1%K, Khi1, phimin(l)) ! get phi at hmin
                h1 = var(i)%h
                h2 = var(i+1)%h
                y1 = var(i)%K*dx(i+1)
                y2 = var(i+1)%K*dx(i)
                ! equate fluxes (K constant) to get initial estimate of h at interface
                hint(l) = (y1*h1+y2*h2+half*gf*(var(i)%K-var(i+1)%K)*dx(i)*dx(i+1))/(y1+y2)
             end if
             hi   = hint(l)
             flag = .true.
             itmp = 0
             ! iterate to get hi at interface for equal fluxes using Newton's method
             ! get dphii1 at interface in upper layer, because of better linearity,
             ! then convert to dhi
             do while (flag)
                itmp = itmp+1
                if (itmp>100) then
                   !write(*,*) "getfluxes: too many iterations finding interface h"
                   !stop
                   call flux(parin(i), var(i), var(i+1), dz(i), q(i), qya(i), qyb(i), qTa(i), qTb(i))
                   goto 111
                end if

                if (hi<parin(i)%he) then
                   vi1%isat = 0
                   call hyofh(hi, parin(i)%lam, parin(i)%eta, parin(i)%Ke, parin(i)%he, vi1%K, Khi1, phii1)
                   vi1%KS = Khi1/vi1%K ! use dK/dphi, not dK/dS
                else
                   vi1%isat = 1
                   vi1%K    = var(i)%Ksat
                   phii1    = var(i)%phie+(hi-parin(i)%he)*var(i)%Ksat
                   vi1%KS   = zero
                end if

                vi1%h    = hi
                vi1%phi  = phii1
                vi1%phiS = one ! use dphi/dphi not dphi/dS
                ! define phiT=0, KT=0 to be consistent with undefined version
                vi1%phiT = zero
                vi1%KT   = zero
                ! macropore_factor was not defined but is used in flux(), set to factor of upper layer
                vi1%macropore_factor = var(i)%macropore_factor
                call flux(parin(i), var(i), vi1, half*dx(i), q(i), qya(i), qyb(i), qTa(i), qTb(i))

                if (hi<parin(i+1)%he) then
                   vi2%isat = 0
                   call hyofh(hi, parin(i+1)%lam, parin(i+1)%eta, parin(i+1)%Ke, parin(i+1)%he, vi2%K, Khi2, vi2%phi)
                   vi2%KS = Khi2/vi2%K ! dK/dphi
                else
                   vi2%isat = 1
                   vi2%K    = var(i+1)%Ksat
                   vi2%phi  = var(i+1)%phie+(hi-parin(i+1)%he)*var(i+1)%Ksat
                end if

                vi2%h    = hi
                vi2%phiS = one ! dphi/dphi
                ! define phiT=0, KT=0 to be consitent with undefined version
                vi2%phiT = zero
                vi2%KT   = zero
                ! macropore_factor was not defined but is used in flux(), set to factor of lower layer
                vi2%macropore_factor = var(i+1)%macropore_factor
                call flux(parin(i+1), vi2, var(i+1), half*dx(i+1), q2, qya2, qyb2, qTa2, qTb2)
                qya2   = qya2*vi2%K/vi1%K ! partial deriv wrt phii1
                ! adjust for equal fluxes
                dphii1 = -(q(i)-q2)/(qyb(i)-qya2)
                limit  = .false.
                if (phii1+dphii1<=phimin(l)) then ! out of range
                   limit  = .true.
                   dphii1 = -half*(phii1-phimin(l))
                end if
                phii1 = phii1+dphii1
                dhi   = dphii1/(vi1%K+half*vi1%KS*dphii1) ! 2nd order Pade approx
                if (-vi1%KS*dphii1 > 1.5_r_2*vi1%K) then ! use 1st order approx for dhi
                   dhi = dphii1/vi1%K
                end if
                hi = hi+dhi

                ! check for convergence - dphi/(mean phi)<=dpmaxr
                if (.not.(limit .or. abs(dphii1/(phii1-half*dphii1))>dpmaxr)) then
                   flag = .false.
                end if
             end do ! while flag

             q(i)    = q(i) + qyb(i)*dphii1
             hint(l) = hi
             ! adjust derivs
             y      = one/(qya2-qyb(i))
             qya(i) = qya(i)*qya2*y
             qyb(i) = -qyb2*qyb(i)*y
          end if
       end if

111    ql(i)  = q(i)
       qTa(i) = qTa(i)+(var(i)%kE+var(i+1)%kE)/thousand/var(1)%lambdav/two/dz(i)
       qTb(i) = qTb(i)-(var(i)%kE+var(i+1)%kE)/thousand/var(1)%lambdav/two/dz(i)
       qvT(i) = (Tsoil(i)-Tsoil(i+1))*(var(i)%kE+var(i+1)%kE)/thousand/var(1)%lambdav/two/dz(i)
       !MC The full description (next two lines) gave problems before ->  third line
       !   Try again original
       !VH Do both formulations come to the same? I found that I could not reproduce
       !   Barnes and Allison semi-analytic solution with original
       !MC Should be re-checked
       ! qvh(i) = ((((Tsoil(i)+Tzero)/Tzero)**1.88+((Tsoil(i+1)+Tzero)/Tzero)**1.88)/two) &
       !      * ((var(i)%cvsat+var(i+1)%cvsat)/two)*(var(i)%phiv-var(i+1)%phiv)/dz(i)
       qvh(i) = (half*(exp(1.88_r_2*log((Tsoil(i)+Tzero)/Tzero))+exp(1.88_r_2*log((Tsoil(i+1)+Tzero)/Tzero)))) &
            * (half*(var(i)%cvsat+var(i+1)%cvsat)) * (var(i)%phiv-var(i+1)%phiv)/dz(i)
       ! qvh(i) = ((var(i)%Dv+var(i+1)%Dv)/two)* ((var(i)%cvsat+var(i+1)%cvsat)/two)*(var(i)%rh-var(i+1)%rh)/dz(i)
       qv(i)  = qvh(i) + qvT(i) ! whole vapour flux has one part from humidity (qvh) and one part from temp diff (qvT)
       q(i)   = qv(i) + ql(i)

       if (var(i)%isat==0) then
          ! qvya(i) = var(i)%phivS/dz(i) *((((Tsoil(i)+Tzero)/Tzero)**1.88_r_2+ &
          !      ((Tsoil(i+1)+Tzero)/Tzero)**1.88_r_2)/two) &
          !      * ((var(i)%cvsat+var(i+1)%cvsat)/two)
          qvya(i) = var(i)%phivS/dz(i) *(half*(exp(1.88_r_2*log((Tsoil(i)+Tzero)/Tzero)) + &
               exp(1.88_r_2*log((Tsoil(i+1)+Tzero)/Tzero)))) &
               * (half*(var(i)%cvsat+var(i+1)%cvsat))
       else
          qvya(i) = zero
       end if

       if (var(i)%isat==0) then
          ! qvyb(i) = -var(i+1)%phivS/dz(i) *((((Tsoil(i)+Tzero)/Tzero)**1.88_r_2+ &
          !      ((Tsoil(i+1)+Tzero)/Tzero)**1.88_r_2)/two) &
          !      * ((var(i)%cvsat+var(i+1)%cvsat)/two)
          qvyb(i) = -var(i+1)%phivS/dz(i) *(half*(exp(1.88_r_2*log((Tsoil(i)+Tzero)/Tzero)) + &
               exp(1.88_r_2*log((Tsoil(i+1)+Tzero)/Tzero)))) &
               * (half*(var(i)%cvsat+var(i+1)%cvsat))
       else
          qvyb(i) = zero
       end if

       qlya(i) = qya(i)
       qlyb(i) = qyb(i)
       qya(i)  = qya(i) + qvya(i)
       qyb(i)  = qyb(i) + qvyb(i)
    end do

    if (iflux==1 .or. var(n)%isat/=0) then ! get bottom flux if required
       if (getqn) then
          call flux(parin(n), var(n), vbot, half*dx(n), q(n), qya(n), qyb(n), qTa(n), qTb(n))
          qvya(n) = zero
          qvyb(n) = zero
          qlya(n) = qya(n)
          qlyb(n) = zero
       else
          qvya(n) = zero
          qvyb(n) = zero
          qlya(n) = qya(n)
          qlyb(n) = zero
       end if
    else
       qvya(n) = zero
       qvyb(n) = zero
       qlya(n) = qya(n)
       qlyb(n) = zero
    end if
    ! otherwise undefined
    ql(n)  = q(n)
    qv(n)  = zero
    qvh(n) = zero
    qvT(n) = zero

    do i=1, n-1
       if (eq(var(i)%Dv, zero) .or. eq(var(i+1)%Dv, zero)) then
          q(i)    = q(i) - qv(i)
          qya(i)  = qya(i) - qvya(i)
          qyb(i)  = qyb(i) - qvyb(i)
          qv(i)   = zero
          !qTa(i)  = zero
          !qTb(i)  = zero
          qvya(i) = zero
          qvyb(i) = zero
       endif
    enddo

  END SUBROUTINE getfluxes_vp_1d

  SUBROUTINE getfluxes_vp_2d(dx, vtop, vbot, parin, var, hint, phimin, i_q, i_qya, i_qyb, i_qTa, i_qTb, &
       i_ql, i_qlya, i_qlyb, i_qv, i_qvT, i_qvh, i_qvya, i_qvyb, iflux, init, getq0, getqn, Tsoil, T0, nsat, nsatlast)

    IMPLICIT NONE

    REAL(r_2),    DIMENSION(:,:), INTENT(IN)    :: dx      ! 1:n
    TYPE(vars),   DIMENSION(:),   INTENT(IN)    :: vtop
    TYPE(vars),   DIMENSION(:),   INTENT(IN)    :: vbot
    TYPE(params), DIMENSION(:,:), INTENT(IN)    :: parin   ! 1:n
    TYPE(vars),   DIMENSION(:,:), INTENT(IN)    :: var     ! 1:n
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: hint    ! 1:n
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: phimin  ! 1:n
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_q       ! 0:n
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_qya     ! 0:n
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_qyb     ! 0:n
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_qTa     ! 0:n
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_qTb     ! 0:n
    REAL(r_2),    DIMENSION(:,:), INTENT(OUT)   :: i_ql      ! 0:n
    REAL(r_2),    DIMENSION(:,:), INTENT(OUT)   :: i_qlya    ! 0:n
    REAL(r_2),    DIMENSION(:,:), INTENT(OUT)   :: i_qlyb    ! 0:n
    REAL(r_2),    DIMENSION(:,:), INTENT(OUT)   :: i_qv      ! 0:n
    REAL(r_2),    DIMENSION(:,:), INTENT(OUT)   :: i_qvT     ! 0:n
    REAL(r_2),    DIMENSION(:,:), INTENT(OUT)   :: i_qvh     ! 0:n
    REAL(r_2),    DIMENSION(:,:), INTENT(OUT)   :: i_qvya    ! 0:n
    REAL(r_2),    DIMENSION(:,:), INTENT(OUT)   :: i_qvyb    ! 0:n
    INTEGER(i_d), DIMENSION(:),   INTENT(IN)    :: iflux
    LOGICAL,      DIMENSION(:),   INTENT(IN)    :: init
    LOGICAL,      DIMENSION(:),   INTENT(IN)    :: getq0
    LOGICAL,      DIMENSION(:),   INTENT(IN)    :: getqn
    REAL(r_2),    DIMENSION(:,:), INTENT(IN)    :: Tsoil   ! 1:n
    REAL(r_2),    DIMENSION(:),   INTENT(IN)    :: T0
    INTEGER(i_d), DIMENSION(:),   INTENT(IN)    :: nsat
    INTEGER(i_d), DIMENSION(:),   INTENT(IN)    :: nsatlast
    ! Gets fluxes q and partial derivs qya, qyb wrt S (if unsat) or phi (if sat).
    ! Fluxes at top and bottom of profile, and fluxes due to plant extraction of
    ! water are included.
    ! Definitions of arguments:
    ! k     - land point
    ! n     - no. of soil layers.
    ! jt(1:n)   - layer soil type nos.
    ! dx(1:n)   - layer thicknesses.
    ! dz(1:n-1)   - distances between layer centres.
    ! vtop    - water vars at soil surface.
    ! vbot    - water vars at bottom of profile.
    ! var(1:n)   - water vars at layer centres.
    ! hint(1:n)   - values of h at interfaces are stored sequentially in hint.
    ! phimin(1:n) - similarly for phi at hmin in layers above interfaces.
    ! q(0:n)   - fluxes; q(i), i=1,...,n-1 is flux from layer i to layer i+1.
    !    q(0) is surface flux and q(n) is flux at bottom of profile.
    ! qya(0:n)   - partial deriv of q(i), i=0,...,n, wrt the variable to be solved
    !    for (S, phi or h) at upper end of flow path.
    ! qyb(0:n)   - ditto for var at lower end.
    ! iflux    - if iflux/=1, get only fluxes involving sat layers.
    ! init    - true if hint and phimin to be initialised.
    ! getq0    - true if q(0) required.
    ! getqn    - true if q(n) required.
    LOGICAL,      DIMENSION(1:size(dx,1))                :: limit, l1, l2, l3
    REAL(r_2),    DIMENSION(1:size(dx,1))                :: dphii1, dhi, h1, h2, hi, Khi1, Khi2, phii1
    REAL(r_2),    DIMENSION(1:size(dx,1))                :: q2, qya2, qyb2, y, y1, y2
    REAL(r_2),    DIMENSION(1:size(dx,1))                :: htmp
    REAL(r_2),    DIMENSION(1:size(dx,1))                :: ztmp1, ztmp2, ztmp3, ztmp4, ztmp5
    TYPE(vars),   DIMENSION(1:size(dx,1))                :: vi1, vi2
    REAL(r_2),    DIMENSION(1:size(dx,1),1:size(dx,2)-1) :: dz
    REAL(r_2),    DIMENSION(1:size(dx,1),0:size(dx,2))   :: q
    REAL(r_2),    DIMENSION(1:size(dx,1),0:size(dx,2))   :: qya
    REAL(r_2),    DIMENSION(1:size(dx,1),0:size(dx,2))   :: qyb
    REAL(r_2),    DIMENSION(1:size(dx,1),0:size(dx,2))   :: qTa
    REAL(r_2),    DIMENSION(1:size(dx,1),0:size(dx,2))   :: qTb
    REAL(r_2),    DIMENSION(1:size(dx,1),0:size(dx,2))   :: ql
    REAL(r_2),    DIMENSION(1:size(dx,1),0:size(dx,2))   :: qlya
    REAL(r_2),    DIMENSION(1:size(dx,1),0:size(dx,2))   :: qlyb
    REAL(r_2),    DIMENSION(1:size(dx,1),0:size(dx,2))   :: qv
    REAL(r_2),    DIMENSION(1:size(dx,1),0:size(dx,2))   :: qvT
    REAL(r_2),    DIMENSION(1:size(dx,1),0:size(dx,2))   :: qvh
    REAL(r_2),    DIMENSION(1:size(dx,1),0:size(dx,2))   :: qvya
    REAL(r_2),    DIMENSION(1:size(dx,1),0:size(dx,2))   :: qvyb
    TYPE(vars)    :: vtmp
    INTEGER(i_d)  :: i, n, mp, itmp

    mp = size(dx,1)
    n  = size(dx,2)
    dz(:,1:n-1) = half*(dx(:,1:n-1)+dx(:,2:n))
    q(:,0:n)    = i_q(:,1:n+1)
    qya(:,0:n)  = i_qya(:,1:n+1)
    qyb(:,0:n)  = i_qyb(:,1:n+1)
    qTa(:,0:n)  = i_qTa(:,1:n+1)
    qTb(:,0:n)  = i_qTb(:,1:n+1)
    ql(:,0:n)   = zero
    qlya(:,0:n) = zero
    qlyb(:,0:n) = zero
    qv(:,0:n)   = zero
    qvT(:,0:n)  = zero
    qvh(:,0:n)  = zero
    qvya(:,0:n) = zero
    qvyb(:,0:n) = zero

    vtmp = zerovars()
    vi1  = spread(vtmp,1,mp)
    vi2  = spread(vtmp,1,mp)
    ztmp1(:) = zero
    ztmp2(:) = zero
    ztmp3(:) = zero
    ztmp4(:) = zero
    ztmp5(:) = zero

    l1(:)   = ((iflux(:)==1) .or. (var(:,1)%isat /= 0)) .and. getq0(:)
    if (any(l1(:))) &
         call flux(parin(:,1), vtop(:), var(:,1), half*dx(:,1), ztmp1(:), ztmp2(:), ztmp3(:), ztmp4(:), ztmp5(:))
    where (l1(:)) ! get top flux if required
       q(:,0)   = ztmp1(:)
       qya(:,0) = ztmp2(:)
       qyb(:,0) = ztmp3(:)
       qTa(:,0) = ztmp4(:)
       qTb(:,0) = ztmp5(:)

       q(:,0)    = q(:,0)+(T0(:)-Tsoil(:,1))*(var(:,1)%kE)/thousand/var(:,1)%lambdav/dx(:,1)*two
       qTa(:,0)  = zero
       qTb(:,0)  = -(var(:,1)%kE)/thousand/var(:,1)%lambdav/dx(:,1)*two
       qv(:,0)   = (vtop(:)%phiv-var(:,1)%phiv)/dx(:,1)*two +(T0(:)-Tsoil(:,1))* &
            (var(:,1)%kE)/thousand/var(:,1)%lambdav/dx(:,1)*two
       ql(:,0)   = q(:,0) - qv(:,0)
       qvya(:,0) = zero
       qvyb(:,0) = zero ! isat==0 below
       qlya(:,0) = qya(:,0) - qvya(:,0)
       qlyb(:,0) = qyb(:,0) - qvyb(:,0)
    endwhere
    where (l1(:) .and. (vtop(:)%isat==0))
       qvyb(:,0) = vtop(:)%phivS/dx(:,1)*two
       qlyb(:,0) = qyb(:,0) - qvyb(:,0)
    endwhere
    ! otherwise undefined
    qvh(:,0) = zero
    qvT(:,0) = zero

    ! get other fluxes
    do i=1, n-1
       l1(:) = (iflux(:)==1 .or. var(:,i)%isat/=0 .or. var(:,i+1)%isat/=0 .or. nsat(:)/=nsatlast(:))
       if (any(l1(:) .and. (parin(:,i)%ishorizon==parin(:,i+1)%ishorizon))) &
            call flux(parin(:,i), var(:,i), var(:,i+1), dz(:,i), ztmp1(:), ztmp2(:), ztmp3(:), ztmp4(:), ztmp5(:))
       where (l1(:) .and. (parin(:,i)%ishorizon==parin(:,i+1)%ishorizon)) ! same soil type, no interface
          q(:,i)   = ztmp1(:)
          qya(:,i) = ztmp2(:)
          qyb(:,i) = ztmp3(:)
          qTa(:,i) = ztmp4(:)
          qTb(:,i) = ztmp5(:)
       endwhere
       l2(:) = l1(:) .and. (parin(:,i)%ishorizon /= parin(:,i+1)%ishorizon) ! interface
       htmp(:) = hmin
       if (any(l2(:) .and. init(:))) &
            call hyofh(htmp(:), parin(:,i)%lam, parin(:,i)%eta, parin(:,i)%Ke, parin(:,i)%he, ztmp1(:), ztmp2(:), ztmp3(:))
       where (l2(:) .and. init(:)) ! interface & get phi at hmin
          vi1(:)%K    = ztmp1(:)
          Khi1(:)     = ztmp2(:)
          phimin(:,i) = ztmp3(:)
          h1(:)       = var(:,i)%h
          h2(:)       = var(:,i+1)%h
          y1(:)       = var(:,i)%K*dx(:,i+1)
          y2(:)       = var(:,i+1)%K*dx(:,i)
          ! equate fluxes (K constant) to get initial estimate of h at interface
          hint(:,i) = (y1(:)*h1(:) + y2(:)*h2(:) + half*gf*(var(:,i)%K-var(:,i+1)%K)*dx(:,i)*dx(:,i+1)) / (y1(:)+y2(:))
       endwhere
       !where ((.not. l2(:)) .and. init(:)) hint(:,i) = zero
       hi(:)   = hint(:,i)
       ! iterate to get hi at interface for equal fluxes using Newton's method
       ! get dphii1 at interface in upper layer, because of better linearity,
       ! then convert to dhi
       if (any(l2(:))) then
          l3(:)    = l2(:)
          limit(:) = .false.
          do itmp=1, 100
             if (any(l3(:) .and. (hi(:)<parin(:,i)%he))) &
                  call hyofh(hi(:), parin(:,i)%lam, parin(:,i)%eta, parin(:,i)%Ke, parin(:,i)%he, ztmp1(:), ztmp2(:), ztmp3(:))
             where (l3(:) .and. (hi(:)<parin(:,i)%he))
                vi1(:)%isat = 0
                vi1(:)%K    = ztmp1(:)
                Khi1(:)     = ztmp2(:)
                phii1(:)    = ztmp3(:)
                vi1(:)%KS   = Khi1(:)/vi1(:)%K ! use dK/dphi, not dK/dS
             endwhere
             where (l3(:) .and. (hi(:)>=parin(:,i)%he))
                vi1(:)%isat = 1
                vi1(:)%K    = var(:,i)%Ksat
                phii1(:)    = var(:,i)%phie+(hi(:)-parin(:,i)%he)*var(:,i)%Ksat
                vi1(:)%KS   = zero
             endwhere

             where (l3(:))
                vi1(:)%h    = hi(:)
                vi1(:)%phi  = phii1(:)
                vi1(:)%phiS = one ! use dphi/dphi not dphi/dS
                !define phiT=0, KT=0 to be consitent with undefined version
                vi1(:)%phiT = zero
                vi1(:)%KT   = zero
                ! macropore_factor was not defined but is used in flux(), set to factor of upper layer
                vi1(:)%macropore_factor = var(:,i)%macropore_factor
             endwhere

             if (any(l3(:))) &
                  call flux(parin(:,i), var(:,i), vi1(:), half*dx(:,i), ztmp1(:), ztmp2(:), ztmp3(:), ztmp4(:), ztmp5(:))
             where (l3(:))
                q(:,i)      = ztmp1(:)
                qya(:,i)    = ztmp2(:)
                qyb(:,i)    = ztmp3(:)
                qTa(:,i)    = ztmp4(:)
                qTb(:,i)    = ztmp5(:)
             endwhere

             if (any(l3(:) .and. (hi(:)<parin(:,i+1)%he))) &
                  call hyofh(hi(:), parin(:,i+1)%lam, parin(:,i+1)%eta, parin(:,i+1)%Ke, parin(:,i+1)%he, &
                  ztmp1(:), ztmp2(:), ztmp3(:))
             where (l3(:) .and. (hi(:)<parin(:,i+1)%he))
                vi2(:)%K    = ztmp1(:)
                Khi2(:)     = ztmp2(:)
                vi2(:)%phi  = ztmp3(:)
                vi2(:)%isat = 0
                vi2(:)%KS   = Khi2(:)/vi2(:)%K ! dK/dphi
             endwhere
             where (l3(:) .and. (hi(:)>=parin(:,i+1)%he))
                vi2(:)%isat = 1
                vi2(:)%K    = var(:,i+1)%Ksat
                vi2(:)%phi  = var(:,i+1)%phie+(hi(:)-parin(:,i+1)%he)*var(:,i+1)%Ksat
             endwhere

             where (l3(:))
                vi2(:)%h    = hi(:)
                vi2(:)%phiS = one ! dphi/dphi
                ! define phiT=0, KT=0 to be consitent with undefined version
                vi2(:)%phiT = zero
                vi2(:)%KT   = zero
                ! macropore_factor was not defined but is used in flux(), set to factor of lower layer
                vi2(:)%macropore_factor = var(:,i+1)%macropore_factor
             endwhere
             if (any(l3(:))) &
                  call flux(parin(:,i+1), vi2(:), var(:,i+1), half*dx(:,i+1), ztmp1(:), ztmp2(:), ztmp3(:), ztmp4(:), ztmp5(:))
             where (l3(:))
                q2(:)     = ztmp1(:)
                qya2(:)   = ztmp2(:)
                qyb2(:)   = ztmp3(:)
                qya2(:)   = qya2(:)*vi2(:)%K/vi1(:)%K ! partial deriv wrt phii1
                ! adjust for equal fluxes
                dphii1(:) = -(q(:,i)-q2(:))/(qyb(:,i)-qya2(:))
                limit(:)  = .false.
             endwhere
             where (l3(:) .and. (phii1(:)+dphii1(:)<=phimin(:,i))) ! out of range
                limit(:)  = .true.
                dphii1(:) = -half*(phii1(:)-phimin(:,i))
             endwhere
             where (l3(:))
                phii1(:) = phii1(:)+dphii1(:)
                dhi(:)   = dphii1(:)/(vi1(:)%K+half*vi1(:)%KS*dphii1(:)) ! 2nd order Pade approx
             endwhere
             where (l3(:) .and. (-vi1%KS*dphii1 > 1.5_r_2*vi1%K)) ! use 1st order approx for dhi
                dhi(:) = dphii1(:)/vi1(:)%K
             endwhere
             where (l3(:))
                hi(:) = hi(:)+dhi(:)
             endwhere

             ! check for convergence - dphi/(mean phi)<=dpmaxr
             where (l3(:) .and. &
                  .not. (limit(:) .or. (abs(dphii1(:)/(phii1(:)-half*dphii1(:)))>dpmaxr))) l3(:) = .false.
             if (.not. any(l3(:))) exit
          end do ! do itmp=1, 100
          if (itmp>=100) then
             !write(*,*) "getfluxes: too many iterations finding interface h"
             !stop
             if (any(l2(:) .and. l3(:))) &
                  call flux(parin(:,i), var(:,i), var(:,i+1), dz(:,i), ztmp1(:), ztmp2(:), ztmp3(:), ztmp4(:), ztmp5(:))
             where (l2(:) .and. l3(:))
                q(:,i)   = ztmp1(:)
                qya(:,i) = ztmp2(:)
                qyb(:,i) = ztmp3(:)
                qTa(:,i) = ztmp4(:)
                qTb(:,i) = ztmp5(:)
             endwhere
          else
             where (l2(:) .and. (.not. l3(:)))
                q(:,i)    = q(:,i) + qyb(:,i)*dphii1(:)
                hint(:,i) = hi(:)
                ! adjust derivs
                y(:)      = one/(qya2(:)-qyb(:,i))
                qya(:,i) = qya(:,i)*qya2(:)*y(:)
                qyb(:,i) = -qyb2(:)*qyb(:,i)*y(:)
             endwhere
          end if
       end if

       ql(:,i)  = q(:,i)
       qTa(:,i) = qTa(:,i)+(var(:,i)%kE+var(:,i+1)%kE)/thousand/var(:,1)%lambdav/two/dz(:,i)
       qTb(:,i) = qTb(:,i)-(var(:,i)%kE+var(:,i+1)%kE)/thousand/var(:,1)%lambdav/two/dz(:,i)
       qvT(:,i) = (Tsoil(:,i)-Tsoil(:,i+1))*(var(:,i)%kE+var(:,i+1)%kE)/thousand/var(:,1)%lambdav/two/dz(:,i)
       !MC The full description (next two lines) gave problems before ->  third line
       !   Try again original
       !VH Do both formulations come to the same? I found that I could not reproduce
       !   Barnes and Allison semi-analytic solution with original
       !MC Should be re-checked
       ! qvh(:,i) = ((((Tsoil(:,i)+Tzero)/Tzero)**1.88+((Tsoil(:,i+1)+Tzero)/Tzero)**1.88)/two) &
       !      * ((var(:,i)%cvsat+var(:,i+1)%cvsat)/two)*(var(:,i)%phiv-var(:,i+1)%phiv)/dz(:,i)
       qvh(:,i) = (half*(exp(1.88_r_2*log((Tsoil(:,i)+Tzero)/Tzero))+exp(1.88_r_2*log((Tsoil(:,i+1)+Tzero)/Tzero)))) &
            * (half*(var(:,i)%cvsat+var(:,i+1)%cvsat)) * (var(:,i)%phiv-var(:,i+1)%phiv)/dz(:,i)
       ! qvh(:,i) = ((var(:,i)%Dv+var(:,i+1)%Dv)/two) * &
       !             ((var(:,i)%cvsat+var(:,i+1)%cvsat)/two)*(var(:,i)%rh-var(:,i+1)%rh)/dz(:,i)
       qv(:,i)  = qvh(:,i) + qvT(:,i) ! whole vapour flux has one part from humidity (qvh) and one part from temp diff (qvT)
       q(:,i)   = qv(:,i) + ql(:,i)

       where (var(:,i)%isat==0)
          ! qvya(:,i) = var(:,i)%phivS/dz(:,i) *((((Tsoil(:,i)+Tzero)/Tzero)**1.88_r_2+ &
          !      ((Tsoil(:,i+1)+Tzero)/Tzero)**1.88_r_2)/two) &
          !      * ((var(:,i)%cvsat+var(:,i+1)%cvsat)/two)
          qvya(:,i) = var(:,i)%phivS/dz(:,i) *(half*(exp(1.88_r_2*log((Tsoil(:,i)+Tzero)/Tzero)) + &
               exp(1.88_r_2*log((Tsoil(:,i+1)+Tzero)/Tzero)))) &
               * (half*(var(:,i)%cvsat+var(:,i+1)%cvsat))
       elsewhere
          qvya(:,i) = zero
       endwhere

       where (var(:,i)%isat==0)
          ! qvyb(:,i) = -var(:,i+1)%phivS/dz(:,i) *((((Tsoil(:,i)+Tzero)/Tzero)**1.88_r_2+ &
          !      ((Tsoil(:,i+1)+Tzero)/Tzero)**1.88_r_2)/two) &
          !      * ((var(:,i)%cvsat+var(:,i+1)%cvsat)/two)
          qvyb(:,i) = -var(:,i+1)%phivS/dz(:,i) *(half*(exp(1.88_r_2*log((Tsoil(:,i)+Tzero)/Tzero)) + &
               exp(1.88_r_2*log((Tsoil(:,i+1)+Tzero)/Tzero)))) &
               * (half*(var(:,i)%cvsat+var(:,i+1)%cvsat))
       elsewhere
          qvyb(:,i) = zero
       endwhere

       qlya(:,i) = qya(:,i)
       qlyb(:,i) = qyb(:,i)
       qya(:,i)  = qya(:,i) + qvya(:,i)
       qyb(:,i)  = qyb(:,i) + qvyb(:,i)
    end do

    l1(:) = (iflux(:)==1) .or. (var(:,n)%isat/=0)
    if (any(l1(:) .and. getqn(:))) &  ! get bottom flux if required
         call flux(parin(:,n), var(:,n), vbot(:), half*dx(:,n), ztmp1(:), ztmp2(:), ztmp3(:), ztmp4(:), ztmp5(:))
    where (l1(:) .and. getqn(:))
       q(:,n)    = ztmp1(:)
       qya(:,n)  = ztmp2(:)
       qyb(:,n)  = ztmp3(:)
       qTa(:,n)  = ztmp4(:)
       qTb(:,n)  = ztmp5(:)
       qvya(:,n) = zero
       qvyb(:,n) = zero
       qlya(:,n) = qya(:,n)
       qlyb(:,n) = zero
    endwhere

    where (l1(:) .and. (.not. getqn(:)))
       qvya(:,n) = zero
       qvyb(:,n) = zero
       qlya(:,n) = qya(:,n)
       qlyb(:,n) = zero
    endwhere
    where (.not. l1(:))
       qvya(:,n) = zero
       qvyb(:,n) = zero
       qlya(:,n) = qya(:,n)
       qlyb(:,n) = zero
    endwhere
    ! otherwise undefined
    ql(:,n)  = q(:,n)
    qv(:,n)  = zero
    qvh(:,n) = zero
    qvT(:,n) = zero

    do i=1, n-1
       where (eq(var(:,i)%Dv, zero) .or. eq(var(:,i+1)%Dv, zero))
          q(:,i)    = q(:,i) - qv(:,i)
          qya(:,i)  = qya(:,i) - qvya(:,i)
          qyb(:,i)  = qyb(:,i) - qvyb(:,i)
          qv(:,i)   = zero
          !qTa(:,i)  = zero
          !qTb(:,i)  = zero
          qvya(:,i) = zero
          qvyb(:,i) = zero
       endwhere
    enddo

    i_q(:,1:n+1)    = q(:,0:n)
    i_qya(:,1:n+1)  = qya(:,0:n)
    i_qyb(:,1:n+1)  = qyb(:,0:n)
    i_qTa(:,1:n+1)  = qTa(:,0:n)
    i_qTb(:,1:n+1)  = qTb(:,0:n)
    i_ql(:,1:n+1)   = ql(:,0:n)
    i_qlya(:,1:n+1) = qlya(:,0:n)
    i_qlyb(:,1:n+1) = qlyb(:,0:n)
    i_qv(:,1:n+1)   = qv(:,0:n)
    i_qvT(:,1:n+1)  = qvT(:,0:n)
    i_qvh(:,1:n+1)  = qvh(:,0:n)
    i_qvya(:,1:n+1) = qvya(:,0:n)
    i_qvyb(:,1:n+1) = qvyb(:,0:n)

  END SUBROUTINE getfluxes_vp_2d

  !**********************************************************************************************************************

  SUBROUTINE getheatfluxes_1d(n, dx, qh, qhya, qhyb, qhTa, qhTb, var, T, &
       q, qya, qyb, qTa, qTb,qadv, qadvya, qadvyb, qadvTa, qadvTb,advection)
    ! modified 25/05/10 to include contribution to heat flux from liquid water flux in the presence of ice
    IMPLICIT NONE

    INTEGER(i_d),               INTENT(IN)    :: n
    REAL(r_2),  DIMENSION(1:n), INTENT(IN)    :: dx
    REAL(r_2),  DIMENSION(0:n), INTENT(INOUT) :: qh, q, qadv
    REAL(r_2),  DIMENSION(0:n), INTENT(INOUT) :: qhya, qya, qadvya
    REAL(r_2),  DIMENSION(0:n), INTENT(INOUT) :: qhyb, qyb, qadvyb
    REAL(r_2),  DIMENSION(0:n), INTENT(INOUT) :: qhTa, qTa, qadvTa
    REAL(r_2),  DIMENSION(0:n), INTENT(INOUT) :: qhTb, qTb, qadvTb
    TYPE(vars), DIMENSION(1:n), INTENT(IN)    :: var
    REAL(r_2),  DIMENSION(1:n), INTENT(IN)    :: T
    INTEGER(i_d),               INTENT(IN)    :: advection
    ! Gets heat fluxes qh and partial derivs qhya, qhyb wrt T and S (if unsat) or phi (if sat).

    INTEGER(i_d)          :: i
    REAL(r_2)             :: rdz, keff ! , w
    REAL(r_2), DIMENSION(1:n-1) :: dz
    REAL(r_2) :: dTqwdTa, dTqwdTb, Tqw

    dz(:) = half*(dx(1:n-1)+dx(2:n))

    do i=1, n-1
       rdz = one/dz(i)
       keff = 2_r_2*(var(i)%kth*var(i+1)%kth)/(var(i)%kth*dx(i)+var(i+1)%kth*dx(i+1))
       ! qh(i) = keff*(T(i)-T(i+1)) +(var(i)%phiv-var(i+1)%phiv)*var(i)%lambdav*thousand*rdz &
       !      * ((((T(i)+Tzero) /Tzero)**1.88_r_2+((T(i+1)+Tzero) /Tzero)**1.88_r_2)/two) &
       !      *((var(i)%cvsat+var(i+1)%cvsat)/two)
       qh(i) = keff*(T(i)-T(i+1)) +(var(i)%phiv-var(i+1)%phiv)*var(i)%lambdav*thousand*rdz &
            * (half*(exp(1.88_r_2*log((T(i)+Tzero) /Tzero))+exp(1.88_r_2*log((T(i+1)+Tzero) /Tzero)))) &
            *(half*(var(i)%cvsat+var(i+1)%cvsat))
       if (var(i)%isat==0) then
          ! qhya(i) = rdz*var(i)%lambdav*thousand*var(i)%phivS*((((T(i)+Tzero) /Tzero)**1.88_r_2+ &
          !      ((T(i+1)+Tzero) /Tzero)**1.88_r_2)/two) &
          !      * ((var(i)%cvsat+var(i+1)%cvsat)/two)
          qhya(i) = rdz*var(i)%lambdav*thousand*var(i)%phivS*(half*(exp(1.88_r_2*log((T(i)+Tzero) /Tzero)) + &
               exp(1.88_r_2*log((T(i+1)+Tzero) /Tzero)))) &
               * (half*(var(i)%cvsat+var(i+1)%cvsat))
       else
          qhya(i) = zero
       end if
       if (var(i+1)%isat==0) then
          ! qhyb(i) = -rdz*var(i)%lambdav*thousand*var(i+1)%phivS*((((T(i)+Tzero) /Tzero)**1.88_r_2+ &
          !      ((T(i+1)+Tzero) /Tzero)**1.88_r_2)/two) &
          !      * ((var(i)%cvsat+var(i+1)%cvsat)/two)
          qhyb(i) = -rdz*var(i)%lambdav*thousand*var(i+1)%phivS*(half*(exp(1.88_r_2*log((T(i)+Tzero) /Tzero)) + &
               exp(1.88_r_2*log((T(i+1)+Tzero) /Tzero)))) &
               * (half*(var(i)%cvsat+var(i+1)%cvsat))
       else
          qhyb(i) = zero
       end if
       qhTa(i) = keff
       qhTb(i) = -keff

       ! add advective terms
       if (advection==1) then
!!$                   if (q(i) > zero) then
!!$                       w = (var(i)%kth/dx(i))/(var(i)%kth/dx(i)+var(i+1)%kth/dx(i+1))
!!$                    else
!!$                       w = (var(i)%kth/dx(i))/(var(i)%kth/dx(i)+var(i+1)%kth/dx(i+1))
!!$                    endif
!!$                    qadv(i) = rhow*cswat*q(i)*(w*(T(i)+zero)+(one-w)*(T(i+1)+zero))
!!$                    qadvya(i) =  rhow*cswat*qya(i)*(w*(T(i)+zero)+(one-w)*(T(i+1)+zero))
!!$                    qadvyb(i) =  rhow*cswat*qyb(i)*(w*(T(i)+zero)+(one-w)*(T(i+1)+zero))
!!$                    qadvTa(i) =  rhow*cswat*q(i)*w
!!$                    qadvTb(i) =  rhow*cswat*q(i)*(one-w)
          Tqw  = merge(T(i), T(i+1), q(i)>zero) +zero

          dTqwdTa = merge(one, zero, (q(i)>zero))

          dTqwdTb = merge(zero,one, q(i)>zero)
          qadv(i) = rhow*cswat*q(i)*Tqw
          qadvya(i) =  rhow*cswat*qya(i)*Tqw
          qadvyb(i) =  rhow*cswat*qyb(i)*Tqw

          qadvTa(i) =  rhow*cswat*q(i)*dTqwdTa + rhow*cswat*Tqw*qTa(i)
          qadvTb(i) =  rhow*cswat*q(i)*dTqwdTb  +  rhow*cswat*Tqw*qTb(i)

          qh(i) = qh(i) + qadv(i)
          qhya(i) = qhya(i) + qadvya(i)
          qhyb(i) = qhyb(i) + qadvyb(i)
          qhTa(i) = qhTa(i) + qadvTa(i)
          qhTb(i) = qhTb(i) + qadvTb(i)
       endif
    enddo

    qh(n)   = zero
    qhya(n) = zero
    qhyb(n) = zero
    qhTa(n) = zero
    qhTb(n) = zero

  END SUBROUTINE getheatfluxes_1d

  SUBROUTINE getheatfluxes_2d(dx, i_qh, i_qhya, i_qhyb, i_qhTa, i_qhTb, var, T, &
       i_q,i_qya,i_qyb,i_qTa,i_qTb, &
       i_qadv,i_qadvya, i_qadvyb, i_qadvTa, i_qadvTb, advection)
    ! modified 25/05/10 to include contribution to heat flux from liquid water flux in the presence of ice
    IMPLICIT NONE

    REAL(r_2),    DIMENSION(:,:), INTENT(IN)    :: dx      ! :,1:n
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_qh    ! :,0:n => :,1:n+1
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_qhya  ! :,0:n => :,1:n+1
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_qhyb  ! :,0:n => :,1:n+1
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_qhTa  ! :,0:n => :,1:n+1
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_qhTb  ! :,0:n => :,1:n+1
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_q    ! :,0:n => :,1:n+1
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_qya  ! :,0:n => :,1:n+1
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_qyb  ! :,0:n => :,1:n+1
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_qTa  ! :,0:n => :,1:n+1
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_qTb  ! :,0:n => :,1:n+1
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_qadv    ! :,0:n => :,1:n+1
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_qadvya  ! :,0:n => :,1:n+1
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_qadvyb  ! :,0:n => :,1:n+1
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_qadvTa  ! :,0:n => :,1:n+1
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_qadvTb  ! :,0:n => :,1:n+1
    TYPE(vars),   DIMENSION(:,:), INTENT(IN)    :: var
    REAL(r_2),    DIMENSION(:,:), INTENT(IN)    :: T       ! :,1:n
    INTEGER(i_d),   INTENT(IN)    :: advection
    ! Gets heat fluxes qh and partial derivs qhya, qhyb wrt T and S (if unsat) or phi (if sat).

    INTEGER(i_d)          :: i, n
    REAL(r_2),  DIMENSION(1:size(dx,1))                :: rdz
    REAL(r_2),  DIMENSION(1:size(dx,1))                :: keff
    REAL(r_2),  DIMENSION(1:size(dx,1),1:size(dx,2)-1) :: dz
    REAL(r_2),  DIMENSION(1:size(dx,1),0:size(dx,2))   :: qh   ! :,0:n
    REAL(r_2),  DIMENSION(1:size(dx,1),0:size(dx,2))   :: qhya ! :,0:n
    REAL(r_2),  DIMENSION(1:size(dx,1),0:size(dx,2))   :: qhyb ! :,0:n
    REAL(r_2),  DIMENSION(1:size(dx,1),0:size(dx,2))   :: qhTa ! :,0:n
    REAL(r_2),  DIMENSION(1:size(dx,1),0:size(dx,2))   :: qhTb ! :,0:n
    REAL(r_2),  DIMENSION(1:size(dx,1),0:size(dx,2))   :: q   ! :,0:n
    REAL(r_2),  DIMENSION(1:size(dx,1),0:size(dx,2))   :: qya ! :,0:n
    REAL(r_2),  DIMENSION(1:size(dx,1),0:size(dx,2))   :: qyb ! :,0:n
    REAL(r_2),  DIMENSION(1:size(dx,1),0:size(dx,2))   :: qTa ! :,0:n
    REAL(r_2),  DIMENSION(1:size(dx,1),0:size(dx,2))   :: qTb ! :,0:n
    REAL(r_2),  DIMENSION(1:size(dx,1),0:size(dx,2))   :: qadv   ! :,0:n
    REAL(r_2),  DIMENSION(1:size(dx,1),0:size(dx,2))   :: qadvya ! :,0:n
    REAL(r_2),  DIMENSION(1:size(dx,1),0:size(dx,2))   :: qadvyb ! :,0:n
    REAL(r_2),  DIMENSION(1:size(dx,1),0:size(dx,2))   :: qadvTa ! :,0:n
    REAL(r_2),  DIMENSION(1:size(dx,1),0:size(dx,2))   :: qadvTb ! :,0:n
    REAL(r_2),  DIMENSION(1:size(dx,1)) :: dTqwdTa, dTqwdTb, Tqw

    n = size(dx,2)
    dz(:,1:n-1) = half*(dx(:,1:n-1)+dx(:,2:n))
    qh(:,0:n)   = i_qh(:,1:n+1)
    qhya(:,0:n) = i_qhya(:,1:n+1)
    qhyb(:,0:n) = i_qhyb(:,1:n+1)
    qhTa(:,0:n) = i_qhTa(:,1:n+1)
    qhTb(:,0:n) = i_qhTb(:,1:n+1)

    q(:,0:n)   = i_q(:,1:n+1)
    qya(:,0:n) = i_qya(:,1:n+1)
    qyb(:,0:n) = i_qyb(:,1:n+1)
    qTa(:,0:n) = i_qTa(:,1:n+1)
    qTb(:,0:n) = i_qTb(:,1:n+1)

    qadv(:,0:n)   = i_qadv(:,1:n+1)
    qadvya(:,0:n) = i_qadvya(:,1:n+1)
    qadvyb(:,0:n) = i_qadvyb(:,1:n+1)
    qadvTa(:,0:n) = i_qadvTa(:,1:n+1)
    qadvTb(:,0:n) = i_qadvTb(:,1:n+1)

    do i=1, n-1
       rdz =  one/dz(:,i)
       keff = 2_r_2*(var(:,i)%kth*var(:,i+1)%kth)/(var(:,i)%kth*dx(:,i)+var(:,i+1)%kth*dx(:,i+1))
       ! qh(:,i) = keff*(T(:,i)-T(:,i+1)) &
       !      + (var(:,i)%phiv-var(:,i+1)%phiv)*var(:,i)%lambdav*thousand*rdz(:) &
       !      * ((((T(:,i)+Tzero) /Tzero)**1.88_r_2+((T(:,i+1)+Tzero) /Tzero)**1.88_r_2)/two) &
       !      *((var(:,i)%cvsat+var(:,i+1)%cvsat)/two)
       qh(:,i) = keff*(T(:,i)-T(:,i+1)) &
            + (var(:,i)%phiv-var(:,i+1)%phiv)*var(:,i)%lambdav*thousand*rdz(:) &
            * (half*exp(1.88_r_2*log(((T(:,i)+Tzero) /Tzero)) + exp(1.88_r_2*log((T(:,i+1)+Tzero) /Tzero)))) &
            *(half*(var(:,i)%cvsat+var(:,i+1)%cvsat))
       where (var(:,i)%isat == 0)
          ! qhya(:,i) = rdz(:)*var(:,i)%lambdav*thousand*var(:,i)%phivS*((((T(:,i)+Tzero) /Tzero)**1.88_r_2+ &
          !      ((T(:,i+1)+Tzero) /Tzero)**1.88_r_2)/two) &
          !      * ((var(:,i)%cvsat+var(:,i+1)%cvsat)/two)
          qhya(:,i) = rdz(:)*var(:,i)%lambdav*thousand*var(:,i)%phivS*(half*(exp(1.88_r_2*log((T(:,i)+Tzero) /Tzero)) + &
               exp(1.88_r_2*log((T(:,i+1)+Tzero) /Tzero)))) &
               * (half*(var(:,i)%cvsat+var(:,i+1)%cvsat))
       elsewhere
          qhya(:,i) = zero
       endwhere
       where (var(:,i+1)%isat == 0)
          ! qhyb(:,i) = -rdz(:)*var(:,i)%lambdav*thousand*var(:,i+1)%phivS*((((T(:,i)+Tzero) /Tzero)**1.88_r_2+ &
          !      ((T(:,i+1)+Tzero) /Tzero)**1.88_r_2)/two) &
          !      * ((var(:,i)%cvsat+var(:,i+1)%cvsat)/two)
          qhyb(:,i) = -rdz(:)*var(:,i)%lambdav*thousand*var(:,i+1)%phivS*(half*exp(1.88_r_2*log(((T(:,i)+Tzero) /Tzero)) + &
               exp(1.88_r_2*log((T(:,i+1)+Tzero) /Tzero)))) &
               * (half*(var(:,i)%cvsat+var(:,i+1)%cvsat))
       elsewhere
          qhyb(:,i) = zero
       endwhere
       qhTa(:,i) = keff
       qhTb(:,i) = -keff

       ! add advective terms
       if (advection==1) then

          Tqw  = merge(T(:,i), T(:,i+1), q(:,i)>zero) +zero
          dTqwdTa = merge(one, zero, q(:,i)>zero)
          dTqwdTb = merge(zero,one, q(:,i)>zero)

          qadv(:,i) = rhow*cswat*q(:,i)*Tqw
          qadvya(:,i) =  rhow*cswat*qya(:,i)*Tqw
          qadvyb(:,i) =  rhow*cswat*qyb(:,i)*Tqw

          qadvTa(:,i) =  rhow*cswat*q(:,i)*dTqwdTa + rhow*cswat*Tqw*qTa(:,i)
          qadvTb(:,i) =  rhow*cswat*q(:,i)*dTqwdTb  +  rhow*cswat*Tqw*qTb(:,i)

          qh(:,i) = qh(:,i) + qadv(:,i)
          qhya(:,i) = qhya(:,i) + qadvya(:,i)
          qhyb(:,i) = qhyb(:,i) + qadvyb(:,i)
          qhTa(:,i) = qhTa(:,i) + qadvTa(:,i)
          qhTb(:,i) = qhTb(:,i) + qadvTb(:,i)
       endif
    enddo

    qh(:,n)   = zero
    qhya(:,n) = zero
    qhyb(:,n) = zero
    qhTa(:,n) = zero
    qhTb(:,n) = zero

    i_qh(:,1:n+1)   = qh(:,0:n)
    i_qhya(:,1:n+1) = qhya(:,0:n)
    i_qhyb(:,1:n+1) = qhyb(:,0:n)
    i_qhTa(:,1:n+1) = qhTa(:,0:n)
    i_qhTb(:,1:n+1) = qhTb(:,0:n)

    if (advection==1) then
       i_qadv(:,1:n+1)   = qadv(:,0:n)
       i_qadvya(:,1:n+1) = qadvya(:,0:n)
       i_qadvyb(:,1:n+1) = qadvyb(:,0:n)
       i_qadvTa(:,1:n+1) = qadvTa(:,0:n)
       i_qadvTb(:,1:n+1) = qadvTb(:,0:n)
    endif

  END SUBROUTINE getheatfluxes_2d

  !**********************************************************************************************************************

  ELEMENTAL PURE SUBROUTINE hyofh(h, lam, eta, Ke, he, K, Kh, phi)

    IMPLICIT NONE

    REAL(r_2),    INTENT(IN)  :: h
    REAL(r_2),    INTENT(IN)  :: lam
    REAL(r_2),    INTENT(IN)  :: eta
    REAL(r_2),    INTENT(IN)  :: Ke
    REAL(r_2),    INTENT(IN)  :: he
    REAL(r_2),    INTENT(OUT) :: K
    REAL(r_2),    INTENT(OUT) :: Kh
    REAL(r_2),    INTENT(OUT) :: phi
    ! Get soil water variables from h.
    ! Definitions of arguments:
    ! h   - matric head.
    ! K   - hydraulic conductivity.
    ! Kh  - derivative dK/dh.
    ! phi - matric flux potential (MFP).
    REAL(r_2) :: a

    a   =  -lam * eta
    K   =  Ke * exp(a*log(h/he))
    Kh  =  a * K / h
    phi =  K * h / (one+a)

  END SUBROUTINE hyofh

  !**********************************************************************************************************************

  ! For debug: remove elemental pure
  ELEMENTAL PURE SUBROUTINE hyofS(S, Tsoil, parin, var)
  !  SUBROUTINE hyofS(S, Tsoil, parin, var)

    IMPLICIT NONE

    REAL(r_2),    INTENT(IN)    :: S
    REAL(r_2),    INTENT(IN)    :: Tsoil
    TYPE(params), INTENT(INOUT)    :: parin
    TYPE(vars),   INTENT(INOUT) :: var
    ! Get soil water variables from S.
    ! Definitions of arguments:
    ! S(1:ms)   - degree of saturation ("effective satn") of layers.
    ! ms        - no. of soil layers.
    ! jt(1:ms)  - layer soil type nos.
    ! var(1:ms) - other water vars of layers.
    REAL(r_2) :: lnS, theta, c, v3, v4
    REAL(r_2) :: dhdS, lambda, crh, int
    REAL(r_2) :: thetal_max
    REAL(r_2) :: Sliq
    REAL(r_2) :: A, B, D, C1
    INTEGER(i_d) :: E
    REAL(r_2) :: F1, F2, F
    ! REAL(r_2) :: macropore_modifier
    REAL(r_2) :: cdry, tmp_thetai
    REAL(r_2), parameter :: tol = 1.e-6_r_2
    theta         = S*(parin%thre) + (parin%the - parin%thre)
    var%lambdav   = rlambda       ! latent heat of vaporisation
    var%lambdav   = 1.91846e6_r_2*((Tsoil+Tzero)/((Tsoil+Tzero)-33.91_r_2))**2  ! Henderson-Sellers, QJRMS, 1984
    var%lambdaf   = lambdaf        ! latent heat of fusion
    !var%Tfrz      = Tfrz(S,parin%he,one/parin%lam)
    var%Tfrz      = Tfrz(S, parin%he, one/(parin%lambc*freezefac)) ! freezefac for test of steep freezing curve

    if ((Tsoil < var%Tfrz-tol) .and. (experiment/=184)) then ! ice
       parin%lam     = parin%lambc * freezefac   ! freezefac>1 -> steeper freezing curve
       parin%eta     = two/parin%lam + two + one ! freezefac>1 -> steeper freezing curve
       thetal_max    = thetalmax(Tsoil,S,parin%he,one/parin%lam,parin%thre,parin%the)
       var%dthetaldT = dthetalmaxdT(Tsoil,S,parin%he,one/parin%lam,parin%thre,parin%the)
       var%iice   = 1
       var%thetai = max((theta - thetal_max),zero) ! volumetric ice content (m3(liq H2O)/m3 soil)
       tmp_thetai = max(min(theta, parin%thre) - thetal_max,zero)
       var%thetal = thetal_max
       ! liquid water content, relative to saturation
       ! Sliq      = (var%thetal - (parin%the-parin%thre))/parin%thre
       if ((parin%thre-tmp_thetai) .le. max(parin%thr,1.e-5_r_2)) then
          Sliq = max(parin%thr,1.e-5_r_2)
       else
         ! Sliq = min((var%thetal-(parin%the-parin%thre))/(parin%thre-var%thetai), one)
          Sliq = min((var%thetal-(parin%the-parin%thre))/(parin%thre-tmp_thetai), one)
       endif
       var%Sliq = Sliq
       ! saturated liquid water content (< 1 for frozen soil)
       ! air entry potential, flux matric potential and hydraulic conductivity for saturated frozen soil
       var%he    = parin%he
       lnS        = log(Sliq)
       v3         = exp(-lnS/parin%lam)
       v4         = exp(parin%eta*lnS)
       var%phie  = parin%phie*v3*v4
       var%Ksat  = parin%Ke*v4

       var%h  = parin%he*v3  ! matric potential
       dhdS   = zero
       var%K  = var%Ksat
       var%KS = zero
       ! var%KT = var%dthetaldT * parin%Ke * parin%eta * exp(lnS*(parin%eta-one))/parin%thre
       var%KT = var%dthetaldT * parin%Ke * parin%eta * exp(lnS*(parin%eta-one))/max(parin%thre-var%thetai,max(parin%thr,1e-5_r_2))
       if (S < one) var%phi = var%phie

       ! var%phiT = parin%phie * exp(lnS*(parin%eta-one/parin%lam-one)) * var%dthetaldT * &
       !      (parin%eta-one/parin%lam)/(parin%thre)
       ! var%phiT = parin%phie * exp(lnS*(parin%eta-one/parin%lam-one)) * var%dthetaldT * &
       !       (parin%eta-one/parin%lam)/(parin%thre-var%thetai)
       ! var%phiT =  -(((-1 +parin%eta*parin%lam)*parin%phie*(theta - parin%thre)* &
       !      (var%thetal/(-theta + var%thetal + parin%thre))**(-1 + parin%eta - 1./parin%lam))/ &
       !      (parin%lam*(-theta + var%thetal + parin%thre)**2))* var%dthetaldT
       if (var%isat==0) then
          var%phiT =  -( ( (-one +parin%eta*parin%lam) * parin%phie * (min(theta,parin%thre) - parin%thre) &
            * exp( (-one + parin%eta - one/parin%lam) &
            * log(var%thetal/(-min(theta,parin%thre) + var%thetal + parin%thre)) ) ) &
            / (parin%lam*(-min(theta,parin%thre) + var%thetal + parin%thre)**2) ) * var%dthetaldT
      ! if (S < one) then
          ! var%phiS = ((parin%eta - 1./parin%lam)*parin%phie* &
          !     (Sliq)**(1 + parin%eta - 1/parin%lam))/var%thetal*parin%thre
          ! var%phiS =((-1. + parin%eta*parin%lam)*parin%phie* &
          !      (var%thetal/(-theta + var%thetal + parin%thre))**(1. + parin%eta - 1./parin%lam))/ &
          !      (parin%lam*var%thetal)*parin%thre
          var%phiS =((-one + parin%eta*parin%lam)*parin%phie* &
               exp((one + parin%eta - one/parin%lam)*log(var%thetal/(-min(theta,parin%thre) + var%thetal + parin%thre)))) / &
               (parin%lam*var%thetal)*parin%thre
       else
          var%phiS = zero
          var%phiT = zero
       endif

       !if ( ((parin%thre-var%thetai)<=max(parin%thr,1e-5_r_2)) .or. (Sliq==one) ) then
       !   var%phiS = zero
       !   var%phiT = zero
       !endif

       var%rh   = max(exp(Mw*gravity*var%h/Rgas/(Tsoil+Tzero)),rhmin)
    else ! no ice
       parin%lam  = parin%lambc               ! freezefac>1 -> steeper freezing curve
       parin%eta  = two/parin%lam + two + one ! freezefac>1 -> steeper freezing curve
       var%he     = parin%he
       var%phie   = parin%phie
       var%Ksat   = parin%Ke
       var%iice   = 0
       var%thetai = zero
       var%thetal = S*(parin%thre) + (parin%the - parin%thre)
       var%dthetaldT = zero
       dhdS   = zero
       Sliq    = S        ! liquid water content, relative to saturation
       lnS     = log(Sliq)
       v3      = exp(-lnS/parin%lam)
       v4      = exp(parin%eta*lnS)
    endif
    if (var%thetal < 1.e-12_r_2) then ! completely frozen
       var%dthetaldT = zero
       var%lambdav   = lambdas ! latent heat of sublimation
       var%KT        = zero
    endif
    if (((Tsoil >= var%Tfrz-tol) .or. (experiment==184)) .and. (S < one)) then ! unsaturated
       var%h    = parin%he*v3  ! matric potential
       ! dhdS     = -parin%he/parin%lam*S**(-one/parin%lam-one)
       dhdS     = -parin%he/parin%lam*exp((-one/parin%lam-one)*log(S))
       var%K    = parin%Ke*v4
       var%phi  = parin%phie*v3*v4
       var%KS   = parin%eta*var%K/S
       var%KT   = zero
       var%phiS = (parin%eta-one/parin%lam)*var%phi/S
       var%phiT = zero
       var%rh   = max(exp(Mw*gravity*var%h/Rgas/(Tsoil+Tzero)),rhmin)
    endif
    if (((Tsoil >= var%Tfrz-tol) .or. (experiment==184)) .and. (S >= one)) then ! saturated
       var%h    = parin%he
       dhdS     = zero
       var%K    = parin%Ke
       var%phi  = parin%phie
       var%KS   = parin%eta*parin%Ke
       var%phiS = (parin%eta-one/parin%lam)*parin%phie
       var%phiT = zero
       var%rh   = one
       var%KT   = zero
    endif

    !  variables required for vapour phase transfer
    theta  =  min(S,one)*(parin%thre) + (parin%the - parin%thre)

    ! if (z.lt.3.and.theta>parin%thfc) then
    !    macropore_modifier = exp(-parin%zeta*(z-3))
    ! ! if (theta > parin%thfc) then
    ! !    macropore_modifier = 1000._r_2*exp(-var%h/(-4._r_2))
    !    var%macropore_factor = macropore_modifier
    ! else
    !    var%macropore_factor = one
    ! endif
    var%macropore_factor = one

    var%sl = slope_esat(Tsoil) * Mw/thousand/Rgas/(Tsoil+Tzero) ! m3 m-3 K-1
    c      = Mw*gravity/Rgas/(Tsoil+Tzero)
    lambda = parin%lam
    if (S < one) then
       var%rh = max(exp(Mw*gravity*var%h/Rgas/(Tsoil+Tzero)),rhmin)
    else
       var%rh = one
    endif
    crh        = c*var%rh
    var%hS     = dhdS
    var%rhS    = crh*dhdS
    var%cvsat  = esat(Tsoil)*Mw/thousand/Rgas/(Tsoil+Tzero) ! m3 m-3
    if (S < one) then
       var%cv     = var%rh*var%cvsat
       var%cvS    = var%rhS *var%cvsat
       var%cvsatT = slope_esat(Tsoil)*Mw/thousand/Rgas/(Tsoil+Tzero) ! m3 m-3 K-1
    else
       var%cv = zero
       var%cvS    = zero
       var%cvsatT = zero
    endif
    ! Penman (1940): tortuosity*theta
    ! var%Dv    = Dva*parin%tortuosity*(parin%the-theta)  * ((Tsoil+Tzero)/Tzero)**1.88_r_2 ! m2 s-1
    var%Dv    = Dva*parin%tortuosity*(parin%the-theta)  * exp(1.88_r_2*log((Tsoil+Tzero)/Tzero)) ! m2 s-1
    ! ! Moldrup et al. (2003), Table 2 repacked soil (from Moldrup et al. (2000)
    ! ! (thetas-theta)**1.5 * (thetas-theta)/thetas
    ! var%Dv    = Dva* (parin%the-theta)**1.5_r_2 * (parin%the-theta)/parin%the * ((Tsoil+Tzero)/Tzero)**1.88_r_2 ! m2 s-1
    ! int       = (-c*parin%he)**lambda * igamma(one-lambda,-c*var%h)
    int       = exp(lambda*log(-c*parin%he)) * igamma(one-lambda,-c*var%h)
    var%phiv  = Dva*parin%tortuosity * (parin%thre*exp(c*var%h) -parin%thre*int)
    ! var%phivS = dhdS * Dva*parin%tortuosity * &
    !      ( (parin%thre)*c*exp(c*var%h) - parin%thre*c*(var%h/parin%he)**(-parin%lam)*exp(c*var%h) )
    var%phivS = dhdS * Dva*parin%tortuosity * &
         ( (parin%thre)*c*exp(c*var%h) - parin%thre*c*exp(-parin%lam*log(var%h/parin%he))*exp(c*var%h) )
    var%kv    = var%Dv * var%cvsat *c * var%rh

    select case (ithermalcond)
    case (0)
       select case (experiment)
       case (11,17) ! Hansson et al. (2004) - special soil properties
          ! Hansson et al. (2004) - Eq. 13b
          ! A=C1, B=C2, C1=C3, D=C4, E=C5
          A  = 0.55_r_2
          B  = 0.8_r_2
          C1 = 3.07_r_2
          D  = 0.13_r_2
          E  = 4
       case default
          ! calculate v%kH as in Campbell (1985) p.32 eq. 4.20
          A  = 0.65_r_2 - 0.78_r_2*parin%rho/thousand + 0.60_r_2*(parin%rho/thousand)**2 ! (4.27)
          B  = 2.8_r_2 * (one-parin%thre) !*theta   ! (4.24)
          if (parin%clay > 0.001_r_2) then ! clay=0.001 -> C1=9.2
             C1 = one + 2.6_r_2/sqrt(parin%clay*100._r_2) ! (4.28)
          else
             C1 = one
          endif
          D  = 0.03_r_2 + 0.7_r_2*(one-parin%thre)**2 ! (4.22)
          E  = 4
       end select
       select case (experiment)
       case (15,17,182:184) ! ice = liquid thermal conductivity
          ! calculate v%kH as in Campbell (1985) p.32 eq. 4.20
          var%kH = A + B*theta-(A-D)*exp(-(C1*theta)**E) ! (4.20)
       case default ! different thermal conductivity in ice and water
          ! Hansson et al. (2004) - Eq. 15
          F1 = 13.05_r_2
          F2 = 1.06_r_2
          if  (Tsoil < var%Tfrz-tol .and.  var%thetai.gt.zero ) then ! ice
             ! F  = one + F1*var%thetai**F2
             F  = one + F1*exp(F2*log(var%thetai))
             if ((C1*(theta+F*var%thetai))**E > 100._r_2) then
                var%kH = A + B*(theta+F*var%thetai)
             else
                var%kH = A + B*(theta+F*var%thetai)-(A-D)*exp(-(C1*(theta+F*var%thetai))**E)
             endif
          else
             var%kH = A + B*(theta)-(A-D)*exp(-(C1*(theta))**E)
          endif
       end select
    case(1) ! van de Griend & O'Neill (1986)
       cdry = 2.0e6_r_2 * (one-parin%thre) ! Sispat Manual (2000), Eq. 2.21
       var%kH = one/(cdry + 4.18e6_r_2*theta) * ( (parin%LambdaS + 2300._r_2*theta - 1890._r_2)/0.654_r_2 )**2
    end select
    var%eta_th = one
    var%kE     = var%Dv*var%rh*var%sl*thousand*var%lambdav*var%eta_th
    var%kth    = var%kE + var%kH ! thermal conductivity of soil (includes contribution from vapour phase)

    var%csoil    = parin%css*parin%rho + rhow*cswat*var%thetal + merge(0._r_2, rhow*csice*var%thetai, experiment==183)
    if (Tsoil < var%Tfrz) then ! increase effective heat capacity due to presence of ice
       var%csoileff = var%csoil + merge(0._r_2, rhow*lambdaf*var%dthetaldT, experiment==183)
    else
       var%csoileff = var%csoil
    endif

  END SUBROUTINE hyofS

  !**********************************************************************************************************************

  SUBROUTINE isosub(iso, c, p, f, fd)

    IMPLICIT NONE

    CHARACTER(LEN=2),       INTENT(IN)    :: iso
    REAL(r_2),               INTENT(IN)    :: c
    REAL(r_2), DIMENSION(:), INTENT(INOUT) :: p
    REAL(r_2),               INTENT(OUT)   :: f, fd
    ! Subroutine to get adsorbed solute (units/g soil) from concn in soil water
    ! according to chosen isotherm code ("Fr" for Freundlich, "La" for Langmuir
    ! and "Ll" for Langmuir-linear).
    ! Definitions of arguments:
    ! iso  - 2 character code.
    ! c    - concn in soil water.
    ! p(:) - isotherm parameters.
    ! f    - adsorbed mass/g soil.
    ! fc   - deriv of f wrt c (slope of isotherm curve).
    REAL(r_2) :: x

    select case (iso)
    case ("Fr")
       if (eq(p(3),zero)) then ! linearise near zero
          p(3) = (0.01_r_2*dsmmax/p(1))**(one/p(2)) ! concn at 0.01*dsmmax
          p(4) = p(1)*p(3)**(p(2)-one) ! slope
       end if
       if (c < p(3)) then
          fd = p(4)
          f  = fd*c
       else
          x  = p(1)*exp((p(2)-one)*log(c))
          f  = x*c
          fd = p(2)*x
       end if
    case ("La")
       x  = one/(one+p(2)*c)
       f  = p(1)*c*x
       fd = p(1)*(x-p(2)*c*x**2)
    case ("Ll")
       x  = one/(one+p(2)*c)
       f  = p(1)*c*x+p(3)*c
       fd = p(1)*(x-p(2)*c*x**2)+p(3)
    case default
       write(*,*) "isosub: illegal isotherm type"
#ifdef __MPI__
       call MPI_Abort(0, 148, ierr) ! Do not know comm nor rank here
#else
       stop 148
#endif
    end select

  END SUBROUTINE isosub

  !**********************************************************************************************************************

  SUBROUTINE massman_sparse_1d(aa, aah, bb, bbh, cc, cch, dd, ddh, ee, eeh, ff, ffh, gg, ggh, dy, dT, condition, err)

    USE cable_def_types_mod, ONLY: r_2, i_d
    USE sli_numbers,       ONLY: zero, one

    IMPLICIT NONE

    ! in/out
    REAL(r_2), DIMENSION(:), INTENT(IN)  :: aa, aah, bb, bbh, ee, eeh, ff, ffh
    REAL(r_2), DIMENSION(:), INTENT(IN)  :: cc, cch, dd, ddh, gg, ggh
    REAL(r_2), DIMENSION(:), INTENT(OUT) :: dy, dT
    INTEGER(i_d),  OPTIONAL, INTENT(IN)  :: condition
    INTEGER(i_d),  OPTIONAL, INTENT(OUT) :: err ! 0: no error; >0: error
    ! local
    INTEGER(i_d)                         :: n, n2
    REAL(r_2), DIMENSION(size(cc),2,2)   :: A, B, C
    REAL(r_2), DIMENSION(size(cc),2)     :: d, x
    ! for conditioning
    REAL(r_2), DIMENSION(2*size(cc))       :: lST, cST
    REAL(r_2), DIMENSION(size(cc))         :: lS, lT, cS, cT
    REAL(r_2), DIMENSION(2*size(cc)*2*size(cc)) :: allvec
    REAL(r_2), DIMENSION(2*size(cc),2*size(cc)) :: allmat
    REAL(r_2)    :: eps
    INTEGER(i_d) :: docond ! 0: no conditioning, 1: columns, 2: lines, 3: both
    INTEGER(i_d) :: iierr ! error code for generic_thomas
    ! CHARACTER(LEN=20) :: form1
    ! integer :: i, nn
    !
    ! check input sizes
    if (.not. all((/size(aa)+1,size(bb)+1,size(dd),size(ee)+1,size(ff)+1,size(gg)/) == size(cc))) then
       write(*,*) 'massman_sparse_1d error1: unequal humidity coeffs.'
#ifdef __MPI__
       call MPI_Abort(0, 149, ierr) ! Do not know comm nor rank here
#else
       stop 149
#endif
    end if
    if (.not. all((/size(aah)+1,size(bbh)+1,size(ddh),size(eeh)+1,size(ffh)+1,size(ggh)/) == size(cch))) then
       write(*,*) 'massman_sparse_1d error2: unequal temperature coeffs.'
#ifdef __MPI__
       call MPI_Abort(0, 150, ierr) ! Do not know comm nor rank here
#else
       stop 150
#endif
    end if
    if (size(cc) /= size(cch)) then
       write(*,*) 'massman_sparse_1d error3: unequal temperature and humidity coeffs.'
#ifdef __MPI__
       call MPI_Abort(0, 151, ierr) ! Do not know comm nor rank here
#else
       stop 151
#endif
    end if
    n = size(cc)
    if (present(condition)) then
       docond = condition
    else
       docond = 0
    endif
    if (present(err)) err = 0
    iierr = 0
    !
    ! Overall matrix
    if (docond >= 1 .and. docond <= 3) then
       eps = epsilon(one)
       n2  = 2*n*2*n
       allvec(:) = zero
       allvec(1:n2:4*n+2)         = cc(1:n)
       allvec(2:n2:4*n+2)         = dd(1:n)
       allvec(3:n2-4*n:4*n+2)     = ee(1:n-1)
       allvec(4:n2-4*n:4*n+2)     = ff(1:n-1)
       allvec(2*n+1:n2:4*n+2)     = cch(1:n)
       allvec(2*n+2:n2:4*n+2)     = ddh(1:n)
       allvec(2*n+3:n2-4*n:4*n+2) = eeh(1:n-1)
       allvec(2*n+4:n2-4*n:4*n+2) = ffh(1:n-1)
       allvec(4*n+1:n2:4*n+2)     = aa(1:n-1)
       allvec(4*n+2:n2:4*n+2)     = bb(1:n-1)
       allvec(6*n+1:n2:4*n+2)     = aah(1:n-1)
       allvec(6*n+1:n2:4*n+2)     = bbh(1:n-1)
       allmat = reshape(allvec,shape=(/2*n,2*n/),order=(/2,1/))
    endif
    ! Get conditioning numbers
    select case (docond)
    case (1)
       cST = maxval(abs(allmat),1)
       where (cST < eps) cST = one
       cST     = one / cST
       cS(1:n) = cST(1:2*n-1:2)
       cT(1:n) = cST(2:2*n:2)
       lS(1:n) = one
       lT(1:n) = one
    case (2)
       lST = maxval(abs(allmat),2)
       where (lST < eps) lST = one
       lST     = one / lST
       lS(1:n) = lST(1:2*n-1:2)
       lT(1:n) = lST(2:2*n:2)
       cS(1:n) = one
       cT(1:n) = one
    case (3)
       lST = maxval(abs(allmat),2)
       where (lST < eps) lST = one
       lST     = one / lST
       lS(1:n) = lST(1:2*n-1:2)
       lT(1:n) = lST(2:2*n:2)
       allmat  = spread(lST,2,2*n)*allmat
       cST = maxval(abs(allmat),1)
       where (cST < eps) cST = one
       cST     = one / cST
       cS(1:n) = cST(1:2*n-1:2)
       cT(1:n) = cST(2:2*n:2)
    case default
       cS(1:n) = one
       cT(1:n) = one
       lS(1:n) = one
       lT(1:n) = one
    end select
    !
    ! fill matrices of generic thomas algorithm
    A(1,1:2,1:2) = zero
    A(2:n,1,1)   = aa(1:n-1)  * lS(2:n)   * cS(1:n-1)
    A(2:n,1,2)   = bb(1:n-1)  * lS(2:n)   * cT(1:n-1)
    A(2:n,2,1)   = aah(1:n-1) * lT(2:n)   * cS(1:n-1)
    A(2:n,2,2)   = bbh(1:n-1) * lT(2:n)   * cT(1:n-1)
    B(1:n,1,1)   = cc(1:n)    * lS(1:n)   * cS(1:n)
    B(1:n,1,2)   = dd(1:n)    * lS(1:n)   * cT(1:n)
    B(1:n,2,1)   = cch(1:n)   * lT(1:n)   * cS(1:n)
    B(1:n,2,2)   = ddh(1:n)   * lT(1:n)   * cT(1:n)
    C(1:n-1,1,1) = ee(1:n-1)  * lS(1:n-1) * cS(2:n)
    C(1:n-1,1,2) = ff(1:n-1)  * lS(1:n-1) * cT(2:n)
    C(1:n-1,2,1) = eeh(1:n-1) * lT(1:n-1) * cS(2:n)
    C(1:n-1,2,2) = ffh(1:n-1) * lT(1:n-1) * cT(2:n)
    C(n,1:2,1:2) = zero
    d(1:n,1)     = gg(1:n)    * lS(1:n)
    d(1:n,2)     = ggh(1:n)   * lT(1:n)
    !
    ! Call Generic Thomas algorithm
    call generic_thomas(n,A,B,C,d,x,iierr)
    if (iierr /= 0) then
       if (present(err)) then
          err = 1
          return
       else
#ifdef __MPI__
          call MPI_Abort(0, 152, ierr) ! Do not know comm nor rank here
#else
          stop 152
#endif
       endif
    endif
    dy(1:n) = x(1:n,1) * cS(1:n)
    dT(1:n) = x(1:n,2) * cT(1:n)
    !
  END SUBROUTINE massman_sparse_1d

  SUBROUTINE massman_sparse_2d(aa, aah, bb, bbh, cc, cch, dd, ddh, ee, eeh, ff, ffh, gg, ggh, dy, dT, condition, err)

    USE cable_def_types_mod, ONLY: r_2, i_d
    USE sli_numbers,       ONLY: zero, one

    IMPLICIT NONE

    ! in/out
    REAL(r_2), DIMENSION(:,:), INTENT(IN)  :: aa, aah, bb, bbh, ee, eeh, ff, ffh
    REAL(r_2), DIMENSION(:,:), INTENT(IN)  :: cc, cch, dd, ddh, gg, ggh
    REAL(r_2), DIMENSION(:,:), INTENT(OUT) :: dy, dT
    INTEGER(i_d),    OPTIONAL, INTENT(IN)  :: condition
    INTEGER(i_d),    OPTIONAL, INTENT(OUT) :: err ! 0: no error; >0: error
    ! local
    INTEGER(i_d)                                        :: n, mp
    INTEGER(i_d)                                        :: n2
    REAL(r_2), DIMENSION(1:size(cc,1),1:size(cc,2),2,2) :: A, B, C
    REAL(r_2), DIMENSION(1:size(cc,1),1:size(cc,2),2)   :: d, x
    ! for conditioning
    REAL(r_2), DIMENSION(1:size(cc,1),2*size(cc,2))     :: lST, cST
    REAL(r_2), DIMENSION(1:size(cc,1),size(cc,2))       :: lS, lT, cS, cT
    REAL(r_2), DIMENSION(1:size(cc,1),2*size(cc,2)*2*size(cc,2)) :: allvec
    REAL(r_2), DIMENSION(1:size(cc,1),2*size(cc,2),2*size(cc,2)) :: allmat
    REAL(r_2)    :: eps
    INTEGER(i_d) :: docond ! 0: no conditioning, 1: columns, 2: lines, 3: both
    INTEGER(i_d) :: iierr ! error code for generic_thomas
    ! CHARACTER(LEN=20) :: form1
    ! integer :: i, k, nn
    !
    ! check input sizes
    if (.not. all((/size(aa,1)+1,size(bb,1)+1,size(dd,1),size(ee,1)+1,size(ff,1)+1,size(gg,1)/) == size(cc,1))) then
       write(*,*) 'massman_sparse_2d error1: unequal humidity coeffs (1st dim).'
#ifdef __MPI__
       call MPI_Abort(0, 153, ierr) ! Do not know comm nor rank here
#else
       stop 153
#endif
    end if
    if (.not. all((/size(aah,1)+1,size(bbh,1)+1,size(ddh,1),size(eeh,1)+1,size(ffh,1)+1,size(ggh,1)/) == size(cch,1))) then
       write(*,*) 'massman_sparse_2d error2: unequal temperature coeffs (1st dim).'
#ifdef __MPI__
       call MPI_Abort(0, 154, ierr) ! Do not know comm nor rank here
#else
       stop 154
#endif
    end if
    if (size(cc,1) /= size(cch,1)) then
       write(*,*) 'massman_sparse_2d error3: unequal temperature and humidity coeffs (1st dim).'
#ifdef __MPI__
       call MPI_Abort(0, 155, ierr) ! Do not know comm nor rank here
#else
       stop 155
#endif
    end if
    if (.not. all((/size(aa,2)+1,size(bb,2)+1,size(dd,2),size(ee,2)+1,size(ff,2)+1,size(gg,2)/) == size(cc,2))) then
       write(*,*) 'massman_sparse_2d error4: unequal humidity coeffs (2nd dim).'
#ifdef __MPI__
       call MPI_Abort(0, 156, ierr) ! Do not know comm nor rank here
#else
       stop 156
#endif
    end if
    if (.not. all((/size(aah,2)+1,size(bbh,2)+1,size(ddh,2),size(eeh,2)+1,size(ffh,2)+1,size(ggh,2)/) == size(cch,2))) then
       write(*,*) 'massman_sparse_2d error5: unequal temperature coeffs (2nd dim).'
#ifdef __MPI__
       call MPI_Abort(0, 157, ierr) ! Do not know comm nor rank here
#else
       stop 157
#endif
    end if
    if (size(cc,2) /= size(cch,2)) then
       write(*,*) 'massman_sparse_2d error6: unequal temperature and humidity coeffs (2nd dim).'
#ifdef __MPI__
       call MPI_Abort(0, 158, ierr) ! Do not know comm nor rank here
#else
       stop 158
#endif
    end if

    mp   = size(cc,1)
    n    = size(cc,2)
    if (present(condition)) then
       docond = condition
    else
       docond = 0
    endif
    if (present(err)) err = 0
    iierr = 0
    !
    ! Overall matrix
    if (docond >= 1 .and. docond <= 3) then
       eps = epsilon(one)
       n2  = 2*n*2*n
       allvec(1:mp,:) = zero
       allvec(1:mp,1:n2:4*n+2)         = cc(1:mp,1:n)
       allvec(1:mp,2:n2:4*n+2)         = dd(1:mp,1:n)
       allvec(1:mp,3:n2-4*n:4*n+2)     = ee(1:mp,1:n-1)
       allvec(1:mp,4:n2-4*n:4*n+2)     = ff(1:mp,1:n-1)
       allvec(1:mp,2*n+1:n2:4*n+2)     = cch(1:mp,1:n)
       allvec(1:mp,2*n+2:n2:4*n+2)     = ddh(1:mp,1:n)
       allvec(1:mp,2*n+3:n2-4*n:4*n+2) = eeh(1:mp,1:n-1)
       allvec(1:mp,2*n+4:n2-4*n:4*n+2) = ffh(1:mp,1:n-1)
       allvec(1:mp,4*n+1:n2:4*n+2)     = aa(1:mp,1:n-1)
       allvec(1:mp,4*n+2:n2:4*n+2)     = bb(1:mp,1:n-1)
       allvec(1:mp,6*n+1:n2:4*n+2)     = aah(1:mp,1:n-1)
       allvec(1:mp,6*n+1:n2:4*n+2)     = bbh(1:mp,1:n-1)
       allmat = reshape(allvec,shape=(/mp,2*n,2*n/),order=(/1,3,2/))
    endif
    ! Get conditioning numbers
    select case (docond)
    case (1)
       cST = maxval(abs(allmat),2)
       where (cST < eps) cST = one
       cST     = one / cST
       cS(1:mp,1:n) = cST(1:mp,1:2*n-1:2)
       cT(1:mp,1:n) = cST(1:mp,2:2*n:2)
       lS(1:mp,1:n) = one
       lT(1:mp,1:n) = one
    case (2)
       lST = maxval(abs(allmat),3)
       where (lST < eps) lST = one
       lST  = one / lST
       lS(1:mp,1:n) = lST(1:mp,1:2*n-1:2)
       lT(1:mp,1:n) = lST(1:mp,2:2*n:2)
       cS(1:mp,1:n) = one
       cT(1:mp,1:n) = one
    case (3)
       lST = maxval(abs(allmat),3)
       where (lST < eps) lST = one
       lST  = one / lST
       lS(1:mp,1:n) = lST(1:mp,1:2*n-1:2)
       lT(1:mp,1:n) = lST(1:mp,2:2*n:2)
       allmat  = spread(lST,2,2*n)*allmat
       cST = maxval(abs(allmat),2)
       where (cST < eps) cST = one
       cST = one / cST
       cS(1:mp,1:n) = cST(1:mp,1:2*n-1:2)
       cT(1:mp,1:n) = cST(1:mp,2:2*n:2)
    case default
       cS(1:mp,1:n) = one
       cT(1:mp,1:n) = one
       lS(1:mp,1:n) = one
       lT(1:mp,1:n) = one
    end select
    !
    ! fill matrices of generic thomas algorithm
    A(1:mp,1,1:2,1:2) = zero
    A(1:mp,2:n,1,1)   = aa(1:mp,1:n-1)  * lS(1:mp,2:n)   * cS(1:mp,1:n-1)
    A(1:mp,2:n,1,2)   = bb(1:mp,1:n-1)  * lS(1:mp,2:n)   * cT(1:mp,1:n-1)
    A(1:mp,2:n,2,1)   = aah(1:mp,1:n-1) * lT(1:mp,2:n)   * cS(1:mp,1:n-1)
    A(1:mp,2:n,2,2)   = bbh(1:mp,1:n-1) * lT(1:mp,2:n)   * cT(1:mp,1:n-1)
    B(1:mp,1:n,1,1)   = cc(1:mp,1:n)    * lS(1:mp,1:n)   * cS(1:mp,1:n)
    B(1:mp,1:n,1,2)   = dd(1:mp,1:n)    * lS(1:mp,1:n)   * cT(1:mp,1:n)
    B(1:mp,1:n,2,1)   = cch(1:mp,1:n)   * lT(1:mp,1:n)   * cS(1:mp,1:n)
    B(1:mp,1:n,2,2)   = ddh(1:mp,1:n)   * lT(1:mp,1:n)   * cT(1:mp,1:n)
    C(1:mp,1:n-1,1,1) = ee(1:mp,1:n-1)  * lS(1:mp,1:n-1) * cS(1:mp,2:n)
    C(1:mp,1:n-1,1,2) = ff(1:mp,1:n-1)  * lS(1:mp,1:n-1) * cT(1:mp,2:n)
    C(1:mp,1:n-1,2,1) = eeh(1:mp,1:n-1) * lT(1:mp,1:n-1) * cS(1:mp,2:n)
    C(1:mp,1:n-1,2,2) = ffh(1:mp,1:n-1) * lT(1:mp,1:n-1) * cT(1:mp,2:n)
    C(1:mp,n,1:2,1:2) = zero
    d(1:mp,1:n,1)     = gg(1:mp,1:n)    * lS(1:mp,1:n)
    d(1:mp,1:n,2)     = ggh(1:mp,1:n)   * lT(1:mp,1:n)
    !
    ! Call Generic Thomas algorithm
    call generic_thomas(mp,n,A,B,C,d,x,iierr)
    if (iierr /= 0) then
       if (present(err)) then
          err = 1
          return
       else
#ifdef __MPI__
          call MPI_Abort(0, 159, ierr) ! Do not know comm nor rank here
#else
          stop 159
#endif
       endif
    endif
    dy(1:mp,1:n) = x(1:mp,1:n,1) * cS(1:mp,1:n)
    dT(1:mp,1:n) = x(1:mp,1:n,2) * cT(1:mp,1:n)
    !
  END SUBROUTINE massman_sparse_2d

  !**********************************************************************************************************************

  ELEMENTAL PURE SUBROUTINE litter_props(S, Tsoil, vlit, plit, h0)

    IMPLICIT NONE

    REAL(r_2),    INTENT(IN)    :: S
    REAL(r_2),    INTENT(IN)    :: Tsoil
    TYPE(vars),   INTENT(INOUT) :: vlit
    TYPE(params), INTENT(IN)    :: plit
    REAL(r_2),    INTENT(IN)    :: h0
    ! Get soil water variables from S.
    ! Definitions of arguments:
    ! S(1:ms)    - degree of saturation ("effective satn") of layers.
    ! ms        - no. of soil layers.
    ! jt(1:ms)    - layer soil type nos.
    ! var(1:ms) - other water vars of layers.
    REAL(r_2), PARAMETER :: u = 1.0 ! wind speed at litter surface
    REAL(r_2) :: theta, c
    REAL(r_2) :: dhdS, sl
    REAL(r_2) :: rhoL, copo, c_w ! params for thermal vapour transfer enhancement factor (Campbell 1985)
    REAL(r_2) :: chi, DT0

    sl      = slope_esat(Tsoil)*Mw/thousand/Rgas/(Tsoil+Tzero) ! m3 m-3 K-1
    vlit%sl = sl
    theta   = S*(plit%thre) + (plit%the - plit%thre)
    c       = Mw*gravity/Rgas/(Tsoil+Tzero)
    rhoL    = plit%rho
    ! dhdS    = -plit%he/plit%lam*S**(-one/plit%lam-one)*(thousand/rhoL*plit%the)**(-one/plit%lam)

    dhdS    = -plit%he/plit%lam*exp((-one/plit%lam-one)*log(S))*exp(-one/plit%lam*log(thousand/rhoL*plit%the))
    !    vlit%hS = dhdS
    ! Mathews (2006), A process-based model of offine fuel moisture,
    !                 International Journal of Wildland Fire 15,155-168
    chi     = 2.08_r_2+u*2.38_r_2 ! (Eq. 45, Tab. 1)
    DT0     = Dva*exp(u*2.6_r_2) ! (Eq. 46, Tab. 1)
    vlit%Dv = DT0*exp(-half*chi) ! heat and vapour diffusivity half-way through layer (Eq. 11)
    !write(*,*) 'litter_props', vlit%Dv
    if (S < one) then
       ! vlit%h = plit%he*(thousand/rhoL*S*plit%the)**(-one/plit%lam)
       vlit%h = plit%he*exp(-one/plit%lam*log(thousand/rhoL*S*plit%the))
       vlit%K = zero
    else
       vlit%K = vlit%Ke
       vlit%h = plit%he
    endif
    vlit%rh     = max(exp(Mw*gravity*vlit%h/Rgas/(Tsoil+Tzero)),rhmin)
    vlit%cvsat  = esat(Tsoil)*Mw/thousand/Rgas/(Tsoil+Tzero) ! m3 m-3
    vlit%cv     = vlit%cvsat*vlit%rh
    vlit%cvS    = vlit%rhS *vlit%cvsat
    vlit%KS     = zero
    vlit%phiv   = vlit%Dv*vlit%cvsat*vlit%rh
    vlit%phivS  = dhdS*vlit%phiv*c
    vlit%phi    = vlit%phiv
    vlit%phiS   = vlit%phivS
    vlit%rhS    = dhdS*c*vlit%rh
    ! sensible heat conductivity of wet soil W m-1 K-1 (Matthews 2006)
    vlit%kH     = 0.2_r_2 + 0.14_r_2*theta*thousand/rhoL
    copo        = 1932._r_2
    c_w         = 4.185e6_r_2
    vlit%csoil  = copo*rhoL + c_w*theta + c_w*h0 !volumetric heat capacity of wet soil
    vlit%eta_th = one !enhancement factor for transport of water vapour due to a temperature gradient
    vlit%kE     = vlit%Dv*vlit%rh*sl*thousand*rlambda
    vlit%kth    = vlit%kE + vlit%kH ! thermal conductivity of soil

  END SUBROUTINE litter_props

  !**********************************************************************************************************************

 ! ELEMENTAL PURE
SUBROUTINE potential_evap(Rn, rbh, rbw, Ta, rha, Tsoil, k, dz,lambdav, &
       Ts, E, H, G, &
       dEdrha, dEdTs, dEdTsoil, dGdTa, dGdTsoil,iice)

    ! Pennman-Monteith equation, with additional account for heat flux into the surface

    IMPLICIT NONE

    REAL(r_2), INTENT(IN)  :: Rn, rbh, rbw, Ta, rha, Tsoil, k, dz, lambdav
    REAL(r_2), INTENT(OUT) :: Ts, E, H, G, dEdrha, dEdTs, dEdTsoil, dGdTa, dGdTsoil
    LOGICAL, INTENT(IN), OPTIONAL :: iice
    REAL(r_2) :: s, es, ea, dEdea, dEdesat, dTsdTa, dEdDa, Da
    REAL(r_2):: rhocp, gamma != 67.0 ! psychrometric constant

    if (present(iice)) then
       if(iice) then
          es = esat_ice(Ta)
          s  = slope_esat_ice(Ta)
       else
          es = esat(Ta)
          s  = slope_esat(Ta)
       endif
    else
       es = esat(Ta)
       s  = slope_esat(Ta)
    endif

    rhocp = rmair*101325._r_2/rgas/(Ta+Tzero)*cpa
    gamma = 101325._r_2*cpa/lambdav/(Mw/rmair)
    ea = es * max(rha, 0.1_r_2)
    Da = ea/max(rha, 0.1_r_2) - ea

    E  = (rhocp*(Da*(k*rbh + dz*rhocp) + rbh*s*(dz*Rn + k*(-Ta + Tsoil)))) / &
         (gamma*rbw*(k*rbh + dz*rhocp) + dz*rbh*rhocp*s)
    Ts = Ta + E*gamma*rbw/s/rhocp - Da/s
    H  = rhocp*(Ts - Ta)/rbh
    G  = k*(Ts-Tsoil)/dz

    dEdDa    = (-(k*rbh*rhocp) - dz*rhocp**2)/(gamma*k*rbh*rbw + dz*gamma*rbw*rhocp + dz*rbh*rhocp*s)
    dEdea    = -dEdDa
    dEdesat  = dEdea
    dEdrha   = dEdea *es
    !dEdTa    = (k*rbh*rhocp*s)/(gamma*k*rbh*rbw + dz*gamma*rbw*rhocp + dz*rbh*rhocp*s) + dEdesat *s
    dEdTsoil = -((k*rbh*rhocp*s)/(gamma*k*rbh*rbw + dz*gamma*rbw*rhocp + dz*rbh*rhocp*s))

    dTsdTa   = (-(dz*gamma*rbw*rhocp) - dz*rbh*rhocp*s)/(gamma*k*rbh*rbw + dz*gamma*rbw*rhocp + dz*rbh*rhocp*s)

    dGdTa    = k/dz * dTsdTa
    dGdTsoil = -k/dz !+k/dz*dEdTsoil*gamma*rbw/s/rhocp
    dEdTs = s*rhocp/gamma/rbw

  END SUBROUTINE potential_evap

  !**********************************************************************************************************************

  SUBROUTINE setlitterpar(mp, veg,index)

    IMPLICIT NONE

    INTEGER(i_d),              INTENT(IN) :: mp
    TYPE(veg_parameter_type), INTENT(IN) :: veg
    integer(i_d), DIMENSION(:),  INTENT(IN) :: index

    allocate(plit(mp))
    allocate(dxL(mp))

    plit%the   = 0.09_r_2
    plit%thre  = 0.09_r_2
    plit%he    = -35.0_r_2
    plit%lam   = one/2.4_r_2
    plit%Ke    = zero
    plit%eta   = zero
    plit%KSe   = zero
    plit%phie  = zero
    plit%phiSe = zero
    plit%rho   = 63.5_r_2
    ! dxL        = zero            ! litter params
    dxL        = real(veg%clitt(index),r_2)*two/plit%rho*0.1_r_2

    plit%ishorizon  = 0
    plit%thw        = zero
    plit%thfc       = zero
    plit%thr        = zero

    plit%eta        = zero
    plit%KSe        = zero
    plit%phie       = zero
    plit%phiSe      = zero
    plit%tortuosity = zero
    plit%clay       = zero
    plit%zeta       = zero
    plit%fsatmax    = zero
    plit%css        = zero
    plit%kd         = zero
    plit%lambc      = zero
    plit%LambdaS    = zero

  END SUBROUTINE setlitterpar

  !**********************************************************************************************************************

  SUBROUTINE setpar(mp, ms, soil, index)

    IMPLICIT NONE

    INTEGER(i_d),                 INTENT(IN) :: mp
    INTEGER(i_d),                 INTENT(IN) :: ms
    TYPE(soil_parameter_type),    INTENT(IN) :: soil
    integer(i_d), DIMENSION(:),   INTENT(IN) :: index

    INTEGER(i_d) :: i

    allocate(par(mp,ms))

    do i=1, ms
       par(:,i)%ishorizon  = soil%ishorizon(index,i)
       par(:,i)%thw        = real(soil%swilt(index),r_2)
       par(:,i)%thfc       = real(soil%sfc(index),r_2)
       par(:,i)%the        = real(soil%ssat(index),r_2)
       par(:,i)%thr        = zero
       par(:,i)%thre       = real(soil%ssat(index),r_2) - par(:,i)%thr
       par(:,i)%he         = real(soil%sucs(index),r_2)
       par(:,i)%Ke         = real(soil%hyds(index),r_2)
       par(:,i)%lam        = one/real(soil%bch(index),r_2)
       par(:,i)%eta        = two/par(:,i)%lam + two + one
       par(:,i)%KSe        = par(:,i)%eta * par(:,i)%Ke    ! dK/dS at he
       par(:,i)%phie       = par(:,i)%Ke * par(:,i)%he / (one - par(:,i)%lam * par(:,i)%eta) ! MFP at he
       par(:,i)%phiSe      = (par(:,i)%eta - one/par(:,i)%lam) * par(:,i)%phie    ! dphi/dS at he
       par(:,i)%kd         = real(soil%cnsd(index),r_2)
       par(:,i)%css        = real(soil%css(index),r_2)
       par(:,i)%rho        = real(soil%rhosoil(index),r_2)
       par(:,i)%tortuosity = 0.67_r_2
       par(:,i)%clay       = real(soil%clay(index),r_2)
       par(:,i)%zeta       = real(soil%zeta(index),r_2)
       par(:,i)%fsatmax    = real(soil%fsatmax(index),r_2)
       par(:,i)%lambc      = par(:,i)%lam
       !par(:,i)%LambdaS    = real(soil%LambdaS(index),r_2)
       par(:,i)%LambdaS    = 2830_r_2 ! Sispat Manual Table 2
    enddo

! special for Cumberland: set soil params to clay:
    if (cable_user%Cumberland_soil) then
       do i=1,ms
          par(:,i)%thw        = 0.286
          par(:,i)%thfc       = 0.367
          par(:,i)%the        = 0.482
          par(:,i)%thre       = real(soil%ssat(index),r_2) - par(:,i)%thr
          par(:,i)%he         = -0.405
          par(:,i)%Ke         = 1e-6
          par(:,i)%lam        = 1./11.4
          par(:,i)%eta        = two/par(:,i)%lam + two + one
          par(:,i)%KSe        = par(:,i)%eta * par(:,i)%Ke    ! dK/dS at he
          par(:,i)%phie       = par(:,i)%Ke * par(:,i)%he / (one - par(:,i)%lam * par(:,i)%eta) ! MFP at he
          par(:,i)%phiSe      = (par(:,i)%eta - one/par(:,i)%lam) * par(:,i)%phie    ! dphi/dS at he
          par(:,i)%clay       = 0.67

!!$          par(:,i)%thw        = 0.178
!!$          par(:,i)%thfc       = 0.367
!!$          par(:,i)%the        = 0.45
!!$          par(:,i)%thre       = real(soil%ssat(index),r_2) - par(:,i)%thr
!!$          par(:,i)%he         = -0.572
!!$          par(:,i)%Ke         = 2.8e-5
!!$          par(:,i)%lam        = 1./8.7
!!$          par(:,i)%eta        = two/par(:,i)%lam + two + one
!!$          par(:,i)%KSe        = par(:,i)%eta * par(:,i)%Ke    ! dK/dS at he
!!$          par(:,i)%phie       = par(:,i)%Ke * par(:,i)%he / (one - par(:,i)%lam * par(:,i)%eta) ! MFP at he
!!$          par(:,i)%phiSe      = (par(:,i)%eta - one/par(:,i)%lam) * par(:,i)%phie    ! dphi/dS at he
!!$          par(:,i)%clay       = 0.2

       enddo
    endif

!!$   ! special for Cumberland: set top 3 layers to sand
!!$   do i=1,3
!!$      par(:,i)%thw        = 0.175
!!$      par(:,i)%thfc       = 0.255
!!$      par(:,i)%the        = 0.420
!!$      par(:,i)%thre       = real(soil%ssat(index),r_2) - par(:,i)%thr
!!$      par(:,i)%he         = -0.299
!!$      par(:,i)%Ke         = 6.e-6
!!$      par(:,i)%lam        = 1./7.12
!!$      par(:,i)%eta        = two/par(:,i)%lam + two + one
!!$      par(:,i)%KSe        = par(:,i)%eta * par(:,i)%Ke    ! dK/dS at he
!!$      par(:,i)%phie       = par(:,i)%Ke * par(:,i)%he / (one - par(:,i)%lam * par(:,i)%eta) ! MFP at he
!!$      par(:,i)%phiSe      = (par(:,i)%eta - one/par(:,i)%lam) * par(:,i)%phie    ! dphi/dS at he
!!$      par(:,i)%clay       = 0.27
!!$
!!$   enddo


  END SUBROUTINE setpar

  !**********************************************************************************************************************

  SUBROUTINE setpar_Loetsch(mp, ms, x2dx)

    IMPLICIT NONE

    INTEGER(i_d),                 INTENT(IN)    :: mp
    INTEGER(i_d),                 INTENT(IN)    :: ms
    REAL(r_2),    DIMENSION(:,:), INTENT(IN)    :: x2dx

    INTEGER(i_d) :: i

    allocate(par(mp,ms))

    do i=1, ms
       where (x2dx(:,i) .lt. 0.15_r_2)
          par(:,i)%ishorizon  = 1
          par(:,i)%thw        = 0.07_r_2
          par(:,i)%thfc       = 0.45_r_2
          par(:,i)%the        = 0.60_r_2
          par(:,i)%thr        = 0.00_r_2
          par(:,i)%he         = -6.76_r_2
          par(:,i)%lam        = 0.70_r_2
          par(:,i)%Ke         = 2.754e-5_r_2
       elsewhere (x2dx(:,i) .gt. 2.0_r_2)
          par(:,i)%ishorizon  = 3
          par(:,i)%thw        = 0.18_r_2
          par(:,i)%thfc       = 0.35_r_2
          par(:,i)%the        = 0.57_r_2
          par(:,i)%thr        = 0.00_r_2
          par(:,i)%he         = -1.23_r_2    ! optimised, was -1.46
          par(:,i)%lam        = 0.67_r_2     ! optimised, was 0.25
          par(:,i)%Ke         = 1.013e-5_r_2 ! optimised, was 2.481e-7
       elsewhere
          par(:,i)%ishorizon  = 2
          par(:,i)%thw        = 0.18_r_2
          par(:,i)%thfc       = 0.35_r_2
          par(:,i)%the        = 0.57_r_2
          par(:,i)%thr        = 0.00_r_2
          par(:,i)%he         = -1.46_r_2
          par(:,i)%lam        = 0.25_r_2
          par(:,i)%Ke         = 2.481e-5_r_2
       end where
    enddo
    par%clay       = 0.06_r_2
    par%tortuosity = 0.67_r_2
    par%thre       = par%the - par%thr
    par%eta        = two/par%lam + two + one
    par%KSe        = par%eta * par%Ke    ! dK/dS at he
    par%phie       = par%Ke * par%he / (one - par%lam * par%eta) ! MFP at he
    par%phiSe      = (par%eta - one/par%lam) * par%phie    ! dphi/dS at he
    par%kd         = zero ! not used
    par%rho        = 2650._r_2 * par%the            ! porosity = soil/stone
    par%css        = 2400._r_2/2650._r_2  * par%rho ! Campbell (1985)
    par%zeta       = 2.0_r_2   ! topmodel exponent baseflow
    par%fsatmax    = 0.30_r_2  ! topmodel maximum saturated fraction
    par%lambc      = par%lam
    par%LambdaS    = 2539._r_2 ! optimised, was 2830 ! Sispat Manual Table 2

    ! ! Cable coarse sand
    ! par%ishorizon  = 1
    ! par%thw        = 0.072_r_2
    ! par%thfc       = 0.143_r_2
    ! par%the        = 0.398_r_2
    ! par%thr        = zero
    ! par%he         = -0.106_r_2
    ! par%Ke         = 166.e-6_r_2
    ! par%lam        = 1.0_r_2/4.2_r_2
    ! par%clay       = 0.09_r_2
    ! par%rho        = 1600._r_2
    ! par%css        = 850._r_2
    ! !
    ! par%tortuosity = 0.67_r_2
    ! par%thre       = par%the - par%thr
    ! par%eta        = two/par%lam + two + one
    ! par%KSe        = par%eta * par%Ke    ! dK/dS at he
    ! par%phie       = par%Ke * par%he / (one - par%lam * par%eta) ! MFP at he
    ! par%phiSe      = (par%eta - one/par%lam) * par%phie    ! dphi/dS at he
    ! par%kd         = zero ! not used
    ! par%zeta       = 2.0_r_2
    ! par%fsatmax    = 0.3_r_2
    ! par%lambc      = par%lam
    ! par%LambdaS    = 2800._r_2 ! Sispat Manual Table 2

  END SUBROUTINE setpar_Loetsch

  !**********************************************************************************************************************

  SUBROUTINE setsol(mp)

    IMPLICIT NONE

    INTEGER(i_d), INTENT(IN) :: mp

    allocate(sol(mp))

    sol%T1      = zero
    sol%Ta      = zero
    sol%cva     = zero
    sol%Rnet    = zero
    sol%hr1     = zero
    sol%hra     = zero
    sol%Dv      = zero
    sol%gv      = zero
    sol%gh      = zero
    sol%Dh      = zero
    sol%dz      = zero
    sol%phie    = zero
    sol%he      = zero
    sol%K1      = zero
    sol%eta     = zero
    sol%lambda  = zero
    sol%Ks      = zero
    sol%lambdav = zero

  END SUBROUTINE setsol

  !**********************************************************************************************************************

  SUBROUTINE setx(mp, ms, soil)

    IMPLICIT NONE

    INTEGER(i_d),              INTENT(IN) :: mp
    INTEGER(i_d),              INTENT(IN) :: ms
    TYPE(soil_parameter_type), INTENT(IN) :: soil

    REAL(r_2), DIMENSION(ms) :: tmp1d
    INTEGER :: i

    allocate(x(mp,ms))
    allocate(dx(mp,ms))

    ! cumulative soil layer depths = bottom of soil layers
    tmp1d(1) = real(soil%zse(1),r_2)
    do i=2, ms
       tmp1d(i) = tmp1d(i-1) + real(soil%zse(i),r_2)
    end do

    ! soil layer depth
    dx(1:mp,1:ms) = spread(real(soil%zse(1:ms),r_2),1,mp)
    ! bottom of each soil layer
    x(1:mp,1:ms)  = spread(tmp1d(1:ms),1,mp)

  END SUBROUTINE setx

  !**********************************************************************************************************************

  SUBROUTINE tri_1d(ns, n, aa, bb, cc, dd, dy)

    IMPLICIT NONE

    INTEGER(i_d),                 INTENT(IN)    :: ns
    INTEGER(i_d),                 INTENT(IN)    :: n
    REAL(r_2),    DIMENSION(0:n), INTENT(IN)    :: aa
    REAL(r_2),    DIMENSION(0:n), INTENT(IN)    :: cc
    REAL(r_2),    DIMENSION(0:n), INTENT(IN)    :: dd
    REAL(r_2),    DIMENSION(0:n), INTENT(INOUT) :: bb
    REAL(r_2),    DIMENSION(0:n), INTENT(INOUT) :: dy
    ! Solves tridiag set of linear eqns. Coeff arrays aa and cc left intact.
    ! Definitions of arguments:
    ! ns      - start index for eqns.
    ! n       - end index.
    ! aa(0:n) - coeffs below diagonal; ns+1:n used.
    ! bb(0:n) - coeffs on diagonal; ns:n used.
    ! cc(0:n) - coeffs above diagonal; ns:n-1 used.
    ! dd(0:n) - rhs coeffs; ns:n used.
    ! ee(0:n) - work space.
    ! dy(0:n) - solution in ns:n.
    REAL(r_2),   DIMENSION(0:n) :: ee
    INTEGER(i_d)                :: i

    dy(ns) = dd(ns) ! decomposition and forward substitution
    do i=ns, n-1
       ee(i)   = cc(i)/bb(i)
       dy(i)   = dy(i)/bb(i)
       bb(i+1) = bb(i+1)-aa(i+1)*ee(i)
       dy(i+1) = dd(i+1)-aa(i+1)*dy(i)
    end do

    dy(n) = dy(n)/bb(n) ! back substitution
    do i=n-1, ns, -1
       dy(i) = dy(i)-ee(i)*dy(i+1)
    end do

  END SUBROUTINE tri_1d

  SUBROUTINE tri_2d(mp, ns, n, aa, bb, cc, dd, dy)

    IMPLICIT NONE

    INTEGER(i_d),                      INTENT(IN)    :: mp
    INTEGER(i_d),                      INTENT(IN)    :: ns
    INTEGER(i_d),                      INTENT(IN)    :: n
    REAL(r_2),    DIMENSION(1:mp,0:n), INTENT(IN)    :: aa
    REAL(r_2),    DIMENSION(1:mp,0:n), INTENT(IN)    :: cc
    REAL(r_2),    DIMENSION(1:mp,0:n), INTENT(IN)    :: dd
    REAL(r_2),    DIMENSION(1:mp,0:n), INTENT(INOUT) :: bb
    REAL(r_2),    DIMENSION(1:mp,0:n), INTENT(INOUT) :: dy
    ! Solves tridiag set of linear eqns. Coeff arrays aa and cc left intact.
    ! Definitions of arguments:
    ! mp      - number of land patches
    ! ns      - start index for eqns.
    ! n       - end index.
    ! aa(0:n) - coeffs below diagonal; ns+1:n used.
    ! bb(0:n) - coeffs on diagonal; ns:n used.
    ! cc(0:n) - coeffs above diagonal; ns:n-1 used.
    ! dd(0:n) - rhs coeffs; ns:n used.
    ! ee(0:n) - work space.
    ! dy(0:n) - solution in ns:n.
    REAL(r_2),   DIMENSION(1:mp,0:n) :: ee
    INTEGER(i_d)                     :: i

    dy(:,ns) = dd(:,ns) ! decomposition and forward substitution
    do i=ns, n-1
       ee(:,i)   = cc(:,i)/bb(:,i)
       dy(:,i)   = dy(:,i)/bb(:,i)
       bb(:,i+1) = bb(:,i+1)-aa(:,i+1)*ee(:,i)
       dy(:,i+1) = dd(:,i+1)-aa(:,i+1)*dy(:,i)
    end do

    dy(:,n) = dy(:,n)/bb(:,n) ! back substitution
    do i=n-1, ns, -1
       dy(:,i) = dy(:,i)-ee(:,i)*dy(:,i+1)
    end do

  END SUBROUTINE tri_2d

  !**********************************************************************************************************************
  ! FUNCTIONS - alphabetical
  !**********************************************************************************************************************

  REAL(r_2) ELEMENTAL PURE FUNCTION csat(T)
    !returns  sat vapour pressure curve in kg m-3
    USE sli_numbers, ONLY: Tzero, Rgas, Mw, esata, esatb, esatc

    IMPLICIT NONE

    real(r_2), intent(in) :: T

    csat = esata * exp(esatb*T/(T+esatc)) * Mw/Rgas/(T+Tzero)

  END FUNCTION csat

  !**********************************************************************************************************************

  REAL(r_2) ELEMENTAL PURE FUNCTION csoil(thetal,thetai,css,rho)
    ! determines heat capacity of soil (Jkg-1K-1)
    USE sli_numbers, ONLY: rhow, cswat, csice

    IMPLICIT NONE

    real(r_2), intent(in) ::  thetal, thetai, css, rho

    csoil = css*rho + rhow*cswat*thetal + rhow*csice*thetai

  END FUNCTION csoil

  !**********************************************************************************************************************

  REAL(r_2) ELEMENTAL PURE FUNCTION dthetalmaxdT(Tin,S,he,b,thre,the)
    ! determines derivative of thetalmax wrt T
    USE sli_numbers, ONLY: gravity, lambdaf, csol, Rgas, Tzero

    IMPLICIT NONE

    real(r_2), intent(in) :: Tin,S,he,b,thre,the
    real(r_2)             :: dthetaldh, h, psi, PI, T, thetalmax

    T = min(Tfrz(min(S,one),he,b),Tin)
    PI  = -csol*Rgas*(T+Tzero)/gravity ! osmotic potential (m)
    psi = lambdaf*T/(gravity*(T+Tzero))! total matric potential in presence of ice
    h   = psi-PI       ! moisture potential in presence of ice

    ! thetalmax = thre*(h/he)**(-one/b) + (the-thre)
    thetalmax = thre*exp(-one/b*log(h/he)) + (the-thre)
    dthetaldh = -thetalmax/b/h

    dthetalmaxdT = dthetaldh*(lambdaf/(T+Tzero)/gravity-lambdaf*T/gravity/(T+Tzero)**2+csol*Rgas/gravity)
    if (S>1) dthetalmaxdT = dthetalmaxdT * S

  END FUNCTION dthetalmaxdT

  !**********************************************************************************************************************

  REAL(r_2) ELEMENTAL PURE FUNCTION dthetalmaxdTh(Tin,S,he,b,thre,the)
    ! determines derivative of thetalmax wrt T
    ! (uses dthetaldh defined in terms of S ...required for crossing freezing point)
    USE sli_numbers, ONLY: gravity, lambdaf, csol, Rgas, Tzero

    IMPLICIT NONE

    real(r_2), INTENT(in) :: Tin,S,he,b,thre,the
    real(r_2)             :: dthetaldh, h, theta, T

    T     = min(Tfrz(S,he,b),Tin)
    ! h     = he * S**(-b)
    h     = he * exp(-b*log(S))
    theta = S*thre + (the-thre)

    dthetaldh = -theta/b/h

    dthetalmaxdTh = dthetaldh*(lambdaf/(T+Tzero)/gravity-lambdaf*T/gravity/(T+Tzero)**2+csol*Rgas/gravity)

  END FUNCTION dthetalmaxdTh

  !*********************************************************************************************************************

  REAL(r_2) ELEMENTAL PURE FUNCTION esat(T)
    !returns sat vapour pressure curve in Pa
    USE sli_numbers, ONLY: esata, esatb, esatc

    IMPLICIT NONE

    real(r_2), intent(in) :: T

    esat = esata * exp(esatb*T/(T+esatc))

  END FUNCTION esat

  !**********************************************************************************************************************

  REAL(r_2) ELEMENTAL PURE FUNCTION gammln(z)
    !  Uses Lanczos-type approximation to ln(gamma) for z > 0.
    !  Reference:
    !       Lanczos, C. 'A precision approximation of the gamma
    !               function', J. SIAM Numer. Anal., B, 1, 86-96, 1964.
    !  Accuracy: About 14 significant digits except for small regions
    !            in the vicinity of 1 and 2.
    !  Programmer: Alan Miller
    !              1 Creswick Street, Brighton, Vic. 3187, Australia
    !  e-mail: amiller @ bigpond.net.au
    !  Latest revision - 14 October 1996
    IMPLICIT NONE

    REAL(r_2), INTENT(IN) :: z
    ! Local variables
    REAL(r_2), PARAMETER :: a(9) = (/ &
         0.9999999999995183_r_2, 676.5203681218835_r_2, -1259.139216722289_r_2, &
         771.3234287757674_r_2, -176.6150291498386_r_2, 12.50734324009056_r_2, &
         -0.1385710331296526_r_2, 0.9934937113930748E-05_r_2, 0.1659470187408462E-06_r_2 /)
    REAL(r_2), PARAMETER :: zero = 0.0_r_2
    REAL(r_2), PARAMETER :: one = 1.0_r_2
    REAL(r_2), PARAMETER :: lnsqrt2pi =  0.9189385332046727_r_2
    REAL(r_2), PARAMETER :: half = 0.5_r_2
    REAL(r_2), PARAMETER :: sixpt5 = 6.5_r_2
    REAL(r_2), PARAMETER :: seven = 7.0_r_2
    REAL(r_2)    :: tmp, tmpgammln
    INTEGER(i_d) :: j

    tmpgammln = zero
    tmp = z + seven
    DO j = 9, 2, -1
       tmpgammln = tmpgammln + a(j)/tmp
       tmp = tmp - one
    END DO
    tmpgammln = tmpgammln + a(1)
    gammln = LOG(tmpgammln) + lnsqrt2pi - (z + sixpt5) + (z - half)*LOG(z + sixpt5)
    RETURN

  END FUNCTION gammln

  !**********************************************************************************************************************

  REAL(r_2) ELEMENTAL PURE FUNCTION gcf(a,x)

    IMPLICIT NONE

    REAL(r_2), INTENT(IN) :: a,x
    INTEGER(i_d), PARAMETER :: ITMAX=100
    REAL(r_2),    PARAMETER :: EPS=epsilon(x)
    REAL(r_2),    PARAMETER :: FPMIN=tiny(x)/EPS
    INTEGER(i_d) :: i
    REAL(r_2)     :: an, b, c, d, del, h

    if (eq(x, 0.0_r_2)) then
       gcf=1.0_r_2
       RETURN
    end if
    b = x + 1.0_r_2 - a
    c = 1.0_r_2/FPMIN
    d = 1.0_r_2/b
    h = d
    do i=1, ITMAX
       an  = -i*(i-a)
       b   = b + 2.0_r_2
       d   = an*d + b
       if (abs(d)<FPMIN) d=FPMIN
       c   = b + an/c
       if (abs(c)<FPMIN) c=FPMIN
       d   = 1.0_r_2/d
       del = d*c
       h   = h*del
       if (abs(del-1.0_r_2) <= EPS) exit
    end do
    gcf = exp(-x + a*log(x) - gammln(a)) * h

  END FUNCTION gcf

  !**********************************************************************************************************************

  REAL(r_2) ELEMENTAL PURE FUNCTION gser(a,x)

    IMPLICIT NONE

    REAL(r_2), INTENT(IN) :: a, x
    INTEGER(i_d), PARAMETER :: ITMAX=100
    REAL(r_2),     PARAMETER :: EPS=epsilon(x)
    INTEGER(i_d) :: n
    REAL(r_2)     :: ap, del, summ

    if (eq(x, 0.0_r_2)) then
       gser = 0.0_r_2
       RETURN
    end if
    ap   = a
    summ = 1.0_r_2/a
    del  = summ
    do n=1, ITMAX
       ap   = ap + 1.0_r_2
       del  = del*x/ap
       summ = summ + del
       if (abs(del) < abs(summ)*EPS) exit
    end do
    gser=summ*exp(-x+a*log(x)-gammln(a))

  END FUNCTION gser

  !**********************************************************************************************************************

  REAL(r_2) ELEMENTAL PURE FUNCTION igamma(a,x)
    USE cable_def_types_mod, ONLY: r_2
    IMPLICIT NONE
    REAL(r_2), INTENT(IN) :: a, x
    REAL(r_2) :: gln
    gln = gammln(a)
    if (x < a+1.0_r_2) then
       igamma = 1.0_r_2 - gser(a,x)
    else
       igamma = gcf(a,x)
    end if
    igamma = igamma * exp(gln)

  END FUNCTION igamma

  !**********************************************************************************************************************

  ! MC this routines has to be adjusted for csol in freezing, probably.
  REAL(r_2) ELEMENTAL PURE FUNCTION phi(hr0, lambda, eta, phie, he, T, Ksat)

    USE sli_numbers, ONLY: zero, one, gravity, lambdaf, Rgas, Tzero ! , csol

    IMPLICIT NONE

    REAL(r_2), INTENT(IN)           :: hr0, lambda, eta, phie, he, T
    REAL(r_2), INTENT(IN), OPTIONAL :: Ksat

    REAL(r_2) :: h, csol1

    ! MC freezing point? csol?
    csol1 = zero ! use zero instead of csol for the moment
    if (present(Ksat)) then
       if (T < zero) then     ! frozen soil  !! need to adjust freezing point for csol and use global csol
          h   = (lambdaf*T/gravity/(T+Tzero)) + csol1*Rgas* (T + Tzero)/gravity
          ! phi = Ksat * he/(one-eta*lambda) * (h/he)**(one-eta*lambda)
          phi = Ksat * he/(one-eta*lambda) * exp((one-eta*lambda)*log(h/he))
       else
          phi = zero
       endif
    else
       if (hr0 < one) then
          ! phi = phie * (Rgas *(T+Tzero) *log(hr0) /(gravity * he * Mw))**(one-eta*lambda)
          phi = phie * exp((one-eta*lambda)*log(Rgas *(T+Tzero) *log(hr0) /(gravity * he * Mw)))
       else
          phi = phie
       endif
    endif

  END FUNCTION phi

  !**********************************************************************************************************************

  REAL(r_2) ELEMENTAL PURE FUNCTION rh0_sol(hr0, solin)

    IMPLICIT NONE

    REAL(r_2),        INTENT(in) :: hr0
    TYPE(solve_type), INTENT(in) :: solin

    !    Rnet = Rnet + rhow*lambdaf*var(1)%iice*((phi(hr1,lambda,eta,phie,he,T1) &
    !         - phi(hr0,lambda,eta,phie,he,T1))/dz - K1)
    rh0_sol = (solin%Dv*(-hr0 + solin%hr1)*csat(solin%T1))/solin%dz - &
         (thousand*(solin%dz*solin%K1 + phi(hr0,solin%lambda,solin%eta,solin%phie,solin%he,solin%T1) &
         - phi(solin%hr1,solin%lambda,solin%eta,solin%phie,solin%he,solin%T1)))/solin%dz - &
         (solin%Dv*solin%hr1*(solin%cva*solin%gv*solin%lambdav + solin%Rnet + solin%gh*(-solin%T1 + solin%Ta) &
         - solin%gv*hr0*solin%lambdav*csat(solin%T1))*slope_csat(solin%T1)) / &
         (solin%DH + solin%dz*solin%gh + solin%dz*solin%gv*hr0*solin%lambdav*slope_csat(solin%T1)) + &
         ((solin%DH + solin%dz*solin%gh)*solin%gv*(solin%cva - hr0*csat(solin%T1)) &
         - solin%dz*solin%gv*hr0*(solin%Rnet + solin%gh*(-solin%T1 + solin%Ta))*slope_csat(solin%T1)) / &
         (solin%DH + solin%dz*solin%gh + solin%dz*solin%gv*hr0*solin%lambdav*slope_csat(solin%T1))

  END FUNCTION rh0_sol

  !**********************************************************************************************************************

  ! Using an elemental subroutine as a function argument is not Fortran90 standard.
  ! This would have been rh0_sol in function rtbis. intel's ifort and nag's f95 accept it but gfortran does not.
  ! Write rh0_sol into rtbis -> rtbis_rh0
  ! It does not check that f(x1) and f(x2) have different signs and not if iteration > MAXIT

  REAL(r_2) ELEMENTAL PURE FUNCTION rtbis_rh0(sol, x1, x2, xacc)

    IMPLICIT NONE

    TYPE(solve_type), INTENT(IN) :: sol
    REAL(r_2),        INTENT(IN) :: x1, x2, xacc

    INTEGER(i_d), PARAMETER :: MAXIT=40
    INTEGER(i_d) :: j
    REAL(r_2)    :: dx, f, fmid, xmid

    fmid = rh0_sol(x2,sol)
    f    = rh0_sol(x1,sol)
    if (f < zero) then
       rtbis_rh0 = x1
       dx        = x2-x1
    else
       rtbis_rh0 = x2
       dx        = x1-x2
    end if
    do j=1, MAXIT
       dx   = dx*half
       xmid = rtbis_rh0+dx
       fmid = rh0_sol(xmid,sol)
       if (fmid <= zero) rtbis_rh0 = xmid
       if (abs(dx) < xacc .or. eq(fmid, zero)) RETURN
    end do

  END FUNCTION rtbis_rh0

  !**********************************************************************************************************************

  REAL(r_2) ELEMENTAL PURE FUNCTION slope_csat(T)
    !returns slope of sat vapour pressure curve in kg m-3 K-1
    USE sli_numbers, ONLY: Tzero, Rgas, Mw, esata, esatb, esatc

    IMPLICIT NONE

    real(r_2), intent(in) :: T
    real(r_2) :: csat

    csat       = esata * exp(esatb*T/(T+esatc)) * Mw/Rgas/(T+Tzero)
    slope_csat = csat * esatb*esatc/(T+esatc)**2

  END FUNCTION slope_csat

  !**********************************************************************************************************************

  REAL(r_2) ELEMENTAL PURE FUNCTION slope_esat(T)
    !returns slope of sat vapour pressure curve in Pa K^-1
    USE sli_numbers, ONLY: esata, esatb, esatc

    IMPLICIT NONE

    real(r_2), intent(in) :: T
    real(r_2) :: esat

    esat       = esata * exp(esatb*T/(T+esatc))
    slope_esat = esat * esatb*esatc/(T+esatc)**2

  END FUNCTION slope_esat

  !*********************************************************************************************************************

  REAL(r_2) ELEMENTAL FUNCTION esat_ice(T)
    !returns sat vapour pressure curve in Pa
    USE sli_numbers, ONLY:  esata_ice, esatb_ice, esatc_ice

    IMPLICIT NONE

    real(r_2), intent(in) :: T

    esat_ice = esata_ice * exp(esatb_ice*T/(T+esatc_ice))

  END FUNCTION esat_ice

  !**********************************************************************************************************************

  REAL(r_2) ELEMENTAL FUNCTION slope_esat_ice(T)
    !returns slope of sat vapour pressure curve in Pa K^-1
    USE sli_numbers, ONLY:  esata_ice, esatb_ice, esatc_ice

    IMPLICIT NONE

    real(r_2), intent(in) :: T
    real(r_2) :: esat_ice

    esat_ice       = esata_ice * exp(esatb_ice*T/(T+esatc_ice))
    slope_esat_ice = esat_ice * esatb_ice*esatc_ice/(T+esatc_ice)**2

  END FUNCTION slope_esat_ice

  !**********************************************************************************************************************

  REAL(r_2) ELEMENTAL PURE FUNCTION Sofh(h,parin)

    IMPLICIT NONE

    REAL(r_2),    INTENT(IN) :: h
    TYPE(params), INTENT(IN) :: parin

    ! Get saturation S from matric head h.
    ! Definitions of arguments:
    ! h   - matric head.

    ! Sofh = (h/parin%he)**(-parin%lam) ! Sofh not used much so ** not an issue
    Sofh = exp(-parin%lam*log(h/parin%he))

  END FUNCTION Sofh

  !**********************************************************************************************************************

  REAL(r_2) ELEMENTAL PURE FUNCTION Tfrz(S,he,b)
    ! determines freezing point temperature for a given moisture content
    USE sli_numbers, ONLY: one, two, four, gravity, lambdaf, csol, Rgas, Tzero

    IMPLICIT NONE

    real(r_2), intent(in) :: S, he, b

    if (csol > e3) then
       ! Tfrz = (gravity*he - (min(S,one))**b * (lambdaf + two*csol*Rgas*Tzero) + &
       !      sqrt(gravity**2 * he**2 - two*gravity*he*lambdaf*(min(S,one))**b + &
       !      lambdaf*(min(S,one))**(two*b)*(lambdaf + four*csol*Rgas*Tzero)))/ &
       !      (two*csol*Rgas*(min(S,one))**b)
       Tfrz = (gravity*he - exp(b*log(min(S,one))) * (lambdaf + two*csol*Rgas*Tzero) + &
            sqrt(gravity**2 * he**2 - two*gravity*he*lambdaf*exp(b*log(min(S,one))) + &
            lambdaf*exp(two*b*log(min(S,one)))*(lambdaf + four*csol*Rgas*Tzero))) / &
            (two*csol*Rgas*exp(b*log(min(S,one))))
    else
       ! Tfrz = (gravity*he*Tzero) / (-(gravity*he) + lambdaf*(S)**b)
       Tfrz = (gravity*he*Tzero) / (-gravity*he + lambdaf*exp(b*log(max(min(S,one),0.01_r_2))))
    endif

  END FUNCTION Tfrz

  !**********************************************************************************************************************

  REAL(r_2) ELEMENTAL PURE FUNCTION Tfrozen(J, dx, theta, thetal, csoil, rhosoil, h0, thetasat)
    ! determines temperature of frozen soil, given total energy content J and liquid water content thetal
    USE sli_numbers, ONLY: csice, cswat, rhow, lambdaf

    IMPLICIT NONE

    real(r_2), intent(in) :: J, dx, theta, thetal, csoil, rhosoil, h0, thetasat

    if (h0 > zero) then
       Tfrozen = (cswat*h0*rhow*theta + h0*lambdaf*rhow*theta - &
            cswat*h0*rhow*thetal - h0*lambdaf*rhow*thetal + &
            J*thetasat - cswat*h0*rhow*thetasat + &
            dx*lambdaf*rhow*theta*thetasat - &
            dx*lambdaf*rhow*thetal*thetasat)/ &
            (csice*h0*rhow*theta - csice*h0*rhow*thetal + &
            csoil*dx*rhosoil*thetasat + csice*dx*rhow*theta*thetasat - &
            csice*dx*rhow*thetal*thetasat +  &
            cswat*dx*rhow*thetal*thetasat)
    else
       Tfrozen = (J + dx*lambdaf*rhow*theta - dx*lambdaf*rhow*thetal)/ &
            (dx*(csoil*rhosoil + csice*rhow*theta - csice*rhow*thetal + cswat*rhow*thetal))
    endif

  END FUNCTION Tfrozen

  !**********************************************************************************************************************

  REAL(r_2) ELEMENTAL PURE FUNCTION GTfrozen(T, J, dx, theta, csoil, rhosoil, h0, thre, the, he, b)
    ! GTfrozen = sensible heat + latent heat - total energy (should be zero)
    USE sli_numbers, ONLY: csice, cswat, rhow, lambdaf

    IMPLICIT NONE

    real(r_2), intent(in) :: T, J, dx, theta, csoil, rhosoil, h0, thre, the, he, b
    real(r_2) :: thetal, S

    S = (theta-(the-thre))/thre

    if (T<Tfrz(min(S,one),he,b)) then
       !thetal = thetalmax(T,min(S,one),he,b,thre,the)
       thetal = thetalmax(T,S,he,b,thre,the)
    else
       thetal = theta
    endif

    GTfrozen = -J + dx*(csoil*rhosoil*T + rhow*(-lambdaf + csice*T)*(theta - thetal) + &
         cswat*rhow*T*thetal) + cswat*T*h0*rhow*(1 - (theta - thetal)/thre) + &
         (h0*rhow*(-lambdaf + csice*T)*(theta - thetal)/thre)

  END FUNCTION GTfrozen

  !**********************************************************************************************************************

  REAL(r_2) ELEMENTAL PURE FUNCTION JSoilLayer(T, dx, theta, csoil, rhosoil, h0, thre, the, he, b)
    ! JSsoilLayer = sensible heat + latent heat (total energy in soil layer J/m2)
    USE sli_numbers, ONLY: one, csice, cswat, rhow, lambdaf

    IMPLICIT NONE

    real(r_2), intent(in) :: T,  dx, theta, csoil, rhosoil, h0, thre, the, he, b
    real(r_2) :: thetal, S

    S = (theta-(the-thre))/thre

    if (T<Tfrz(S,he,b)) then
       !thetal = thetalmax(T,min(S,one),he,b,thre,the)
       thetal = thetalmax(T,S,he,b,thre,the)
    else
       thetal = theta
    endif

    JSoilLayer =  dx*(csoil*rhosoil*T + rhow*(-lambdaf + csice*T)*(theta - thetal) + &
         cswat*rhow*T*thetal) + cswat*T*h0*rhow*(one - (theta - thetal)/thre) + &
         (h0*rhow*(-lambdaf + csice*T)*(theta - thetal)/thre)

  END FUNCTION JSoilLayer


  !**********************************************************************************************************************

  ! Using an elemental subroutine as a function argument is not Fortran90 standard.
  ! This would have been rh0_sol in function rtbis. intel's ifort and nag's f95 accept it but gfortran does not.
  ! Write rh0_sol into rtbis -> rtbis_rh0
  ! It does not check that f(x1) and f(x2) have different signs and not if iteration > MAXIT

  REAL(r_2) ELEMENTAL PURE FUNCTION rtbis_Tfrozen(J, dxsoil, theta,csoil, rhosoil, h0, thre, the, he, b, x1, x2)

    IMPLICIT NONE

    real(r_2), intent(in) :: J, dxsoil, theta, csoil, rhosoil, h0, thre, the, he, b
    REAL(r_2), INTENT(IN) :: x1, x2

    INTEGER(i_d), PARAMETER :: MAXIT=80
    INTEGER(i_d) :: k
    REAL(r_2)    :: dx, f, fmid, xmid

    fmid = GTfrozen(x2, J, dxsoil, theta, csoil, rhosoil, h0, thre, the, he, b)
    f    = GTfrozen(x1, J, dxsoil, theta, csoil, rhosoil, h0, thre, the, he, b)
    if (f < zero) then
       rtbis_Tfrozen = x1
       dx            = x2-x1
    else
       rtbis_Tfrozen = x2
       dx            = x1-x2
    end if
    do k=1, MAXIT
       dx   = dx*half
       xmid = rtbis_Tfrozen+dx
       fmid = GTfrozen(xmid, J, dxsoil, theta, csoil, rhosoil, h0, thre, the, he, b)
       if (fmid <= zero) rtbis_Tfrozen = xmid
       ! if (abs(dx) < xacc .or. fmid == zero) RETURN
       if (abs(fmid) < one) then
          rtbis_Tfrozen = xmid
          return
       endif
    end do

  END FUNCTION rtbis_Tfrozen

  !*********************************************************************************************************************

  REAL(r_2) ELEMENTAL PURE FUNCTION thetalmax(Tin,S,he,b,thre,the)
    ! determines maximum liquid water content, given T
    USE sli_numbers, ONLY: gravity, lambdaf, csol, Rgas, Tzero

    IMPLICIT NONE

    real(r_2), intent(in) :: Tin, S, he, b, thre, the
    real(r_2)             :: opsi, psi, h, T

    T    = min(Tfrz(min(S,one),he,b),Tin)
    opsi = -csol *Rgas *(T+Tzero)/gravity ! osmotic potential (m)
    psi  = lambdaf*T/(gravity*(T+Tzero))  ! matric potential in presence of ice
    h    = psi-opsi                       ! moisture potential in presence of ice

    ! thetalmax = thre*(h/he)**(-1/b) + (the-thre)
    thetalmax = thre*exp(-one/b*log(h/he)) + (the-thre)
    if (S > one) thetalmax = S * thetalmax

  END FUNCTION thetalmax

  !**********************************************************************************************************************

  REAL(r_2) ELEMENTAL PURE FUNCTION Tthetalmax(thetal,Tin,S,he,b,thre,the)
    ! determines T, given maximum liquid water content
    USE sli_numbers, ONLY: gravity, lambdaf, csol, Rgas, Tzero

    IMPLICIT NONE

    real(r_2), intent(in) :: thetal,S,he,b,thre,the,Tin
    real(r_2)             :: opsi, psi, h, T, Tfreezing

    !integer(i_d) :: k

    T = Tin
    !do k=1,20
    Tfreezing  = Tfrz(S,he,b)
    T          = min(Tfreezing,T)
    opsi       = -csol *Rgas *(T+Tzero)/gravity                ! osmotic potential (m)
    ! h          = he*(max((thetal-(the-thre)),0.01_r_2)/thre)**(-b) ! moisture potential in presence of ice
    h          = he*exp(-b*log(max((thetal-(the-thre)),0.01_r_2)/thre)) ! moisture potential in presence of ice
    psi        = h + opsi
    Tthetalmax = psi*(gravity*(T+Tzero)) / lambdaf
    !T = T + 0.1*(Tthetalmax-T)
    !enddo

  END FUNCTION Tthetalmax

  !**********************************************************************************************************************

  REAL(r_2) ELEMENTAL PURE FUNCTION weight(parin, h, K, phi, dz)

    IMPLICIT NONE

    TYPE(params), INTENT(IN) :: parin
    REAL(r_2),    INTENT(IN) :: h
    REAL(r_2),    INTENT(IN) :: K
    REAL(r_2),    INTENT(IN) :: phi
    REAL(r_2),    INTENT(IN) :: dz
    ! Get conductivity weighting for gravity flux calculations.
    ! Definitions of arguments:
    ! l   - land point
    ! j   - soil type no.
    ! h   - matric head.
    ! K   - conductivity.
    ! phi - MFP.
    ! dz  - flow path length.
    LOGICAL   :: done
    REAL(r_2) :: a, hz, Khz, Kz, phiz, w, x

    done = .false.
    hz   = h-gf*dz ! gf is gravity fac in direction of dz
    if (h<parin%he) then
       a = parin%lam*parin%eta
       x = -gf*dz/h
       if (a <= 3.0_r_2 .or. x*(a-3.0_r_2) <= 4.0_r_2) then ! use predetermined approx.
          w    = (60.0_r_2+x*(70.0_r_2+10.0_r_2*a+x*(16.0_r_2+a*(5.0_r_2+a))))/ &
               (120.0_r_2+x*(120.0_r_2+x*(22.0_r_2+2.0_r_2*a**2)))
          done = .true.
       end if
    end if
    if (.not. done) then
       call hyofh(hz, parin%lam, parin%eta, parin%Ke, parin%he, Kz, Khz, phiz) ! accurate but slower
       w = -((phiz-phi)/(gf*dz)+K)/(Kz-K)
    end if
    weight = min(max(w,zero),one)

  END FUNCTION weight

  !**********************************************************************************************************************

  FUNCTION zerovars()

    ! Sets all fields of type vars to zero

    IMPLICIT NONE

    TYPE(vars) :: zerovars

    zerovars%isat      = 0
    zerovars%h         = zero
    zerovars%phi       = zero
    zerovars%phiS      = zero
    zerovars%K         = zero
    zerovars%KS        = zero
    zerovars%Dv        = zero
    zerovars%cvsat     = zero
    zerovars%rh        = zero
    zerovars%phiv      = zero
    zerovars%phivS     = zero
    zerovars%kH        = zero
    zerovars%kE        = zero
    zerovars%kth       = zero
    zerovars%csoil     = zero
    zerovars%eta_th    = zero
    zerovars%hS        = zero
    zerovars%rhS       = zero
    zerovars%sl        = zero
    zerovars%cv        = zero
    zerovars%cvsatT    = zero
    zerovars%cvS       = zero
    zerovars%kv        = zero
    zerovars%iice      = 0
    zerovars%thetai    = zero
    zerovars%thetal    = zero
    zerovars%phiT      = zero
    zerovars%KT        = zero
    zerovars%lambdav   = zero
    zerovars%lambdaf   = zero
    zerovars%he        = zero
    zerovars%phie      = zero
    zerovars%Ksat      = zero
    zerovars%dthetaldT = zero
    zerovars%Tfrz      = zero
    zerovars%csoileff  = zero
    zerovars%zsat      = zero
    zerovars%macropore_factor = zero

  END FUNCTION zerovars

  !**********************************************************************************************************************

  subroutine bracket(x, xval, left, right)

    !*****************************************************************************80
    !
    !! BRACKET searches a sorted vector for successive brackets of a value.
    !
    !  Discussion:
    !
    !    A vector is an array of double precision real values.
    !    If the values in the vector are thought of as defining intervals
    !    on the real line, then this routine searches for the interval
    !    nearest to or containing the given value.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    06 April 1999
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real(r_2) X(N), an array sorted into ascending order.
    !
    !    Input, real(r_2) XVAL(M), values to be bracketed.
    !
    !    Output, integer(i_d) LEFT(M), RIGHT(M), the results of the search.
    !    Either:
    !      XVAL < X(1), when LEFT = 1, RIGHT = 2;
    !      X(N) < XVAL, when LEFT = N-1, RIGHT = N;
    !    or
    !      X(LEFT) <= XVAL <= X(RIGHT).
    !
    implicit none

    real(r_2), dimension(:), intent(in) :: x
    real(r_2), dimension(:), intent(in) :: xval
    integer(i_d), dimension(size(xval)), intent(out) :: left
    integer(i_d), dimension(size(xval)), intent(out) :: right

    integer(i_d) :: n, m
    integer(i_d) :: i, j
    logical :: istheend

    n = size(x)
    m = size(xval)
    do j=1, m
       istheend = .true.
       do i = 2, n - 1
          if ( xval(j) < x(i) ) then
             left(j)  = i - 1
             right(j) = i
             istheend = .false.
             exit
          end if
       end do
       if (istheend) then
          left(j)  = n - 1
          right(j) = n
       endif
    end do

    return

  end subroutine bracket

  function spline_b(tval, tdata, ydata)

    !*****************************************************************************80
    !
    !! SPLINE_B_VAL evaluates a cubic B spline approximant.
    !
    !  Discussion:
    !
    !    The cubic B spline will approximate the data, but is not
    !    designed to interpolate it.
    !
    !    In effect, two "phantom" data values are appended to the data,
    !    so that the spline will interpolate the first and last data values.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    11 February 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Carl deBoor,
    !    A Practical Guide to Splines,
    !    Springer, 2001,
    !    ISBN: 0387953663.
    !
    !  Parameters:
    !
    !    Input, real(r_2) TVAL, the points at which the spline is to be evaluated.
    !
    !    Input, real(r_2) TDATA(NDATA), the abscissas of the data.
    !
    !    Input, real(r_2) YDATA(NDATA), the data values.
    !
    !    Output, real(r_2) SPLINE_B, the values of the function at TVAL.
    !
    implicit none

    real(r_2), dimension(:), intent(in) :: tval
    real(r_2), dimension(:), intent(in) :: tdata
    real(r_2), dimension(:), intent(in) :: ydata
    real(r_2), dimension(size(tval))    :: spline_b

    integer(i_d) :: ndata
    real(r_2),    dimension(size(tval)) :: bval
    integer(i_d), dimension(size(tval)) :: left
    integer(i_d), dimension(size(tval)) :: right
    real(r_2),    dimension(size(tval)) :: u

    ndata = size(tdata)
    !
    !  Find the nearest interval [ TDATA(LEFT), TDATA(RIGHT) ] to TVAL.
    !
    call bracket(tdata, tval, left, right)
    !
    !  Evaluate the 5 nonzero B spline basis functions in the interval,
    !  weighted by their corresponding data values.
    !
    u = ( tval - tdata(left) ) / ( tdata(right) - tdata(left) )
    spline_b = 0.0_r_2
    !
    !  B function associated with node LEFT - 1, (or "phantom node"),
    !  evaluated in its 4th interval.
    !
    bval = ( ( (     - 1.0_r_2   &
         * u + 3.0_r_2 ) &
         * u - 3.0_r_2 ) &
         * u + 1.0_r_2 ) / 6.0_r_2

    where ( 0 < left-1 )
       spline_b = spline_b + ydata(left-1) * bval
    elsewhere
       spline_b = spline_b + ( 2.0_r_2 * ydata(1) - ydata(2) ) * bval
    end where
    !
    !  B function associated with node LEFT,
    !  evaluated in its third interval.
    !
    bval = ( ( (       3.0_r_2   &
         * u - 6.0_r_2 ) &
         * u + 0.0_r_2 ) &
         * u + 4.0_r_2 ) / 6.0_r_2

    spline_b = spline_b + ydata(left) * bval
    !
    !  B function associated with node RIGHT,
    !  evaluated in its second interval.
    !
    bval = ( ( (     - 3.0_r_2   &
         * u + 3.0_r_2 ) &
         * u + 3.0_r_2 ) &
         * u + 1.0_r_2 ) / 6.0_r_2

    spline_b = spline_b + ydata(right) * bval
    !
    !  B function associated with node RIGHT+1, (or "phantom node"),
    !  evaluated in its first interval.
    !
    bval = u**3 / 6.0_r_2

    where ( right+1 <= ndata )
       spline_b = spline_b + ydata(right+1) * bval
    elsewhere
       spline_b = spline_b + ( 2.0_r_2 * ydata(ndata) - ydata(ndata-1) ) * bval
    end where

    return

  end function spline_b

  !**********************************************************************************************************************

  FUNCTION mean_1d(dat, mask)

    IMPLICIT NONE

    REAL(r_2), DIMENSION(:),           INTENT(IN)  :: dat
    LOGICAL,   DIMENSION(:), OPTIONAL, INTENT(IN)  :: mask
    REAL(r_2)                                      :: mean_1d

    REAL(r_2) :: n

    LOGICAL, DIMENSION(size(dat)) :: maske

    if (present(mask)) then
       if (size(mask) /= size(dat)) then
          write(*,*) 'Error mean_1d: size(mask) /= size(dat)'
#ifdef __MPI__
          call MPI_Abort(0, 160, ierr) ! Do not know comm nor rank here
#else
          stop 160
#endif
       endif
       maske = mask
       n = real(count(maske),r_2)
    else
       maske(:) = .true.
       n = real(size(dat),r_2)
    endif
    if (n <= (1.0_r_2+tiny(1.0_r_2))) then
       write(*,*) 'mean_1d: n must be at least 2'
#ifdef __MPI__
       call MPI_Abort(0, 161, ierr) ! Do not know comm nor rank here
#else
       stop 161
#endif
    endif

    ! Mean
    mean_1d  = sum(dat(:), mask=maske)/n

  END FUNCTION mean_1d

  FUNCTION mean_2d(dat, mask)

    IMPLICIT NONE

    REAL(r_2), DIMENSION(:,:),           INTENT(IN)  :: dat
    LOGICAL,   DIMENSION(:,:), OPTIONAL, INTENT(IN)  :: mask
    REAL(r_2)                                        :: mean_2d

    REAL(r_2) :: n

    LOGICAL, DIMENSION(size(dat,1),size(dat,2)) :: maske

    if (present(mask)) then
       if (size(mask) /= size(dat)) then
          write(*,*) 'Error mean_2d: size(mask) /= size(dat)'
#ifdef __MPI__
          call MPI_Abort(0, 162, ierr) ! Do not know comm nor rank here
#else
          stop 162
#endif
       endif
       maske = mask
       n = real(count(maske),r_2)
    else
       maske = .true.
       n = real(size(dat),r_2)
    endif
    if (n <= (1.0_r_2+tiny(1.0_r_2))) then
       write(*,*) 'mean_2d: n must be at least 2'
#ifdef __MPI__
       call MPI_Abort(0, 163, ierr) ! Do not know comm nor rank here
#else
       stop 163
#endif
    endif

    ! Mean
    mean_2d  = sum(dat, mask=maske)/n

  END FUNCTION mean_2d

  !**********************************************************************************************************************

  ! Nash-Sutcliffe Efficiency: 1-sum(obs-model)/sum(obs-mean(obs))
  FUNCTION nse_1d(x, y, mask)

    IMPLICIT NONE

    REAL(r_2), DIMENSION(:),           INTENT(IN) :: x, y
    LOGICAL,   DIMENSION(:), OPTIONAL, INTENT(IN) :: mask
    REAL(r_2)                                     :: nse_1d

    REAL(r_2)                                     :: xmean
    REAL(r_2), DIMENSION(size(x))                 :: v1, v2
    LOGICAL,   DIMENSION(size(x))                 :: maske

    if (present(mask)) then
       maske = mask
    else
       maske = .true.
    endif

    ! mean of x
    xmean = mean(x, mask=maske)
    ! nse
    v1 = merge(y - x    , 0.0_r_2, maske)
    v2 = merge(x - xmean, 0.0_r_2, maske)
    nse_1d = 1.0_r_2 - dot_product(v1,v1) / dot_product(v2,v2)

  END FUNCTION nse_1d


  FUNCTION nse_2d(x, y, mask)

    IMPLICIT NONE

    REAL(r_2), DIMENSION(:,:),           INTENT(IN) :: x, y
    LOGICAL,   DIMENSION(:,:), OPTIONAL, INTENT(IN) :: mask
    REAL(r_2)                                       :: nse_2d

    REAL(r_2)                                       :: xmean
    REAL(r_2), DIMENSION(size(x,1),size(x,2))       :: v1, v2
    LOGICAL,   DIMENSION(size(x,1),size(x,2))       :: maske

    INTEGER(i_d) :: i

    if (present(mask)) then
       maske = mask
    else
       maske = .true.
    endif

    ! nse
    v1 = merge(y - x, 0.0_r_2, maske)
    do i=1, size(x,1)
       xmean   = mean(x(i,:), mask=maske(i,:))
       v2(i,:) = merge(x(i,:) - xmean, 0.0_r_2, maske(i,:))
    enddo
    nse_2d = 1.0_r_2 - sum(v1*v1) / sum(v2*v2)

  END FUNCTION nse_2d

END MODULE sli_utils
