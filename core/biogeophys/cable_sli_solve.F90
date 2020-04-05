! Original code by P.J. Ross 2001-2007  See: Ross, P.J. (2003)
!     Modeling soil water and solute transport - fast, simplified numerical solutions. Agron. J. 95:1352-1361.
!
! Modified by V. Haverd 2008: matrix expanded by factor of 2 to allow solution of coupled heat/moisture fluxes
! Explicit calculation of heat/moisture fluxes at surface
! Switchable option for litter
! Isotope subroutine (V. Haverd and M. Cuntz, September 2008)
!
! Frozen Soil (V. Haverd, M. Cuntz Sep 2010 - Jan 2011)
! Pond lumped with top soil layer (V. Haverd Jan 2011)
! Snow included in layers -2:0 (V. Haverd Jan 2011)
! Convergence test for dthetaldT (V.Haverd Feb 2011)
! Include heat advection by liquid water flux (V.Haverd Feb 2011)

! Code cleaned and isotope routine updated (V. Haverd, M. Cuntz Aug-Sep 2012)
! lots of changes, especially relating to ponding and freezing together
! Outstanding problems with surface runoff and advection (set h0max large and turn off advection for now)
!
! Sep 17 2012 (V.Haverd)
! remove separate solution for maxpond (h0 > h0max)
! Instead, convert excess pond to runoff at t=tfin, and correct energy stores for loss of pond
! Advection and surface runoff now functioning
!
! Dec 31 2012 (V.Haverd)
! bug fixes to improve energy conservation associated with change in ice,
! especially correction at time of pond disappearance (removal of negative pond)
!
! Jan 5 2013 (V.Haverd)
! Remove iteration over frozen soil within "(while iok==0)" loop
! For frozen soil, at the time of updating variables, evaluate new T and ice content,
! consistent with J0 + deltaJ, where deltaJ is the change in energy obtained from the matrix solution (=LHS_h*dt)
!
! Jan 7-8 2013 (V.Haverd)
! Calls to soil_snow albedo and density
! Improve iteration convergence at time of updating new T and ice content in frozen soil
! Bug fix relating to moisture conservation in snow
!
! Jan 14-15 2013 (V. Haverd)
! Modified definition of Sliq for frozen soil (in hyofS): now defined relative to (theta_sat - theta_r - theta_ice)
! Allow for snow-pack initalisation in the absence of pond
!
! Feb 7 2013 (V.Haverd)
! Revised formulation of dphidT for frozen soil (see hyofS): fixes a few remaining -ve T spikes in frozen soil
!
! Feb 9 2013 (V.Haverd)
! Moved snow layer adjustments to subroutine
!
! Feb 24 2013 (V. Haverd)
! Snow structure expanded to accommodate 3 layers
!
! March 30 (V. Haverd)
!
! Jan 7 2014 (V. Haverd)
! Adjust snow surface bcs such that snow surface temperature can't be more than 0
!
! Jan 8 2014 (V Haverd)
! enable 2nd snow layer so that if snowpack exceeds 3 cm , the excess goes into the layer below.
!
! March 2014 (V. Haverd)
! Extract SEB calcs to subroutine in utils, and trial Force-Restore method
!
! October 2016 (S Wales)
! Restructured code of solve into several subroutines

MODULE sli_solve

  USE cable_def_types_mod, ONLY: r_2, i_d
  USE sli_numbers,         ONLY: &
       experiment, &
       zero, one, two, half, thousand, e5, &  ! numbers
       Tzero, rlambda, lambdaf, lambdas, Dva, rhocp, rhow, gf, hmin, & ! parameters
       csice , cswat, &
       snmin, nsnow_max, &
       params, vars_aquifer, vars_met, vars, vars_snow, & ! types
       dSfac, h0min, Smax, h0max, dSmax, dSmaxr, dtmax, dSmax, dSmaxr, & ! numerical limits
       dtmax, dtmin, dTsoilmax, dTLmax, nsteps_ice_max, nsteps_max,tol_dthetaldT, &
       hbot, botbc, &
       Mw, Rgas, & ! boundary condition
       freezefac
  USE sli_utils,           ONLY: &
       Sofh, hyofh, hyofS, litter_props, massman_sparse, tri, &
       Tfrz, thetalmax, dthetalmaxdTh, &
       getfluxes_vp, getheatfluxes, flux, sol, &
       tri, setsol, zerovars, &
       esat_ice, slope_esat_ice, rtbis_Tfrozen, GTFrozen, &
       JSoilLayer, SEB
  USE cable_IO_vars_module, ONLY: wlogn
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: solve ! solution routine

  INTEGER(i_d), DIMENSION(:), ALLOCATABLE :: nless, n_noconverge ! global counters

  ! Definitions of public entities and private parameters (see above for default
  ! values):
  ! botbc    - bottom boundary condn for water; "constant head", "free drainage",
  !    "seepage", or "zero flux". Constant head means that matric head h
  !    is specified. Free drainage means zero gradient of matric head,
  !    i.e. unit hydraulic gradient. Seepage means zero flux when the
  !    matric head is below zero and an upper limit of zero for the head.
  ! h0max    - max pond depth allowed before runoff.
  ! hbot    - matric head at bottom of profile when botbc set to "constant head".
  ! dSmax    - max change in S (the "effective saturation") of any unsaturated
  !    layer to aim for each time step; controls time step size.
  ! dSmaxr   - maximum negative relative change in S each time step. This
  !    parameter helps avoid very small or negative S.
  ! dtmax    - max time step allowed.
  ! dsmmax   - max solute change per time step (see dSmax); user should set this
  !    according to solute units used. Units for different solutes can be
  !    scaled by the user (e.g. to an expected max of around 1.0).
  ! dSfac    - a change in S of up to dSfac*dSmax is accepted.
  ! dpmaxr   - relative change in matric flux potential (MFP) phi that is
  !    accepted for convergence when finding head h at soil interfaces.
  ! h0min    - min (negative) value for surface pond when it empties.
  ! Smax    - max value for layer saturation to allow some overshoot.
  ! solve    - sub to call to solve RE
  !

CONTAINS

  SUBROUTINE get_fluxes_and_derivs( &
       mp, qprec, qprec_snow, n, dx, h0, &
       Tsoil, &
       qh, vmet, vlit, vsnow, var, T0, Tsurface, &
       SL, Tl, &
       plit, par, &
       qali, &
       qvh, &
       qevap, &
       again, getq0,getqn,init, ns, nsat, &
       nsatlast, &
       phip, qpme, hint, phimin, &
       qexd, &
       q, qya, qyb, qTa, qTb,qhya, qhyb, qhTa, qhTb, qadv, qadvya, qadvyb, &
       qadvTa, qadvTb, vtmp, qliq, qv, qvT, qlya, qlyb, &
       qvya, qvyb, qlTb, qvTa, qvTb, &
       vtop, vbot, &
       lE0, G0, &
       Epot, surface_case, &
       iflux, j, kk, &
       advection, dTqwdTa, dTqwdTb, Tqw, keff, &
       hice &
       )
    
    IMPLICIT NONE
    
    INTEGER(i_d)                                           :: mp
    REAL(r_2),      DIMENSION(1:mp)                        :: qprec
    REAL(r_2),      DIMENSION(1:mp)                        :: qprec_snow
    INTEGER(i_d)                                           :: n
    REAL(r_2),      DIMENSION(1:n) :: dx
    REAL(r_2),      DIMENSION(1:mp)                        :: h0
    REAL(r_2),      DIMENSION(1:n) :: Tsoil
    REAL(r_2),      DIMENSION(-nsnow_max:n) :: qh
    TYPE(vars_met), DIMENSION(1:mp)                        :: vmet
    TYPE(vars),     DIMENSION(1:mp)                        :: vlit
    TYPE(vars_snow), DIMENSION(1:mp)                       :: vsnow
    TYPE(vars),     DIMENSION(1:n) :: var
    REAL(r_2),      DIMENSION(1:mp)                        :: T0, Tsurface
    REAL(r_2),      DIMENSION(1:mp)                        :: SL, Tl
    TYPE(params),   DIMENSION(1:mp)                        :: plit
    TYPE(params),   DIMENSION(1:n) :: par
    REAL(r_2),      DIMENSION(1:mp), OPTIONAL                        :: qali
    REAL(r_2),      DIMENSION(-nsnow_max:n) :: qvh
    REAL(r_2),    DIMENSION(1:mp)                          :: qevap
    LOGICAL,      DIMENSION(1:mp)                          :: again, getq0,getqn,init
    INTEGER(i_d), DIMENSION(1:mp)                          :: ns, nsat, nsatlast
    REAL(r_2),    DIMENSION(1:mp)                          :: phip
    REAL(r_2),    DIMENSION(1:mp)                          :: qpme
    REAL(r_2),    DIMENSION(1:n) :: hint, phimin, qexd
    REAL(r_2),    DIMENSION(-nsnow_max:n) :: q, qya, qyb, qTa, qTb,qhya, qhyb, qhTa, qhTb
    REAL(r_2),    DIMENSION(-nsnow_max:n) :: qadv, qadvya, qadvyb, qadvTa, qadvTb
    TYPE(vars)                                             :: vtmp
    REAL(r_2),    DIMENSION(-nsnow_max:n) :: qliq, qv, qvT, qlya, qlyb, qvya, qvyb, qlTb, qvTa, qvTb
    TYPE(vars),         DIMENSION(1:mp)                    :: vtop, vbot
    REAL(r_2),          DIMENSION(1:mp)                    :: lE0, G0, Epot
    INTEGER(i_d),       DIMENSION(1:mp)                    :: surface_case
    INTEGER(i_d),       DIMENSION(1:mp)                    :: iflux
    INTEGER(i_d)                                           :: j, kk
    INTEGER(i_d)                                           :: advection ! switches
    REAL(r_2)                                              :: dTqwdTa, dTqwdTb, Tqw, keff
    REAL(r_2),          DIMENSION(1:mp)                    :: hice


    !----- get fluxes and derivs
    ! get surface condition
    ! ns==1 if no pond or full pond, i.e. do not solve for change in pond height
    ! ns==0 then change for pond height
    if ((var(1)%phi <= phip(kk) .and. h0(kk) <= zero .and. nsat(kk) < n) .or. &
         (var(1)%isat==0)) then ! no ponding
       ns(kk)    = 1 ! start index for eqns
    else ! ponding
       ns(kk)    = 0
       var(1)%phi = (one+e5)*var(1)%phie+(h0(kk)-hice(kk))*var(1)%Ksat
       TL(kk)    = vmet(kk)%Ta ! initialise pond T
       vtop(kk)%isat    = 1
       vtop(kk)%h       = h0(kk)
       vtop(kk)%phi     = max((var(1)%phie -var(1)%he*var(1)%Ksat), &
            (one+e5)*var(1)%phie)+(h0(kk)-hice(kk))*var(1)%Ksat
       vtop(kk)%K       = var(1)%Ksat
       vtop(kk)%rh      = one
       vtop(kk)%lambdav = rlambda
       vtop(kk)%lambdaf = lambdaf
       ! var(1)%phi = var(1)%phie
       ! calculates phi1,eff (pond + top soil layer)   !!vh!! does this get used???
       call flux(par(1), vtop(kk), var(1), half*dx(1), &
            q(0), qya(0), qyb(0), qTa(0), qTb(0))


    endif

    if (Sl(kk) >= one) then
       vlit(kk)%isat = 1
       vlit(kk)%h    = plit(kk)%he
    endif

    ! get bottom boundary condn
    if (botbc=="seepage") then
       vtmp = zerovars()
       vbot = spread(vtmp,1,mp)
       if (var(n)%h > -half*gf*dx(n)) then
          getqn(kk) = .true.
          vbot(kk)%isat    = 1
          vbot(kk)%phi     = (zero-var(n)%he)*var(n)%Ksat+var(n)%phie ! (special for frozen soil)
          vbot(kk)%K       = var(n)%Ksat
          vbot(kk)%rh      = one
          vbot(kk)%lambdav = rlambda
          vbot(kk)%lambdaf = lambdaf
       else
          getqn(kk) = .false.
       endif
    end if

    ! Fluxes and derivatives at the air/surface interface
    if (vsnow(kk)%nsnow > 0) then
       surface_case(kk) = 2
    else
       surface_case(kk) = 1
    endif
    select case (surface_case(kk))
    case (1) ! no snow
       CALL SEB(n, par(:), vmet(kk), vsnow(kk), var(:), qprec(kk), qprec_snow(kk), dx(:), &
            h0(kk), Tsoil(:), &
            Tsurface(kk), G0(kk), lE0(kk),Epot(kk),  &
            q(0), qevap(kk), qliq(0), qv(0), &
            qyb(0), qTb(0), qlyb(0), qvyb(0), qlTb(0), qvTb(0), qh(0), &
            qadv(0), qhyb(0), qhTb(0), qadvyb(0), qadvTb(0))
       qya(0)    = zero
       qTa(0)    = zero
       qlya(0)   = zero
       qvya(0)   = zero
       qvTa(0)   = zero
       qhya(0)   = zero
       qhTa(0)   = zero
       qadvya(0) = zero
       qadvTa(0) = zero

    case (2) ! snow
       CALL SEB(n, par(:), vmet(kk), vsnow(kk), var(:), qprec(kk), qprec_snow(kk), dx(:), &
            h0(kk), Tsoil(:), &
            Tsurface(kk), G0(kk), lE0(kk), Epot(kk),  &
            q(-vsnow(kk)%nsnow), qevap(kk), qliq(-vsnow(kk)%nsnow), qv(-vsnow(kk)%nsnow), &
            qyb(-vsnow(kk)%nsnow), qTb(-vsnow(kk)%nsnow), qlyb(-vsnow(kk)%nsnow), &
            qvyb(-vsnow(kk)%nsnow), qlTb(-vsnow(kk)%nsnow), qvTb(-vsnow(kk)%nsnow), &
            qh(-vsnow(kk)%nsnow), qadv(-vsnow(kk)%nsnow), qhyb(-vsnow(kk)%nsnow), &
            qhTb(-vsnow(kk)%nsnow), qadvyb(-vsnow(kk)%nsnow), qadvTb(-vsnow(kk)%nsnow))
       qya(-vsnow(kk)%nsnow)    = zero
       qTa(-vsnow(kk)%nsnow)    = zero
       qlya(-vsnow(kk)%nsnow)   = zero
       qvya(-vsnow(kk)%nsnow)   = zero
       qvTa(-vsnow(kk)%nsnow)   = zero
       qhya(-vsnow(kk)%nsnow)   = zero
       qhTa(-vsnow(kk)%nsnow)   = zero
       qadvya(-vsnow(kk)%nsnow) = zero
       qadvTa(-vsnow(kk)%nsnow) = zero

    case default
       write(*,*) "solve: illegal surface case."
       stop 2
    end select ! surface_case
    qpme(kk) = q(0) ! water flux into top of soil column
    ! finished all the surfaces

    ! get moisture fluxes and derivatives (at time t=0, i.e. q0 etc.)

    call getfluxes_vp(n, dx(1:n), vtop(kk), vbot(kk), par(1:n), var(1:n), & ! moisture fluxes
         hint(1:n), phimin(1:n), q(0:n), qya(0:n), qyb(0:n), qTa(0:n), qTb(0:n), &
         qliq(0:n), qlya(0:n), qlyb(0:n), qv(0:n), qvT(0:n), qvh(0:n), qvya(0:n), &
         qvyb(0:n), &
         iflux(kk), init(kk), getq0(kk), getqn(kk), Tsoil(1:n), T0(kk), nsat(kk), nsatlast(kk))
    qTa(n) = zero
    qTb(n) = zero
    qvTa(1:n) = qTa(1:n)
    qvTb(1:n) = qTb(1:n)
    qlTb(1:n) = zero


    ! get  fluxes heat and derivatives (at time t=0, i.e. q0 etc.)
    call getheatfluxes(n, dx(1:n), &
         qh(0:n), qhya(0:n), qhyb(0:n), qhTa(0:n), qhTb(0:n), &
         var(1:n), Tsoil(1:n), &
         q(0:n), qya(0:n), qyb(0:n), qTa(0:n), qTb(0:n), &
         qadv(0:n),qadvya(0:n), qadvyb(0:n), qadvTa(0:n), qadvTb(0:n), &
         advection) ! heat fluxes

    ! get heat and vapour fluxes and derivatives in snow-pack
    if (vsnow(kk)%nsnow>1) then
       ! vapour flux at interfaces between snow layers
       do j=1, vsnow(kk)%nsnow-1
          keff = 2_r_2*((vsnow(kk)%kE(j)/(thousand*lambdaf))*(vsnow(kk)%kE(j+1)/(thousand*lambdaf))/ &
               ((vsnow(kk)%kE(j)/(thousand*lambdaf))*vsnow(kk)%depth(j+1)+(vsnow(kk)%kE(j+1)/ &
               (thousand*lambdaf))*vsnow(kk)%depth(j)) )
          q(j-vsnow(kk)%nsnow) = keff*(vsnow(kk)%tsn(j)-vsnow(kk)%tsn(j+1))
          qTa(j-vsnow(kk)%nsnow) = merge(keff,zero,vsnow(kk)%hliq(j)<=zero)
          qTb(j-vsnow(kk)%nsnow) = merge(-keff,zero,vsnow(kk)%hliq(j+1)<=zero)
          qya(j-vsnow(kk)%nsnow) = zero
          qyb(j-vsnow(kk)%nsnow) = zero
          ! conductive heat flux at interface between snow layers
          keff = 2_r_2*(vsnow(kk)%kth(j+1)*vsnow(kk)%kth(j))/ &
               (vsnow(kk)%kth(j+1)*vsnow(kk)%depth(j)+vsnow(kk)%kth(j)*vsnow(kk)%depth(j+1))  ! check this!
          qh(j-vsnow(kk)%nsnow) = keff*(vsnow(kk)%tsn(j)-vsnow(kk)%tsn(j+1))
          qhTa(j-vsnow(kk)%nsnow) = merge(zero,keff,vsnow(kk)%hliq(j)>zero)
          qhTb(j-vsnow(kk)%nsnow) = merge(zero,-keff,vsnow(kk)%hliq(j+1)>zero)

          ! advective heat flux at interface between snow layers
          Tqw  = merge(vsnow(kk)%tsn(j), vsnow(kk)%tsn(j+1), q(j-vsnow(kk)%nsnow)>zero)
          dTqwdTb = merge(zero,one, q(j-vsnow(kk)%nsnow)>zero)
          dTqwdTa = merge(one,zero, q(j-vsnow(kk)%nsnow)>zero)
          if (vsnow(kk)%hliq(j)>zero) then
             qadv(j-vsnow(kk)%nsnow) = rhow*q(j-vsnow(kk)%nsnow)*cswat*Tqw
             qadvTa(j-vsnow(kk)%nsnow) = zero
          else
             qadv(j-vsnow(kk)%nsnow) = rhow*q(j-vsnow(kk)%nsnow)*cswat*Tqw
             qadvTa(j-vsnow(kk)%nsnow) = rhow*cswat*q(j-vsnow(kk)%nsnow)*dTqwdTa  + &
                  rhow*cswat*Tqw*qTa(j-vsnow(kk)%nsnow)
          endif
          qadvTb(0) = rhow*cswat*q(j-vsnow(kk)%nsnow)*dTqwdTb + rhow*cswat*Tqw*qTb(j-vsnow(kk)%nsnow)

          qh(j-vsnow(kk)%nsnow) = qh(j-vsnow(kk)%nsnow) + qadv(j-vsnow(kk)%nsnow)
          qhTa(j-vsnow(kk)%nsnow) = qhTa(j-vsnow(kk)%nsnow) +  qadvTa(j-vsnow(kk)%nsnow)
          qhTb(j-vsnow(kk)%nsnow) = qhTb(j-vsnow(kk)%nsnow) +  qadvTb(j-vsnow(kk)%nsnow)
       enddo
    endif ! end fluxes at snow/snow interfaces

    if (vsnow(kk)%nsnow.ge.1) then
       ! vapour flux at soil/snow interface
       if (var(1)%isat==1) then
          keff = vsnow(kk)%kE(vsnow(kk)%nsnow)/(thousand*lambdaf)/vsnow(kk)%depth(vsnow(kk)%nsnow)/2_r_2
       endif
       if (var(1)%isat==0) then
          keff = 2_r_2*((vsnow(kk)%kE(vsnow(kk)%nsnow)/(thousand*lambdaf))*(var(1)%kE/thousand/var(1)%lambdav))/ &
               ((vsnow(kk)%kE(vsnow(kk)%nsnow)/(thousand*lambdaf))*dx(1)+(var(1)%kE/thousand/var(1)%lambdav)* &
               vsnow(kk)%depth(vsnow(kk)%nsnow))
       endif
       q(0) = keff*(vsnow(kk)%tsn(vsnow(kk)%nsnow)-Tsoil(1))
       qTa(0) = keff
       qTb(0) =-keff

       qya(0) = zero
       qyb(0) = zero
       qv(0)   = q(0)
       qvyb(0) = qyb(0)
       qvTb(0) = qTb(0)
       qvya(0) = qya(0)
       qvTa(0) = qTa(0)
       qliq(0) = zero
       qlyb(0) = zero
       qlTb(0) = zero
       qlya(0) = zero
       qya(0) = zero;

       if (vsnow(kk)%hliq(vsnow(kk)%nsnow)>zero) then
          qhTa(0) = zero
          qTa(0) = zero
       endif

       ! conductive heat flux at snow/soil interface  ! check this!
       keff = 2._r_2*(vsnow(kk)%kth(vsnow(kk)%nsnow)*var(1)%kth)/ &
            (vsnow(kk)%kth(vsnow(kk)%nsnow)*dx(1)+var(1)%kth*vsnow(kk)%depth(vsnow(kk)%nsnow))
       qh(0) = keff*(vsnow(kk)%tsn(vsnow(kk)%nsnow)-Tsoil(1))
       if (vsnow(kk)%hliq(1)>zero) then
          qhTa(0) = zero
       else
          qhTa(0) = keff
       endif
       qhTb(0) = -keff

       ! advective heat flux at snow/soil interface
       Tqw  = merge(vsnow(kk)%tsn(vsnow(kk)%nsnow), Tsoil(1), q(0)>zero)
       dTqwdTb = merge(zero,one, q(0)>zero)
       dTqwdTa = merge(one,zero, q(0)>zero)
       if (vsnow(kk)%hliq(vsnow(kk)%nsnow)>zero) then
          qadv(0) = rhow*q(0)*cswat*Tqw
          qadvTa(0) = zero
       else
          qadv(0) = rhow*q(0)*cswat*Tqw
          qadvTa(0) = rhow*cswat*q(0)*dTqwdTa  +  rhow*cswat*Tqw*qTa(0)
       endif
       qadvTb(0) = rhow*cswat*q(0)*dTqwdTb + rhow*cswat*Tqw*qTb(0)
    endif ! end of heat and vapour fluxes at soil/snow interface

    if (ns(kk)==0) then ! pond included in top soil layer
       ! change qya(1) from dq/dphi (returned by getfluxes) to dq/dh
       qya(1) = var(1)%Ksat*qya(1)
       if (advection==1) then
          qhya(1) = qhya(1) - qadvya(1)
          Tqw  = merge(Tsoil(1), Tsoil(2), q(1)>zero)
          qadvya(1) =  rhow*cswat*qya(1)*Tqw  ! apply corrected qya(1) to qadvya(1)
          qhya(1) = qhya(1) + qadvya(1)
       endif
    endif

    ! adjust for bottom boundary condition
    if (botbc=="zero flux") then
       qliq(n) = zero
       qv(n)   = zero
       q(n)    = zero
       qya(n)  = zero
       qlya(n) = zero
       qvya(n) = zero
    endif

    ! specify mositure flux at bottom of soil column (heat flux set to zero)
    if (botbc /= "constant head") then
       select case (botbc)
       case ("zero flux")
          q(n)   = zero
          qya(n) = zero
       case ("free drainage")
          q(n) = gf*var(n)%K
          if (var(n)%isat == 0) then
             qya(n) = gf*var(n)%KS
          else
             qya(n) = zero
          end if
       case ("seepage")
          if (var(n)%h <= -half*gf*dx(n)) then
             q(n)   = zero
             qya(n) = zero
          end if
       case default
          write(*,*) "solve: illegal bottom boundary condition."
          stop 2
       end select
    end if
    if (present(qali)) then
       if (qali(kk)>zero) then
          q(n)   = -qali(kk)
          qya(n) = zero
       end if
    endif
    if (experiment==7 .or. experiment==8) then
       q(n)   = q(0)
       qya(n) = zero
    endif
    if (experiment==8) qh(n) = G0(kk)

    ! adjust lower heat flux for advection
    if (advection==1) then
       qadv(n) = rhow*cswat*(Tsoil(n))*q(n)
       qadvya(n) = rhow*cswat*(Tsoil(n))*qya(n)
       qadvTa(n)= rhow*cswat*q(n)
       qh(n) = qh(n) + qadv(n)
       qhya(n) = qhya(n) + qadvya(n)
       qhTa(n) = qhTa(n) + qadvTa(n)
    else
       qadv(:)   = zero
       qadvya(:) = zero
       qadvyb(:) = zero
       qadvTa(:) = zero
       qadvTb(:) = zero
    endif

    qexd(1:n) = zero ! time derivative for root extraction (assumed fixed at input value)
    again(kk)  = .false. ! flag for recalcn of fluxes (default=false)
    !----- end get fluxes and derivs
  END SUBROUTINE get_fluxes_and_derivs
  

  SUBROUTINE estimate_timestep( &
       tfin, mp, n, dx, h0, &
       qh, &
       nsteps, vlit, vsnow, var, &
       dxL, &
       plit, par, &
       deltaTa, &
       qL, qhL, &
       again, ns, nsat, &
       nsatlast, nsteps0, dmax, dt, &
       qpme, t, &
       q, qadv, &
       tmp2d1, tmp2d2, &
       tmp1d1, tmp1d3, &
       iflux, litter, kk, &
       advection, &
       iqex &
       )
    
    IMPLICIT NONE
    
    REAL(r_2)                                              :: tfin
    INTEGER(i_d)                                           :: mp
    INTEGER(i_d)                                           :: n
    REAL(r_2),      DIMENSION(1:n) :: dx
    REAL(r_2),      DIMENSION(1:mp)                        :: h0
    REAL(r_2),      DIMENSION(-nsnow_max:n) :: qh
    INTEGER(i_d),   DIMENSION(1:mp)                        :: nsteps
    TYPE(vars),     DIMENSION(1:mp)                        :: vlit
    TYPE(vars_snow), DIMENSION(1:mp)                       :: vsnow
    TYPE(vars),     DIMENSION(1:n) :: var
    REAL(r_2),      DIMENSION(1:mp)                        :: dxL
    TYPE(params),   DIMENSION(1:mp)                        :: plit
    TYPE(params),   DIMENSION(1:n) :: par
    REAL(r_2),      DIMENSION(1:mp)                        :: deltaTa
    REAL(r_2),    DIMENSION(1:mp)                          :: qL, qhL
    LOGICAL,      DIMENSION(1:mp)                          :: again
    INTEGER(i_d), DIMENSION(1:mp)                          :: ns, nsat, nsatlast, nsteps0
    REAL(r_2),    DIMENSION(1:mp)                          :: dmax, dt
    REAL(r_2),    DIMENSION(1:mp)                          :: qpme, t
    REAL(r_2),    DIMENSION(-nsnow_max:n) :: q
    REAL(r_2),    DIMENSION(-nsnow_max:n) :: qadv
    REAL(r_2),    DIMENSION(0:n) :: tmp2d1, tmp2d2,  deltaTmax
    REAL(r_2),          DIMENSION(1:mp)                    :: tmp1d1, tmp1d3
    INTEGER(i_d),       DIMENSION(1:mp)                    :: iflux
    LOGICAL                                                :: litter
    INTEGER(i_d)                                           :: kk
    INTEGER(i_d)                                           :: advection ! switches
    REAL(r_2),          DIMENSION(1:n) :: iqex

    !----- first estimate of time step dt before the calculation
    !      gets revised after the calculation
    dmax(kk)  = zero
    tmp2d1(:) = zero ! 1st order estim of rate of change of moisture storage
    tmp2d2(:) = zero ! 1st order estim of rate of change of temperature
    ! estimate rate of change of moisture storage [m/s]
    where (var(1:n)%isat==0.and.var(1:n)%iice==0.) tmp2d1(1:n) = &
         abs(q(1:n)-q(0:n-1)-iqex(1:n))/(par(1:n)%thre*dx(1:n))
    where (var(1:n)%iice==1)  tmp2d1(1:n) =  tmp2d1(1:n)/2._r_2
    ! estimate rate of change of temperature [K/s]
    tmp2d2(1:n) = abs(qh(1:n)-qh(0:n-1))/(var(1:n)%csoileff*dx(1:n))
    if (advection==1) then
       ! first order estimate of rate of temperature change
       tmp2d2(1:n) = abs((qh(1:n)-qadv(1:n))-(qh(0:n-1)-qadv(0:n-1)))/(var(1:n)%csoileff*dx(1:n))
       deltaTmax(1:n) = tmp2d2(1:n)*(tfin-t(kk)) ! maximum  T change
    endif
    if (litter .and. ns(kk)==1 ) then ! litter , no pond
       if (vlit(kk)%isat==0) then ! estimate rate of change of moisture storage [m/s]
          write(*,*) 'Should not be here - QL 01 ', qL(kk)
          tmp2d1(0) = abs(q(0) - qL(kk))/(plit(kk)%thre*dxL(kk))
       endif
       write(*,*) 'Should not be here - QHL 01 ', qhL(kk)
       tmp2d2(0)  = abs(qh(0) - qhL(kk))/(vlit(kk)%csoileff*dxL(kk)) ! estimate rate of change of heat storage [K/s]
       tmp1d3(kk) = (dTLmax-abs(deltaTa(kk))) / tmp2d2(0)
    else
       tmp2d1(0)  = zero
       tmp2d2(0)  = zero
       tmp1d3(kk) = dtmax
    endif

    dmax(kk)   = maxval(tmp2d1(1:n),1) ! max derivative |dS/dt|
    tmp1d1(kk) = maxval(tmp2d2(1:n),1) ! max derivative |dTsoil/dt|
    if (dmax(kk) > zero) then
       dt(kk) = min(dSmax/dmax(kk), dTLmax/tmp1d1(kk), tmp1d3(kk)) ! constrained either by moisture or temp
    else ! steady state flow
       if (qpme(kk)>=q(n)) then ! if saturated soil columnn and more precip then drainige -> finish
          dt(kk) = tfin-t(kk) ! step to finish
       else ! otherwise adjust dt because change of pond height
          dt(kk) = -(h0(kk)-half*h0min)/(qpme(kk)-q(n))
       end if
       dt(kk) = min(dt(kk), dTLmax/tmp1d1(kk), tmp1d3(kk)) ! constrained by  temp
    end if

    ! check that time step is short enough to prevent melting a top snow layer
    if (vsnow(kk)%nsnow.gt.0 ) then
       ! energy required to melt ice in snow-pack
       tmp1d1(kk) = -rhow*(vsnow(kk)%hsnow(1)-vsnow(kk)%hliq(1))*(csice*vsnow(kk)%tsn(1) - lambdaf)
       if (qh(-vsnow(kk)%nsnow)*dt(kk).gt.tmp1d1(kk)) then
          dt(kk) = 0.9_r_2*tmp1d1(kk)/qh(-vsnow(kk)%nsnow)
       endif
    endif

    if (dt(kk)>dtmax) dt(kk) = dtmax ! user's limit

    ! if initial step, improve phi where S>=1
    ! might be that you get better derivatives, especially at sat/non-sat interfaces
    if (nsteps(kk)==nsteps0(kk) .and. nsat(kk)>0 .and. iflux(kk)==1) then
       again(kk) = .true.
       dt(kk)    = dtmin
    end if
    ! if fully saturated but was not fully saturated before, adjust even within each time step iteration
    if (nsat(kk)==n .and. nsatlast(kk)<n .and. iflux(kk)==1) then
       again(kk) = .true. ! profile has just become saturated so adjust phi values
       dt(kk)    = dtmin
    end if
    ! sprint to the end
    if (t(kk)+1.1_r_2*dt(kk)>tfin) then ! step to finish
       dt(kk) = tfin-t(kk)
       t(kk)  = tfin
    else
       t(kk) = t(kk)+dt(kk) ! tentative update
       if (again(kk)) t(kk) = t(kk)-dt(kk)
    end if

    !----- end estimate time step dt


  END SUBROUTINE estimate_timestep

  SUBROUTINE iflux_loop( &
       tfin, irec, mp, qprec, qprec_snow, n, dx, h0, S, thetai, &
       Jsensible, Tsoil, evap, infil, drainage, discharge, &
       qh, nsteps, vmet, vlit, vsnow, var, T0, Tsurface, Hcum, lEcum, &
       Gcum, Qadvcum, Jcol_sensible, &
       Jcol_latent_S, Jcol_latent_T, csoil, kth, phi, dxL, zdelta, SL, Tl, &
       plit, par, wex, ciso_snow, cisoice_snow, &
       qali, &
       qvsig, qlsig, qvTsig, qvh, deltaTa, &
       precip, qevap, qL, qhL, qybL, qTbL, qhTbL, qhybL, rexcol, wcol, &
       again, getq0,getqn,init, again_ice, ih0, iok, itmp, ns, nsat, &
       nsatlast, nsteps0, accel, dmax, dt, dwinfil, dwoff, fac, &
       phip, qpme, rsig, rsigdt, sig, t, hint, phimin, &
       qexd, aa, bb, cc, dd, ee, ff, gg, dy, aah, bbh, cch, ddh, eeh, ffh, ggh, &
       de, q, qya, qyb, qTa, qTb,qhya, qhyb, qhTa, qhTb, qadv, qadvya, qadvyb, &
       qadvTa, qadvTb, vtmp, qsig, qhsig, qadvsig, qliq, qv, qvT, qlya, qlyb, &
       qvya, qvyb, qlTb, qvTa, qvTb, deltaS, dTsoil, tmp2d1, tmp2d2, &
       cv0, &
       cisoliqice_snow, &
       dthetaldT, thetal, isave, nsteps_ice, imelt, vtop, vbot, v_aquifer, &
       dwcol, dwdrainage, drn,inlit, dwinlit, drexcol, dwdischarge, &
       dJcol_latent_S, dJcol_latent_T, dJcol_sensible, deltaJ_latent_S, &
       deltaJ_latent_T, deltaJ_sensible_S, deltaJ_sensible_T, qevapsig, &
       qrunoff, tmp1d1, tmp1d2, tmp1d3,  tmp1d4, deltah0, deltaSL, &
       lE0, G0, &
       Epot, Tfreezing, dtdT, LHS, RHS, LHS_h, RHS_h, surface_case, nns, &
       iflux, litter, i, j, k, kk, condition, littercase, isotopologue, &
       advection, c2, theta, dTqwdTa, dTqwdTb, Tqw, keff, cp, &
       cpeff, hice, h0_0, hice_0, h0_tmp, hice_tmp, qmelt, &
       qtransfer, delta_snowcol, delta_snowT, &
       delta_snowliq, thetai_0, J0, tmp1, tmp2, iqex, &
       nfac1, nfac2, nfac3, nfac4, nfac5, nfac6, nfac7, nfac8, nfac9, &
       nfac10, nfac11, nfac12, J0snow, wcol0snow, h_ex, wpi, err &
       )
    IMPLICIT NONE
    REAL(r_2)                                              :: tfin
    INTEGER(i_d)                                           :: irec, mp
    REAL(r_2),      DIMENSION(1:mp)                        :: qprec
    REAL(r_2),      DIMENSION(1:mp)                        :: qprec_snow
    INTEGER(i_d)                                           :: n
    REAL(r_2),      DIMENSION(1:n) :: dx
    REAL(r_2),      DIMENSION(1:mp)                        :: h0
    REAL(r_2),      DIMENSION(1:n) :: S
    REAL(r_2),      DIMENSION(1:n) :: thetai
    REAL(r_2),      DIMENSION(1:n) :: Jsensible
    REAL(r_2),      DIMENSION(1:n) :: Tsoil
    REAL(r_2),      DIMENSION(1:mp)                        :: evap
    REAL(r_2),      DIMENSION(1:mp)                        :: infil
    REAL(r_2),      DIMENSION(1:mp)                        :: drainage, discharge
    REAL(r_2),      DIMENSION(1:mp,-nsnow_max:n)           :: qh
    INTEGER(i_d),   DIMENSION(1:mp)                        :: nsteps
    TYPE(vars_met), DIMENSION(1:mp)                        :: vmet
    TYPE(vars),     DIMENSION(1:mp)                        :: vlit
    TYPE(vars_snow), DIMENSION(1:mp)                       :: vsnow
    TYPE(vars),     DIMENSION(1:n) :: var
    REAL(r_2),      DIMENSION(1:mp)                        :: T0, Tsurface
    REAL(r_2),      DIMENSION(1:mp)                        :: Hcum, lEcum
    REAL(r_2),      DIMENSION(1:mp)                        :: Gcum, Qadvcum
    REAL(r_2),      DIMENSION(1:mp)                        :: Jcol_sensible, Jcol_latent_S, Jcol_latent_T
    REAL(r_2),      DIMENSION(1:n) :: csoil, kth
    REAL(r_2),      DIMENSION(1:n) :: phi
    REAL(r_2),      DIMENSION(1:mp)                        :: dxL
    REAL(r_2),      DIMENSION(1:mp)                        :: zdelta, SL, Tl
    TYPE(params),   DIMENSION(1:mp)                        :: plit
    TYPE(params),   DIMENSION(1:n) :: par
    REAL(r_2),      DIMENSION(1:mp,1:n), OPTIONAL                    :: wex
    REAL(r_2),      DIMENSION(1:mp,1:nsnow_max)            :: ciso_snow
    REAL(r_2),      DIMENSION(1:mp,1:nsnow_max)            :: cisoice_snow
    REAL(r_2),      DIMENSION(1:mp), OPTIONAL                        :: qali
    REAL(r_2),      DIMENSION(1:mp,-nsnow_max:n)           :: qvsig, qlsig, qvTsig, qvh
    REAL(r_2),      DIMENSION(1:mp)                        :: deltaTa
    REAL(r_2),    DIMENSION(1:mp)                          :: precip, qevap
    REAL(r_2),    DIMENSION(1:mp)                          :: qL, qhL, qybL, qTbL, qhTbL, qhybL, rexcol, wcol
    LOGICAL,      DIMENSION(1:mp)                          :: again, getq0,getqn,init
    LOGICAL,      DIMENSION(1:mp,1:n)                      :: again_ice
    INTEGER(i_d), DIMENSION(1:mp)                          :: ih0, iok, itmp, ns, nsat, nsatlast, nsteps0
    REAL(r_2),    DIMENSION(1:mp)                          :: accel, dmax, dt, dwinfil, dwoff, fac, phip
    REAL(r_2),    DIMENSION(1:mp)                          :: qpme, rsig, rsigdt, sig, t
    REAL(r_2),    DIMENSION(1:n) :: hint, phimin, qexd
    REAL(r_2),    DIMENSION(-nsnow_max+1:n) :: aa, bb, cc, dd, ee, ff, gg, dy
    REAL(r_2),    DIMENSION(-nsnow_max+1:n) :: aah, bbh, cch, ddh, eeh, ffh, ggh, de
    REAL(r_2),    DIMENSION(-nsnow_max:n) :: q, qya, qyb, qTa, qTb,qhya, qhyb, qhTa, qhTb
    REAL(r_2),    DIMENSION(-nsnow_max:n) :: qadv, qadvya, qadvyb, qadvTa, qadvTb
    TYPE(vars)                                             :: vtmp
    REAL(r_2),    DIMENSION(-nsnow_max:n) :: qsig, qhsig, qadvsig
    REAL(r_2),    DIMENSION(-nsnow_max:n) :: qliq, qv, qvT, qlya, qlyb, qvya, qvyb, qlTb, qvTa, qvTb
    REAL(r_2),    DIMENSION(1:n) :: deltaS, dTsoil
    REAL(r_2),    DIMENSION(0:n) :: tmp2d1, tmp2d2
    REAL(r_2),    DIMENSION(1:n) :: cv0
    REAL(r_2),      DIMENSION(1:nsnow_max) :: cisoliqice_snow
    REAL(r_2),    DIMENSION(1:n) :: dthetaldT, thetal
    INTEGER(i_d), DIMENSION(1:n) :: isave, nsteps_ice, imelt
    TYPE(vars),         DIMENSION(1:mp)                    :: vtop, vbot
    TYPE(vars_aquifer), DIMENSION(1:mp)                    :: v_aquifer
    REAL(r_2),          DIMENSION(1:mp)                    :: dwcol, dwdrainage, drn,inlit, dwinlit, drexcol, dwdischarge
    REAL(r_2),          DIMENSION(1:mp)                    :: dJcol_latent_S, dJcol_latent_T, dJcol_sensible
    REAL(r_2),          DIMENSION(1:n) :: deltaJ_latent_S, deltaJ_latent_T, deltaJ_sensible_S, deltaJ_sensible_T
    REAL(r_2),          DIMENSION(1:mp)                    :: qevapsig
    REAL(r_2),          DIMENSION(1:mp)                    :: qrunoff
    REAL(r_2),          DIMENSION(1:mp)                    :: tmp1d1, tmp1d2, tmp1d3,  tmp1d4
    REAL(r_2),          DIMENSION(1:mp)                    :: deltah0
    REAL(r_2),          DIMENSION(1:mp)                    :: deltaSL
    REAL(r_2),          DIMENSION(1:mp)                    :: lE0, G0, Epot
    REAL(r_2),          DIMENSION(1:mp)                    :: Tfreezing
    REAL(r_2),          DIMENSION(1:mp)                    :: dtdT
    REAL(r_2),          DIMENSION(-nsnow_max+1:n) :: LHS, RHS, LHS_h, RHS_h
    INTEGER(i_d),       DIMENSION(1:mp)                    :: surface_case
    INTEGER(i_d),       DIMENSION(1:mp)                    :: nns, iflux
    LOGICAL                                                :: litter
    INTEGER(i_d)                                           :: i, j, k, kk, condition
    INTEGER(i_d)                                           :: littercase, isotopologue, advection ! switches
    REAL(r_2)                                              :: c2, theta
    REAL(r_2)                                              :: dTqwdTa, dTqwdTb, Tqw, keff
    REAL(r_2),          DIMENSION(1:mp)                    :: cp, cpeff, hice, h0_0, hice_0, h0_tmp, hice_tmp
    REAL(r_2),          DIMENSION(nsnow_max)          :: qmelt, hsnow
    REAL(r_2),  DIMENSION(1:mp)                            :: qtransfer
    REAL(r_2),          DIMENSION(nsnow_max) :: delta_snowcol, delta_snowT, delta_snowliq
    REAL(r_2),      DIMENSION(1:n) :: thetai_0, J0
    REAL(r_2)                                              :: tmp1, tmp2
    REAL(r_2),          DIMENSION(1:n) :: iqex
    INTEGER(i_d),       DIMENSION(1:mp)                    :: nfac1, nfac2, nfac3, nfac4, nfac5, &
         nfac6, nfac7, nfac8, nfac9, nfac10, nfac11, nfac12
    REAL(r_2),          DIMENSION(1:mp)                    :: J0snow, wcol0snow
    REAL(r_2), DIMENSION(1:n)                              :: h_ex
    REAL(r_2)                                              :: wpi
    ! Error flag if nstep of SLI > nsteps_max: err=0 -> no error; err/=0 -> error
    INTEGER(i_d),        INTENT(INOUT), OPTIONAL           :: err

    hsnow(:) = zero

    do while (again(kk)) ! sometimes need twice to adjust phi at satn

       nsatlast(kk) = nsat(kk) ! for detecting onset of profile saturation
       nsat(kk)     = sum(var(:)%isat,1) ! no. of sat layers
       sig(kk)      = half
       if (nsat(kk) /= 0) sig(kk) = one ! time weighting sigma
       rsig(kk)     = one/sig(kk)

       ! update variables
       if (iflux(kk)==1) then
          ! Calc flux matric potentials (and derivatives) from S
          ! this set var-structure
          isave(:) = var(:)%isat
          ! Debug for mp=1: remove elemental from hyofS and do loop instead of next line
          call hyofS(S(1:n), Tsoil(1:n), par(1:n), var(1:n))
          !do i=1, n
          !   call hyofS(S(i), Tsoil(i), par(i), var(i))
          !end do

          ! End debug hyofS
          cp(kk) = real(1-var(1)%iice,r_2)*cswat*rhow & ! heat capacity of pond
               + real(var(1)%iice,r_2)*rhow* &
               ((one-var(1)%thetai/par(1)%thre)*cswat + (var(1)%thetai/par(1)%thre)*csice)
          cpeff(kk) = cp(kk) + rhow*lambdaf*var(1)%dthetaldT/par(1)%thre
          ! phip(kk) = max(var(1)%phie-var(1)%he*var(1)%Ksat, (one+e5)*var(1)%phie) !at onset of ponding
          phip(kk) = (one+e5)*var(1)%phie !at onset of ponding
          var(:)%isat  = isave(:)
          thetai(:)    = var(:)%thetai ! ice profile
          thetal(:)    = var(:)%thetal ! liq water profile
          thetai_0(:)  = thetai(:) ! initial ice profile
          dthetaldT(:) = var(:)%dthetaldT
          hice(kk)   = h0(kk)*var(1)%thetai/par(1)%thre
          hice_0(kk) = hice(kk)
          h0_0(kk)   = h0(kk)
          ! sensible heat stored in top layer + pond
          Jsensible(1)   = (var(1)%csoil* dx(1)+h0(kk)*cp(kk))*(Tsoil(1))
          ! sensible heat stored in soil column (2:n)
          Jsensible(2:n) = var(2:n)%csoil*(Tsoil(2:n))* dx(2:n)


          ! calculate this here in case snow pack disappears: ciso of transferred water needs to be retained
          if ((isotopologue /= 0).and.(vsnow(kk)%hsnow(1).gt.zero)) then
             cisoliqice_snow(1:nsnow_max) = (ciso_snow(kk,1:nsnow_max)*vsnow(kk)%hliq(1:nsnow_max) + &
                  cisoice_snow(kk,1:nsnow_max)*(vsnow(kk)%hsnow(1:nsnow_max)-vsnow(kk)%hliq(1:nsnow_max))) / &
                  vsnow(kk)%hsnow(1:nsnow_max)

          endif

          CALL snow_adjust(irec, mp, n, kk, ns, h0, hice, thetai(:), dx(:), vsnow, var, par, S(:), Tsoil(:), &
               Jcol_latent_S, Jcol_latent_T, Jcol_sensible, deltaJ_sensible_S(:), qmelt(:), qtransfer, j0snow)
          thetai(1) = var(1)%thetai  ! this is the value of thetaice prior to matrix call

          call hyofS(S(1), Tsoil(1), par(1), var(1))

       endif ! iflux==1

       ! initialise litter vars
       if (iflux(kk)==1 .and. (littercase==1 .or. littercase == 2)) then
          ! call litter_props(Sl(kk), Tl(kk), vlit(kk), plit(kk), h0(kk))
          call litter_props(Sl(kk), Tsoil(1), vlit(kk), plit(kk), h0(kk))
       endif

       ! ! phi is solution var at satn, so h calc from phi where S>=1 - done in hyofS above for S<1
       ! where (S(:) >= one) & ! (special for frozen soil)
       !      var(:)%h = var(:)%he + (var(:)%phi-var(:)%phie)/var(:)%Ksat


       CALL get_fluxes_and_derivs( &
            mp, qprec, qprec_snow, n, dx(:), h0, &
            Tsoil(:), &
            qh(kk,:), vmet, vlit, vsnow, var, T0, Tsurface, &
            SL, Tl, &
            plit, par, &
            qali, &
            qvh(kk,:), &
            qevap, &
            again, getq0,getqn,init, ns, nsat, &
            nsatlast, &
            phip, qpme, hint(:), phimin(:), &
            qexd(:), &
            q(:), qya(:), qyb(:), qTa(:), qTb(:),qhya(:), qhyb(:), qhTa(:), qhTb(:), qadv(:), qadvya(:), qadvyb(:), &
            qadvTa(:), qadvTb(:), vtmp, qliq(:), qv(:), qvT(:), qlya(:), qlyb(:), &
            qvya(:), qvyb(:), qlTb(:), qvTa(:), qvTb(:), &
            vtop, vbot, &
            lE0, G0, &
            Epot, surface_case, &
            iflux, j, kk, &
            advection, dTqwdTa, dTqwdTb, Tqw, keff, &
            hice &
            )


       CALL estimate_timestep( &
            tfin, mp, n, dx(:), h0, &
            qh(kk,:), &
            nsteps, vlit, vsnow, var, &
            dxL, &
            plit, par, &
            deltaTa, &
            qL, qhL, &
            again, ns, nsat, &
            nsatlast, nsteps0, dmax, dt, &
            qpme, t, &
            q(:), qadv(:), &
            tmp2d1(:), tmp2d2(:), &
            tmp1d1, tmp1d3, &
            iflux, litter, kk, &
            advection, &
            iqex(:) &
            )

       !----- get and solve eqns
       rsigdt(kk) = one/(sig(kk)*dt(kk))
       if (.not. again(kk))  then
          dwoff(kk)   = max(h0(kk)*(one-var(1)%thetai/par(1)%thre)-h0max,zero)
          qrunoff(kk) = dwoff(kk)/dt(kk)
       else
          qrunoff(kk) = zero
       endif

       ! aa, bb, cc and dd hold coeffs and rhs of linear equation eqn set
       CALL get_and_solve_eqn( &
            irec, mp, qprec,  n, dx(:), h0, S(:), thetai(:), &
            Tsoil(:),       &
            qh(kk,:), nsteps,  vlit, vsnow, var,     &
            dxL,    &
            plit, par,        &
            deltaTa, &
            qevap, qL, qhL, qybL, qTbL, qhTbL, qhybL,    &
            again,  again_ice,  iok, itmp, ns,  &
            accel,  dt,   fac,   &
            phip,  rsig, rsigdt, sig, t,      &
            qexd(:), aa(:), bb(:), cc(:), dd(:), ee(:), ff(:), gg(:), dy(:), &
            aah(:), bbh(:), cch(:), ddh(:), eeh(:), ffh(:), ggh(:), &
            de(:), q(:), qya(:), qyb(:), qTa(:), qTb(:),qhya(:), qhyb(:), qhTa(:), qhTb(:), qadv(:), qadvya(:), qadvyb(:), &
            qadvTa(:), qadvTb(:),  qsig(:), qhsig(:), qadvsig(:),      &
            dTsoil(:),   &
            dthetaldT(:),   nsteps_ice, imelt,     &
            deltaJ_latent_T(:), deltaJ_sensible_S(:),   &
            qrunoff, tmp1d1, tmp1d2, tmp1d3,  tmp1d4,     &
            G0, &
            Tfreezing,   LHS(:), RHS(:), LHS_h(:), RHS_h(:),  nns, &
            iflux, litter, i,   kk, condition,   &
            advection,  c2, theta,     cp, &
            cpeff, hice,    h0_tmp,   hsnow(:), &
            delta_snowT(:), &
            delta_snowliq(:),    J0(:),   iqex(:),  &
            nfac1, nfac2, nfac3, nfac4, nfac5, nfac6, nfac7, nfac8, nfac9, &
            nfac10, nfac11, nfac12, err     &
            )
       if (err /= 0) return

       CALL update_unknowns( &
            irec, mp, qprec, qprec_snow, n, dx(:), h0, S(:),  &
            Tsoil(:), evap,   infil, drainage, discharge, &
            qh(kk,:), nsteps, vmet, vlit, vsnow, var,   Hcum, lEcum, &
            Gcum, Qadvcum,  &
            csoil(:), kth(:), phi(:),  zdelta, SL, Tl, &
            par,  wex,      &
            qvsig(kk,:), qlsig(kk,:), qvTsig(kk,:),   &
            precip, qevap,       rexcol, wcol,  &
            again,   ih0,   ns,  &
            dt, dwinfil, dwoff,    &
            sig,       &
            qexd(:),        dy(:),        &
            de(:), q(:), qya(:), qyb(:), qTa(:), qTb(:),qhya(:), qhyb(:), qhTa(:), qhTb(:),    &
            qsig(:), qhsig(:), qadvsig(:), qliq(:), qv(:), qvT(:),  qlyb(:), &
            qvya(:), qvyb(:), qlTb(:), qvTa(:), qvTb(:),   dTsoil(:),    &
            v_aquifer,  &
            dwcol, dwdrainage, drn,inlit, dwinlit, drexcol, dwdischarge, &
            dJcol_latent_S, dJcol_latent_T, dJcol_sensible, deltaJ_latent_S(:), &
            deltaJ_latent_T(:), deltaJ_sensible_S(:), deltaJ_sensible_T(:), qevapsig, &
            qrunoff, tmp1d1, tmp1d2, tmp1d3,   deltah0,  deltaSL,  &
            LHS_h(:),  surface_case,  &
            litter, i, j,  kk,    &
            advection,    theta,   Tqw,  cp, &
            hsnow(:), &
            delta_snowcol(:), delta_snowT(:), &
            delta_snowliq(:),    J0(:),   iqex(:),  &
            J0snow, wcol0snow   &
            )

       CALL update_s_t( &
            mp,   n, dx(:), h0, S(:), thetai(:), &
            Tsoil(:),    infil,   &
            nsteps,    var,     &
            par,        &
            qlsig(kk,:),    &
            again,   ih0,   ns,  &
            dt,      &
            sig,       &
            dy(:),        &
            qhsig(:),       &
            deltaS(:), dTsoil(:),    &
            dJcol_latent_S, dJcol_latent_T, dJcol_sensible, deltaJ_latent_S(:), &
            deltaJ_latent_T(:), deltaJ_sensible_S(:), deltaJ_sensible_T(:),  &
            qrunoff, tmp1d1, tmp1d2, tmp1d3,  tmp1d4, deltah0,    &
            Tfreezing,  dtdT,   LHS_h(:),    &
            i, j, k, kk,    &
            theta,     cp, &
            hice,     hice_tmp,   &
            J0(:), tmp1, tmp2,   &
            h_ex, wpi &
            )

       if (.not. again(kk)) then
          cv0(1:n)   = var(1:n)%cv  ! save cv before updating
          isave(1:n) = var(1:n)%iice ! save isat before updating
          ! update variables (particularly thetaice)
          ! Debug for mp=1: remove elemental from hyofS and do loop instead of next line
          call hyofS(S(1:n), Tsoil(1:n), par(1:n), var(1:n))
          !do i=1, n
          !   call hyofS(S(i), Tsoil(i), par(i), var(i))
          !end do
          ! End debug hyofS
          var(1:n)%iice = isave(1:n)
       endif ! if .not.again

       iflux(kk) = iflux(kk) + 1

    end do ! do while (again(kk)) ! iflux loop
  END SUBROUTINE iflux_loop

  SUBROUTINE update_s_t( &
       mp,   n, dx, h0, S, thetai, &
       Tsoil,    infil,   &
       nsteps,    var,     &
       par,        &
       qlsig,    &
       again,   ih0,   ns,  &
       dt,      &
       sig,       &
       dy,        &
       qhsig,       &
       deltaS, dTsoil,    &
       dJcol_latent_S, dJcol_latent_T, dJcol_sensible, deltaJ_latent_S, &
       deltaJ_latent_T, deltaJ_sensible_S, deltaJ_sensible_T,  &
       qrunoff, tmp1d1, tmp1d2, tmp1d3,  tmp1d4, deltah0,    &
       Tfreezing,  dtdT,   LHS_h,    &
       i, j, k, kk,    &
       theta,     cp, &
       hice,     hice_tmp,   &
       J0, tmp1, tmp2,   &
       h_ex, wpi &
       )
    IMPLICIT NONE
    INTEGER(i_d)                                           :: mp
    INTEGER(i_d)                                           :: n
    REAL(r_2),      DIMENSION(1:n) :: dx
    REAL(r_2),      DIMENSION(1:mp)                        :: h0
    REAL(r_2),      DIMENSION(1:n) :: S
    REAL(r_2),      DIMENSION(1:n) :: thetai
    REAL(r_2),      DIMENSION(1:n) :: Tsoil
    REAL(r_2),      DIMENSION(1:mp)                        :: infil
    INTEGER(i_d),   DIMENSION(1:mp)                        :: nsteps
    TYPE(vars),     DIMENSION(1:n) :: var
    TYPE(params),   DIMENSION(1:n) :: par
    REAL(r_2),      DIMENSION(-nsnow_max:n) :: qlsig
    LOGICAL,      DIMENSION(1:mp)                          :: again
    INTEGER(i_d), DIMENSION(1:mp)                          :: ih0, ns
    REAL(r_2),    DIMENSION(1:mp)                          :: dt
    REAL(r_2),    DIMENSION(1:mp)                          :: sig
    REAL(r_2),    DIMENSION(-nsnow_max+1:n) :: dy
    REAL(r_2),    DIMENSION(-nsnow_max:n) :: qhsig
    REAL(r_2),    DIMENSION(1:n) :: deltaS, dTsoil
    REAL(r_2),          DIMENSION(1:mp)                    :: dJcol_latent_S, dJcol_latent_T, dJcol_sensible
    REAL(r_2),          DIMENSION(1:n) :: deltaJ_latent_S, deltaJ_latent_T, deltaJ_sensible_S, deltaJ_sensible_T
    REAL(r_2),          DIMENSION(1:mp)                    :: qrunoff
    REAL(r_2),          DIMENSION(1:mp)                    :: tmp1d1, tmp1d2, tmp1d3,  tmp1d4
    REAL(r_2),          DIMENSION(1:mp)                    :: deltah0
    REAL(r_2),          DIMENSION(1:mp)                    :: Tfreezing
    REAL(r_2),          DIMENSION(1:mp)                    :: dtdT
    REAL(r_2),          DIMENSION(-nsnow_max+1:n) :: LHS_h
    INTEGER(i_d)                                           :: i, j, k, kk
    REAL(r_2)                                              :: theta
    REAL(r_2),          DIMENSION(1:mp)                    :: cp, hice, hice_tmp
    REAL(r_2),      DIMENSION(1:n) :: J0
    REAL(r_2)                                              :: tmp1, tmp2

    REAL(r_2), DIMENSION(1:n)                              :: h_ex
    REAL(r_2)                                              :: wpi
    ! update variables (S,T) to end of time step
    ! update variables to sig for use in isotope routine
    do i=1, n

       if (.not.again(kk)) Tsoil(i)  = Tsoil(i) + dTsoil(i)
       if (var(i)%isat==0) then
          if (.not.again(kk)) then
             deltaS(i) = dy(i)  !required for isotope subroutine
             S(i)      = S(i)+dy(i)
             if (S(i)>one .and. dy(i)>zero) then ! onset of saturation of layer
                var(i)%isat = 1
                var(i)%K    = var(i)%Ksat
                var(i)%phi  = var(i)%phie
             endif
          endif
       else
          deltaS(i)  = zero    !required for isotope subroutine
          if (i==1) then
             var(i)%phi = var(i)%phi + dy(i)*real(ns(kk),r_2) ! pond included in top soil layer
          else
             var(i)%phi = var(i)%phi + dy(i)
          endif

          ! var(i)%phi = zero from Ross
          ! VH thinks it is o.k. because hyofS is called later again
          ! MC think that we should probably get rid of -he*Ksat in phip etc.
          !   because we merged the pond with the first layer.
          if (i==1 .and. ih0(kk)/=0 .and. var(i)%phi>=var(i)%phie) var(i)%phi = zero ! pond gone
          if (var(i)%phi<var(i)%phie .and. var(i)%iice==0) then ! desaturation of layer
             var(i)%isat = 0
             var(i)%K    = var(i)%Ksat
             var(i)%phi  = var(i)%phie
             var(i)%KS   = par(i)%KSe
             var(i)%phiS = par(i)%phiSe
          elseif (var(i)%phi<var(i)%phie .and. var(i)%iice==1) then
             var(i)%isat = 0
             var(i)%K    = var(i)%Ksat
             var(i)%phi  = var(i)%phie
             var(i)%KS   = var(i)%KS
             var(i)%phiS = var(i)%phiS
          endif
       end if

       if (.not. again(kk)) then
!!$                   if (h0(kk)<zero .and. var(1)%isat==0 .and. (i==1)) then ! start negative pond correction
!!$                      infil(kk)     = infil(kk)+h0(kk)
!!$                      ! negative pond converted to loss of soil moisture in top two layers
!!$                      S(1) = S(1) + h0(kk)/(dx(1)*par(1)%thre)* dx(1)/(dx(1)+dx(2))
!!$                      deltaS(1) = h0(kk)/(dx(1)*par(1)%thre)* dx(1)/(dx(1)+dx(2))
!!$                      ! adjust flux between layer 1 and layer 2
!!$                      qlsig(1) = qlsig(1) + h0(kk)* dx(2)/(dx(1)+dx(2))/dt(kk)
!!$                      deltaS(2) =   h0(kk)/(dx(2)*par(2)%thre)* dx(2)/(dx(1)+dx(2))
!!$                      if (var(2)%isat.eq.0) then
!!$                         dy(2) = dy(2) + deltaS(2)
!!$                      else
!!$                         dy(2) = deltaS(2)
!!$                         var(2)%isat = 0
!!$                      endif
!!$
!!$                      deltah0(kk) = -(h0(kk)-deltah0(kk)) ! whole pond lost (new change = - original pond height)
!!$                      h0(kk) = zero  ! zero pond remaining
!!$                   endif

          ! water in profile initially
          wpi = sum((par(:)%thr + (par(:)%the-par(:)%thr)*S(:))*dx(:))
          if (h0(kk)<zero .and. var(1)%isat==0 .and. (i==1)) then ! start negative pond correction
             wpi = sum((par(:)%thr + (par(:)%the-par(:)%thr)*S(:))*dx(:))

             infil(kk)     = infil(kk)+h0(kk)
             do j=1,n
                h_ex(j) = -h0(kk)* (par(j)%the-par(j)%thr)*S(j)* dx(j)/wpi
             enddo

             do j=1,n-1
                qlsig(j) = qlsig(j) +sum(h_ex(j:n))
             enddo

             do j = 1,n
                if (j<n .and. ( S(j)*(dx(j)*par(j)%thre) + h_ex(j)).lt.zero) then
                   h_ex(j+1) = h_ex(j+1) +(h_ex(j)- S(j)*(dx(j)*par(j)%thre))
                   h_ex(j) = S(j)*(dx(j)*par(j)%thre)
                endif
                deltaS(j)= -min(h_ex(j)/(dx(j)*par(j)%thre),S(j))
                if (var(j)%isat.eq.0) then
                   dy(j) = dy(j) + deltaS(j)
                else
                   dy(j) = deltaS(j)
                   var(j)%isat = 0
                endif
                if (j==1) S(j)=S(j)+deltaS(j)

             enddo
             deltah0(kk) = -(h0(kk)-deltah0(kk)) ! whole pond lost (new change = - original pond height)
             h0(kk) = zero  ! zero pond remaining

          endif


          if (S(i)<zero .and. i<n) then ! start correction for negative soil moisture in any layer
             wpi = sum((par(i+1:n)%thr + (par(i+1:n)%the-par(i+1:n)%thr)*S(i+1:n))*dx(i+1:n))

             do j=i+1,n
                h_ex(j) = -deltaS(i)*dx(i)*par(i)%thre * &
                     (par(j)%thr + (par(j)%the-par(j)%thr)*S(j))*dx(j) /wpi
             enddo

             do j=i+1,n
                qlsig(j) = qlsig(j) +sum(h_ex(j:n))
             enddo

             do j = i+1,n
                if (j<n .and. ( S(j)*(dx(j)*par(j)%thre) - h_ex(j)).lt.zero) then
                   h_ex(j+1) = h_ex(j+1) +(h_ex(j)- S(j)*(dx(j)*par(j)%thre))
                   h_ex(j) = S(j)*(dx(j)*par(j)%thre)
                endif
                deltaS(j)= -min(h_ex(j)/(dx(j)*par(j)%thre),S(j))
                if (var(j)%isat.eq.0) then
                   dy(j) = dy(j) + deltaS(j)
                else
                   dy(j) = deltaS(j)
                   var(j)%isat = 0
                endif
                if (j==i) S(j)=S(j)+deltaS(j)

             enddo

             S(i)=S(i)-deltaS(i)
             deltaS(i) = zero

          endif

          ! corrections for onset of freezing and melting
          ! calculate freezing point temperature at new S
          Tfreezing(kk) = Tfrz(S(i), par(i)%he, one/(par(i)%lambc*freezefac))
          ! check for onset of freezing and adjust (increase) temperature to account for latent heat
          ! release by ice formation

          if (experiment /= 184) then
             if (Tsoil(i) < Tfreezing(kk) .and. var(i)%iice==0) then  ! start correction for onset of freezing
                dtdT(kk)      = dthetalmaxdTh(Tfreezing(kk), S(i), par(i)%he, one/(par(i)%lambc*freezefac), &
                     par(i)%thre,par(i)%the)
                tmp1d3(kk) = 1000._r_2
                tmp1d2(kk) = Tsoil(i) - dTsoil(i)
                k = 1
                if ((i==1).and.(h0(kk)- deltah0(kk))>zero) then
                   h0(kk) = h0(kk) - deltah0(kk)  ! calculate stored heat using h0 at t-dt

                   do while ((k<nsteps_ice_max).and.(tmp1d3(kk)>tol_dthetaldT))
                      tmp1d1(kk)  = LHS_h(i)*dt(kk)-deltaJ_sensible_S(i) + tmp1d2(kk)* &
                           (var(i)%csoil*dx(i)+cp(kk)*h0(kk)) + &
                           Tfreezing(kk)*rhow*(lambdaf+(tmp1d2(kk))*(cswat-csice)) * &
                           dtdT(kk)*(dx(i)+h0(kk)/par(i)%thre)
                      tmp1d1(kk)  = tmp1d1(kk) /( (var(i)%csoil*dx(i)+cp(kk)*h0(kk)) + &
                           rhow*(lambdaf+(tmp1d2(kk))*(cswat-csice)) * &
                           dtdT(kk)*(dx(i)+h0(kk)/par(i)%thre))
                      dtdT(kk) = (thetalmax(tmp1d1(kk), S(i), par(i)%he, one/(par(i)%lambc*freezefac), &
                           par(i)%thre, par(i)%the) &
                           - (S(i)*par(i)%thre + (par(i)%the-par(i)%thre)))/ &
                           (tmp1d1(kk) - Tfreezing(kk))
                      tmp1d3(kk) = abs((tmp1d1(kk)- Tfreezing(kk))*dtdT(kk))
                      k = k + 1
                   enddo  ! end of do while
                   dTsoil(i) = tmp1d1(kk) - tmp1d2(kk)
                   Tsoil(i) = tmp1d1(kk)
                   deltaJ_sensible_T(i) = (var(i)%csoil*dx(i)+ cp(kk)*h0(kk))*dTsoil(i) + &
                        (Tsoil(i)-Tfreezing(kk))*rhow*dtdT(kk)*(dx(1)+h0(kk)/par(1)%thre)* &
                        ((cswat-csice)*(tmp1d2(kk)))
                   deltaJ_latent_T(i) = (Tsoil(i)-Tfreezing(kk))*rhow*dtdT(kk)*lambdaf* &
                        (dx(1)+h0(kk)/par(1)%thre)
                   h0(kk) = h0(kk) + deltah0(kk)

                   ! change in heat stored in soil column
                   dJcol_latent_S(kk) = sum(deltaJ_latent_S(1:n))
                   dJcol_latent_T(kk) = sum(deltaJ_latent_T(1:n))
                   dJcol_sensible(kk) = sum(deltaJ_sensible_T(1:n)) + sum(deltaJ_sensible_S(1:n))
                   ! tmp1d1(kk) = dJcol_latent_S(kk) + dJcol_latent_T(kk) + &
                   !      dJcol_sensible(kk) - (qhsig(0)-qhsig(6))*dt(kk) + &
                   !      qrunoff(kk)*cswat*rhow*(Tsoil(1) + sig(kk)*dTsoil(1))*dt(kk)
                   tmp1d1(kk) = dJcol_latent_S(kk) + dJcol_latent_T(kk) + &
                        dJcol_sensible(kk) - (qhsig(0)-qhsig(n))*dt(kk) + &
                        qrunoff(kk)*cswat*rhow*(Tsoil(1) + sig(kk)*dTsoil(1))*dt(kk)
                   ! if (abs(tmp1d1(kk)).gt.10.) then
                   !    write(*,*) 'E imbalance check 2'
                   ! endif
                else
                   do while ((k<nsteps_ice_max).and.(tmp1d3(kk)>tol_dthetaldT))
                      ! zeroth estimate of corrected temperature
                      tmp1d1(kk)  = LHS_h(i)*dt(kk)-deltaJ_sensible_S(i) + tmp1d2(kk)* &
                           (var(i)%csoil*dx(i)) + Tfreezing(kk)*rhow*(lambdaf+(tmp1d2(kk))*(cswat-csice)) * &
                           dtdT(kk)*dx(i)
                      tmp1d1(kk)  = tmp1d1(kk) /( (var(i)%csoil*dx(i)) + &
                           rhow*(lambdaf+(tmp1d2(kk))*(cswat-csice)) * dtdT(kk)*dx(i))
                      dtdT(kk) =  (thetalmax(tmp1d1(kk), S(i), par(i)%he, &
                           one/(par(i)%lambc*freezefac), par(i)%thre, par(i)%the) &
                           - (S(i)*par(i)%thre + (par(i)%the-par(i)%thre)))/ &
                           (tmp1d1(kk) - Tfreezing(kk))
                      tmp1d3(kk) = abs((tmp1d1(kk)- Tfreezing(kk))*dtdT(kk))
                      k = k + 1
                   enddo  ! end of do while
                   dTsoil(i) = tmp1d1(kk)- (Tsoil(i)-dTsoil(i))
                   Tsoil(i) = tmp1d1(kk)
                   deltaJ_sensible_T(i) = (var(i)%csoil*dx(i))*dTsoil(i) + &
                        (Tsoil(i)-Tfreezing(kk))*rhow*dtdT(kk)*dx(i)* &
                        ((cswat-csice)*(tmp1d2(kk)))
                   deltaJ_latent_T(i) = (Tsoil(i)-Tfreezing(kk))*rhow*dtdT(kk)*lambdaf*dx(i)
                endif
             endif  ! end correction for onset of freezing
          endif

          if (var(i)%iice==1) then
             theta         = S(i)*(par(i)%thre) + (par(i)%the - par(i)%thre)
             tmp1d1(kk) = Tsoil(i) - dTsoil(i)
             ! energy for complete melting
             if ((i==1).and.(h0(kk)>zero)) then
                tmp1d3(kk) = thetai(i)*dx(i)*rhow*lambdaf + &
                     thetai(i)*(h0(kk)- deltah0(kk))/par(i)%thre*rhow*lambdaf

                tmp1d2(kk) = (Tfreezing(kk))*(rhow*cswat*(theta*dx(1)+h0(kk))+ &
                     par(1)%rho*par(1)%css*dx(1)) - &
                     (tmp1d1(kk))*(rhow*cswat*(var(i)%thetal*dx(1)+(h0(kk)- deltah0(kk))* &
                     (one-thetai(1)/par(1)%thre))+par(1)%rho*par(1)%css*dx(1)) - &
                     (tmp1d1(kk))*(rhow*csice*(thetai(i)*dx(1)+(h0(kk)- deltah0(kk))* &
                     thetai(1)/par(1)%thre))
             else
                tmp1d3(kk) = thetai(i)*dx(i)*rhow*lambdaf

                tmp1d2(kk) = (Tfreezing(kk))*(rhow*cswat*(theta*dx(i))+par(i)%rho*par(i)%css*dx(i)) - &
                     (tmp1d1(kk))*(rhow*cswat*(var(i)%thetal*dx(i))+par(i)%rho*par(i)%css*dx(i)) - &
                     (tmp1d1(kk))*(rhow*csice*(thetai(i)*dx(i)))
             endif
          endif
          ! check for onset of thawing and decrease temperature to account for latent heat required to melt ice

          if ((var(i)%iice==1).and.((J0(i) + LHS_h(i)*dt(kk)).ge.JSoilLayer(Tfreezing(kk), &
               dx(i), theta,par(i)%css, par(i)%rho, &
               merge(h0(kk),zero,i==1), par(i)%thre, par(i)%the, &
               par(i)%he, one/(par(i)%lambc*freezefac)))) then
             ! start correction for onset of thawing
             ! Correct for "overthawing". This extra energy is used to heat the soil even more
             ! The correction comes from equating the energy balances before and after correction.
             tmp1d1(kk) = Tsoil(i) - dTsoil(i)

             if ((i==1) .and. (h0(kk)>zero)) then
                deltaJ_latent_T(i)=      tmp1d3(kk)
                deltaJ_latent_S(i) = zero
                h0(kk) = h0(kk) - deltah0(kk)

                deltaJ_sensible_S(i) = (tmp1d1(kk))*(rhow*cswat*(theta*dx(1)+(h0(kk)+deltah0(kk)))+ &
                     par(1)%rho*par(1)%css*dx(1)) - &
                     (tmp1d1(kk))*(rhow*cswat*(var(i)%thetal*dx(1)+(h0(kk))* &
                     (one-thetai(1)/par(1)%thre))+par(1)%rho*par(1)%css*dx(1)) - &
                     (tmp1d1(kk))*(rhow*csice*(thetai(i)*dx(1)+(h0(kk))* &
                     thetai(1)/par(1)%thre))

                dTsoil(i) = (LHS_h(i)*dt(kk) - (deltaJ_latent_T(i) + deltaJ_sensible_S(i)))/ &
                     (rhow*cswat*(theta*dx(1)+(h0(kk)+deltah0(kk)))+par(1)%rho*par(1)%css*dx(1))

                deltaJ_sensible_T(i) = dTsoil(i)* &
                     (rhow*cswat*(theta*dx(1)+(h0(kk)+deltah0(kk)))+par(1)%rho*par(1)%css*dx(1))

                h0(kk) = h0(kk) + deltah0(kk)
             else
                deltaJ_latent_T(i) = thetai(i)*dx(i)*rhow*lambdaf
                deltaJ_latent_S(i) = zero

                theta         = S(i)*(par(i)%thre) + (par(i)%the - par(i)%thre)
                var(i)%csoil = cswat*theta*rhow + par(i)%css*par(i)%rho

                deltaJ_sensible_S(i) = (tmp1d1(kk))*(rhow*cswat*(theta*dx(i))+ &
                     par(i)%rho*par(i)%css*dx(i)) - &
                     (tmp1d1(kk))*(rhow*cswat*(var(i)%thetal*dx(i))+par(i)%rho*par(i)%css*dx(i)) - &
                     (tmp1d1(kk))*(rhow*csice*(thetai(i)*dx(i)))
                dTsoil(i) = (LHS_h(i)*dt(kk) - (deltaJ_latent_T(i) + deltaJ_sensible_S(i)))/ &
                     (dx(i)*var(i)%csoil)
                deltaJ_sensible_T(i) = var(i)%csoil*dTsoil(i)*dx(i)
             endif
             Tsoil(i) = tmp1d1(kk) + dTsoil(i)
             ! correct thetai and T in frozen soil
          elseif ((var(i)%iice==1).and.((J0(i) + LHS_h(i)*dt(kk))<JSoilLayer(Tfreezing(kk), &
               dx(i), theta,par(i)%css, par(i)%rho, &
               merge(h0(kk),zero,i==1), par(i)%thre, par(i)%the, &
               par(i)%he, one/(par(i)%lambc*freezefac)))) then

             tmp1d2(kk) = J0(i) + LHS_h(i)*dt(kk) ! total energy in  soil layer
             if (Tsoil(i) .lt. -40.0_r_2) then
                ! assume no lquid below min temp threshhold
                var(i)%thetal = 0.0_r_2
                var(i)%thetai = theta
                if (i.eq.1) hice(kk) = h0(kk)
                tmp1d3(kk) = (tmp1d2(kk) + rhow*lambdaf*(theta*dx(i) +  merge(h0(kk),zero,i==1)))/ &
                     (dx(i)*par(i)%css*par(i)%rho + rhow*csice*(theta*dx(i) + &
                     merge(h0(kk),zero,i==1)))
             else
                !check there is a zero
                tmp1 = GTfrozen(Tsoil(i)-dTsoil(i)-50._r_2, tmp1d2(kk), dx(i), theta,par(i)%css, par(i)%rho, &
                     merge(h0(kk),zero,i==1), par(i)%thre, par(i)%the, par(i)%he, one/(par(i)%lambc*freezefac))

                tmp2 = GTFrozen(Tfreezing(kk), tmp1d2(kk), dx(i), theta,par(i)%css, par(i)%rho, &
                     merge(h0(kk),zero,i==1), par(i)%thre, par(i)%the, par(i)%he, one/(par(i)%lambc*freezefac))
                ! there is a zero in between
                if ((tmp1*tmp2) < zero) then
                   tmp1d3(kk) = rtbis_Tfrozen(tmp1d2(kk), dx(i), theta,par(i)%css, par(i)%rho, &
                        merge(h0(kk),zero,i==1), par(i)%thre, par(i)%the, &
                        par(i)%he, one/(par(i)%lambc*freezefac), &
                        Tsoil(i)-dTsoil(i)-50._r_2, Tfreezing(kk))
                   tmp1d4(kk) = thetalmax(tmp1d3(kk), S(i), par(i)%he, one/(par(i)%lambc*freezefac), &
                        par(i)%thre, par(i)%the) ! liquid content at solution for Tsoil
                else
                   write(wlogn,*) "Found no solution for Tfrozen 1. ", kk, i
                   write(wlogn,*) "Assume soil is totally frozen"
                   var(i)%thetal = 0.0_r_2
                   var(i)%thetai = theta
                   if (i.eq.1) hice(kk) = h0(kk)
                   tmp1d3(kk) = (tmp1d2(kk) + rhow*lambdaf*(theta*dx(i) +  merge(h0(kk),zero,i==1)))/ &
                        (dx(i)*par(i)%css*par(i)%rho + rhow*csice*(theta*dx(i) + &
                        merge(h0(kk),zero,i==1)))
                   write(wlogn,*) "frozen soil temperature: ", tmp1d3(kk)
                   write(wlogn,*) nsteps(kk), S(i), Tsoil(i), dTsoil(i), h0(kk), tmp1, tmp2, tmp1d2(kk), theta, &
                        JSoilLayer(Tfreezing(kk), &
                        dx(i), theta,par(i)%css, par(i)%rho, &
                        merge(h0(kk),zero,i==1), par(i)%thre, par(i)%the, &
                        par(i)%he, one/(par(i)%lambc*freezefac)), J0(i) + LHS_h(i)*dt(kk), Tfreezing(kk)

                endif
                var(i)%thetal = tmp1d4(kk)
                var(i)%thetai = theta - tmp1d4(kk)
             endif ! Tsoil <-40
             if (i==1) then
                hice_tmp(kk) = hice(kk)
                hice(kk) = h0(kk)*var(1)%thetai/par(1)%thre

                ! correct total energy stored in pond + soil
                deltaJ_latent_S(i) = - rhow*lambdaf*((hice(kk)-hice_tmp(kk)) + &
                     dx(1)*(var(1)%thetai-thetai(1)))
             else
                ! correct total energy stored in  soil
                deltaJ_latent_S(i) = - rhow*lambdaf*( + &
                     dx(i)*(var(i)%thetai-thetai(i)))


             endif

             deltaJ_latent_T(i) = zero

             deltaJ_sensible_T(i) = JSoilLayer(tmp1d3(kk), &
                  dx(i), theta,par(i)%css, par(i)%rho, &
                  merge(h0(kk),zero,i==1), par(i)%thre, par(i)%the, &
                  par(i)%he, one/(par(i)%lambc*freezefac)) - &
                  J0(i)- &
                  deltaJ_latent_S(i)

             deltaJ_sensible_S(i) = 0._r_2

             thetai(i) = var(i)%thetai
             Tsoil(i) = tmp1d3(kk)

          endif ! soil remains frozen

       endif ! if .not.again

    end do ! i=1, n => update variables S,T
  END SUBROUTINE update_s_t


  SUBROUTINE update_unknowns( &
       irec, mp, qprec, qprec_snow, n, dx, h0, S,  &
       Tsoil, evap,   infil, drainage, discharge, &
       qh, nsteps, vmet, vlit, vsnow, var,   Hcum, lEcum, &
       Gcum, Qadvcum,  &
       csoil, kth, phi,  zdelta, SL, Tl, &
       par,  wex,      &
       qvsig, qlsig, qvTsig,   &
       precip, qevap,       rexcol, wcol,  &
       again,   ih0,   ns,  &
       dt, dwinfil, dwoff,    &
       sig,       &
       qexd,        dy,        &
       de, q, qya, qyb, qTa, qTb,qhya, qhyb, qhTa, qhTb,    &
       qsig, qhsig, qadvsig, qliq, qv, qvT,  qlyb, &
       qvya, qvyb, qlTb, qvTa, qvTb,   dTsoil,    &
       v_aquifer,  &
       dwcol, dwdrainage, drn,inlit, dwinlit, drexcol, dwdischarge, &
       dJcol_latent_S, dJcol_latent_T, dJcol_sensible, deltaJ_latent_S, &
       deltaJ_latent_T, deltaJ_sensible_S, deltaJ_sensible_T, qevapsig, &
       qrunoff, tmp1d1, tmp1d2, tmp1d3,   deltah0,  deltaSL,  &
       LHS_h,  surface_case,  &
       litter, i, j,  kk,    &
       advection,    theta,   Tqw,  cp, &
       hsnow, &
       delta_snowcol, delta_snowT, &
       delta_snowliq,    J0,   iqex,  &
       J0snow, wcol0snow   &
       )
    IMPLICIT NONE
    INTEGER(i_d)                                           :: irec, mp
    REAL(r_2),      DIMENSION(1:mp)                        :: qprec
    REAL(r_2),      DIMENSION(1:mp)                        :: qprec_snow
    INTEGER(i_d)                                           :: n
    REAL(r_2),      DIMENSION(1:n) :: dx
    REAL(r_2),      DIMENSION(1:mp)                        :: h0
    REAL(r_2),      DIMENSION(1:n) :: S
    REAL(r_2),      DIMENSION(1:n) :: Tsoil
    REAL(r_2),      DIMENSION(1:mp)                        :: evap
    REAL(r_2),      DIMENSION(1:mp)                        :: infil
    REAL(r_2),      DIMENSION(1:mp)                        :: drainage, discharge
    REAL(r_2),      DIMENSION(-nsnow_max:n) :: qh
    INTEGER(i_d),   DIMENSION(1:mp)                        :: nsteps
    TYPE(vars_met), DIMENSION(1:mp)                        :: vmet
    TYPE(vars),     DIMENSION(1:mp)                        :: vlit
    TYPE(vars_snow), DIMENSION(1:mp)                       :: vsnow
    TYPE(vars),     DIMENSION(1:n) :: var
    REAL(r_2),      DIMENSION(1:mp)                        :: Hcum, lEcum
    REAL(r_2),      DIMENSION(1:mp)                        :: Gcum, Qadvcum
    REAL(r_2),      DIMENSION(1:n) :: csoil, kth
    REAL(r_2),      DIMENSION(1:n) :: phi
    REAL(r_2),      DIMENSION(1:mp)                        :: zdelta, SL, Tl
    TYPE(params),   DIMENSION(1:n) :: par
    REAL(r_2),      DIMENSION(1:mp,1:n), OPTIONAL                    :: wex
    REAL(r_2),      DIMENSION(-nsnow_max:n) :: qvsig, qlsig, qvTsig
    REAL(r_2),    DIMENSION(1:mp)                          :: precip, qevap
    REAL(r_2),    DIMENSION(1:mp)                          :: rexcol, wcol
    LOGICAL,      DIMENSION(1:mp)                          :: again
    INTEGER(i_d), DIMENSION(1:mp)                          :: ih0, ns
    REAL(r_2),    DIMENSION(1:mp)                          :: dt, dwinfil, dwoff
    REAL(r_2),    DIMENSION(1:mp)                          :: sig
    REAL(r_2),    DIMENSION(1:n) :: qexd
    REAL(r_2),    DIMENSION(-nsnow_max+1:n) :: dy
    REAL(r_2),    DIMENSION(-nsnow_max+1:n) :: de
    REAL(r_2),    DIMENSION(-nsnow_max:n) :: q, qya, qyb, qTa, qTb,qhya, qhyb, qhTa, qhTb
    REAL(r_2),    DIMENSION(-nsnow_max:n) :: qsig, qhsig, qadvsig
    REAL(r_2),    DIMENSION(-nsnow_max:n) :: qliq, qv, qvT, qlyb, qvya, qvyb, qlTb, qvTa, qvTb
    REAL(r_2),    DIMENSION(1:n) :: dTsoil
    TYPE(vars_aquifer), DIMENSION(1:mp)                    :: v_aquifer
    REAL(r_2),          DIMENSION(1:mp)                    :: dwcol, dwdrainage, drn,inlit, dwinlit, drexcol, dwdischarge
    REAL(r_2),          DIMENSION(1:mp)                    :: dJcol_latent_S, dJcol_latent_T, dJcol_sensible
    REAL(r_2),          DIMENSION(1:n) :: deltaJ_latent_S, deltaJ_latent_T, deltaJ_sensible_S, deltaJ_sensible_T
    REAL(r_2),          DIMENSION(1:mp)                    :: qevapsig
    REAL(r_2),          DIMENSION(1:mp)                    :: qrunoff
    REAL(r_2),          DIMENSION(1:mp)                    :: tmp1d1, tmp1d2, tmp1d3
    REAL(r_2),          DIMENSION(1:mp)                    :: deltah0
    REAL(r_2),          DIMENSION(1:mp)                    :: deltaSL
    REAL(r_2),          DIMENSION(-nsnow_max+1:n) :: LHS_h
    INTEGER(i_d),       DIMENSION(1:mp)                    :: surface_case
    LOGICAL                                                :: litter
    INTEGER(i_d)                                           :: i, j, kk
    INTEGER(i_d)                                           :: advection ! switches
    REAL(r_2)                                              :: theta
    REAL(r_2)                                              :: Tqw
    REAL(r_2),          DIMENSION(1:mp)                    :: cp
    REAL(r_2),          DIMENSION(nsnow_max) :: hsnow
    REAL(r_2),          DIMENSION(nsnow_max) :: delta_snowcol, delta_snowT, delta_snowliq
    REAL(r_2),      DIMENSION(1:n) :: J0
    REAL(r_2),          DIMENSION(1:n) :: iqex

    REAL(r_2),          DIMENSION(1:mp)                    :: J0snow, wcol0snow
    !----- update unknowns
    ! i.e. update state variables to end of time step
    !      cumulate some surface fluxes
    ih0(kk) = 0
    if (.not. again(kk)) then
       ! get surface fluxes at sigma of the time step
       dwoff(kk)  = zero
       deltah0(kk) = zero
       precip(kk) = precip(kk) + (qprec(kk)+qprec_snow(kk))*dt(kk)
       dwcol(kk) = sum(par(1:n)%thre*dx(1:n)*dy(1:n)*(abs(var(1:n)%isat-1)),1) + &
            real(1-ns(kk),r_2)*dy(1)

       if (vsnow(kk)%nsnow==1) dwcol(kk) = dwcol(kk) + dy(0)

       ! absolute heat stored in soil column before update
       do j=1, n
          theta  =  S(j)*(par(j)%thre) + (par(j)%the - par(j)%thre)
          J0(j) = JSoilLayer(Tsoil(j), &
               dx(j), theta,par(j)%css, par(j)%rho, &
               merge(h0(kk),zero,j.eq.1), par(j)%thre, par(j)%the, &
               par(j)%he, one/(par(j)%lambc*freezefac))
       end do

       ! change in heat stored in soil column
       do j=1, n
          if (j==1) then
             ! pond included in top soil layer
             deltaJ_latent_S(j) = -dx(j)*dy(j)*par(j)%thre*real(var(j)%iice,r_2)* &
                  real(1-var(j)%isat,r_2)*rhow*lambdaf - &
                  real(1-ns(kk),r_2)*dy(1)*real(var(j)%iice,r_2)*rhow*lambdaf*(var(1)%thetai/par(1)%thre)

             deltaJ_latent_T(j) = var(j)%dthetaldT*dx(j)*dTsoil(j)*rhow*lambdaf* &
                  real(var(j)%iice,r_2) + &
                  h0(kk)/par(1)%thre*var(j)%dthetaldT*dTsoil(j)*rhow*lambdaf*real(var(j)%iice,r_2)

             deltaJ_sensible_T(j) = var(j)%csoil*dx(j)*dTsoil(j) + &
                  cp(kk)*h0(kk)*dTsoil(j) + &
                  var(1)%iice*var(j)%dthetaldT*(Tsoil(j))*rhow*(cswat-csice)*dx(j)*dTsoil(j) + &
                  var(j)%iice*var(j)%dthetaldT*(Tsoil(j))*rhow*(cswat-csice)* &
                  h0(kk)/par(j)%thre*dTsoil(j)

             if (advection==1) then
                deltaJ_sensible_S(j) = rhow*cswat*(Tsoil(j))*dx(j)*par(j)%thre*dy(j)* &
                     real(1-var(j)%isat,r_2)*real(1-var(j)%iice,r_2)  + &
                     rhow*csice*(Tsoil(j))*dx(j)*par(j)%thre*dy(j)* &
                     real(1-var(j)%isat,r_2)*real(var(j)%iice,r_2) + &
                     real(1-ns(kk),r_2)*real(1-var(1)%iice,r_2)* &
                     dy(1)*rhow*cswat*(Tsoil(1))  + &
                     real(1-ns(kk),r_2)*real(var(1)%iice,r_2)* &
                     dy(1)*rhow*csice*(Tsoil(1))*(var(1)%thetai/par(1)%thre) + &
                     real(1-ns(kk),r_2)*real(var(1)%iice,r_2)* &
                     dy(1)*rhow*cswat*(Tsoil(1))*(one-(var(1)%thetai/par(1)%thre))
             endif
          else ! (j>1)
             deltaJ_latent_S(j) = -dx(j)*dy(j)*par(j)%thre*real(var(j)%iice,r_2)* &
                  real(1-var(j)%isat,r_2)*rhow*lambdaf
             deltaJ_sensible_T(j) = var(j)%csoil*dx(j)*dTsoil(j) + &
                  var(j)%iice*var(j)%dthetaldT*(Tsoil(j))*rhow*(cswat-csice)*dx(j)*dTsoil(j)


             if (advection==1) then
                deltaJ_sensible_S(j) = rhow*cswat*(Tsoil(j))*dx(j)*dy(j)*par(j)%thre* &
                     real(1-var(j)%isat,r_2)*real(1-var(j)%iice,r_2) &
                     + rhow*csice*(Tsoil(j))*dx(j)*dy(j)*par(j)%thre &
                     *real(1-var(j)%isat,r_2)*real(var(j)%iice,r_2)
             endif

             deltaJ_latent_T(j) = var(j)%dthetaldT*dx(j)*dTsoil(j)*rhow*lambdaf*real(var(j)%iice,r_2)
          endif
       end do

       ! change in heat stored in soil column
       dJcol_latent_S(kk) = sum(deltaJ_latent_S(1:n))
       dJcol_latent_T(kk) = sum(deltaJ_latent_T(1:n))
       dJcol_sensible(kk) = sum(deltaJ_sensible_T(1:n)) + sum(deltaJ_sensible_S(1:n))
       ! tmp1d1(kk) = dJcol_latent_S(kk) + dJcol_latent_T(kk) + dJcol_sensible(kk) - (qhsig(0)-qhsig(6))*dt(kk) &
       !      + qrunoff(kk)*cswat*rhow*(Tsoil(1) + sig(kk)*dTsoil(1))*dt(kk)
       tmp1d1(kk) = dJcol_latent_S(kk) + dJcol_latent_T(kk) + dJcol_sensible(kk) - (qhsig(0)-qhsig(n))*dt(kk) &
            + qrunoff(kk)*cswat*rhow*(Tsoil(1) + sig(kk)*dTsoil(1))*dt(kk)

       if (ns(kk)<1) then ! change in pond height
          h0(kk) = h0(kk) + dy(1)
          deltah0(kk) = dy(1)
          if (h0(kk)<zero .and. dy(1)<zero) then
             ih0(kk)=1 ! pond gone
          endif
       endif

       ! cumulate evaporation from top of soil column or top of litter/pond or top of snow pack
       select case (surface_case(kk))
       case(1)  ! no snow
          evap(kk)     = evap(kk) +qevap(kk)*dt(kk) - sig(kk)*(qyb(0)*dy(1)+qTb(0)*dTsoil(1))*dt(kk)
          qevapsig(kk) = qevap(kk) - sig(kk)*(qyb(0)*dy(1)+qTb(0)*dTsoil(1))
          Gcum(kk) = Gcum(kk)+(qhsig(0)-qadvsig(0))*dt(kk)
          Qadvcum(kk) = Qadvcum(kk) + qadvsig(0)*dt(kk) - qadvsig(n)*dt(kk)
          lEcum(kk)    = lEcum(kk) + qevapsig(kk)*thousand*var(1)%lambdav*dt(kk)
          Hcum(kk) = Hcum(kk) + (vmet(kk)%Rn*dt(kk)-(qhsig(0)-qadvsig(0))*dt(kk)- &
               qevapsig(kk)*thousand*var(1)%lambdav*dt(kk))
          dwinfil(kk) = (qprec(kk)+qprec_snow(kk)-qevapsig(kk))*dt(kk)
       case(2) ! dedicated snow layer
          evap(kk)     = evap(kk) +qevap(kk)*dt(kk) - sig(kk)*(qyb(-vsnow(kk)%nsnow)*dy(1-vsnow(kk)%nsnow)+ &
               qTb(-vsnow(kk)%nsnow)*de(1-vsnow(kk)%nsnow))*dt(kk)
          qevapsig(kk) = qevap(kk) - sig(kk)*(qyb(-vsnow(kk)%nsnow)*dy(1-vsnow(kk)%nsnow)+ &
               qTb(-vsnow(kk)%nsnow)*de(1-vsnow(kk)%nsnow))
          Gcum(kk) = Gcum(kk)+(qhsig(-vsnow(kk)%nsnow)-qadvsig(-vsnow(kk)%nsnow))*dt(kk)
          Qadvcum(kk) = Qadvcum(kk) + qadvsig(-vsnow(kk)%nsnow)*dt(kk) - qadvsig(n)*dt(kk)
          if (vsnow(kk)%hliq(1).gt.zero) then
             lEcum(kk)    = lEcum(kk) + qevapsig(kk)*thousand*lambdaf*dt(kk)
             Hcum(kk) = Hcum(kk) + (vmet(kk)%Rn*dt(kk)-(qhsig(-vsnow(kk)%nsnow)-qadvsig(-vsnow(kk)%nsnow))*dt(kk)- &
                  qevapsig(kk)*thousand*lambdaf*dt(kk))
          else
             lEcum(kk)    = lEcum(kk) + qevapsig(kk)*thousand*lambdas*dt(kk)
             Hcum(kk) = Hcum(kk) + (vmet(kk)%Rn*dt(kk)-(qhsig(-vsnow(kk)%nsnow)-qadvsig(-vsnow(kk)%nsnow))*dt(kk)- &
                  qevapsig(kk)*thousand*lambdas*dt(kk))
          endif
          dwinfil(kk) = (qprec(kk)+qprec_snow(kk)-qevapsig(kk))*dt(kk)
       end select ! surface_case

       dwdrainage(kk)     = q(n)*dt(kk) +sig(kk)*dt(kk)*(qya(n)*dy(n)+qTa(n)*dTsoil(n))
       if (botbc=="aquifer" .and. v_aquifer(kk)%isat==0) then
          dwdischarge(kk) = v_aquifer(kk)%discharge*dt(kk)
       else
          dwdischarge(kk) = dwdrainage(kk)
       endif

       drexcol(kk) = sum(iqex(:),1)*dt(kk)
       ! cumulative fluxes
       infil(kk)     = infil(kk) + dwinfil(kk)
       inlit(kk)     = inlit(kk) + dwinlit(kk)
       wcol(kk)      = wcol(kk) + dwcol(kk)
       drainage(kk)  = drainage(kk) + dwdrainage(kk)
       discharge(kk) = discharge(kk) + qrunoff(kk)*dt(kk)
       rexcol(kk)    = rexcol(kk) + drexcol(kk)
       Qadvcum(kk) = Qadvcum(kk) - sum(iqex(:)*(Tsoil(:)),1)*dt(kk)*rhow*cswat - &
            qrunoff(kk)*(Tsoil(1))*dt(kk)*rhow*cswat
       ! if (irec.eq.109) then
       !     write(*,"(1i8,16e16.6)") irec, infil(kk)-(wcol(kk)+discharge(kk)+drainage(kk))-rexcol(kk), &
       !         wcol(kk), deltah0(kk), h0(kk), discharge(kk), evap(kk), &
       !         drainage(kk), rexcol(kk)
       !     write(*,"(1i8,16e16.6)") irec, dwinfil(kk)-dwcol(kk)-dwdrainage(kk)-drexcol(kk)-qrunoff(kk)*dt(kk)
       ! endif

       ! evaluate soil fluxes at sigma of time step for use in isotope routine

       select case (surface_case(kk))
       case(1)  ! no snow
          qsig(0)  = q(0)  + sig(kk)*(qyb(0)*dy(1)  + qya(0)*dy(0)  + qTb(0)*de(1) &
               + qTa(0)*de(0))
          qhsig(0) = qh(0) + sig(kk)*(qhyb(0)*dy(1) + qhya(0)*dy(0) + qhTb(0)*de(1) &
               + qhTa(0)*de(0))
          qvsig(0) = qv(0)+sig(kk)*qvyb(0)*dy(1)+sig(kk)*qvTb(0)*de(1)
          qlsig(0) = qliq(0)+sig(kk)*qlyb(0)*dy(1)+sig(kk)*qlTb(0)*de(1)

          qsig(1:n-1)   = q(1:n-1) + sig(kk)*(qya(1:n-1)*dy(1:n-1) + qyb(1:n-1)*dy(2:n) &
               + qTa(1:n-1)*de(1:n-1) + qTb(1:n-1)*de(2:n))
          qsig(n)       = q(n) + sig(kk)*(qya(n)*dy(n) + qTa(n)*de(n))
          qhsig(1:n-1)  = qh(1:n-1) + sig(kk)*(qhya(1:n-1)*dy(1:n-1) + qhyb(1:n-1)*dy(2:n) &
               + qhTa(1:n-1)*de(1:n-1) + qhTb(1:n-1)*de(2:n))
          qhsig(n)      = qh(n) + sig(kk)*(qhya(n)*dy(n) + qhTa(n)*de(n))
          qvsig(1:n-1)  = qv(1:n-1) + sig(kk)*(qvya(1:n-1)*dy(1:n-1) + qvyb(1:n-1)*dy(2:n) &
               + qvTa(1:n-1)*de(1:n-1) + qvTb(1:n-1)*de(2:n))
          qvsig(n)      = zero

          qvTsig(0)     = qvT(0) + sig(kk)*qvTb(0)*de(1)
          qvTsig(1:n-1) = qvT(1:n-1) + sig(kk)*(qvTa(1:n-1)*de(1:n-1) + qvTb(1:n-1)*de(2:n))
          qvTsig(n)     = zero
          qlsig(1:n-1)  = qsig(1:n-1) - qvsig(1:n-1)
          qlsig(n)      = qsig(n)

       case(2) ! dedicated snow layer
          qsig(-vsnow(kk)%nsnow)  = q(-vsnow(kk)%nsnow)  + &
               sig(kk)*(qyb(-vsnow(kk)%nsnow)*dy(-vsnow(kk)%nsnow+1)  + &
                                ! qya(-vsnow(kk)%nsnow)*dy(-vsnow(kk)%nsnow)  + &
               qTb(-vsnow(kk)%nsnow)*de(-vsnow(kk)%nsnow+1)) ! &
          !  + qTa(-vsnow(kk)%nsnow)*de(-vsnow(kk)%nsnow))
          qhsig(-vsnow(kk)%nsnow) = qh(-vsnow(kk)%nsnow) + &
               sig(kk)*(qhyb(-vsnow(kk)%nsnow)*dy(-vsnow(kk)%nsnow+1) + &
                                ! qhya(-vsnow(kk)%nsnow)*dy(-vsnow(kk)%nsnow) + &
               qhTb(-vsnow(kk)%nsnow)*de(-vsnow(kk)%nsnow+1)) !&
          !  + qhTa(-vsnow(kk)%nsnow)*de(-vsnow(kk)%nsnow))
          qvsig(-vsnow(kk)%nsnow) = qv(-vsnow(kk)%nsnow) + &
               sig(kk)*qvyb(-vsnow(kk)%nsnow)*dy(-vsnow(kk)%nsnow+1) + &
               sig(kk)*qvTb(-vsnow(kk)%nsnow)*de(-vsnow(kk)%nsnow+1)

          qlsig(-vsnow(kk)%nsnow:-1) = zero
          qlsig(0) = qliq(0)+sig(kk)*qlyb(0)*dy(1)+sig(kk)*qlTb(0)*de(1)

          qsig(-vsnow(kk)%nsnow+1:n-1)   = q(-vsnow(kk)%nsnow+1:n-1) + &
               sig(kk)*(qya(-vsnow(kk)%nsnow+1:n-1)*dy(-vsnow(kk)%nsnow+1:n-1) + &
               qyb(-vsnow(kk)%nsnow+1:n-1)*dy(-vsnow(kk)%nsnow+2:n) &
               + qTa(-vsnow(kk)%nsnow+1:n-1)*de(-vsnow(kk)%nsnow+1:n-1) + &
               qTb(-vsnow(kk)%nsnow+1:n-1)*de(-vsnow(kk)%nsnow+2:n))
          qsig(n)       = q(n) + sig(kk)*(qya(n)*dy(n) + qTa(n)*de(n))
          qhsig(-vsnow(kk)%nsnow+1:n-1)  = qh(-vsnow(kk)%nsnow+1:n-1) + &
               sig(kk)*(qhya(-vsnow(kk)%nsnow+1:n-1)*dy(-vsnow(kk)%nsnow+1:n-1) + &
               qhyb(-vsnow(kk)%nsnow+1:n-1)*dy(-vsnow(kk)%nsnow+2:n) &
               + qhTa(-vsnow(kk)%nsnow+1:n-1)*de(-vsnow(kk)%nsnow+1:n-1) + &
               qhTb(-vsnow(kk)%nsnow+1:n-1)*de(-vsnow(kk)%nsnow+2:n))
          qhsig(n)      = qh(n) + sig(kk)*(qhya(n)*dy(n) + qhTa(n)*de(n))
          qvsig(-vsnow(kk)%nsnow+1:n-1)  = qv(-vsnow(kk)%nsnow+1:n-1) + &
               sig(kk)*(qvya(-vsnow(kk)%nsnow+1:n-1)*dy(-vsnow(kk)%nsnow+1:n-1) + &
               qvyb(-vsnow(kk)%nsnow+1:n-1)*dy(-vsnow(kk)%nsnow+2:n) &
               + qvTa(-vsnow(kk)%nsnow+1:n-1)*de(-vsnow(kk)%nsnow+1:n-1) + &
               qvTb(-vsnow(kk)%nsnow+1:n-1)*de(-vsnow(kk)%nsnow+2:n))
          qvsig(n)      = zero

          qvTsig(0)     = qvT(0) + sig(kk)*qvTb(0)*de(1)
          qvTsig(-vsnow(kk)%nsnow+1:n-1) = qvT(-vsnow(kk)%nsnow+1:n-1) + &
               sig(kk)*(qvTa(-vsnow(kk)%nsnow+1:n-1)*de(-vsnow(kk)%nsnow+1:n-1) + &
               qvTb(-vsnow(kk)%nsnow+1:n-1)*de(-vsnow(kk)%nsnow+2:n))
          qvTsig(n)     = zero
          qlsig(-vsnow(kk)%nsnow+1:n-1)  = qsig(-vsnow(kk)%nsnow+1:n-1) - qvsig(-vsnow(kk)%nsnow+1:n-1)
          qlsig(n)      = qsig(n)
       end select ! surface_case

       if (botbc=="aquifer") then ! update aquifer props
          if (v_aquifer(kk)%isat==0) then
             if (v_aquifer(kk)%WA+(dwdrainage(kk)-dwdischarge(kk))*dt(kk) > &
                  (v_aquifer(kk)%zzero-v_aquifer(kk)%zsoil)*v_aquifer(kk)%Sy) then
                v_aquifer(kk)%isat = 1  ! aquifer saturated
                S(n) = S(n) + (-(v_aquifer(kk)%WA+(dwdrainage(kk)-dwdischarge(kk))*dt(kk)) &
                     +(v_aquifer(kk)%zzero-v_aquifer(kk)%zsoil)*v_aquifer(kk)%Sy)/(dx(n)*par(n)%thre)
             endif
             v_aquifer(kk)%WA        = min(v_aquifer(kk)%WA+(dwdrainage(kk)-dwdischarge(kk))*dt(kk), &
                  (v_aquifer(kk)%zzero-v_aquifer(kk)%zsoil)*v_aquifer(kk)%Sy)
             v_aquifer(kk)%zdelta    = v_aquifer(kk)%zzero - v_aquifer(kk)%Wa/v_aquifer(kk)%Sy
             ! new discharge rate
             v_aquifer(kk)%discharge = v_aquifer(kk)%Rsmax*exp(-v_aquifer(kk)%f*v_aquifer(kk)%zdelta)
          elseif (v_aquifer(kk)%isat==1) then
             ! check for desat of aquifer
             if (dwdrainage(kk) < v_aquifer(kk)%Rsmax*exp(-v_aquifer(kk)%f*v_aquifer(kk)%zsoil)) then
                v_aquifer(kk)%isat      = 0
                v_aquifer(kk)%zdelta    = v_aquifer(kk)%zsoil
                ! new discharge rate
                v_aquifer(kk)%discharge = v_aquifer(kk)%Rsmax*exp(-v_aquifer(kk)%f*v_aquifer(kk)%zdelta)
             endif
          endif
          zdelta(kk) = v_aquifer(kk)%zdelta
       endif       ! end update aquifer props

       if (botbc=="constant head") then
          drn(kk) = drn(kk)+(q(n)+sig(kk)*qya(n)*dy(n))*dt(kk)
       else
          drn(kk) = drn(kk)+(q(n)+sig(kk)*qya(n)*dy(n))*dt(kk)
       end if

       if (present(wex)) then
          wex(kk,1:n) = wex(kk,1:n)+(iqex(1:n)+spread(sig(kk),1,n)*qexd(1:n)*dy(1:n))*spread(dt(kk),1,n)
       end if

       if (litter .and. ns(kk)==1 ) then ! adjust litter moisture content if no ponding
          SL(kk)      = SL(kk) + dy(0)
          deltaSL(kk) = dy(0)
       else
          deltaSL(kk) = zero
       endif

       TL(kk) = TL(kk) + de(0)   ! Tlitter assumed the same as pond T
       if (SL(kk) >= one) vlit(kk)%isat = 1   ! check for litter saturation

       csoil(:) = var(1:n)%csoil
       kth(:) = var(1:n)%kth
       phi(:) = var(1:n)%phi

       ! update snow pack water content and temperature (or liquid water)
       if (vsnow(kk)%nsnow>0) then
          vsnow(kk)%Qadv_snow = vsnow(kk)%Qadv_snow + rhow*(qprec_snow(kk))* &
               (csice*(min(vmet(kk)%Ta,zero))-lambdaf)*dt(kk)
          vsnow(kk)%Qadv_rain = vsnow(kk)%Qadv_rain + rhow*(qprec(kk))*cswat*(max(vmet(kk)%Ta,zero))*dt(kk)

          ! update heat flux components
          Tqw  = merge(vmet(kk)%Ta, vsnow(kk)%tsn(vsnow(kk)%nsnow), -qevap(kk)>zero)

          !!vsnow(kk)%Qadv_vap = vsnow(kk)%Qadv_vap + &
          !   rhow*(-qevapsig(kk))*cswat*Tqw*dt(kk) - &
          !    qadvsig(0)*dt(kk)

          vsnow(kk)%Qadv_vap =vsnow(kk)%Qadv_vap + (qadvsig(-vsnow(kk)%nsnow)-qadvsig(0))*dt(kk) - &
               rhow*(qprec_snow(kk))*(csice*(min(vmet(kk)%Ta,zero))-lambdaf)*dt(kk) - &
               rhow*(qprec(kk))*cswat*(vmet(kk)%Ta)*dt(kk)

          vsnow(kk)%Qcond_net = vsnow(kk)%Qcond_net + (qhsig(-vsnow(kk)%nsnow) -  qhsig(0))*dt(kk) - &
               (qadvsig(-vsnow(kk)%nsnow) -  qadvsig(0))*dt(kk)

          do i=1, vsnow(kk)%nsnow
             if (vsnow(kk)%hliq(i)>zero) then

                if ((vsnow(kk)%hliq(i) + de(i-vsnow(kk)%nsnow))>(dy(i-vsnow(kk)%nsnow) + vsnow(kk)%hsnow(i))) then
                   ! account here for new hliq exceeding new hsnow
                   tmp1d1(kk) = vsnow(kk)%hliq(i)*cswat*rhow*zero + &
                        (vsnow(kk)%hsnow(i)-vsnow(kk)%hliq(i))*rhow*(zero*csice-lambdaf)
                   tmp1d2(kk) = ((tmp1d1(kk) + LHS_h(i-vsnow(kk)%nsnow)*dt(kk)))/ &
                        ((dy(i-vsnow(kk)%nsnow) + vsnow(kk)%hsnow(i))*cswat*rhow)
                   vsnow(kk)%tsn(i) = tmp1d2(kk)

                   if (vsnow(kk)%tsn(i)>zero) then
                      write(*,*) 'T>zero ', irec, nsteps(kk), i, vsnow(kk)%tsn(i)
                   endif

                   vsnow(kk)%deltaJlatent(i) = vsnow(kk)%deltaJlatent(i) + &
                        lambdaf*(vsnow(kk)%hsnow(i)-vsnow(kk)%hliq(i))*rhow
                   vsnow(kk)%deltaJsensible(i) = vsnow(kk)%deltaJsensible(i) + LHS_h(i-vsnow(kk)%nsnow)*dt(kk) - &
                        lambdaf*(vsnow(kk)%hsnow(i)-vsnow(kk)%hliq(i))*rhow
                   vsnow(kk)%hliq(i) = vsnow(kk)%hsnow(i) + dy(i-vsnow(kk)%nsnow)
                elseif ((vsnow(kk)%hliq(i)+de(i-vsnow(kk)%nsnow))<zero) then
                   ! account here for new hliq less than zero
                   tmp1d1(kk) = (csice*(vsnow(kk)%tsn(i))-lambdaf)*rhow*(vsnow(kk)%hsnow(i)-vsnow(kk)%hliq(i)) + &
                        cswat*(vsnow(kk)%tsn(i))*rhow*vsnow(kk)%hliq(i) + LHS_h(i-vsnow(kk)%nsnow)*dt(kk)
                   delta_snowT(i)  = (tmp1d1(kk)/rhow/(dy(i-vsnow(kk)%nsnow) + vsnow(kk)%hsnow(i))+lambdaf)/csice
                   vsnow(kk)%tsn(i) = delta_snowT(i)

                   vsnow(kk)%deltaJlatent(i) = vsnow(kk)%deltaJlatent(i) - lambdaf*dy(i-vsnow(kk)%nsnow)*rhow + &
                        lambdaf*vsnow(kk)%hliq(i)*rhow
                   vsnow(kk)%deltaJsensible(i) = vsnow(kk)%deltaJsensible(i) + LHS_h(i-vsnow(kk)%nsnow)*dt(kk) - &
                        (- lambdaf*dy(i-vsnow(kk)%nsnow)*rhow + lambdaf*vsnow(kk)%hliq(i)*rhow)
                   vsnow(kk)%hliq(i) = zero
                else
                   vsnow(kk)%hliq(i) = vsnow(kk)%hliq(i) + de(i-vsnow(kk)%nsnow)
                   delta_snowliq(i) = de(i-vsnow(kk)%nsnow)
                   delta_snowT(i) = zero
                   vsnow(kk)%deltaJlatent(i) = vsnow(kk)%deltaJlatent(i) - &
                        lambdaf*dy(i-vsnow(kk)%nsnow)*rhow + lambdaf*de(i-vsnow(kk)%nsnow)*rhow
                   vsnow(kk)%deltaJsensible(i) = vsnow(kk)%deltaJsensible(i) + &
                        dy(i-vsnow(kk)%nsnow)*rhow*csice*(vsnow(kk)%tsn(i)) + &
                        de(i-vsnow(kk)%nsnow)*rhow*(cswat-csice)*(vsnow(kk)%tsn(i))
                endif
             else !hliq = zero
                if ((vsnow(kk)%tsn(i) + de(i-vsnow(kk)%nsnow))<zero) then
                   vsnow(kk)%deltaJlatent(i) = vsnow(kk)%deltaJlatent(i) - lambdaf*dy(i-vsnow(kk)%nsnow)*rhow
                   vsnow(kk)%deltaJsensible(i) = vsnow(kk)%deltaJsensible(i) + &
                        dy(i-vsnow(kk)%nsnow)*rhow*csice*(vsnow(kk)%tsn(i)) + &
                        de(i-vsnow(kk)%nsnow)*csice*rhow*hsnow(i)
                   vsnow(kk)%tsn(i) = vsnow(kk)%tsn(i) + de(i-vsnow(kk)%nsnow)
                   delta_snowT(i) = de(i-vsnow(kk)%nsnow)
                   delta_snowliq(i) = zero
                   vsnow(kk)%hliq(i) = zero
                else ! vsnow(kk)%tsn + de(0))>zero
                   delta_snowT(i) = zero - vsnow(kk)%tsn(i)
                   delta_snowliq(i) = (LHS_h(i-vsnow(kk)%nsnow)*dt(kk) - &
                        delta_snowT(i)*csice*hsnow(i) *rhow - &
                        dy(i-vsnow(kk)%nsnow)*rhow*(csice*(vsnow(kk)%tsn(i))-lambdaf))/ &
                        (rhow*(zero*(cswat-csice)+lambdaf))

                   vsnow(kk)%hliq(i) = vsnow(kk)%hliq(i) + delta_snowliq(i)
                   vsnow(kk)%deltaJlatent(i) = vsnow(kk)%deltaJlatent(i) &
                        - lambdaf*dy(i-vsnow(kk)%nsnow)*rhow + &
                        lambdaf*rhow*delta_snowliq(i)
                   vsnow(kk)%deltaJsensible(i) = vsnow(kk)%deltaJsensible(i)  &
                        + dy(i-vsnow(kk)%nsnow)*rhow*csice*(vsnow(kk)%tsn(i)) + &
                        delta_snowT(i)*csice*rhow*hsnow(i)
                   vsnow(kk)%tsn(i) = vsnow(kk)%tsn(i) + delta_snowT(i)
                endif ! vsnow(kk)%tsn + de(0))>zero
             endif ! hliq = zero
             vsnow(kk)%hsnow(i) = dy(i-vsnow(kk)%nsnow) + vsnow(kk)%hsnow(i)
             vsnow(kk)%Jsensible(i) = (vsnow(kk)%hsnow(i)-vsnow(kk)%hliq(i))*rhow*csice*(vsnow(kk)%Tsn(i))
             vsnow(kk)%Jlatent(i) = (vsnow(kk)%hsnow(i)-vsnow(kk)%hliq(i))*rhow*(-lambdaf)

             delta_snowcol(i) = dy(i-vsnow(kk)%nsnow)
          enddo ! loop over snow layers
          ! total internal energy and water content of snowpack
          vsnow(kk)%wcol = sum(vsnow(kk)%hsnow(1:vsnow(kk)%nsnow))
          vsnow(kk)%deltawcol = vsnow(kk)%wcol - wcol0snow(kk)
          vsnow(kk)%Qprec = vsnow(kk)%Qprec + qprec_snow(kk)*dt(kk) + qprec(kk)*dt(kk)
          vsnow(kk)%Qevap = vsnow(kk)%Qevap + qevapsig(kk)*dt(kk)
          vsnow(kk)%Qvap = vsnow(kk)%Qvap + qsig(0)*dt(kk)

          vsnow(kk)%J = sum(vsnow(kk)%Jsensible(1:vsnow(kk)%nsnow)+vsnow(kk)%Jlatent(1:vsnow(kk)%nsnow))

          vsnow(kk)%deltaJ = vsnow(kk)%J - J0snow(kk)

          if (vsnow(kk)%nsnow.lt.nsnow_max) then
             vsnow(kk)%hsnow(vsnow(kk)%nsnow+1:nsnow_max) = zero
          endif

          ! update densities and depths of snow layers

          ! effect of new snowfall

          ! snowfall density (tmp1d1), LaChapelle 1969
          if (vmet(kk)%Ta > 2.0_r_2) then
             tmp1d1(kk) = 189.0_r_2
          elseif ((vmet(kk)%Ta > -15.0_r_2) .and. (vmet(kk)%Ta <= 2.0_r_2)) then
             tmp1d1(kk) = 50.0_r_2 + 1.7_r_2*(vmet(kk)%Ta+15.0_r_2)**1.5_r_2
          else
             tmp1d1(kk) = 50.0_r_2
          endif

          if ((vsnow(kk)%hsnow(1)-qprec_snow(kk)*dt(kk)) > 1.e-3_r_2) then
             vsnow(kk)%dens(1) = (vsnow(kk)%dens(1)*(vsnow(kk)%hsnow(1)-qprec_snow(kk)*dt(kk)) &
                  + tmp1d1(kk)*qprec_snow(kk)*dt(kk))/vsnow(kk)%hsnow(1)
          endif

          do i=1, vsnow(kk)%nsnow
             vsnow(kk)%depth(i) = vsnow(kk)%hsnow(i)/(vsnow(kk)%dens(i)/rhow)

             ! adjust depth for compaction
             if (vsnow(kk)%hliq(i) .gt. zero) then
                tmp1d1(kk) = merge(one, exp(-0.06_r_2*(vsnow(kk)%dens(i)-150._r_2)), vsnow(kk)%dens(i)<150._r_2)
                tmp1d2(kk) = 2
             else
                tmp1d1(kk) = merge(one, exp(-0.046_r_2*(vsnow(kk)%dens(i)-150._r_2)), vsnow(kk)%dens(i)<150._r_2)
                tmp1d2(kk) = 1
             endif
             tmp1d3(kk) =  vsnow(kk)%hsnow(i)*rhow/2._r_2    ! overburden kg/m2
             if (i.gt.1) then
                tmp1d3(kk) = tmp1d3(kk) + sum(vsnow(kk)%hsnow(1:i-1)*rhow)
             endif

             if (vsnow(kk)%depth(i)*dt(kk)* &
                  (2.8e-6_r_2*tmp1d1(kk)*tmp1d2(kk)*exp(-vsnow(kk)%tsn(i)/25.0_r_2) + &
                  tmp1d3(kk)/3.6e6_r_2*exp(-0.08_r_2*vsnow(kk)%tsn(i))*exp(-0.023_r_2*vsnow(kk)%dens(i))) &
                  .lt.0.5_r_2*vsnow(kk)%depth(i)) then
                vsnow(kk)%depth(i) = vsnow(kk)%depth(i) - vsnow(kk)%depth(i)*dt(kk)* &
                     (2.8e-6_r_2*tmp1d1(kk)*tmp1d2(kk)*exp(-vsnow(kk)%tsn(i)/25.0_r_2) + & ! metamorphism
                     tmp1d3(kk)/3.6e6_r_2*exp(-0.08_r_2*vsnow(kk)%tsn(i))*exp(-0.023_r_2*vsnow(kk)%dens(i))) ! overburden
             endif
             vsnow(kk)%dens(i) = vsnow(kk)%hsnow(i)*rhow/vsnow(kk)%depth(i)
          enddo

       endif ! end update snowlayers
       vsnow(kk)%MoistureFluxDivergence = vsnow(kk)%Qprec - vsnow(kk)%Qevap - vsnow(kk)%Qvap &
            - vsnow(kk)%Qmelt +  vsnow(kk)%Qtransfer
       vsnow(kk)%FluxDivergence = vsnow(kk)%Qadv_rain + vsnow(kk)%Qadv_snow + vsnow(kk)%Qadv_vap + &
            vsnow(kk)%Qadv_melt + vsnow(kk)%Qcond_net  +  vsnow(kk)%Qadv_transfer

    endif ! .not. again

  END SUBROUTINE update_unknowns

  SUBROUTINE get_and_solve_eqn( &
       irec, mp, qprec,  n, dx, h0, S, thetai, &
       Tsoil,       &
       qh, nsteps,  vlit, vsnow, var,     &
       dxL,    &
       plit, par,        &
       deltaTa, &
       qevap, qL, qhL, qybL, qTbL, qhTbL, qhybL,    &
       again,  again_ice,  iok, itmp, ns,  &
       accel,  dt,   fac,   &
       phip,  rsig, rsigdt, sig, t,      &
       qexd, aa, bb, cc, dd, ee, ff, gg, dy, aah, bbh, cch, ddh, eeh, ffh, ggh, &
       de, q, qya, qyb, qTa, qTb,qhya, qhyb, qhTa, qhTb, qadv, qadvya, qadvyb, &
       qadvTa, qadvTb,  qsig, qhsig, qadvsig,      &
       dTsoil,   &
       dthetaldT,   nsteps_ice, imelt,     &
       deltaJ_latent_T, deltaJ_sensible_S,   &
       qrunoff, tmp1d1, tmp1d2, tmp1d3,  tmp1d4,     &
       G0, &
       Tfreezing,   LHS, RHS, LHS_h, RHS_h,  nns, &
       iflux, litter, i,   kk, condition,   &
       advection,  c2, theta,     cp, &
       cpeff, hice,    h0_tmp,   hsnow, &
       delta_snowT, &
       delta_snowliq,    J0,   iqex,  &
       nfac1, nfac2, nfac3, nfac4, nfac5, nfac6, nfac7, nfac8, nfac9, &
       nfac10, nfac11, nfac12, err  &
       )
    IMPLICIT NONE
    INTEGER(i_d)                                           :: irec, mp
    REAL(r_2),      DIMENSION(1:mp)                        :: qprec
    INTEGER(i_d)                                           :: n
    REAL(r_2),      DIMENSION(1:n) :: dx
    REAL(r_2),      DIMENSION(1:mp)                        :: h0
    REAL(r_2),      DIMENSION(1:n) :: S
    REAL(r_2),      DIMENSION(1:n) :: thetai
    REAL(r_2),      DIMENSION(1:n) :: Tsoil
    REAL(r_2),      DIMENSION(-nsnow_max:n) :: qh
    INTEGER(i_d),   DIMENSION(1:mp)                        :: nsteps
    TYPE(vars),     DIMENSION(1:mp)                        :: vlit
    TYPE(vars_snow), DIMENSION(1:mp)                       :: vsnow
    TYPE(vars),     DIMENSION(1:n) :: var
    REAL(r_2),      DIMENSION(1:mp)                        :: dxL
    TYPE(params),   DIMENSION(1:mp)                        :: plit
    TYPE(params),   DIMENSION(1:n) :: par
    REAL(r_2),      DIMENSION(1:mp)                        :: deltaTa
    REAL(r_2),    DIMENSION(1:mp)                          :: qevap
    REAL(r_2),    DIMENSION(1:mp)                          :: qL, qhL, qybL, qTbL, qhTbL, qhybL
    LOGICAL,      DIMENSION(1:mp)                          :: again
    LOGICAL,      DIMENSION(1:mp,1:n)                      :: again_ice
    INTEGER(i_d), DIMENSION(1:mp)                          :: iok, itmp, ns
    REAL(r_2),    DIMENSION(1:mp)                          :: accel, dt, fac, phip
    REAL(r_2),    DIMENSION(1:mp)                          :: rsig, rsigdt, sig, t
    REAL(r_2),    DIMENSION(1:n) :: qexd
    REAL(r_2),    DIMENSION(-nsnow_max+1:n) :: aa, bb, cc, dd, ee, ff, gg, dy
    REAL(r_2),    DIMENSION(-nsnow_max+1:n) :: aah, bbh, cch, ddh, eeh, ffh, ggh, de
    REAL(r_2),    DIMENSION(-nsnow_max:n) :: q, qya, qyb, qTa, qTb,qhya, qhyb, qhTa, qhTb
    REAL(r_2),    DIMENSION(-nsnow_max:n) :: qadv, qadvya, qadvyb, qadvTa, qadvTb
    REAL(r_2),    DIMENSION(-nsnow_max:n) :: qsig, qhsig, qadvsig
    REAL(r_2),    DIMENSION(1:n) :: dTsoil
    REAL(r_2),    DIMENSION(1:n) :: dthetaldT
    INTEGER(i_d), DIMENSION(1:n) :: nsteps_ice, imelt
    REAL(r_2),          DIMENSION(1:n) :: deltaJ_latent_T, deltaJ_sensible_S
    REAL(r_2),          DIMENSION(1:mp)                    :: qrunoff
    REAL(r_2),          DIMENSION(1:mp)                    :: tmp1d1, tmp1d2, tmp1d3,  tmp1d4
    REAL(r_2),          DIMENSION(1:mp)                    :: G0
    REAL(r_2),          DIMENSION(1:mp)                    :: Tfreezing
    REAL(r_2),          DIMENSION(-nsnow_max+1:n) :: LHS, RHS, LHS_h, RHS_h
    INTEGER(i_d),       DIMENSION(1:mp)                    :: nns, iflux
    LOGICAL                                                :: litter
    INTEGER(i_d)                                           :: i, kk, condition
    INTEGER(i_d)                                           :: advection ! switches
    REAL(r_2)                                              :: c2, theta
    REAL(r_2),          DIMENSION(1:mp)                    :: cp, cpeff, hice, h0_tmp
    REAL(r_2),          DIMENSION(nsnow_max) :: hsnow
    REAL(r_2),          DIMENSION(nsnow_max) :: delta_snowT, delta_snowliq
    REAL(r_2),      DIMENSION(1:n) :: J0
    REAL(r_2),          DIMENSION(1:n) :: iqex
    INTEGER(i_d),       DIMENSION(1:mp)                    :: nfac1, nfac2, nfac3, nfac4, nfac5, &
         nfac6, nfac7, nfac8, nfac9, nfac10, nfac11, nfac12
    ! Error flag if nstep of SLI > nsteps_max: err=0 -> no error; err/=0 -> error
    INTEGER(i_d),       INTENT(INOUT), OPTIONAL            :: err

    iok(kk)         = 0 ! flag for time step test
    itmp(kk)        = 0 ! counter to abort if not getting solution
    again_ice(kk,:) = .false.
    nsteps_ice(1:n) = 0
    h0_tmp(kk) = h0(kk)
    do while (iok(kk)==0) ! keep reducing time step until all ok
       itmp(kk)  = itmp(kk) + 1
       accel(kk) = one - 0.05_r_2*real(min(10,max(0,itmp(kk)-4)),r_2) ! acceleration [0.5,1], start with 1
       if (itmp(kk) > 1000) then
          write(wlogn,*) "Solve: too many iterations of equation solution"
          write(wlogn,*) " irec, kk, S"
          write(wlogn,*) irec, kk, S(:)
          write(wlogn,*) " irec, kk, Tsoil"
          write(wlogn,*) irec, kk, Tsoil(:)
          write(wlogn,*) " irec, kk, qex"
          write(wlogn,*) irec, kk, iqex(:)
          write(wlogn,*) " irec, kk, h0, hsnow, hsnowliq"
          write(wlogn,*) irec, kk, h0(kk), vsnow(kk)%hsnow, vsnow(kk)%hliq
          write(wlogn,*) nfac1(kk), nfac2(kk), nfac3(kk), nfac4(kk), nfac5(kk), &
               nfac6(kk), nfac7(kk), nfac8(kk), nfac9(kk), nfac10(kk), nfac11(kk), nfac12(kk)
          err = 1
          return
       end if

       ! prelim estimate of new top snow layer depth for use in energy cons eq'n             !
       if ((vsnow(kk)%nsnow>0))  then
          hsnow(1:vsnow(kk)%nsnow) = vsnow(kk)%hsnow(1:vsnow(kk)%nsnow) + &
               (-q(1-vsnow(kk)%nsnow:0) + q(-vsnow(kk)%nsnow:-1))*dt(kk)
       endif

       aa(1:n)   =  qya(0:n-1)
       ee(0:n-1) = -qyb(0:n-1)
       bb(1:n)   =  qTa(0:n-1)
       ff(0:n-1) = -qTb(0:n-1)
       gg(1:n) = -(q(0:n-1)-q(1:n)-iqex(1:n))*rsig(kk)
       gg(1) = gg(1)+qrunoff(kk)*rsig(kk)

       aah(1:n)   =  qhya(0:n-1)
       eeh(0:n-1) = -qhyb(0:n-1)
       bbh(1:n)   =  qhTa(0:n-1)
       ffh(0:n-1) = -qhTb(0:n-1)
       ggh(1:n) =  -(qh(0:n-1)-qh(1:n))*rsig(kk)

       if (advection==1) then
          ggh(1:n) = ggh(1:n) + iqex(1:n)*rsig(kk)*(Tsoil(1:n))*cswat*rhow
          ggh(1) = ggh(1) + qrunoff(kk)*rsig(kk)*(Tsoil(1))*cswat*rhow
       endif

       if (litter) then ! full litter model: litter in zeroth layer
          ! only use zeroth layer for litter (pond included in layer 1)
          write(*,*) 'Should not be here - QYBL 01 ', qybl(kk)
          cc(0) = -qya(0) - rsigdt(kk)*plit(kk)%thre*dxL(kk) + qybL(kk)
          gg(0)  = -(qprec(kk)-qevap(kk)-q(0))*rsig(kk)
          ggh(0) = -(G0(kk)-qh(0))*rsig(kk)
          dd(0)  = -qTa(0)
       endif
       aa(0)  = qya(0)
       aah(0) = zero
       bbh(0) = zero

       where (var(1:n)%isat==0) ! unsaturated layers
          cc(1:n) = qyb(0:n-1) - qya(1:n) - par(1:n)%thre*dx(1:n)*rsigdt(kk) - qexd(1:n)
       elsewhere ! saturated layers
          cc(1:n) = qyb(0:n-1) - qya(1:n) - qexd(1:n)
       endwhere

       if (ns(kk)<1) then ! pond included in top soil layer, solving for change in pond height
          cc(1) = -qya(1)-rsigdt(kk) -qexd(1)
       endif

       cch(1:n) = qhyb(0:n-1)-qhya(1:n) +   &
            real(var(1:n)%iice,r_2)*real(1-var(1:n)%isat,r_2)*rhow*lambdaf*par(1:n)%thre*dx(1:n)*rsigdt(kk)

       if (ns(kk)==0) then ! change in pond height (top layer saturated)
          cch(1) = cch(1) +  &
               real(var(1)%iice,r_2)*rhow*lambdaf*rsigdt(kk)*(var(1)%thetai/par(1)%thre)
       endif

       ! modification to cch for advection
       if (advection==1) then
          if (ns(kk)==0) then  ! changing pond included in top soil layer
             cch(1) = cch(1) -rhow*cswat*(Tsoil(1))*rsigdt(kk)*real(1-var(1)%iice,r_2) &
                  -rhow*csice*(Tsoil(1))*rsigdt(kk)*real(var(1)%iice,r_2) &
                  *(var(1)%thetai/par(1)%thre) &
                  -rhow*cswat*(Tsoil(1))*rsigdt(kk)*real(var(1)%iice,r_2)* &
                  (one-(var(1)%thetai/par(1)%thre))
          else ! no pond
             cch(1) = cch(1) -rhow*(Tsoil(1))*par(1)%thre*dx(1)*rsigdt(kk)* &
                  real(1-var(1)%isat,r_2)*(cswat*real(1-var(1)%iice,r_2) &
                  +csice*real(var(1)%iice,r_2))
          endif

          cch(2:n) = cch(2:n) -rhow*(Tsoil(2:n))*par(2:n)%thre*dx(2:n)*rsigdt(kk)* &
               real(1-var(2:n)%isat,r_2)*(cswat*real(1-var(2:n)%iice,r_2) &
               +csice*real(var(2:n)%iice,r_2))
       endif

       dd(1:n)  = qTb(0:n-1)-qTa(1:n)

       ddh(1:n) = qhTb(0:n-1)-qhTa(1:n) - &
            ! Only apply latent heat component of heat capacity to total deltaT if soil remains frozen
            var(1:n)%csoileff*dx(1:n)*rsigdt(kk)- &
            (cswat-csice)*dx(1:n)*var(1:n)%dthetaldt*rhow*(Tsoil(1:n))* &
            rsigdt(kk)*real(var(1:n)%iice,r_2)
       ! modification to ddh(1) for pond
       ddh(1) = ddh(1) - cpeff(kk)*h0_tmp(kk)*rsigdt(kk) - &
            (cswat-csice)*h0(kk)/par(1)%thre*var(1)%dthetaldt*rhow*(Tsoil(1))* &
            rsigdt(kk)*real(var(1)%iice,r_2)

       if (advection==1) then
          ddh(1:n) = ddh(1:n) - iqex(1:n)*cswat*rhow
          ddh(1) = ddh(1) - qrunoff(kk)*cswat*rhow
       endif

       ! modification of matrix to incorporate single snow layer
       if (vsnow(kk)%nsnow==1) then
          cc(0) = qyb(-1)-qya(0)-rsigdt(kk)
          dd(0) =  qTb(-1)-qTa(0)
          ee(0) =   -qyb(0)
          ff(0) =   -qTb(0)
          gg(0) = -(q(-1)-q(0))*rsig(kk)
          if (vsnow(kk)%hliq(1)>zero) then ! liquid phase present, solve for change in liq content
             ddh(0) =  - rhow*((vsnow(kk)%tsn(1))*(cswat-csice)+lambdaf)*rsigdt(kk)
          else ! solid phase only; solve for change in snow t
             ddh(0) = qhTb(-1) - qhTa(0) - rhow*csice*hsnow(1)*rsigdt(kk)
          endif
          cch(0) = qhyb(-1)-qhya(0) - rhow*(csice*(vsnow(kk)%tsn(1))-lambdaf)*rsigdt(kk)
          eeh(0) = -qhyb(0)
          ffh(0) = -qhTb(0)
          ggh(0) = -(qh(-1)-qh(0))*rsig(kk)
       endif
       ! modification of matrix to incorporate more than one snow layer
       if (vsnow(kk)%nsnow>1) then
          aa(1-vsnow(kk)%nsnow:0)   =  qya(-vsnow(kk)%nsnow:-1)
          bb(1-vsnow(kk)%nsnow:0)   =  qTa(-vsnow(kk)%nsnow:-1)
          cc(1-vsnow(kk)%nsnow:0) = qyb(-vsnow(kk)%nsnow:-1)-qya(1-vsnow(kk)%nsnow:0)-rsigdt(kk)
          dd(1-vsnow(kk)%nsnow:0) =  qTb(-vsnow(kk)%nsnow:-1)-qTa(1-vsnow(kk)%nsnow:0)
          ee(1-vsnow(kk)%nsnow:0) =   -qyb(1-vsnow(kk)%nsnow:0)
          ff(1-vsnow(kk)%nsnow:0) =   -qTb(1-vsnow(kk)%nsnow:0)
          gg(1-vsnow(kk)%nsnow:0) = -(q(-vsnow(kk)%nsnow:-1)-q(1-vsnow(kk)%nsnow:0))*rsig(kk)
          ddh(1-vsnow(kk)%nsnow:0) = merge( & ! liquid phase present, solve for change in liq content
               -rhow*((vsnow(kk)%tsn(1:vsnow(kk)%nsnow))*(cswat-csice)+lambdaf)*rsigdt(kk), &
               qhTb(-vsnow(kk)%nsnow:-1) - qhTa(1-vsnow(kk)%nsnow:0) - &
               rhow*csice*hsnow(1:vsnow(kk)%nsnow)*rsigdt(kk), & ! solid phase only; solve for change in snow t
               vsnow(kk)%hliq(1:vsnow(kk)%nsnow) > zero)
          aah(1-vsnow(kk)%nsnow:0)   =  qhya(-vsnow(kk)%nsnow:-1)
          bbh(1-vsnow(kk)%nsnow:0)   =  qhTa(-vsnow(kk)%nsnow:-1)

          cch(1-vsnow(kk)%nsnow:0) = qhyb(-vsnow(kk)%nsnow:-1)-qhya(1-vsnow(kk)%nsnow:0) - &
               rhow*(csice*(vsnow(kk)%tsn(1:vsnow(kk)%nsnow))-lambdaf)*rsigdt(kk)
          eeh(1-vsnow(kk)%nsnow:0) = -qhyb(1-vsnow(kk)%nsnow:0)
          ffh(1-vsnow(kk)%nsnow:0) = -qhTb(1-vsnow(kk)%nsnow:0)
          ggh(1-vsnow(kk)%nsnow:0) = -(qh(-vsnow(kk)%nsnow:-1)-qh(1-vsnow(kk)%nsnow:0))*rsig(kk)
       endif
       if (litter .and. ns(kk)==1) then ! litter and no pond
          ! watch for deltaTa
          write(*,*) 'Should not be here - QYBL 02 ', qybl(kk)
          cc(0)  = qybL(kk) -qya(0) -rsigdt(kk)*plit(kk)%thre*dxL(kk)
          dd(0)  = qTbL(kk) - qTa(0)
          ee(0)  = -qyb(0)
          ff(0)  = -qTb(0)
          gg(0)  = (q(0) - qL(kk))*rsig(kk) + deltaTa(kk)*(qTbL(kk)-qTa(0))
          gg(1)  = -(q(0)-q(1)-iqex(1))*rsig(kk) + deltaTa(kk)*qTa(0)
          cch(0) = qhybL(kk) - qhya(0)
          ddh(0) = -qhTa(0)-vlit(kk)%csoil*dxL(kk)*rsigdt(kk)+qhTbL(kk)
          eeh(0) = -qhyb(0)
          ffh(0) = -qhTb(0)
          ggh(0) = (qh(0)-qhL(kk))*rsig(kk) + deltaTa(kk)*(qhTbL(kk) - qhTa(0))
          ggh(1) = -(qh(0)-qh(1))*rsig(kk) + deltaTa(kk)*qhTa(0)
       endif

       ! litter and pond !!vh!! need to check this now that pond is lumped with top soil layer
       if (litter .and. ns(kk)==0) then
          ddh(0) = -qhTa(0)-cswat*rhow*h0(kk)*rsigdt(kk) -vlit(kk)%csoil*dxL(kk)*rsigdt(kk)
       endif

       nns(kk) = 1  ! pond included in top soil layer
       if (vsnow(kk)%nsnow>0) then
          nns(kk) = 1-vsnow(kk)%nsnow
       endif

       call massman_sparse(aa(nns(kk)+1:n), aah(nns(kk)+1:n), bb(nns(kk)+1:n), bbh(nns(kk)+1:n), &
            cc(nns(kk):n), cch(nns(kk):n), dd(nns(kk):n), ddh(nns(kk):n), ee(nns(kk):n-1), &
            eeh(nns(kk):n-1), &
            ff(nns(kk):n-1), ffh(nns(kk):n-1), gg(nns(kk):n), ggh(nns(kk):n), &
            dy(nns(kk):n), de(nns(kk):n), condition=condition, err=err)
       if (err /= 0) then
          write(wlogn,*) "Sparse matrix solution failed ", irec, kk
          write(wlogn,*) Tsoil(1), S(1)
          return
       endif

       dTsoil(1:n) = de(1:n)


       where (vsnow(kk)%hliq(:)>zero)
          delta_snowT(:) = zero
       elsewhere
          delta_snowT(:) = de(lbound(de,1):0)
       endwhere

       ! evaluate soil fluxes at sigma of time step
       qsig(0)  = q(0)+sig(kk)*qyb(0)*dy(1) + sig(kk)*qya(0)*dy(0) + &
            sig(kk)*qTb(0)*dTsoil(1) + sig(kk)*qTa(0)*de(0)

       qhsig(0) = qh(0) + sig(kk)*qhyb(0)*dy(1) + sig(kk)*qhya(0)*dy(0) + &
            sig(kk)*qhTb(0)*dTsoil(1) + sig(kk)*qhTa(0)*de(0)

       qadvsig(0) = qadv(0) + sig(kk)*qadvyb(0)*dy(1) &
            + sig(kk)*qadvya(0)*dy(0) + sig(kk)*qadvTb(0)*dTsoil(1) + sig(kk)*qadvTa(0)*de(0)

       qsig(1:n-1) = q(1:n-1) + sig(kk)*(qya(1:n-1)*dy(1:n-1)+qyb(1:n-1)*dy(2:n) &
            +qTa(1:n-1)*dTsoil(1:n-1)+qTb(1:n-1)*dTsoil(2:n))
       qsig(n)     = q(n) + sig(kk)*(qya(n)*dy(n)+qTa(n)*dTsoil(n))

       qhsig(1:n-1) = qh(1:n-1) + sig(kk)*(qhya(1:n-1)*dy(1:n-1)+qhyb(1:n-1)*dy(2:n) &
            +qhTa(1:n-1)*dTsoil(1:n-1)+qhTb(1:n-1)*dTsoil(2:n))
       qhsig(n) = qh(n) + sig(kk)*(qhya(n)*dy(n)+qhTa(n)*dTsoil(n))

       qadvsig(1:n-1) = qadv(1:n-1) + sig(kk)*(qadvya(1:n-1)*dy(1:n-1)+qadvyb(1:n-1)*dy(2:n) &
            +qadvTa(1:n-1)*dTsoil(1:n-1)+qadvTb(1:n-1)*dTsoil(2:n))
       qadvsig(n) = qadv(n) + sig(kk)*(qadvya(n)*dy(n)+qadvTa(n)*dTsoil(n))

       LHS_h(1) = (dx(1)*var(1)%csoileff+h0_tmp(kk)*cpeff(kk))*dTsoil(1)/dt(kk) - &
            dy(1)*par(1)%thre*dx(1)*rhow*lambdaf/dt(kk)*var(1)%iice* &
            real(1-var(1)%isat,r_2) - &
            dy(1)*(var(1)%thetai/par(1)%thre)*rhow*lambdaf/dt(kk)*real(1-ns(kk),r_2)* &
            var(1)%iice  + &
            var(1)%dthetaldT*(Tsoil(1))*rhow*(cswat-csice)*dx(1)*dTsoil(1)/dt(kk)* &
            var(1)%iice + &
            var(1)%dthetaldT*(Tsoil(1))*rhow*(cswat-csice)*h0(kk)/par(1)%thre* &
            dTsoil(1)/dt(kk)*var(1)%iice

       LHS_h(2:n) = (dx(2:n)*var(2:n)%csoileff)*dTsoil(2:n)/dt(kk)  - &
            dy(2:n)*par(2:n)%thre*dx(2:n)*rhow*lambdaf/dt(kk)*var(2:n)%iice* &
            real(1-var(2:n)%isat,r_2) + &
            var(2:n)%dthetaldT*(Tsoil(2:n))*rhow*(cswat-csice)*dx(2:n)* &
            dTsoil(2:n)/dt(kk)*var(2:n)%iice

       RHS_h(1:n) = qhsig(0:n-1) -qhsig(1:n)

       if (advection==1) then
          LHS_h(1:n) = LHS_h(1:n) &
               + real(1-var(1:n)%iice,r_2) * real(1-var(1:n)%isat,r_2) * &
               dx(1:n) * dy(1:n)/dt(kk) * par(1:n)%thre * (Tsoil(1:n)) * rhow * cswat &
               + real(var(1:n)%iice,r_2) * real(1-var(1:n)%isat,r_2) * &
               dx(1:n) * dy(1:n)/dt(kk) * par(1:n)%thre * rhow * csice * (Tsoil(1:n))

          LHS_h(1) = LHS_h(1) &
               + real(1-ns(kk),r_2) * real(1-var(1)%iice,r_2) * &
               dy(1)/dt(kk) * rhow * cswat * (Tsoil(1)) &
               + real(1-ns(kk),r_2) * real(var(1)%iice,r_2)* &
               dy(1)/dt(kk) * rhow * csice * (Tsoil(1)) * (var(1)%thetai/par(1)%thre) &
               + real(1-ns(kk),r_2) * real(var(1)%iice,r_2)  * &
               dy(1)/dt(kk) * rhow * cswat * (Tsoil(1)) * (one-(var(1)%thetai/par(1)%thre))

          RHS_h(1:n) = RHS_h(1:n) - iqex(1:n)*cswat*rhow*(Tsoil(1:n)+ sig(kk)*dTsoil(1:n))
          RHS_h(1) = RHS_h(1) - qrunoff(kk)*cswat*rhow*(Tsoil(1) + sig(kk)*dTsoil(1))
       endif

       ! check mass balance on top soil layer
       if (ns(kk)==0) then  ! pond included in top soil layer
          LHS(1) = dy(1)/dt(kk)
       else
          LHS(1) = dy(1)*par(1)%thre*dx(1)*real(-var(1)%isat+1)/dt(kk)
       endif
       RHS(1) = qsig(0) - qsig(1)- iqex(1) - qrunoff(kk)
       LHS(2:n) = dy(2:n)*par(2:n)%thre*dx(2:n)*real(-var(2:n)%isat+1)/dt(kk)
       RHS(2:n) = qsig(1:n-1) - qsig(2:n)- iqex(2:n)
       ! snow pack
       if (vsnow(kk)%nsnow>0) then

          qsig(-vsnow(kk)%nsnow:-1)  = q(-vsnow(kk)%nsnow:-1)+sig(kk)*qyb(-vsnow(kk)%nsnow:-1)* &
               dy(1-vsnow(kk)%nsnow:0) +  sig(kk)*qTb(-vsnow(kk)%nsnow:-1)*de(1-vsnow(kk)%nsnow:0)

          qhsig(-vsnow(kk)%nsnow:-1) = qh(-vsnow(kk)%nsnow:-1) + sig(kk)*qhyb(-vsnow(kk)%nsnow:-1)* &
               dy(1-vsnow(kk)%nsnow:0) + sig(kk)*qhTb(-vsnow(kk)%nsnow:-1)*de(1-vsnow(kk)%nsnow:0)

          qadvsig(-vsnow(kk)%nsnow:-1) = qadv(-vsnow(kk)%nsnow:-1) + &
               sig(kk)*qadvyb(-vsnow(kk)%nsnow:-1)*dy(1-vsnow(kk)%nsnow:0) + &
               sig(kk)*qadvTb(-vsnow(kk)%nsnow:-1)*de(1-vsnow(kk)%nsnow:0)

          if (vsnow(kk)%nsnow>1) then
             qsig(1-vsnow(kk)%nsnow:-1) = qsig(1-vsnow(kk)%nsnow:-1) + &
                  sig(kk)*qya(1-vsnow(kk)%nsnow:-1)*dy(1-vsnow(kk)%nsnow:-1) + &
                  sig(kk)*qTa(1-vsnow(kk)%nsnow:-1)*de(1-vsnow(kk)%nsnow:-1)

             qhsig(1-vsnow(kk)%nsnow:-1) = qhsig(1-vsnow(kk)%nsnow:-1)  + &
                  sig(kk)*qhya(1-vsnow(kk)%nsnow:-1)*dy(1-vsnow(kk)%nsnow:-1) + &
                  sig(kk)*qhTa(1-vsnow(kk)%nsnow:-1)*de(1-vsnow(kk)%nsnow:-1)

             qadvsig(1-vsnow(kk)%nsnow:-1) = qadvsig(1-vsnow(kk)%nsnow:-1) + &
                  sig(kk)*qadvya(1-vsnow(kk)%nsnow:-1)*dy(1-vsnow(kk)%nsnow:-1) + &
                  sig(kk)*qadvTa(1-vsnow(kk)%nsnow:-1)*de(1-vsnow(kk)%nsnow:-1)

          endif

          ! rate of change of snow water col
          RHS(1-vsnow(kk)%nsnow:0) = qsig(-vsnow(kk)%nsnow:-1) - qsig(1-vsnow(kk)%nsnow:0)
          LHS(1-vsnow(kk)%nsnow:0) = dy(1-vsnow(kk)%nsnow:0)/dt(kk)
          RHS_h(1-vsnow(kk)%nsnow:0) = qhsig(-vsnow(kk)%nsnow:-1) - qhsig(1-vsnow(kk)%nsnow:0)
          LHS_h(1-vsnow(kk)%nsnow:0) = merge( &  ! liquid phase present, solve for change in liq content
               (dy(1-vsnow(kk)%nsnow:0)* rhow*(csice*(vsnow(kk)%tsn(1:vsnow(kk)%nsnow))-lambdaf) + &
               de(1-vsnow(kk)%nsnow:0)*rhow*((cswat-csice)*(vsnow(kk)%tsn(1:vsnow(kk)%nsnow))+lambdaf))/dt(kk), &
               (dy(1-vsnow(kk)%nsnow:0)*rhow* (csice*(vsnow(kk)%tsn(1:vsnow(kk)%nsnow))-lambdaf) + &
               de(1-vsnow(kk)%nsnow:0)*csice*rhow*hsnow(1:vsnow(kk)%nsnow))/dt(kk), &
               vsnow(kk)%hliq(1:vsnow(kk)%nsnow)>zero)

       endif ! snow pack

       tmp1d1 = nless
       ! dy contains dS or, for sat layers, dphi values
       iok(kk) = 1
       fac(kk) = one

       if (.not. again(kk)) then

          ! check if time step ok, if not then set fac to make it less
          iok(kk) = 1
          if (vsnow(kk)%nsnow>0) then
             do i=1, vsnow(kk)%nsnow
                if ((vsnow(kk)%hsnow(i) + dy(i-vsnow(kk)%nsnow))<zero) then
                   fac(kk) = 0.99_r_2*vsnow(kk)%hsnow(i)/abs(dy(i-vsnow(kk)%nsnow))
                   nfac6(kk) = nfac6(kk)+1
                   iok(kk) = 0
                   exit
                endif

                if (vsnow(kk)%hliq(i)>zero) then
                   if ((vsnow(kk)%hliq(i) + de(i-vsnow(kk)%nsnow))>(dy(i-vsnow(kk)%nsnow) + &
                        vsnow(kk)%hsnow(i))) then

                      ! adjust fac such that the residual snow layer will be removed at the next call to snow_adjust
                      fac(kk) = 0.99_r_2*(vsnow(kk)%hsnow(i) - vsnow(kk)%hliq(i))/ &
                           (de(i-vsnow(kk)%nsnow)-dy(i-vsnow(kk)%nsnow))
                      nfac7(kk) = nfac7(kk)+1
                      iok(kk) =0

                      exit
                   endif

                   if ((vsnow(kk)%hliq(i) + de(i-vsnow(kk)%nsnow)) < -0.01_r_2*vsnow(kk)%hsnow(i)) then
                      fac(kk) = (-0.009_r_2*vsnow(kk)%hsnow(i)-vsnow(kk)%hliq(i))/de(i-vsnow(kk)%nsnow)
                      nfac8(kk) = nfac8(kk)+1
                      iok(kk) =0
                      exit
                   endif
                elseif ((vsnow(kk)%tsn(i) +  de(i-vsnow(kk)%nsnow)).gt.zero) then
                   delta_snowT(i) = zero - vsnow(kk)%tsn(i)
                   delta_snowliq(i) = (LHS_h(i-vsnow(kk)%nsnow)*dt(kk) - &
                        delta_snowT(i)*csice*hsnow(i) *rhow - &
                        dy(i-vsnow(kk)%nsnow)*rhow*(csice*(vsnow(kk)%tsn(i))-lambdaf))/ &
                        (rhow*(zero*(cswat-csice)+lambdaf))
                   if (delta_snowliq(i).gt.(vsnow(kk)%hsnow(i) + dy(i-vsnow(kk)%nsnow))) then
                      fac(kk) = 0.9_r_2*vsnow(kk)%hsnow(i)/(delta_snowliq(i)- dy(i-vsnow(kk)%nsnow))
                      nfac12(kk) = nfac12(kk)+1
                      iok(kk) = 0
                   endif
                endif

             enddo
          endif

          do i=1, n
             if (var(i)%isat==0) then ! check change in S in initially unsaturated layers
                if (abs(dy(i)) > dSfac*dSmax) then
                   fac(kk) = max(half,accel(kk)*abs(dSmax/dy(i)))
                   nfac1(kk) = nfac1(kk)+1
                   iok(kk) = 0
                   exit
                end if
                if (-dy(i) > dSmaxr*S(i)) then ! Check relative moisture change
                   fac(kk) = max(half,accel(kk)*dSmaxr*S(i)/(-dSfac*dy(i)))
                   nfac2(kk) = nfac2(kk)+1
                   iok(kk) = 0
                   exit
                end if
                if (S(i) < one .and. S(i)+dy(i) > Smax) then ! Check for oversaturating,

                   fac(kk) = accel(kk)*(half*(one+Smax)-S(i))/dy(i)
                   nfac3(kk) = nfac3(kk)+1
                   !write(*,*) 'incrementing nfac3 due to layer', i
                   iok(kk) = 0
                   exit
                end if
                if (S(i) >= one .and. dy(i) > half*(Smax-one)) then ! Check for limit at oversaturation
                   fac(kk) = 0.25_r_2*(Smax-one)/dy(i)
                   nfac4(kk) = nfac4(kk)+1
                   iok(kk) = 0
                   exit
                end if
             end if
          end do

          ! Check absolute soil temperature change in frozen soil layers where updated Tsoil exceeds freezing point
          if (iok(kk)==1) then
             do i=1, n
                if(var(i)%iice==1) then
                   Tfreezing(kk) = Tfrz(S(i)+dy(i)*real(1-var(i)%isat), &
                        par(i)%he, one/(par(i)%lambc*freezefac))
                   if ((Tsoil(i)+dTsoil(i))>Tfreezing(kk)) then

                      theta = (S(i)+dy(i)*real(1-var(i)%isat))*(par(i)%thre) + &
                           (par(i)%the - par(i)%thre)
                      deltaJ_latent_T(i) = thetai(i)*dx(i)*rhow*lambdaf
                      tmp1d1(kk) = Tsoil(i)
                      deltaJ_sensible_S(i) = (tmp1d1(kk))*(rhow*cswat*(theta*dx(i)) + &
                           par(i)%rho*par(i)%css*dx(i)) - &
                           (tmp1d1(kk))*(rhow*cswat*(var(i)%thetal*dx(i)) + &
                           par(i)%rho*par(i)%css*dx(i)) - &
                           (tmp1d1(kk))*(rhow*csice*(thetai(i)*dx(i)))

                      if ((i==1).and.(h0(kk)>zero)) then
                         deltaJ_latent_T(i) = deltaJ_latent_T(i) + &
                              thetai(i)*h0(kk)/par(i)%thre*rhow*lambdaf
                         deltaJ_sensible_S(i) = deltaJ_sensible_S(i) + &
                              (tmp1d1(kk))*(rhow*cswat*(h0(kk)+dy(1)*real(1-ns(kk)))) - &
                              (tmp1d1(kk))*(rhow*cswat*(h0(kk)* &
                              (one-thetai(1)/par(1)%thre))) - &
                              (tmp1d1(kk))*(rhow*csice*(h0(kk)* &
                              thetai(1)/par(1)%thre))

                         tmp1d1(kk) = (LHS_h(i)*dt(kk) - (deltaJ_latent_T(i) + deltaJ_sensible_S(i)))/ &
                              (rhow*cswat*(theta*dx(1)+(h0(kk)+dy(1) * &
                              real(1-ns(kk))))+par(1)%rho*par(1)%css*dx(1))
                      else
                         tmp1d1(kk) = (LHS_h(i)*dt(kk) - (deltaJ_latent_T(i) + deltaJ_sensible_S(i)))/ &
                              (dx(i)*(cswat*theta*rhow + par(i)%css*par(i)%rho))
                      endif
                      if (tmp1d1(kk) > dTsoilmax) then ! Check absolute soil temperature change
                         fac(kk) = min(fac(kk),max(half,accel(kk)*abs(dTsoilmax/tmp1d1(kk))))
                         nfac11(kk) = nfac11(kk)+1
                         iok(kk) = 0
                         exit
                      end if
                   endif
                endif
             enddo
          endif

          do i=1, n
             if (abs(dTsoil(i)) > dTsoilmax.and.(var(i)%iice==0)) then ! Check absolute soil temperature change
                fac(kk) = min(fac(kk),max(half,accel(kk)*abs(dTsoilmax/dTsoil(i))))
                nfac5(kk) = nfac5(kk)+1

                iok(kk) = 0
                exit
             end if
          enddo

          if (litter .and. ns(kk)==0) then ! litter and pond
             if (abs(de(0)) > dTLmax) then ! Check absolute litter temperature change
                fac(kk) = min(fac(kk),max(half,accel(kk)*abs(dTLmax/de(0))))
                iok(kk) = 0
             end if
          endif

          ! pond decreasing or runoff starts
!!$          if (iok(kk)==1 .and. ns(kk)<1 .and. h0(kk)+dy(1)<h0min) then ! pond going
!!$             fac(kk) = -(h0(kk)-half*h0min)/dy(1)
!!$             nfac9(kk) = nfac9(kk) +1
!!$             iok(kk) = 0
!!$          end if

!!$           ! pond decreasing or runoff starts
!!$          if (iok(kk)==1 .and. ns(kk)<1 .and.S(1)*par(1)%thre*dx(1) + h0(kk)+dy(1)<0.1*dx(1)) then ! pond going
!!$             fac(kk) = -(S(1)*par(1)%thre*dx(1)+h0(kk)-half*0.1*dx(1))/dy(1)
!!$             nfac9(kk) = nfac9(kk) +1
!!$             iok(kk) = 0
!!$          end if

          ! J0+deltaJ<J(thetai=theta) i.e. too much energy extracted
          J0(1) = rhow*cswat*(Tsoil(1))*dx(1)*var(1)%thetal + &
               rhow*dx(1)*thetai(1)*(csice*(Tsoil(1))-lambdaf) + &
               par(1)%rho*par(1)%css*dx(1)*(Tsoil(1)) + &
               (h0(kk)-hice(kk))*cswat*rhow*(Tsoil(1)) + &
               hice(kk)*rhow*(csice*(Tsoil(1))-lambdaf)
          tmp1d2(kk) = J0(1) + LHS_h(1)*dt(kk)  ! available energy
          tmp1d3(kk) = Tsoil(1)+dTsoil(1)        ! projected temperature
          tmp1d4(kk) = h0(kk) + dy(1)*real(one-ns(kk)) ! projected pond height
          theta      = (S(1)+ dy(1)*real(one-var(1)%isat))*(par(1)%thre) + (par(1)%the - par(1)%thre) ! projected moisture
          c2 = rhow*((tmp1d3(kk))*csice-lambdaf)*(dx(1)*theta+tmp1d4(kk)*theta/par(1)%thre) + & ! projected energy content
               rhow*(tmp1d3(kk))*cswat*tmp1d4(kk)*(one-theta/par(1)%thre) + &
               dx(1)*(tmp1d3(kk))*par(1)%rho*par(1)%css
!!$          if(((tmp1d2(kk))-c2)<zero) then
!!$             fac(kk) = 0.5_r_2
!!$             iok(kk) = 0
!!$             nfac10(kk) = nfac10(kk)+1
!!$          endif

          if (fac(kk) < one) then
             again_ice(kk,1:n) = .false.  ! reset all again_ice if calc is to be repeated with smaller time-step
             imelt(:) = 0
             var(:)%thetai = thetai(:)
             var(:)%dthetaldT = dthetaldT(:)
             do i=1,n
                theta         = S(i)*(par(i)%thre) + (par(i)%the - par(i)%thre)
                var(i)%csoil = par(i)%css*par(i)%rho + rhow*cswat*(theta-var(i)%thetai) + &
                     rhow*csice*var(i)%thetai
                if ((i==1) .and. (ns(kk)==0)) then
                   cp(kk) = real(1-var(1)%iice,r_2)*cswat*rhow & ! heat capacity of pond
                        + real(var(1)%iice,r_2)*rhow* &
                        ((one-var(1)%thetai/par(1)%thre)*cswat + (var(1)%thetai/par(1)%thre)*csice)
                endif
             enddo
             var(:)%csoileff = var(:)%csoil + rhow*lambdaf*var(:)%dthetaldT*real(var(:)%iice,r_2)
             cpeff(kk)= cp(kk) + rhow*lambdaf*var(1)%dthetaldT/par(1)%thre*real(var(1)%iice,r_2)
             h0_tmp(kk) = h0(kk)

          endif

          ! reduce time step
          if (iok(kk)==0) then
             t(kk)      = t(kk)-dt(kk)
             ! dt(kk)     = max(fac(kk)*dt(kk), dtmin)
             dt(kk)     = fac(kk)*dt(kk)
             t(kk)      = t(kk)+dt(kk)
             rsigdt(kk) = one/(sig(kk)*dt(kk))
             nless(kk)  = nless(kk) + 1 ! count step size reductions
          end if
          if (var(1)%isat/=0 .and. iflux(kk)==1 .and. var(1)%phi<phip(kk) .and. &
               var(1)%phi+dy(1)>phip(kk)) then
             ! incipient (onset of) ponding - adjust state of saturated regions
             t(kk)      = t(kk)-dt(kk)
             dt(kk)     = dtmin
             rsigdt(kk) = one/(sig(kk)*dt(kk))
             again(kk)  = .true.
             iok(kk)    = 0
          end if
       end if  ! (.not. again(kk))
       nsteps(kk) = nsteps(kk) + 1

!!$                if ((irec.eq.8992).and.(kk.eq.1) ) then
!!$                   !if ((irec.eq.5).and.(kk.eq.1626)  .and. wlogn == 1011) then
!!$                    write(*,*) 'writing diags', again(kk), nsteps(kk)
!!$
!!$                    ! if (.not. again(kk)) then
!!$ !write(345,"(13i8,1500e16.6)")
!!$                    write(346,"(13i8,1500e16.6)") nsteps(kk), nfac1(kk), nfac2(kk), nfac3(kk), &
!!$                         nfac4(kk), nfac5(kk), nfac6(kk), nfac7(kk), nfac8(kk), nfac9(kk), nfac10(kk), &
!!$                         nfac11(kk), nfac12(kk) , q(:), qsig(:), qH(:), qhsig(:), &
!!$                         dy(0:n), de(0:n), dTsoil(:), S(:),thetai(:), Tsoil(:), &
!!$                         real(var(1)%iice), real(var(1)%isat), &
!!$                         h0(kk), real(iok(kk)), var(1)%phie, var(1)%phi, phip(kk),var(2)%phi, &
!!$                         vsnow(kk)%wcol, &
!!$                         qadv(:), qadvsig(:), qhya(:), qhyb(:), qhTa(:), qhTb(:), &
!!$                         qya(:), qyb(:), qTa(:), qTb(:), &
!!$                         var(1:n)%kH, LHS_h(1:n)*dt(kk), &
!!$                         RHS(1:n)*dt(kk), LHS(1:n)*dt(kk), par(1:n)%thre,dx(1:n), &
!!$                         real(-var(1:n)%isat), dt(kk), real(ns(kk)), vsnow(kk)%tsn(1), vsnow(kk)%hsnow(1)
!!$                  !  endif
!!$                    if (nsteps(kk).gt.1000) STOP
!!$                 endif

       if (nsteps(kk) > nsteps_max) then
          write(wlogn,*) "nsteps > nsteps_max ", irec, kk
          write(wlogn,*) Tsoil(1), S(1)
          write(wlogn,*) nfac1(kk), nfac2(kk), nfac3(kk), nfac4(kk), nfac5(kk), &
               nfac6(kk), nfac7(kk), nfac8(kk), nfac9(kk), nfac10(kk), nfac11(kk), nfac12(kk)
          err = 1
          return
       endif

    end do ! while (iok==0) ----- end get and solve eqns

  END SUBROUTINE get_and_solve_eqn

  !*********************************************************************************************************************

  SUBROUTINE timestep_loop( &
       tfin, irec, mp, qprec, qprec_snow, n, dx, h0, S, thetai, Jsensible, Tsoil, evap, runoff, infil, &
       drainage, discharge, qh, nsteps, vmet, vlit, vsnow, var, T0, Tsurface, Hcum, lEcum, deltaice_cum_T, &
       deltaice_cum_S, Gcum, &
       Qadvcum, Jcol_sensible, Jcol_latent_S, Jcol_latent_T, csoil, kth, phi, dxL, zdelta, SL, Tl, plit, par, wex, ciso, &
       cisoice, ciso_snow, cisoice_snow, cisos, cprec, cprec_snow, qali, qiso_in, qiso_out, qiso_evap_cum, qiso_trans_cum, &
       qiso_liq_adv, qiso_vap_adv, qiso_liq_diff, qiso_vap_diff, qvsig, qlsig, qvTsig, qvh, deltaTa, &
       precip, qevap, qL, qhL, qybL, qTbL, qhTbL, qhybL, rexcol, wcol, ql0, &
       qv0, again, getq0,getqn,init, again_ice, ih0, iok, itmp, ns, nsat, nsatlast, nsteps0, accel, dmax, dt, dwinfil, &
       dwoff, fac, &
       phip, qpme, rsig, rsigdt, sig, t, hint, phimin, qexd, &
       vtmp, deltaS, dTsoil, tmp2d1, &
       tmp2d2, vtop, &
       vbot, v_aquifer, dwcol, dwdrainage, drn,inlit, dwinlit, drexcol, dwdischarge, dJcol_latent_S, dJcol_latent_T, &
       dJcol_sensible, deltaJ_latent_S, deltaJ_latent_T, deltaJ_sensible_S, deltaJ_sensible_T, qevapsig, qrunoff, &
       tmp1d1, tmp1d2, &
       tmp1d3, tmp1d4, deltah0, SL0, deltaSL, cvL0, SLliq0, deltacvL, SLliq, deltaSLliq, qiso_evap, qiso_trans, &
       lE0, G0, Epot, &
       Tfreezing, dtdT, LHS, RHS, LHS_h, RHS_h, surface_case, nns, iflux, litter, i, j, k, kk, condition, littercase, &
       isotopologue, advection, c2, theta, dTqwdTa, dTqwdTb, Tqw, keff, cp, cpeff, hice, deltahice, h0_0, hice_0, &
       h0_tmp, hice_tmp, qtransfer, qmelt_ss , qprec_ss, cprec_ss, delta_snowcol, delta_snowT, delta_snowliq, &
       thetai_0, J0, tmp1, tmp2, iqex, icali, nfac1, nfac2, nfac3, nfac4, nfac5, nfac6, nfac7, nfac8, nfac9, &
       nfac10, nfac11, nfac12, J0snow, wcol0snow, h_ex, wpi, err &
       )
    IMPLICIT NONE
    REAL(r_2)              :: tfin
    INTEGER(i_d)              :: irec, mp
    REAL(r_2),      DIMENSION(1:mp)              :: qprec
    REAL(r_2),      DIMENSION(1:mp)              :: qprec_snow
    INTEGER(i_d)              :: n
    REAL(r_2),      DIMENSION(1:n) :: dx
    REAL(r_2),      DIMENSION(1:mp)           :: h0
    REAL(r_2),      DIMENSION(1:n) :: S
    REAL(r_2),      DIMENSION(1:n) :: thetai
    REAL(r_2),      DIMENSION(1:n) :: Jsensible
    REAL(r_2),      DIMENSION(1:n) :: Tsoil
    REAL(r_2),      DIMENSION(1:mp)           :: evap
    REAL(r_2),      DIMENSION(1:mp)             :: runoff, infil
    REAL(r_2),      DIMENSION(1:mp)             :: drainage, discharge
    REAL(r_2),      DIMENSION(1:mp,-nsnow_max:n)      :: qh
    INTEGER(i_d),   DIMENSION(1:mp)             :: nsteps
    TYPE(vars_met), DIMENSION(1:mp)           :: vmet
    TYPE(vars),     DIMENSION(1:mp)           :: vlit
    TYPE(vars_snow), DIMENSION(1:mp)           :: vsnow
    TYPE(vars),     DIMENSION(1:n) :: var
    REAL(r_2),      DIMENSION(1:mp)           :: T0, Tsurface
    REAL(r_2),      DIMENSION(1:mp)             :: Hcum, lEcum,  deltaice_cum_T, deltaice_cum_S
    REAL(r_2),      DIMENSION(1:mp)             :: Gcum, Qadvcum
    REAL(r_2),      DIMENSION(1:mp)             :: Jcol_sensible, Jcol_latent_S, Jcol_latent_T
    REAL(r_2),      DIMENSION(1:n) :: csoil, kth
    REAL(r_2),      DIMENSION(1:n) :: phi
    REAL(r_2),      DIMENSION(1:mp)              :: dxL
    REAL(r_2),      DIMENSION(1:mp)           :: zdelta, SL, Tl
    TYPE(params),   DIMENSION(1:mp)              :: plit
    TYPE(params),   DIMENSION(1:n) :: par
    REAL(r_2),      DIMENSION(1:mp,1:n), OPTIONAL :: wex
    REAL(r_2),      DIMENSION(1:mp,1:n), OPTIONAL :: ciso
    REAL(r_2),      DIMENSION(1:mp,1:n), OPTIONAL :: cisoice
    REAL(r_2),      DIMENSION(1:mp,1:nsnow_max), OPTIONAL :: ciso_snow
    REAL(r_2),      DIMENSION(1:mp,1:nsnow_max), OPTIONAL :: cisoice_snow
    REAL(r_2),      DIMENSION(1:mp), OPTIONAL :: cisos
    REAL(r_2),      DIMENSION(1:mp),    OPTIONAL :: cprec
    REAL(r_2),      DIMENSION(1:mp),    OPTIONAL :: cprec_snow
    REAL(r_2),      DIMENSION(1:mp),    OPTIONAL :: qali
    REAL(r_2),      DIMENSION(1:mp),   OPTIONAL :: qiso_in, qiso_out
    REAL(r_2),      DIMENSION(1:mp),   OPTIONAL :: qiso_evap_cum, qiso_trans_cum
    REAL(r_2),      DIMENSION(1:mp,-nsnow_max+1:n),   OPTIONAL :: qiso_liq_adv, qiso_vap_adv
    REAL(r_2),      DIMENSION(1:mp,-nsnow_max+1:n-1),   OPTIONAL :: qiso_liq_diff, qiso_vap_diff
    REAL(r_2),      DIMENSION(1:mp,-nsnow_max:n),   OPTIONAL :: qvsig, qlsig, qvTsig, qvh
    REAL(r_2),      DIMENSION(1:mp), OPTIONAL :: deltaTa
    ! Error flag if nstep of SLI > nsteps_max: err=0 -> no error; err/=0 -> error
    INTEGER(i_d),     INTENT(INOUT), OPTIONAL :: err

    ! Solves the RE and, optionally, the ADE from time ts to tfin.
    ! Definitions of arguments:
    ! Required args:
    ! ts   - start time (h).
    ! tfin   - finish time.
    ! qprec   - precipitation (or water input) rate (fluxes are in cm/h).
    ! qevap   - potl evaporation rate from soil surface.
    ! n    - no. of soil layers.
    ! nsol   - no. of solutes.
    ! dx(1:n) - layer thicknesses.
    ! h0   - surface head, equal to depth of surface pond.
    ! S(1:n)  - degree of saturation ("effective satn") of layers.
    ! evap   - cumulative evaporation from soil surface (cm, not initialised).
    ! runoff  - cumulative runoff.
    ! infil   - cumulative net infiltration (time integral of flux across surface).
    ! drn   - cumulative net drainage (time integral of flux across bottom).
    ! nsteps  - cumulative no. of time steps for RE soln.
    ! Optional args:
    ! heads(1:n)   - matric heads h of layers at finish.
    ! qexsub    - subroutine to get layer water extraction rates (cm/h) by
    !     plants. Note that there is no solute extraction and osmotic
    !     effects due to solute are ignored. Arguments:
    !     qex(1:n) - layer extraction rates; qexh(1:n) - partial
    !     derivs of qex wrt h.
    ! wex(1:n)    - cumulative water extraction from layers.
    ! cin(1:nsol)   - solute concns in water input (user's units/cc).
    ! c0(1:nsol)   - solute concns in surface pond.
    ! sm(1:n,1:nsol)  - solute (mass) concns in layers.
    ! soff(1:nsol)   - cumulative solute runoff (user's units).
    ! sinfil(1:nsol)  - cumulative solute infiltration.
    ! sdrn(1:nsol)   - cumulative solute drainage.
    ! nssteps(1:nsol) - cumulative no. of time steps for ADE soln.
    ! isosub    - subroutine to get adsorbed solute (units/g soil) from concn
    !     in soil water according to chosen isotherm code.
    !     Arguments: iso - 2 character code; c - concn in soil water;
    !     p(:) - isotherm parameters; f - adsorbed mass/g soil;
    !     fc - deriv of f wrt c (slope of isotherm curve). Note that
    !     linear adsorption does not require a sub, and other types
    !     are available in sub isosub.

    REAL(r_2),    DIMENSION(1:mp)       :: precip, qevap
    REAL(r_2),    DIMENSION(1:mp)       :: qL, qhL, qybL, qTbL, qhTbL, qhybL, rexcol, wcol, ql0, qv0
    LOGICAL,      DIMENSION(1:mp)       :: again, getq0,getqn,init
    LOGICAL,      DIMENSION(1:mp,1:n)   :: again_ice
    INTEGER(i_d), DIMENSION(1:mp)       :: ih0, iok, itmp, ns, nsat, nsatlast, nsteps0
    REAL(r_2),    DIMENSION(1:mp)       :: accel, dmax, dt, dwinfil, dwoff, fac, phip
    REAL(r_2),    DIMENSION(1:mp)       :: qpme, rsig, rsigdt, sig, t
    REAL(r_2),    DIMENSION(1:n) :: hint, phimin, qexd
    REAL(r_2),    DIMENSION(-nsnow_max+1:n) :: aa, bb, cc, dd, ee, ff, gg, dy
    REAL(r_2),    DIMENSION(-nsnow_max+1:n) :: aah, bbh, cch, ddh, eeh, ffh, ggh, de
    REAL(r_2),    DIMENSION(-nsnow_max:n) :: q, qya, qyb, qTa, qTb,qhya, qhyb, qhTa, qhTb
    REAL(r_2),    DIMENSION(-nsnow_max:n) :: qadv, qadvya, qadvyb, qadvTa, qadvTb


    TYPE(vars)                          :: vtmp
    !TYPE(vars),   DIMENSION(1:n) :: var
    REAL(r_2),    DIMENSION(-nsnow_max:n) :: qsig, qhsig, qadvsig
    REAL(r_2),    DIMENSION(-nsnow_max:n) :: qliq, qv, qvT, qlya, qlyb, qvya, qvyb, qlTb, qvTa, qvTb
    REAL(r_2),    DIMENSION(1:n) :: deltaS, dTsoil
    REAL(r_2),    DIMENSION(0:n) :: tmp2d1, tmp2d2
    REAL(r_2),    DIMENSION(1:n) :: S0, Sliq0, Sliq, deltaSliq, cv0, deltacv
    REAL(r_2),    DIMENSION(1:n) :: Sliqice0, Sliqice, deltaSliqice
    REAL(r_2),    DIMENSION(1:n) :: Sice0, Sice, deltaSice

    REAL(r_2),    DIMENSION(-nsnow_max+1:n) :: Sliq0_ss, Sliq_ss, deltaSliq_ss
    REAL(r_2),    DIMENSION(-nsnow_max+1:n) :: Sliqice0_ss, Sliqice_ss, deltaSliqice_ss
    REAL(r_2),    DIMENSION(-nsnow_max+1:n) :: Sice0_ss, Sice_ss, deltaSice_ss
    REAL(r_2),    DIMENSION(-nsnow_max+1:n) :: S0_ss, S_ss, Tsoil_ss, dTsoil_ss
    REAL(r_2),    DIMENSION(-nsnow_max+1:n) :: cv0_ss, cv_ss, Dv_ss, deltacv_ss, dx_ss
    REAL(r_2),    DIMENSION(-nsnow_max+1:n-1) :: dz_ss
    REAL(r_2),      DIMENSION(1:nsnow_max) :: cisoliqice_snow
    INTEGER(i_d) :: itop ! integer corresponding to top of soil-snow column
    INTEGER(i_d) :: nsnow ! number of dedicated snow layers
    REAL(r_2),    DIMENSION(1:n) :: tmp_thetasat, tmp_thetar
    REAL(r_2),    DIMENSION(-nsnow_max+1:n) :: thetasat_ss, thetar_ss
    REAL(r_2),    DIMENSION(-nsnow_max+1:n) :: tmp_tortuosity
    REAL(r_2),    DIMENSION(-nsnow_max+1:n) :: ciso_ss, cisoice_ss
    REAL(r_2),    DIMENSION(1:n)   :: delthetai, dthetaldT, thetal
    INTEGER(i_d), DIMENSION(1:n) :: isave, nsteps_ice, imelt


    TYPE(vars),         DIMENSION(1:mp) :: vtop, vbot
    TYPE(vars_aquifer), DIMENSION(1:mp) :: v_aquifer
    REAL(r_2),          DIMENSION(1:mp) :: dwcol, dwdrainage, drn,inlit, dwinlit, drexcol, dwdischarge
    REAL(r_2),          DIMENSION(1:mp) :: dJcol_latent_S, dJcol_latent_T, dJcol_sensible
    REAL(r_2),          DIMENSION(1:n) :: deltaJ_latent_S, deltaJ_latent_T, deltaJ_sensible_S, deltaJ_sensible_T
    REAL(r_2),          DIMENSION(1:mp) :: qevapsig
    REAL(r_2),          DIMENSION(1:mp) :: qrunoff
    REAL(r_2),          DIMENSION(1:mp) :: tmp1d1, tmp1d2, tmp1d3,  tmp1d4
    REAL(r_2),          DIMENSION(1:mp) :: deltah0
    REAL(r_2),          DIMENSION(1:mp) :: SL0, deltaSL, cvL0, SLliq0, deltacvL, SLliq, deltaSLliq
    REAL(r_2),          DIMENSION(1:mp) :: qiso_evap, qiso_trans
    REAL(r_2),          DIMENSION(1:mp) :: lE0, G0, Epot
    REAL(r_2),          DIMENSION(1:mp) :: Tfreezing
    REAL(r_2),          DIMENSION(1:mp) :: dtdT
    REAL(r_2),          DIMENSION(-nsnow_max+1:n) :: LHS, RHS, LHS_h, RHS_h
    INTEGER(i_d),       DIMENSION(1:mp) :: surface_case

    INTEGER(i_d),       DIMENSION(1:mp) :: nns, iflux
    LOGICAL      :: litter
    INTEGER(i_d) :: i, j, k, kk, condition
    INTEGER(i_d) :: littercase, isotopologue, advection ! switches
    REAL(r_2)    :: c2, theta
    REAL(r_2)    :: dTqwdTa, dTqwdTb, Tqw, keff
    REAL(r_2),          DIMENSION(1:mp) :: cp, cpeff, hice, deltahice, h0_0, hice_0, h0_tmp, hice_tmp
    REAL(r_2),          DIMENSION(nsnow_max) :: qmelt
    ! tranfer of water  from soil to   snow (+ve) or soil to snow (-ve) at init or termination of snowpack
    REAL(r_2),  DIMENSION(1:mp) :: qtransfer
    REAL(r_2),  DIMENSION(1:mp) :: qmelt_ss ! melt water from bottom snow layer to soil
    REAL(r_2),  DIMENSION(1:mp) :: qprec_ss
    REAL(r_2),  DIMENSION(1:mp) :: cprec_ss
    REAL(r_2),          DIMENSION(nsnow_max) :: delta_snowcol, delta_snowT, delta_snowliq
    REAL(r_2),      DIMENSION(1:n) :: thetai_0, J0
    REAL(r_2) :: tmp1, tmp2
    REAL(r_2),          DIMENSION(1:n) :: iqex
    REAL(r_2),          DIMENSION(1:mp)     :: icali
    INTEGER(i_d),       DIMENSION(1:mp) :: nfac1, nfac2, nfac3, nfac4, nfac5, &
         nfac6, nfac7, nfac8, nfac9, nfac10, nfac11, nfac12
    REAL(r_2),          DIMENSION(1:mp)     :: J0snow, wcol0snow
    REAL(r_2), DIMENSION(1:n) :: h_ex
    REAL(r_2) :: wpi


    aa(:)     = zero
    aah(:)    = zero
    bb(:)     = zero
    bbh(:)    = zero
    cc(:)     = zero
    cch(:)    = zero
    dd(:)     = zero
    ddh(:)    = zero
    ee(:)     = zero
    eeh(:)    = zero
    ff(:)     = zero
    ffh(:)    = zero
    gg(:)     = zero
    ggh(:)    = zero
    dy(:)     = zero
    q(:)      = zero
    qya(:)    = zero
    qyb(:)    = zero
    qTa(:)    = zero
    qTb(:)    = zero
    qhya(:)   = zero
    qhyb(:)   = zero
    qhTa(:)   = zero
    qhTb(:)   = zero
    qadvyb(:) = zero
    qadvya(:) = zero
    de(:)            = zero
    ! initialise diagnostic vars for input to isotope_vap
    Sliq0_ss  = zero
    Sliq_ss= zero
    deltaSliq_ss= zero
    Sliqice0_ss= zero
    Sliqice_ss= zero
    deltaSliqice_ss= zero
    Sice0_ss= zero
    Sice_ss= zero
    deltaSice_ss= zero
    S0_ss= zero
    S_ss= zero
    cv0_ss= zero
    cv_ss= zero
    Dv_ss = zero
    deltacv_ss= zero

    do while (t(kk) < tfin)
       !----- take next time step
       iflux(kk)=1
       again(kk)  = .true. ! flag for recalcn of fluxes (default=false)
       imelt = 0 ! initialise imelt (==1 at onset of melting)
       qmelt(:) = zero
       vsnow(kk)%nsnow_last = vsnow(kk)%nsnow ! for detecting onset and disappearance of dedicated snow layer
       S0(1:n) = S(1:n)
       Sice0(1:n)    = var(1:n)%thetai/par(1:n)%thre
       Sliqice0(1:n) = (S(1:n) - var(1:n)%cv)/(one-var(1:n)%cv)
       Sliq0(1:n)    = Sliqice0(1:n) - Sice0(1:n)
       Sice0(1) = Sice0(1)  + hice_0(kk)/(dx(1)*par(1)%thre)
       Sliq0(1) = Sliq0(1) + (h0_0(kk)-hice_0(kk))/(dx(1)*par(1)%thre) ! add pond component to Sliq(1)

       if (nsnow_max.gt.0) S0_ss(-nsnow_max+1:0) = vsnow(kk)%hsnow      ! divide by zero after changes to layer depth
       if (nsnow_max.gt.0) Sliq0_ss(-nsnow_max+1:0) = vsnow(kk)%hliq ! divide by zero after changes to layer depth
       ! divide by zero after changes to layer depth
       if (nsnow_max.gt.0) Sice0_ss(-nsnow_max+1:0) = vsnow(kk)%hsnow - vsnow(kk)%hliq
       if (nsnow_max.gt.0) cv0_ss(-nsnow_max+1:0) = vsnow(kk)%cv
       if (nsnow_max.gt.0) dx_ss(-nsnow_max+1:0) = vsnow(kk)%depth

       S0_ss(1:n) = S(1:n)
       Sliqice0_ss(1:n) = Sliqice0(1:n)
       Sliq0_ss(1:n) = Sliq0(1:n)
       cv0_ss(1:n) =  var(1:n)%cv
       dx_ss(1:n) = dx(1:n)

       CALL iflux_loop(&
            tfin, irec, mp, qprec, qprec_snow, n, dx(:), h0, S(:), thetai(:), &
            Jsensible(:), Tsoil(:), evap, infil, drainage, discharge, &
            qh, nsteps, vmet, vlit, vsnow, var(:), T0, Tsurface, Hcum, lEcum, &
            Gcum, Qadvcum, Jcol_sensible, &
            Jcol_latent_S, Jcol_latent_T, csoil(:), kth(:), phi(:), dxL, zdelta, SL, Tl, &
            plit, par(:), wex, ciso_snow, cisoice_snow, &
            qali, &
            qvsig, qlsig, qvTsig, qvh, deltaTa, &
            precip, qevap, qL, qhL, qybL, qTbL, qhTbL, qhybL, rexcol, wcol, &
            again, getq0,getqn,init, again_ice, ih0, iok, itmp, ns, nsat, &
            nsatlast, nsteps0, accel, dmax, dt, dwinfil, dwoff, fac, &
            phip, qpme, rsig, rsigdt, sig, t, hint(:), phimin(:), &
            qexd(:), aa(:), bb(:), cc(:), dd(:), ee(:), ff(:), gg(:), dy(:), &
            aah(:), bbh(:), cch(:), ddh(:), eeh(:), ffh(:), ggh(:), &
            de(:), q(:), qya(:), qyb(:), qTa(:), qTb(:),qhya(:), qhyb(:), qhTa(:), qhTb(:), qadv(:), qadvya(:), qadvyb(:), &
            qadvTa(:), qadvTb(:), vtmp, qsig(:), qhsig(:), qadvsig(:), qliq(:), qv(:), qvT(:), qlya(:), qlyb(:), &
            qvya(:), qvyb(:), qlTb(:), qvTa(:), qvTb(:), deltaS(:), dTsoil(:), tmp2d1(:), tmp2d2(:), &
            cv0(:), &
            cisoliqice_snow, &
            dthetaldT(:), thetal(:), isave, nsteps_ice, imelt, vtop, vbot, v_aquifer, &
            dwcol, dwdrainage, drn,inlit, dwinlit, drexcol, dwdischarge, &
            dJcol_latent_S, dJcol_latent_T, dJcol_sensible, deltaJ_latent_S(:), &
            deltaJ_latent_T(:), deltaJ_sensible_S(:), deltaJ_sensible_T(:), qevapsig, &
            qrunoff, tmp1d1, tmp1d2, tmp1d3,  tmp1d4, deltah0, deltaSL, &
            lE0, G0, &
            Epot, Tfreezing, dtdT, LHS, RHS, LHS_h, RHS_h, surface_case, nns, &
            iflux, litter, i, j, k, kk, condition, littercase, isotopologue, &
            advection, c2, theta, dTqwdTa, dTqwdTb, Tqw, keff, cp, &
            cpeff, hice, h0_0, hice_0, h0_tmp, hice_tmp, qmelt, &
            qtransfer, delta_snowcol, delta_snowT, &
            delta_snowliq, thetai_0(:), J0(:), tmp1, tmp2, iqex(:), &
            nfac1, nfac2, nfac3, nfac4, nfac5, nfac6, nfac7, nfac8, nfac9, &
            nfac10, nfac11, nfac12, J0snow, wcol0snow, h_ex, wpi, err &
            )
       if (err /= 0) return

       runoff(kk) = runoff(kk) + qrunoff(kk)*dt(kk)

       where ((S0(:).ge.one).or.(S(:).ge.one))
          var(1:n)%cv = zero
          cv0(1:n) = zero
       endwhere
       if (h0(kk).gt.zero.or.h0_0(kk).gt.zero) then
          var(1)%cv = zero
          cv0(1) = zero
       endif

       Sliqice0(1:n) = (S0(1:n) - cv0(1:n))/(one-cv0(1:n))
       deltacv(1:n)   = var(1:n)%cv - cv0(1:n)
       deltah0(kk) = h0(kk) - h0_0(kk)
       Sliqice(1:n)      = (S(1:n) - var(1:n)%cv)/(one-var(1:n)%cv)
       deltaSliqice(1:n) = Sliqice(1:n) - Sliqice0(1:n)
       ! add pond component to Sliqice(1) and deltaSliqice(1)
       Sliqice(1) = Sliqice(1) + h0(kk)/(dx(1)*par(1)%thre) ! add pond component to Sliq(1)
       deltaSliqice(1) =  deltaSliqice(1) + deltah0(kk)/(dx(1)*par(1)%thre)

       delthetai(1:n) = (var(1:n)%thetai-thetai_0(1:n))

       hice(kk) = h0(kk)*var(1)%thetai/par(1)%thre
       deltahice(kk) =  hice(kk) - hice_0(kk)
       deltahice(kk) =  hice(kk) - h0_0(kk)*thetai_0(1)/par(1)%thre

       ! change in snow pack (lumped with top layer)
       if (vsnow(kk)%nsnow==0) then
          if (var(1)%iice==1) then
             vsnow(kk)%wcol = vsnow(kk)%wcol + &
                  min(qprec_snow(kk)*dt(kk), max(zero,(deltahice(kk)+dx(1)*delthetai(1)))) + & ! accumulation
                  max(-vsnow(kk)%wcol,min(zero,(deltahice(kk)+dx(1)*delthetai(1)))) ! melting

             if (h0(kk)>zero) then
                vsnow(kk)%wcol = min(vsnow(kk)%wcol, hice(kk))
             else
                vsnow(kk)%wcol = min(vsnow(kk)%wcol, dx(1)*var(1)%thetai)
             endif
          else
             vsnow(kk)%wcol = zero
          endif

          vsnow(kk)%tsn = merge(Tsoil(1),zero,vsnow(kk)%wcol>zero)
          if (vsnow(kk)%wcol>zero) then
             vsnow(kk)%depth(1) = vsnow(kk)%wcol/(vsnow(kk)%dens(1)/rhow)
          else
             vsnow(kk)%depth(1) = zero
          endif
       endif

       Sice(1:n) = var(1:n)%thetai/par(1:n)%thre
       deltaSice(1:n) = delthetai(1:n)/par(1:n)%thre
       deltahice(kk) = hice(kk)-hice_0(kk)
       Sice(1) = Sice(1)  + hice(kk)/(dx(1)*par(1)%thre)
       deltaSice(1) = deltaSice(1) + deltahice(kk)/(dx(1)*par(1)%thre)

       Sliq(1:n) = Sliqice(1:n) - Sice(1:n)
       deltaSliq(1:n) = deltaSliqice(1:n) - deltaSice(1:n)


       thetai_0(1:n)    = var(1:n)%thetai
       thetai(1:n) = var(1:n)%thetai
       dthetaldt(1:n) = var(1:n)%dthetaldT
       ! h0_0(kk) = h0(kk)
       ! hice_0(kk) = hice(kk)
       ! Sliq0(1:n) = Sliq(1:n)
       ! Sice0(1:n)  = Sice(1:n)
       ! Sliqice0(1:n) = Sliqice(1:n)


       ! extension of diagnostics for isotope_vap to snow

       S_ss(1:n) = S(1:n)
       S_ss(1) = S_ss(1) + h0(kk)/(dx(1)*par(1)%thre)
       S0_ss(1) = S0_ss(1) + (h0(kk)-deltah0(kk))/(dx(1)*par(1)%thre)
       Sliqice_ss(1:n) = Sliqice(1:n)
       Sliq_ss(1:n) = Sliq(1:n)
       cv_ss(1:n) =  var(1:n)%cv
       Dv_ss(1:n) =  var(1:n)%Dv
       Sice_ss(1:n) =  Sice(1:n)

       deltaSice_ss(1:n) = deltaSice(1:n)
       deltaSliq_ss(1:n) =  deltaSliq(1:n)
       deltaSliqice_ss(1:n) = deltaSliqice(1:n)
       deltacv_ss(1:n) =  deltacv(1:n)

       Tsoil_ss(1:n) = Tsoil(1:n)
       dTsoil_ss(1:n) = dTsoil(1:n)
       if (isotopologue /= 0) then
          ciso_ss(1:n) = ciso(kk,1:n)
          cisoice_ss(1:n) = cisoice(kk,1:n)
       endif


       if (nsnow_max.gt.0) then

          dx_ss(-nsnow_max+1:0) = vsnow(kk)%depth       ! fixed snow layer depth for isotope calc
          tmp1d1(kk) = S0_ss(0) ! initial SWE
          where ( dx_ss(-nsnow_max+1:0).gt.zero)


             cv_ss(-nsnow_max+1:0)   = vsnow(kk)%cv
             !   vsnow(kk)%cv = zero  ! test vh
             !   cv_ss(-nsnow_max+1:0)   = zero
             !   cv0_ss(-nsnow_max+1:0)   = zero
             deltacv_ss(-nsnow_max+1:0)   = cv_ss(-nsnow_max+1:0) - cv0_ss(-nsnow_max+1:0)
             ! deltaSliqice_ss(-nsnow_max+1:0) =   (vsnow(kk)%hsnow - S0_ss(-nsnow_max+1:0)- &
             !      dx_ss(-nsnow_max+1:0)*deltacv_ss(-nsnow_max+1:0))/ dx_ss(-nsnow_max+1:0)

             Sliqice0_ss(-nsnow_max+1:0) = (S0_ss(-nsnow_max+1:0)/dx_ss(-nsnow_max+1:0) - &
                  cv0_ss(-nsnow_max+1:0))/(one-cv0_ss(-nsnow_max+1:0))
             Sliqice_ss(-nsnow_max+1:0) = (vsnow(kk)%hsnow(1:nsnow_max)/dx_ss(-nsnow_max+1:0) - &
                  cv_ss(-nsnow_max+1:0))/(one-cv_ss(-nsnow_max+1:0))
             deltaSliqice_ss(-nsnow_max+1:0) = Sliqice_ss(-nsnow_max+1:0) - Sliqice0_ss(-nsnow_max+1:0)


             S0_ss(-nsnow_max+1:0)    = S0_ss(-nsnow_max+1:0)/dx_ss(-nsnow_max+1:0)
             S_ss(-nsnow_max+1:0) =  vsnow(kk)%hsnow/dx_ss(-nsnow_max+1:0)


             Sliq0_ss(-nsnow_max+1:0)    = Sliq0_ss(-nsnow_max+1:0)/dx_ss(-nsnow_max+1:0)
             Sliq_ss(-nsnow_max+1:0) =  vsnow(kk)%hliq/dx_ss(-nsnow_max+1:0)


             Sice0_ss(-nsnow_max+1:0) = Sliqice0_ss(-nsnow_max+1:0) - Sliq0_ss(-nsnow_max+1:0)
             Sice_ss(-nsnow_max+1:0) = Sliqice_ss(-nsnow_max+1:0) - Sliq_ss(-nsnow_max+1:0)


             deltaSliq_ss(-nsnow_max+1:0) = Sliq_ss(-nsnow_max+1:0) - Sliq0_ss(-nsnow_max+1:0)
             deltaSice_ss(-nsnow_max+1:0) = Sice_ss(-nsnow_max+1:0) - Sice0_ss(-nsnow_max+1:0)

             cv_ss(-nsnow_max+1:0)   = vsnow(kk)%cv
             Dv_ss(-nsnow_max+1:0)   = vsnow(kk)%Dv

             Tsoil_ss(-nsnow_max+1:0) =  vsnow(kk)%tsn(1:nsnow_max)
             dTsoil_ss(-nsnow_max+1:0) = delta_snowT(1:nsnow_max)
          elsewhere
             S_ss(-nsnow_max+1:0) = zero
             Sliq_ss(-nsnow_max+1:0) = zero
             Sliqice_ss(-nsnow_max+1:0) = zero
             Sice_ss(-nsnow_max+1:0) = zero
             deltaSliq_ss(-nsnow_max+1:0) = zero
             deltaSliqice_ss(-nsnow_max+1:0) =  zero
             deltaSice_ss(-nsnow_max+1:0) =zero
             deltacv_ss(-nsnow_max+1:0) = zero
             cv_ss(-nsnow_max+1:0)   = zero
             Dv_ss(-nsnow_max+1:0)   = zero
             Tsoil_ss(-nsnow_max+1:0) = zero
             dTsoil_ss(-nsnow_max+1:0) = zero
          endwhere

          if (isotopologue /= 0) then
             where ( dx_ss(-nsnow_max+1:0) .gt. zero)
                where (vsnow(kk)%hliq(1:nsnow_max).gt.zero)
                   ciso_ss(-nsnow_max+1:0) = ciso_snow(kk,1:nsnow_max)
                elsewhere
                   ciso_ss(-nsnow_max+1:0) = cisoice_snow(kk,1:nsnow_max)
                endwhere
                cisoice_ss(-nsnow_max+1:0) =  cisoice_snow(kk,1:nsnow_max)

             elsewhere
                ciso_ss(-nsnow_max+1:0) = zero
                cisoice_ss(-nsnow_max+1:0) = zero
             endwhere
             if (vsnow(kk)%hsnow(1).gt.zero) then

                cisoliqice_snow(1:nsnow_max) = (ciso_snow(kk,1:nsnow_max)*vsnow(kk)%hliq(1:nsnow_max) + &
                     cisoice_snow(kk,1:nsnow_max)*(vsnow(kk)%hsnow(1:nsnow_max)-vsnow(kk)%hliq(1:nsnow_max))) / &
                     vsnow(kk)%hsnow(1:nsnow_max)

                cisos(kk) =  cisoliqice_snow(1)

             endif
          endif

       endif

       if (littercase==1) then  !!!! vh needs attention !!!!
          SL0(kk)    = SL(kk) - deltaSL(kk)
          cvL0(kk)   = vlit(kk)%cv
          SLliq0(kk) = (SL0(kk) - cvL0(kk))/(one-cvL0(kk))
          call litter_props(Sl(kk), Tl(kk), vlit(kk), plit(kk), h0(kk))
          deltacvL(kk)   = vlit(kk)%cv - cvL0(kk)
          SLliq(kk)      = (SL(kk) - vlit(kk)%cv)/(one-vlit(kk)%cv)
          deltaSLliq(kk) = SLliq(kk) - SLliq0(kk)
       endif

       Jsensible(1) = (var(1)%csoil* dx(1)+h0(kk)*cswat)*(Tsoil(1))
       Jsensible(2:n) = var(2:n)%csoil*(Tsoil(2:n))* dx(2:n)

       ! change in heat stored in soil column
       dJcol_latent_S(kk) = sum(deltaJ_latent_S(1:n))
       Jcol_latent_S(kk) = Jcol_latent_S(kk) + dJcol_latent_S(kk)
       dJcol_latent_T(kk) = sum(deltaJ_latent_T(1:n))
       Jcol_latent_T(kk) = Jcol_latent_T(kk) + dJcol_latent_T(kk)
       dJcol_sensible(kk) = sum(deltaJ_sensible_T(1:n)) + sum(deltaJ_sensible_S(1:n))
       Jcol_sensible(kk) = Jcol_sensible(kk) + dJcol_sensible(kk)

       deltaice_cum_S(kk) = -Jcol_latent_S(kk)/rhow/lambdaf
       deltaice_cum_T(kk) = -Jcol_latent_T(kk)/rhow/lambdaf

       deltaTa(kk)       = zero
       init(kk)          = .false.

       if (isotopologue /= 0) then

          tmp_thetasat(1:n)   = par(1:n)%thre
          tmp_tortuosity(1:n) = par(1:n)%tortuosity
          tmp_thetar(1:n)     = par(1:n)%the - par(1:n)%thre

          nns(kk) = 1  ! pond included in top soil layer
          if (vsnow(kk)%nsnow>0) then
             nns(kk) = 1-vsnow(kk)%nsnow
          endif

          itop = -nsnow_max+1
          ql0(kk) = qlsig(kk,-vsnow(kk)%nsnow)
          qv0(kk) = qvsig(kk,-vsnow(kk)%nsnow)
          nsnow = vsnow(kk)%nsnow

          dz_ss(itop:n-1) =   half*(dx_ss(itop:n-1)+dx_ss(itop+1:n)) ! flow paths
          thetasat_ss(itop:0) = 1.0_r_2
          thetar_ss(itop:0) = 1.0_r_2
          thetasat_ss(1:n) = tmp_thetasat(1:n)
          thetar_ss(1:n) = tmp_thetar(1:n)
          qprec_ss(kk) = qprec(kk)
          cprec_ss(kk) = cprec(kk)
          qmelt_ss(kk) = zero
          if (vsnow(kk)%nsnow.gt.0) then
             qmelt_ss(kk) = qmelt(vsnow(kk)%nsnow)/dt(kk)
             tmp_tortuosity(itop:0) = 0.8_r_2
             thetasat_ss(itop:0) = 1.0_r_2
             thetar_ss(itop:0) = 0.0_r_2
          elseif (vsnow(kk)%nsnow_last.gt.0.and.vsnow(kk)%nsnow.eq.0) then
             !qmelt_ss(kk) = -qtransfer(kk)/dt(kk)   ! transfer from terminated snowpack to soil column
             if (abs(qprec(kk)-qtransfer(kk)/dt(kk)).gt.zero) then
                cprec_ss(kk) = (cprec(kk)*qprec(kk)-cisoliqice_snow(1)*qtransfer(kk)/dt(kk))/(qprec(kk)-qtransfer(kk)/dt(kk))
                qprec_ss(kk) = qprec(kk)  -qtransfer(kk)/dt(kk)
             endif
             qtransfer(kk) = zero ! otherwise qtransfer retains value for transfer from soil to new snowpack
          endif

          call isotope_vap(irec,isotopologue, n, nsnow, vsnow(kk)%nsnow_last, &
               itop, dx_ss(itop:n), dz_ss(itop:n-1), sig(kk), dt(kk), &
               Tsoil_ss(itop:n), dTsoil_ss(itop:n), Sliqice_ss(itop:n), deltaSliqice_ss(itop:n), &
               Sliq_ss(itop:n), deltaSliq_ss(itop:n),Sice_ss(itop:n), deltaSice_ss(itop:n), &
               Tsurface(kk), vmet(kk)%Ta, &
               qsig(itop-1:n), qlsig(kk,itop-1:n), qvsig(kk,itop-1:n), &
               qmelt_ss(kk),qtransfer(kk)/dt(kk), &  ! melt water to soil and water transfer from soil to snow (+ve)
               qprec_ss(kk),qprec_snow(kk), qevapsig(kk), qrunoff(kk), iqex(1:n), &
               cv_ss(itop:n),Dv_ss(itop:n), thetasat_ss(itop:n), thetar_ss(itop:n), tmp_tortuosity(itop:n), &
               deltacv_ss(itop:n), vmet(kk)%rbw, vmet(kk)%cva, vmet(kk)%civa, &
               cprec_ss(kk),cprec_snow(kk), icali(kk), &
               ql0(kk), qv0(kk), &
               ciso_ss(itop:n), cisoice_ss(itop:n), cisos(kk), &
               qiso_in(kk), qiso_out(kk), qiso_evap(kk), qiso_trans(kk), &
               qiso_liq_adv(kk,itop:n), qiso_vap_adv(kk,itop:n), qiso_liq_diff(kk,itop:n-1), qiso_vap_diff(kk,itop:n-1))



          qiso_evap_cum(kk)  = qiso_evap_cum(kk)  + qiso_evap(kk)*dt(kk)
          qiso_trans_cum(kk) = qiso_trans_cum(kk) + qiso_trans(kk)*dt(kk)
          ciso(kk,1:n) = ciso_ss(1:n)
          cisoice(kk,1:n) = cisoice_ss(1:n)
          ciso_snow(kk,1:nsnow_max) = ciso_ss(itop:0)
          cisoice_snow(kk,1:nsnow_max) = cisoice_ss(itop:0)
          if (vsnow(kk)%nsnow.gt.1) then
             cisoliqice_snow(1:nsnow_max) = (ciso_snow(kk,1:nsnow_max)*vsnow(kk)%hliq(1:nsnow_max) + &
                  cisoice_snow(kk,1:nsnow_max)*(vsnow(kk)%hsnow(1:nsnow_max)-vsnow(kk)%hliq(1:nsnow_max))) / &
                  vsnow(kk)%hsnow(1:nsnow_max)
          endif

       endif ! isotopologue/=0

    end do ! while (t<tfin)
  END SUBROUTINE timestep_loop

  SUBROUTINE solve(ts, tfin, irec, mp, qprec, qprec_snow, n, dx, h0, S, thetai, Jsensible, Tsoil, evap, evap_pot, runoff, &
       infil, drainage, discharge, qh, nsteps, vmet, vlit, vsnow, var, csoil, kth, phi, T0, Tsurface, Hcum, lEcum, &
       Gcum, Qadvcum, Jcol_sensible, Jcol_latent_S, Jcol_latent_T, deltaice_cum_T, deltaice_cum_S, dxL, zdelta, &
       SL, TL, plit, par, qex, wex, heads,  &
       ciso, cisoice, ciso_snow, cisoice_snow, cisos, cprec, cprec_snow, cali, &
       qali, qiso_in, qiso_out, qiso_evap_cum, qiso_trans_cum, qiso_liq_adv, &
       qiso_vap_adv, qiso_liq_diff, qiso_vap_diff, qvsig, qlsig, qvTsig, qvh, deltaTa, lE_old, &
       dolitter, doisotopologue, docondition, doadvection, err)

    IMPLICIT NONE

    REAL(r_2),                             INTENT(IN)              :: ts, tfin
    INTEGER(i_d),                          INTENT(IN)              :: irec, mp
    REAL(r_2),      DIMENSION(1:mp),       INTENT(IN)              :: qprec
    REAL(r_2),      DIMENSION(1:mp),       INTENT(INOUT)           :: qprec_snow
    INTEGER(i_d),                          INTENT(IN)              :: n
    REAL(r_2),      DIMENSION(1:mp,1:n),   INTENT(IN)              :: dx
    REAL(r_2),      DIMENSION(1:mp),       INTENT(INOUT)           :: h0
    REAL(r_2),      DIMENSION(1:mp,1:n),   INTENT(INOUT)           :: S
    REAL(r_2),      DIMENSION(1:mp,1:n),   INTENT(OUT)             :: thetai
    REAL(r_2),      DIMENSION(1:mp,1:n),   INTENT(OUT)             :: Jsensible
    REAL(r_2),      DIMENSION(1:mp,1:n),   INTENT(INOUT)           :: Tsoil
    REAL(r_2),      DIMENSION(1:mp),       INTENT(INOUT)           :: evap
    REAL(r_2),      DIMENSION(1:mp),       INTENT(OUT)             :: evap_pot, runoff, infil
    REAL(r_2),      DIMENSION(1:mp),       INTENT(OUT)             :: drainage, discharge
    REAL(r_2),      DIMENSION(1:mp,-nsnow_max:n), INTENT(OUT)      :: qh
    INTEGER(i_d),   DIMENSION(1:mp),       INTENT(OUT)             :: nsteps
    TYPE(vars_met), DIMENSION(1:mp),       INTENT(INOUT)           :: vmet
    TYPE(vars),     DIMENSION(1:mp),       INTENT(INOUT)           :: vlit
    TYPE(vars_snow), DIMENSION(1:mp),      INTENT(INOUT)           :: vsnow
    TYPE(vars),     DIMENSION(1:mp,1:n),   INTENT(INOUT)           :: var
    REAL(r_2),      DIMENSION(1:mp),       INTENT(INOUT)           :: T0, Tsurface
    REAL(r_2),      DIMENSION(1:mp),       INTENT(OUT)             :: Hcum, lEcum,  deltaice_cum_T, deltaice_cum_S
    REAL(r_2),      DIMENSION(1:mp),       INTENT(OUT)             :: Gcum, Qadvcum
    REAL(r_2),      DIMENSION(1:mp),       INTENT(OUT)             :: Jcol_sensible, Jcol_latent_S, Jcol_latent_T
    REAL(r_2),      DIMENSION(1:mp,1:n),   INTENT(OUT)             :: csoil, kth
    REAL(r_2),      DIMENSION(1:mp,1:n),   INTENT(INOUT)           :: phi
    REAL(r_2),      DIMENSION(1:mp),       INTENT(IN)              :: dxL
    REAL(r_2),      DIMENSION(1:mp),       INTENT(INOUT)           :: zdelta, SL, Tl
    TYPE(params),   DIMENSION(1:mp),       INTENT(IN)              :: plit
    TYPE(params),   DIMENSION(1:mp,1:n),   INTENT(INOUT)           :: par
    REAL(r_2),       DIMENSION(1:mp,1:n),  INTENT(IN), OPTIONAL    :: qex
    REAL(r_2),      DIMENSION(1:mp,1:n),   INTENT(INOUT), OPTIONAL :: wex
    REAL(r_2),      DIMENSION(1:mp,1:n),   INTENT(OUT),   OPTIONAL :: heads
    REAL(r_2),      DIMENSION(1:mp,1:n),   INTENT(INOUT), OPTIONAL :: ciso
    REAL(r_2),      DIMENSION(1:mp,1:n),   INTENT(INOUT), OPTIONAL :: cisoice
    REAL(r_2),      DIMENSION(1:mp,1:nsnow_max),   INTENT(INOUT), OPTIONAL :: ciso_snow
    REAL(r_2),      DIMENSION(1:mp,1:nsnow_max),   INTENT(INOUT), OPTIONAL :: cisoice_snow
    REAL(r_2),      DIMENSION(1:mp),       INTENT(INOUT), OPTIONAL :: cisos
    REAL(r_2),      DIMENSION(1:mp),       INTENT(IN),    OPTIONAL :: cprec
    REAL(r_2),      DIMENSION(1:mp),       INTENT(IN),    OPTIONAL :: cprec_snow
    REAL(r_2),      DIMENSION(1:mp),       INTENT(IN),    OPTIONAL :: cali
    REAL(r_2),      DIMENSION(1:mp),       INTENT(IN),    OPTIONAL :: qali
    REAL(r_2),      DIMENSION(1:mp),       INTENT(OUT),   OPTIONAL :: qiso_in, qiso_out
    REAL(r_2),      DIMENSION(1:mp),       INTENT(OUT),   OPTIONAL :: qiso_evap_cum, qiso_trans_cum
    REAL(r_2),      DIMENSION(1:mp,-nsnow_max+1:n),   INTENT(OUT),   OPTIONAL :: qiso_liq_adv, qiso_vap_adv
    REAL(r_2),      DIMENSION(1:mp,-nsnow_max+1:n-1), INTENT(OUT),   OPTIONAL :: qiso_liq_diff, qiso_vap_diff
    REAL(r_2),      DIMENSION(1:mp,-nsnow_max:n),   INTENT(OUT),   OPTIONAL :: qvsig, qlsig, qvTsig, qvh
    REAL(r_2),      DIMENSION(1:mp),       INTENT(INOUT), OPTIONAL :: deltaTa
    REAL(r_2),      DIMENSION(1:mp),       INTENT(IN),    OPTIONAL :: lE_old
    INTEGER(i_d),                          INTENT(IN),    OPTIONAL :: dolitter       ! 0: no; 1: normal; 2: resistance
    INTEGER(i_d),                          INTENT(IN),    OPTIONAL :: doisotopologue ! 0: no isotope; 1: HDO; 2: H218O
    INTEGER(i_d),                          INTENT(IN),    OPTIONAL :: docondition    ! 0: no cond., 1: columns, 2: lines, 3: both
    INTEGER(i_d),                          INTENT(IN),    OPTIONAL :: doadvection    ! 0: off; 1: onn
    ! Error flag if nstep of SLI > nsteps_max: err=0 -> no error; err/=0 -> error
    INTEGER(i_d),  DIMENSION(1:mp),        INTENT(INOUT), OPTIONAL :: err

    ! Solves the RE and, optionally, the ADE from time ts to tfin.
    ! Definitions of arguments:
    ! Required args:
    ! ts   - start time (h).
    ! tfin   - finish time.
    ! qprec   - precipitation (or water input) rate (fluxes are in cm/h).
    ! qevap   - potl evaporation rate from soil surface.
    ! n    - no. of soil layers.
    ! nsol   - no. of solutes.
    ! dx(1:n) - layer thicknesses.
    ! h0   - surface head, equal to depth of surface pond.
    ! S(1:n)  - degree of saturation ("effective satn") of layers.
    ! evap   - cumulative evaporation from soil surface (cm, not initialised).
    ! runoff  - cumulative runoff.
    ! infil   - cumulative net infiltration (time integral of flux across surface).
    ! drn   - cumulative net drainage (time integral of flux across bottom).
    ! nsteps  - cumulative no. of time steps for RE soln.
    ! Optional args:
    ! heads(1:n)   - matric heads h of layers at finish.
    ! qexsub    - subroutine to get layer water extraction rates (cm/h) by
    !     plants. Note that there is no solute extraction and osmotic
    !     effects due to solute are ignored. Arguments:
    !     qex(1:n) - layer extraction rates; qexh(1:n) - partial
    !     derivs of qex wrt h.
    ! wex(1:n)    - cumulative water extraction from layers.
    ! cin(1:nsol)   - solute concns in water input (user's units/cc).
    ! c0(1:nsol)   - solute concns in surface pond.
    ! sm(1:n,1:nsol)  - solute (mass) concns in layers.
    ! soff(1:nsol)   - cumulative solute runoff (user's units).
    ! sinfil(1:nsol)  - cumulative solute infiltration.
    ! sdrn(1:nsol)   - cumulative solute drainage.
    ! nssteps(1:nsol) - cumulative no. of time steps for ADE soln.
    ! isosub    - subroutine to get adsorbed solute (units/g soil) from concn
    !     in soil water according to chosen isotherm code.
    !     Arguments: iso - 2 character code; c - concn in soil water;
    !     p(:) - isotherm parameters; f - adsorbed mass/g soil;
    !     fc - deriv of f wrt c (slope of isotherm curve). Note that
    !     linear adsorption does not require a sub, and other types
    !     are available in sub isosub.

    REAL(r_2),    DIMENSION(1:mp)       :: precip, qevap
    REAL(r_2),    DIMENSION(1:mp)       :: qL, qhL, qybL, qTbL, qhTbL, qhybL, rexcol, wcol, ql0, qv0
    LOGICAL,      DIMENSION(1:mp)       :: again, getq0,getqn,init
    LOGICAL,      DIMENSION(1:mp,1:n)   :: again_ice
    INTEGER(i_d), DIMENSION(1:mp)       :: ih0, iok, itmp, ns, nsat, nsatlast, nsteps0
    REAL(r_2),    DIMENSION(1:mp)       :: accel, dmax, dt, dwinfil, dwoff, fac, Khmin1, Kmin1, phimin1, phip
    REAL(r_2),    DIMENSION(1:mp)       :: qpme, rsig, rsigdt, sig, t
    REAL(r_2),    DIMENSION(1:mp,1:n)   :: Sbot, Tbot
    REAL(r_2),    DIMENSION(1:mp,1:n-1) :: dz
    REAL(r_2),    DIMENSION(1:mp,1:n)   :: hint, phimin, qexd
    !   REAL(r_2),    DIMENSION(1:mp,-nsnow_max+1:n)   :: aa, bb, cc, dd, ee, ff, gg, dy
    !   REAL(r_2),    DIMENSION(1:mp,-nsnow_max+1:n)   :: aah, bbh, cch, ddh, eeh, ffh, ggh, de
    !   REAL(r_2),    DIMENSION(1:mp,-nsnow_max:n)   :: q, qya, qyb, qTa, qTb,qhya, qhyb, qhTa, qhTb
    !   REAL(r_2),    DIMENSION(1:mp,-nsnow_max:n)   :: qadv, qadvya, qadvyb, qadvTa, qadvTb

    TYPE(vars)                          :: vtmp
    !TYPE(vars),   DIMENSION(1:mp,1:n)   :: var
    !    REAL(r_2),    DIMENSION(1:mp,-nsnow_max:n)   :: qsig, qhsig, qadvsig
    !    REAL(r_2),    DIMENSION(1:mp,-nsnow_max:n)   :: qliq, qv, qvT, qlya, qlyb, qvya, qvyb, qlTb, qvTa, qvTb
    TYPE(vars),   DIMENSION(1:mp,1:n)   :: vcall
    REAL(r_2),    DIMENSION(1:mp,1:n)   :: deltaS, dTsoil
    REAL(r_2),    DIMENSION(1:mp,0:n)   :: tmp2d1, tmp2d2
    !    REAL(r_2),    DIMENSION(1:mp,1:n)   :: S0, Sliq0, Sliq, deltaSliq, cv0, deltacv
    !    REAL(r_2),    DIMENSION(1:mp,1:n)   :: Sliqice0, Sliqice, deltaSliqice
    !    REAL(r_2),    DIMENSION(1:mp,1:n)   :: Sice0, Sice, deltaSice
    !
    !    REAL(r_2),    DIMENSION(1:mp,-nsnow_max+1:n)   :: Sliq0_ss, Sliq_ss, deltaSliq_ss
    !    REAL(r_2),    DIMENSION(1:mp,-nsnow_max+1:n)   :: Sliqice0_ss, Sliqice_ss, deltaSliqice_ss
    !    REAL(r_2),    DIMENSION(1:mp,-nsnow_max+1:n)   :: Sice0_ss, Sice_ss, deltaSice_ss
    !    REAL(r_2),    DIMENSION(1:mp,-nsnow_max+1:n)   :: S0_ss, S_ss, deltaS_ss, Tsoil_ss, dTsoil_ss
    !    REAL(r_2),    DIMENSION(1:mp,-nsnow_max+1:n)   :: cv0_ss, cv_ss, Dv_ss, deltacv_ss, dx_ss
    !    REAL(r_2),    DIMENSION(1:mp,-nsnow_max+1:n-1)   :: dz_ss
    !    REAL(r_2),      DIMENSION(1:mp,1:nsnow_max) :: cisoliqice_snow
    !    INTEGER(i_d) :: itop ! integer corresponding to top of soil-snow column
    !    INTEGER(i_d) :: nsnow ! number of dedicated snow layers
    !    REAL(r_2),    DIMENSION(1:mp,1:n)   :: tmp_thetasat, tmp_thetar
    !    REAL(r_2),    DIMENSION(1:mp,-nsnow_max+1:n)   :: thetasat_ss, thetar_ss
    !    REAL(r_2),    DIMENSION(1:mp,-nsnow_max+1:n)   ::  tmp_tortuosity
    !    REAL(r_2),    DIMENSION(1:mp,-nsnow_max+1:n)   ::  ciso_ss, cisoice_ss
    !    REAL(r_2),    DIMENSION(1:mp,1:n)   :: delthetai, dthetaldT, thetal
    INTEGER(i_d), DIMENSION(1:mp,1:n)   :: isave !, nsteps_ice, imelt

    TYPE(vars),         DIMENSION(1:mp) :: vtop, vbot
    TYPE(vars_aquifer), DIMENSION(1:mp) :: v_aquifer
    REAL(r_2),          DIMENSION(1:mp) :: qd, dwcol, dwdrainage, drn,inlit, dwinlit, drexcol, dwdischarge
    REAL(r_2),          DIMENSION(1:mp) :: dJcol_latent_S, dJcol_latent_T, dJcol_sensible
    REAL(r_2),          DIMENSION(1:mp,1:n):: deltaJ_latent_S, deltaJ_latent_T, deltaJ_sensible_S, deltaJ_sensible_T
    REAL(r_2),          DIMENSION(1:mp) :: qevapsig
    REAL(r_2),          DIMENSION(1:mp) :: qrunoff
    REAL(r_2),          DIMENSION(1:mp) :: tmp1d1, tmp1d2, tmp1d3,  tmp1d4
    REAL(r_2),          DIMENSION(1:mp) :: deltah0
    REAL(r_2),          DIMENSION(1:mp) :: SL0, deltaSL, cvL0, SLliq0, deltacvL, SLliq, deltaSLliq
    REAL(r_2),          DIMENSION(1:mp) :: qiso_evap, qiso_trans
    REAL(r_2),          DIMENSION(1:mp) :: lE0, G0, Epot
    REAL(r_2),          DIMENSION(1:mp) :: Tfreezing, dT0
    REAL(r_2),          DIMENSION(1:mp) :: dtdT
    REAL(r_2),          DIMENSION(1:mp,-nsnow_max+1:n) :: LHS, RHS, LHS_h, RHS_h
    INTEGER(i_d),       DIMENSION(1:mp) :: surface_case

    INTEGER(i_d),       DIMENSION(1:mp) :: nns, iflux
    LOGICAL      :: litter
    INTEGER(i_d) :: i, j, k, kk, condition
    INTEGER(i_d) :: littercase, isotopologue, advection ! switches
    REAL(r_2)    :: ztmp, c2, theta
    REAL(r_2)    :: dTqwdTa, dTqwdTb, Tqw, keff
    REAL(r_2),          DIMENSION(1:mp) :: cp, cpeff, hice, deltahice, h0_0, hice_0, h0_tmp, hice_tmp
    !REAL(r_2),          DIMENSION(1:mp,nsnow_max) :: qmelt, hsnow
    ! tranfer of water  from soil to   snow (+ve) or soil to snow (-ve) at init or termination of snowpack
    REAL(r_2),  DIMENSION(1:mp) :: qtransfer
    REAL(r_2),  DIMENSION(1:mp) :: qmelt_ss ! melt water from bottom snow layer to soil
    REAL(r_2),  DIMENSION(1:mp) :: qprec_ss
    REAL(r_2),  DIMENSION(1:mp) :: cprec_ss
    REAL(r_2),          DIMENSION(1:mp,nsnow_max) :: delta_snowcol, delta_snowT, delta_snowliq, dTsnow
    REAL(r_2),          DIMENSION(1:mp) :: melt ! cumulative loss of snow pack as melt water
    REAL(r_2),      DIMENSION(1:mp,1:n)       :: thetai_0, J0
    REAL(r_2) :: tmp1, tmp2
    REAL(r_2),          DIMENSION(1:mp,1:n) :: iqex, thetal_max
    REAL(r_2),          DIMENSION(1:mp)     :: icali
    INTEGER(i_d),       DIMENSION(1:mp) :: nfac1, nfac2, nfac3, nfac4, nfac5, &
         nfac6, nfac7, nfac8, nfac9, nfac10, nfac11, nfac12
    REAL(r_2),          DIMENSION(1:mp)     :: J0snow, wcol0snow
    REAL(r_2), DIMENSION(1:n) :: h_ex
    REAL(r_2) :: wpi

    !open (unit=7, file="Test.out", status="replace", position="rewind")
    ! The derived types params and vars hold soil water parameters and variables.
    ! Parameter names often end in e, which loosely denotes "air entry", i.e.,
    ! values at h=he. While values of water content th and hydraulic conductivity K
    ! at h=he are equal to those at saturation, the derivs wrt S are nonzero. The
    ! MFP phi for h>he is given by phi=phie+Ke*(h-he). The saturation status of a
    ! layer is stored as 0 or 1 in isat since S may be >1 (because of previous
    ! overshoot) when a layer desaturates. Fluxes at the beginning of a time step
    ! and their partial derivs wrt S or phi of upper and lower layers or boundaries
    ! are stored in q, qya and qyb.


    ! set switches
    if (present(dolitter)) then
       littercase = dolitter
    else
       littercase = 0
    endif
    if (present(doadvection)) then
       advection = doadvection
    else
       advection = 0
    endif
    if (littercase > 2) then
       write(*,*) 'dolitter not in [0-2]: ', littercase
       stop 2
    endif

    if (present(doisotopologue)) then
       isotopologue = doisotopologue
    else
       isotopologue = 0
    endif
    if (isotopologue > 2) then
       write(*,*) 'doisotopologue not in [0-2]: ', isotopologue
       stop 2
    endif
    if (isotopologue /= 0 .and. (.not. present(ciso))) then
       write(*,*) 'doisotopologue /= 0 but no ciso present.'
       stop 2
    endif

    if (present(docondition)) then
       condition = docondition
    else
       condition = 0
    endif
    if (condition < 0 .or. condition > 3) then
       write(*,*) 'docondition not in [0-3]: ', condition
       stop 2
    endif

    if (present(qex)) then
       iqex = qex
    else
       iqex = zero
    endif

    if (present(qali) .and. present(cali)) then
       where (qali>zero)
          icali = cali
       elsewhere
          icali = zero
       endwhere
    else
       icali = zero
    endif

    ! global counters
    if (.not. allocated(nless)) allocate(nless(mp))
    nless(:) = 0
    if (.not. allocated(n_noconverge)) allocate(n_noconverge(mp))
    n_noconverge(:) = 0

    ! set solve_type for numerical derivatives
    if (.not. allocated(sol)) call setsol(mp)

    ! initialise cumulative variables
    wcol(:)      = zero
    Jcol_sensible(:) = zero
    Jcol_latent_S(:) = zero
    Jcol_latent_T(:) = zero
    deltaice_cum_T(:) = zero
    deltaice_cum_S(:) = zero
    deltaJ_sensible_S(:,:) = zero
    deltaJ_sensible_T(:,:) = zero
    deltaJ_latent_S(:,:) = zero
    deltaJ_latent_T(:,:) = zero
    drainage(:)  = zero
    discharge(:) = zero
    infil(:)     = zero
    inlit(:)     = zero
    dwinlit(:)   = zero
    evap(:)      = zero
    evap_pot(:)  = zero
    runoff(:)    = zero
    melt(:) = zero
    rexcol(:)    = zero
    Hcum(:)      = zero
    Gcum(:)      = zero
    lEcum(:)     = zero
    Qadvcum(:)  = zero
    wex(:,:)          = zero
    precip(:)         = zero
    drn(:)            = zero
    if (isotopologue /= 0) then
       qiso_evap_cum(:)  = zero
       qiso_trans_cum(:) = zero
    endif
    deltah0(:) = zero
    ! zero var-structure that contains all the hydrological variables
    vtmp = zerovars()
    vtmp%h       = one
    vtmp%lambdav = rlambda
    vtmp%lambdaf = lambdaf
    ! zero vars at the bottom and top of the soil column
    vtop = spread(vtmp,1,mp)
    vbot = spread(vtmp,1,mp)
    vlit = spread(vtmp,1,mp)
    ! Vanessa: try this with limited stacksize, otherwise the double-loop
    !var = spread(vtop,2,n)
    ! do i=1, mp
    !    do k=1, n
    !       var(i,k) = vtmp
    !    end do
    ! end do
    hint(:,:)   = zero
    phimin(:,:) = zero
    dTsoil(:,:) = zero
    keff               = zero
    dTqwdTa            = zero
    dTqwdTb            = zero
    Tqw                = zero
    ztmp               = zero
    c2                 = zero
    theta              = zero
    cp                 = zero
    cpeff              = zero
    hice               = zero
    deltahice          = zero
    h0_0               = zero
    hice_0             = zero
    h0_tmp             = zero
    hice_tmp           = zero
    !qmelt(:,:)         = zero
    qmelt_ss(:)        = zero
    cprec_ss(:) = zero
    qprec_ss(:) = zero
    qtransfer(:)       = zero
    !hsnow(:,:)         = zero
    delta_snowcol(:,:) = zero
    delta_snowT(:,:)   = zero
    delta_snowliq(:,:) = zero
    dTsnow(:,:)        = zero
    melt(:)            = zero
    thetai_0(:,:)      = zero
    J0(:,:)            = zero
    dT0(:)             = zero
    thetal_max         = zero
    nsteps             = 0
    qvsig(:,:)         = zero
    qlsig(:,:)         = zero

    ! initialise snow Ebal diagnostics
    vsnow(:)%Qadv_rain = zero
    vsnow(:)%Qadv_snow = zero
    vsnow(:)%Qadv_vap  = zero
    vsnow(:)%Qcond_net = zero

    do k=1, nsnow_max
       vsnow(:)%deltaJlatent(k)   = zero
       vsnow(:)%deltaJsensible(k) = zero
    enddo
    vsnow(:)%Qadv_transfer  = zero
    vsnow(:)%Qadv_melt      = zero
    vsnow(:)%FluxDivergence = zero

    vsnow(:)%MoistureFluxDivergence = zero
    vsnow(:)%Qprec = zero
    vsnow(:)%Qvap = zero
    vsnow(:)%Qevap = zero
    vsnow(:)%Qtransfer = zero
    vsnow(:)%Qmelt = zero

    litter = .false.
    if (littercase == 1) litter=.true. ! full litter model

    qexd(:,:) = zero
    phip(:)   = zero !max(par(:,1)%phie-par(:,1)%he*par(:,1)%Ke, 1.00001_r_2*par(:,1)%phie) ! phi at h=0

    ! get K, Kh and phi at hmin (hmin is smallest h, stored in hy-props)
    do k=1, mp
       call hyofh(hmin, par(k,1)%lam, par(k,1)%eta, par(k,1)%Ke, par(k,1)%he, Kmin1(k), Khmin1(k), phimin1(k))
    end do

    dz(:,:) = half*(dx(:,1:n-1)+dx(:,2:n)) ! flow paths

    !----- set up for boundary conditions
    getq0(:) = .false.
    getqn(:) = .false.
    if (botbc == "constant head") then ! h at bottom bdry specified
       getqn(:)  = .true.
       tmp1d1(:)  = hbot
       ! for hbot < he
       Sbot(:,:) = spread(Sofh(tmp1d1,par(:,n)),1,n)

       Tbot(:,:) = spread(Tsoil(:,n),1,n)
       ! Debug for mp=1: remove elemental from hyofS and do loop instead of next line
       call hyofS(Sbot, Tbot, par, vcall)
       !do i=1, n
       !   call hyofS(Sbot(1,i), Tbot(1,i), par(1,i), vcall(1,i))
       !end do
       ! End debug hyofS
       ! for hbot >= he
       vtmp = zerovars()
       vtmp%isat    = 1
       vtmp%h       = hbot
       vtmp%rh      = one
       vtmp%lambdav = rlambda
       vtmp%lambdaf = lambdaf
       vbot = spread(vtmp,1,mp)
       vbot(:)%phi = (hbot-par(:,n)%he)*par(:,n)%Ke+var(:,n)%phie
       vbot(:)%K   = par(:,n)%Ke
       where (par(:,n)%he > hbot)
          vbot(:)      = vcall(:,n)
          vbot(:)%isat = 0
       endwhere
    end if
    !----- end set up for boundary conditions

    !----- initialise
    nfac1  = 0
    nfac2  = 0
    nfac3  = 0
    nfac4  = 0
    nfac5  = 0
    nfac6  = 0
    nfac7  = 0
    nfac8  = 0
    nfac9  = 0
    nfac10 = 0
    nfac11 = 0
    nfac12=0
    t(:)       = ts
    nsteps0(:) = nsteps
    nsat(:)    = 0
    qd(:)      = zero
    ! initialise saturated regions
    var(:,:)%isat  = 0
    where (S(:,:) >= one)
       var(:,:)%K    = par(:,:)%Ke
       var(:,:)%isat = 1
       var(:,:)%phi = phi
    endwhere

    vlit(:)%isat = 0
    where (SL(:) >= one) vlit(:)%isat = 1

    ! initialise litter
    if (littercase == 1 .or. littercase == 2) then
       call litter_props(Sl(:), Tl(:), vlit(:), plit(:), h0(:))
    endif
    ! Add resistance through litter for simple litter model
    if (littercase == 2) then
       ztmp        = one/rhocp
       where (vsnow(:)%nsnow == 0)
          vmet(:)%rbw = vmet(:)%rbw + dxL(:)/vlit(:)%Dv
          vmet(:)%rbh = vmet(:)%rbh + dxL(:)/(vlit(:)%kH*ztmp)
          vmet(:)%rrc = vmet(:)%rrc + dxL(:)/(vlit(:)%kH*ztmp)
       endwhere
    endif

    lE0(:) = lE_old(:) ! used for initial guess of litter temperature
    !----- end initialise

    !----- solve until tfin
    init(:) = .true. ! flag to initialise h at soil interfaces
    do kk=1, mp
       !  rewind(3337)
       !  write(3337,*) irec, kk, Tsoil(kk,1)

       !  CALL snow_augment(irec, mp, n, kk, ns, qprec_snow, vmet%Ta, tfin, h0, hice, thetai, dx, vsnow, var, par, S, Tsoil, &
       !                 Jcol_latent_S, Jcol_latent_T, Jcol_sensible, deltaJ_sensible_S, qmelt, qtransfer, j0snow)
       J0snow(kk) = vsnow(kk)%J ! for tracking change in internal energy of snowpack


       wcol0snow(kk) = sum(vsnow(kk)%hsnow(1:nsnow_max)) ! for tracking change in water content of snowpack

       CALL timestep_loop( &
            tfin, irec, mp, qprec, qprec_snow, n, dx(kk,:), h0, S(kk,:), thetai(kk,:), Jsensible(kk,:), &
            Tsoil(kk,:), evap, runoff, &
            infil, drainage, discharge, qh, nsteps, vmet, vlit, vsnow, var(kk,:), T0, Tsurface, Hcum, lEcum, deltaice_cum_T, &
            deltaice_cum_S, Gcum, Qadvcum, Jcol_sensible, Jcol_latent_S, Jcol_latent_T, csoil(kk,:), kth(kk,:), &
            phi(kk,:), dxL, zdelta, &
            SL, Tl, plit, par(kk,:), wex, ciso, cisoice, ciso_snow, cisoice_snow, cisos, cprec, cprec_snow, qali, &
            qiso_in, qiso_out, &
            qiso_evap_cum, qiso_trans_cum, qiso_liq_adv, qiso_vap_adv, qiso_liq_diff, qiso_vap_diff, qvsig, qlsig, &
            qvTsig, qvh, deltaTa, &
            precip, qevap, qL, qhL, qybL, qTbL, qhTbL, qhybL, rexcol, wcol, ql0, qv0, again, getq0,getqn,init, &
            again_ice, ih0, iok, itmp, &
            ns, nsat, nsatlast, nsteps0, accel, dmax, dt, dwinfil, dwoff, fac, phip, qpme, rsig, rsigdt, sig, t, hint(kk,:), &
            phimin(kk,:), qexd(kk,:), &
            vtmp, deltaS(kk,:), dTsoil(kk,:), tmp2d1(kk,:), tmp2d2(kk,:), vtop, vbot, v_aquifer, dwcol, &
            dwdrainage, drn,inlit, &
            dwinlit, drexcol, dwdischarge, dJcol_latent_S, dJcol_latent_T, dJcol_sensible, deltaJ_latent_S(kk,:), &
            deltaJ_latent_T(kk,:), &
            deltaJ_sensible_S(kk,:), deltaJ_sensible_T(kk,:), qevapsig, qrunoff, tmp1d1, tmp1d2, tmp1d3, tmp1d4, &
            deltah0, SL0, deltaSL, &
            cvL0, SLliq0, deltacvL, SLliq, deltaSLliq, qiso_evap, qiso_trans, lE0, G0, Epot, Tfreezing, dtdT, &
            LHS(kk,:), RHS(kk,:), LHS_h(kk,:), RHS_h(kk,:), &
            surface_case, nns, iflux, litter, i, j, k, kk, condition, littercase, isotopologue, advection, &
            c2, theta, dTqwdTa, &
            dTqwdTb, Tqw, keff, cp, cpeff, hice, deltahice, h0_0, hice_0, h0_tmp, hice_tmp, qtransfer, qmelt_ss, qprec_ss, &
            cprec_ss, delta_snowcol(kk,:), delta_snowT(kk,:), delta_snowliq(kk,:), thetai_0(kk,:), J0(kk,:), &
            tmp1, tmp2, iqex(kk,:), icali, nfac1, nfac2, &
            nfac3, nfac4, nfac5, nfac6, nfac7, nfac8, nfac9, nfac10, nfac11, nfac12, J0snow, wcol0snow, h_ex, wpi,  &
            err=err(kk))
       ! only if melt_transfer=.false.: runoff(kk) = runoff(kk) + vsnow(kk)%Qmelt
       !runoff(kk) = runoff(kk) + vsnow(kk)%Qmelt
    end do ! kk=1, mp

    ! get heads if required
    if (present(heads)) then
       isave    = var%isat
       var%isat = isave
       heads    = var%h
       where (S(:,:) >= one) heads(:,:) = par(:,:)%he + (var(:,:)%phi-par(:,:)%phie)/par(:,:)%Ke
    end if

  END SUBROUTINE solve

  !*********************************************************************************************************************

  ! SUBROUTINE solute(ti,tf,thi,thf,win,cin,n,ns,dx,jt,dsmmax,sm,sdrn,nssteps,c, &
  !      isosub)

  !   USE sli_utils, ONLY: dis, isotype, bd, isopar

  !   IMPLICIT NONE

  !   INTEGER(i_d),INTENT(IN)::n,ns,jt(n)
  !   REAL(r_2),INTENT(IN)::ti,tf,thi(n),thf(n),win,cin(ns),dx(n),dsmmax
  !   INTEGER(i_d),INTENT(INOUT)::nssteps(ns)
  !   REAL(r_2),INTENT(INOUT)::sm(n,ns),sdrn(ns),c(n,ns)
  !   OPTIONAL::isosub
  !   INTERFACE
  !      SUBROUTINE isosub(iso,c,p,f,fc)
  !        USE sli_numbers, ONLY: r_2
  !        CHARACTER(LEN=2),INTENT(IN)::iso
  !        REAL(r_2),INTENT(IN)::c
  !        REAL(r_2),DIMENSION(:),INTENT(INOUT)::p
  !        REAL(r_2),INTENT(OUT)::f, fc
  !      END SUBROUTINE isosub
  !   END INTERFACE
  !   ! Solves the ADE from time ti to tf. Diffusion of solute ignored - dispersion
  !   ! coeff = dispersivity * abs(pore water velocity).
  !   ! Definitions of arguments:
  !   ! Required args:
  !   ! ti   - start time (h).
  !   ! tf   - finish time.
  !   ! thi(1:n)  - initial layer water contents.
  !   ! thf(1:n)  - final layer water contents.
  !   ! win   - water in at top of profile.
  !   ! cin(1:ns)  - solute concn in win.
  !   ! n    - no. of soil layers.
  !   ! ns   - no. of solutes.
  !   ! dx(1:n)  - layer thicknesses.
  !   ! jt(1:n)  - layer soil type nos.
  !   ! dsmmax(1:ns) - max change in sm of any layer to aim for each time step;
  !   !      controls time step size.
  !   ! sm(1:n,1:ns) - layer masses of solute per cc.
  !   ! sdrn(1:ns) - cumulative solute drainage.
  !   ! nssteps(1:ns) - cumulative no. of time steps for ADE soln.
  !   ! Optional args:
  !   ! isosub  - subroutine to get adsorbed solute (units/g soil) from concn
  !   !      in soil water according to chosen isotherm code.
  !   !      Arguments: iso - 2 character code; c - concn in soil water;
  !   !      p(:) - isotherm parameters; f - adsorbed mass/g soil;
  !   !      fc - deriv of f wrt c (slope of isotherm curve).
  !   INTEGER(i_d),PARAMETER::itmax=20 ! max iterations for finding c from sm
  !   REAL(r_2),PARAMETER::eps=0.00001 ! for stopping
  !   INTEGER(i_d)::i,it,j,k
  !   REAL(r_2)::dc,dm,dmax,dt,dz(n-1),f,fc,r,rsig,rsigdt,sig,sigdt,t,tfin,th,v1,v2
  !   REAL(r_2),DIMENSION(n-1)::coef1,coef2
  !   REAL(r_2),DIMENSION(n)::csm,tht
  !   REAL(r_2),DIMENSION(0:n)::aa,bb,cc,dd,dy,q,qw,qya,qyb
  !   INTEGER(i_d) :: info

  !   sig=half
  !   rsig=one/sig
  !   tfin=tf
  !   dz=half*(dx(1:n-1)+dx(2:n))
  !   !get average water fluxes
  !   r=one/(tf-ti)
  !   qw(0)=r*win
  !   tht=r*(thf-thi)
  !   do i=1,n
  !      qw(i)=qw(i-1)-dx(i)*tht(i)
  !   end do
  !   !get constant coefficients
  !   do i=1,n-1
  !      v1=half*qw(i)
  !      v2=half*(dis(jt(i))+dis(jt(i+1)))*abs(qw(i))/dz(i)
  !      coef1(i)=v1+v2
  !      coef2(i)=v1-v2
  !   end do
  !   do j=1,ns
  !      t=ti
  !      if (qw(0)>zero) then
  !         q(0)=qw(0)*cin(j)
  !      else
  !         q(0)=zero
  !      end if
  !      qyb(0)=zero
  !      do while (t<tfin)
  !         ! get fluxes
  !         do i=1,n
  !            ! get c and csm=dc/dsm (with theta constant)
  !            k=jt(i)
  !            th=thi(i)+(t-ti)*tht(i)
  !            if (isotype(k,j)=="no" .or. sm(i,j)<zero) then ! handle sm<0 here
  !               csm(i)=one/th
  !               c(i,j)=csm(i)*sm(i,j)
  !            else if (isotype(k,j)=="li") then
  !               csm(i)=one/(th+bd(k)*isopar(k,j)%p(1))
  !               c(i,j)=csm(i)*sm(i,j)
  !            else
  !               do it=1,itmax ! get c from sm using Newton's method and bisection
  !                  if (c(i,j)<zero) c(i,j)=zero ! c and sm are >=0
  !                  call isosub(isotype(k,j), c(i,j), isopar(k,j)%p(:), f, fc)
  !                  csm(i)=one/(th+bd(k)*fc)
  !                  dm=sm(i,j)-(bd(k)*f+th*c(i,j))
  !                  dc=dm*csm(i)
  !                  if (sm(i,j)>=zero .and. c(i,j)+dc<zero) then
  !                     c(i,j)=half*c(i,j)
  !                  else
  !                     c(i,j)=c(i,j)+dc
  !                  end if
  !                  if (abs(dm)<eps*(sm(i,j)+10.0_r_2*dsmmax)) exit
  !                  if (it==itmax) then
  !                     write(*,*) "solute: too many iterations getting c"
  !                     stop
  !                  end if
  !               end do
  !            end if
  !         end do
  !         q(1:n-1)=coef1*c(1:n-1,j)+coef2*c(2:n,j)
  !         qya(1:n-1)=coef1*csm(1:n-1)
  !         qyb(1:n-1)=coef2*csm(2:n)
  !         q(n)=qw(n)*c(n,j)
  !         qya(n)=qw(n)*csm(n)
  !         ! get time step
  !         dmax=maxval(abs(q(1:n)-q(0:n-1))/dx)
  !         if (dmax==zero) then
  !            dt=tfin-t
  !         elseif (dmax<zero) then
  !            write(*,*) "solute: errors in fluxes prevent continuation"
  !            stop
  !         else
  !            dt=dsmmax/dmax
  !         end if
  !         if (t+1.1_r_2*dt>tfin) then
  !            dt=tfin-t
  !            t=tfin
  !         else
  !            t=t+dt
  !         end if
  !         sigdt=sig*dt
  !         rsigdt=one/sigdt
  !         ! adjust q for change in theta
  !         q(1:n-1)=q(1:n-1)-sigdt*(qya(1:n-1)*tht(1:n-1)*c(1:n-1,j)+ &
  !              qyb(1:n-1)*tht(2:n)*c(2:n,j))
  !         q(n)=q(n)-sigdt*qya(n)*tht(n)*c(n,j)
  !         ! get and solve eqns
  !         aa(2:n)=qya(1:n-1)
  !         cc(1:n-1)=-qyb(1:n-1)
  !         bb(1:n)=qyb(0:n-1)-qya(1:n)-dx*rsigdt
  !         dd(1:n)=-(q(0:n-1)-q(1:n))*rsig
  !         call dgtsv(n, 1, aa(2:n), bb(1:n), cc(1:n-1), dd(1:n), n, info)
  !         if (info > 0) then
  !            write(*,*) 'solute: singular matrix (01).'
  !            stop
  !         endif
  !         dy(1:n) = dd(1:n)

  !         ! update unknowns
  !         sdrn(j)=sdrn(j)+(q(n)+sig*qya(n)*dy(n))*dt
  !         sm(:,j)=sm(:,j)+dy(1:n)
  !         nssteps(j)=nssteps(j)+1
  !      end do
  !   end do
  ! END SUBROUTINE solute
  !*********************************************************************************************************************
  
  ! SUBROUTINE snow_augment( mp, kk,  qprec_snow, Ta, tfin,     vsnow     &
  !      )

  !   INTEGER(i_d),                               INTENT(IN)    :: mp    ! # of grid-cells
  !   INTEGER(i_d),                               INTENT(IN)    :: kk    ! grid-cell reference
  !   REAL(r_2),    DIMENSION(mp),             INTENT(INOUT) :: qprec_snow    ! snowfall ms-1
  !   REAL(r_2),    DIMENSION(mp),             INTENT(IN) :: Ta    ! air temp
  !   REAL(r_2),                 INTENT(IN) :: tfin    ! time
  !   TYPE(vars_snow), DIMENSION(1:mp),           INTENT(INOUT) :: vsnow
  !   REAL(r_2),       DIMENSION(1:mp) :: tmp1d1, tmp1d2

  !   if (qprec_snow(kk).gt.0) then
  !      ! total energy of augmented snow pack
  !      tmp1d1(kk) = rhow*(vsnow(kk)%hsnow(1)-vsnow(kk)%hliq(1))*(csice*vsnow(kk)%tsn(1) - lambdaf) + &
  !           rhow*(qprec_snow(kk)*tfin)* &
  !           (csice*min(Ta(kk),zero) - lambdaf)
  !      ! total water content of augmented snow pack
  !      vsnow(kk)%hsnow(1) = vsnow(kk)%hsnow(1) + qprec_snow(kk)*tfin
  !      vsnow(kk)%wcol = vsnow(kk)%wcol+qprec_snow(kk)*tfin

  !      !  write(*,*) 'augment_snow', tfin,kk, qprec_snow(kk), vsnow(kk)%hsnow(1)
  !      ! temperature and liquid moisture content of augmented snow pack
  !      if (tmp1d1(kk).lt. -vsnow(kk)%hsnow(1)*rhow*lambdaf) then
  !         vsnow(kk)%hliq(1)= zero
  !         vsnow(kk)%tsn(1)= (tmp1d1(kk)/(rhow*vsnow(kk)%hsnow(1))+lambdaf)/csice
  !      else
  !         vsnow(kk)%tsn(1)= zero
  !         vsnow(kk)%hliq(1)= vsnow(kk)%hsnow(1)+ tmp1d1(kk)/(rhow*lambdaf)
  !      endif

  !      ! density of augmented snow pack
  !      ! snowfall density (tmp1d1), LaChapelle 1969
  !      if (Ta(kk) > 2.0_r_2) then
  !         tmp1d2(kk) = 189.0_r_2
  !      elseif ((Ta(kk) > -15.0_r_2) .and. (Ta(kk) <= 2.0_r_2)) then
  !         tmp1d2(kk) = 50.0_r_2 + 1.7_r_2*(Ta(kk)+15.0_r_2)**1.5_r_2
  !      else
  !         tmp1d2(kk) = 50.0_r_2
  !      endif

  !      vsnow(kk)%dens(1) = (vsnow(kk)%dens(1)*(vsnow(kk)%hsnow(1)-qprec_snow(kk)*tfin) &
  !           + tmp1d2(kk)*qprec_snow(kk)*tfin)/vsnow(kk)%hsnow(1)

  !      vsnow(kk)%depth(1) = vsnow(kk)%hsnow(1)/(vsnow(kk)%dens(1)/rhow)

  !      vsnow(kk)%Qprec = qprec_snow(kk)*tfin
  !      vsnow(kk)%Qadv_snow=rhow*(qprec_snow(kk))* &
  !           (csice*(min(Ta(kk),zero))-lambdaf)*tfin
  !      qprec_snow(kk)=0
  !      if (vsnow(kk)%nsnow==0) then
  !         vsnow(kk)%nsnow=1
  !      endif
  !   endif

  ! END SUBROUTINE snow_augment

  !*********************************************************************************************************************

  SUBROUTINE snow_adjust(irec, mp, n, kk, ns, h0, hice, thetai, dx, vsnow, var, par, S, Tsoil, &
       Jcol_latent_S, Jcol_latent_T, Jcol_sensible, deltaJ_sensible_S, qmelt, qtransfer, j0snow)

    INTEGER(i_d),                               INTENT(IN)    :: irec  ! # of grid-cells
    INTEGER(i_d),                               INTENT(IN)    :: mp    ! # of grid-cells
    INTEGER(i_d),                               INTENT(IN)    :: n     ! # of soil layers
    INTEGER(i_d),                               INTENT(IN)    :: kk    ! grid-cell reference
    INTEGER(i_d),    DIMENSION(mp),             INTENT(INOUT) :: ns    ! pond (0), np ond (1)
    REAL(r_2),       DIMENSION(1:n),         INTENT(IN)    :: dx    ! soil depths
    REAL(r_2),       DIMENSION(1:n),         INTENT(INOUT) :: Tsoil ! soil temperatures soil
    REAL(r_2),       DIMENSION(1:n),         INTENT(INOUT) :: S     ! soil temperatures soil
    REAL(r_2),       DIMENSION(mp),             INTENT(INOUT) :: h0, hice, j0snow ! pond
    REAL(r_2),       DIMENSION(1:n),       INTENT(INOUT) :: thetai
    REAL(r_2),       DIMENSION(1:mp),           INTENT(INOUT) :: Jcol_latent_S, Jcol_latent_T, Jcol_sensible
    REAL(r_2),       DIMENSION(1:n),       INTENT(INOUT) :: deltaJ_sensible_S
    REAL(r_2),       DIMENSION(nsnow_max), INTENT(INOUT) :: qmelt
    REAL(r_2),       DIMENSION(1:mp), INTENT(OUT) :: qtransfer
    TYPE(vars),      DIMENSION(1:n),       INTENT(INOUT) :: var
    TYPE(params),    DIMENSION(1:n),       INTENT(IN) :: par
    TYPE(vars_snow), DIMENSION(1:mp),           INTENT(INOUT) :: vsnow
    REAL(r_2),       DIMENSION(1:mp) :: tmp1d1, tmp1d2, tmp1d3,  tmp1d4
    REAL(r_2),       DIMENSION(1:mp) :: h0_tmp, hice_tmp
    REAL(r_2) :: theta, tmp1, tmp2 ,Tfreezing(1:mp), Jsoil, theta_tmp
    INTEGER(i_d) :: i,j ! counters
    LOGICAL :: melt_transfer=.true.

    tmp1d1(kk) = h0(kk)+dx(1)*(var(1)%thetai+var(1)%thetal) ! total moisture content of top soil layer + pond
    ! tmp1d1(kk) = hice(kk)+dx(1)*(var(1)%thetai) ! total ice  moisture content of top soil layer + pond
    vsnow(kk)%melt = zero
    qmelt = zero
    qtransfer = zero
    ! no dedicated snow pack if solid part of cumulated snow is less than min thresshold
    ! also, don't initialise dedicated snowpack if this would deplete water in top soil layer to less than 1 mm
    if (((vsnow(kk)%wcol-sum(vsnow(kk)%hliq(:)))<snmin*(vsnow(kk)%dens(1)/rhow)).or. &
         ((vsnow(kk)%nsnow==0).and.(tmp1d1(kk)-vsnow(kk)%wcol)<0.001_r_2)) then
       vsnow(kk)%nsnow = 0
       theta         = S(1)*(par(1)%thre) + (par(1)%the - par(1)%thre)
       if (vsnow(kk)%hsnow(1)>zero)  then ! termination of dedicated snow layer
          ! total energy in old snow layer
          tmp1d1(kk) = (vsnow(kk)%tsn(1))*rhow*vsnow(kk)%hliq(1)*cswat + &
               (vsnow(kk)%hsnow(1)-vsnow(kk)%hliq(1))*rhow*((vsnow(kk)%tsn(1))*csice-lambdaf)
          vsnow(kk)%Qadv_transfer = vsnow(kk)%Qadv_transfer - tmp1d1(kk) ! transfer of energy from soil to snow
          vsnow(kk)%Qtransfer = vsnow(kk)%Qtransfer -vsnow(kk)%hsnow(1) ! transfer of water from soil to snow
          qtransfer(kk) = -vsnow(kk)%hsnow(1) ! transfer of water from soil to snow
          vsnow(kk)%deltaJlatent(1) = vsnow(kk)%deltaJlatent(1) +lambdaf*(vsnow(kk)%hsnow(1)-vsnow(kk)%hliq(1))*rhow
          vsnow(kk)%deltaJsensible(1) = vsnow(kk)%deltaJsensible(1) -(vsnow(kk)%tsn(1))*rhow*vsnow(kk)%hliq(1)*cswat- &
               (vsnow(kk)%hsnow(1)-vsnow(kk)%hliq(1))*rhow*(vsnow(kk)%tsn(1))*csice

          ! total energy in old top soil layer
          tmp1d2(kk) = JSoilLayer(Tsoil(1), &
               dx(1), theta,par(1)%css, par(1)%rho, &
               h0(kk), par(1)%thre, par(1)%the, &
               par(1)%he, one/(par(1)%lambc*freezefac))

          tmp1d3(kk) = Tsoil(1)
          !calculate new thetal, consistent with total energy and new pond height
          h0_tmp(kk) = h0(kk)
          if (h0(kk)>zero) then
             h0(kk) = h0(kk) + vsnow(kk)%hsnow(1)
          elseif ((theta+vsnow(kk)%hsnow(1)/dx(1))<par(1)%the) then
             h0(kk)=zero
             theta = theta+vsnow(kk)%hsnow(1)/dx(1)
             S(1) = (theta - (par(1)%the - par(1)%thre) )/(par(1)%thre)
          elseif ((theta+vsnow(kk)%hsnow(1)/dx(1))>par(1)%the) then
             h0(kk) = ((theta+vsnow(kk)%hsnow(1)/dx(1))-par(1)%the)*dx(1)
             theta = par(1)%the
             S(1) = one
             var(1)%isat = 1
             ns(kk) = 0
          endif

          Tfreezing(kk) = Tfrz(S(1), par(1)%he, one/(par(1)%lambc*freezefac))

          ! check if total energy in old snow layer and top soil (+pond) is enough for complete melting
          tmp1d3(kk) = var(1)%thetai*dx(1)*rhow*lambdaf + &
               var(1)%thetai*h0_tmp(kk)/par(1)%thre*rhow*lambdaf
          tmp1d4(kk) = (Tfreezing(kk))*(rhow*cswat*(theta*dx(1)+h0(kk))+par(1)%rho*par(1)%css*dx(1))

          if ((tmp1d1(kk)+tmp1d2(kk)-tmp1d4(kk))>zero) then !  complete melting

             if (Tsoil(1).gt.var(1)%Tfrz) then
                Jcol_latent_T(kk) =     Jcol_latent_T(kk) + tmp1d3(kk)
                deltaJ_sensible_S(1) = zero
                Tsoil(1) = Tsoil(1) + tmp1d1(kk)/(rhow*cswat*(theta*dx(1)+h0(kk))+par(1)%rho*par(1)%css*dx(1))
             else

                Jcol_latent_T(kk) =     Jcol_latent_T(kk) + tmp1d3(kk)
                deltaJ_sensible_S(1) = zero
                Tsoil(1) = var(1)%Tfrz + (tmp1d1(kk)+tmp1d2(kk) -tmp1d4(kk))/ &
                     (rhow*cswat*(theta*dx(1)+h0(kk))+par(1)%rho*par(1)%css*dx(1))
             endif
             var(1)%iice = 0
             var(1)%thetai = zero
             var(1)%thetal = theta
             hice(kk) = zero
             vsnow(kk)%hsnow(1) = zero
             vsnow(kk)%hliq(1)= zero
          else  ! soil remains frozen

             Jsoil = tmp1d1(kk)+tmp1d2(kk)! total energy in  soil layer
             !check there is a zero
             tmp1 = GTfrozen(Tsoil(1)-50._r_2, Jsoil, dx(1), theta,par(1)%css, par(1)%rho, &
                  h0(kk), par(1)%thre, par(1)%the, par(1)%he, one/(par(1)%lambc*freezefac))

             tmp2 = GTFrozen(Tfreezing(kk), Jsoil, dx(1), theta,par(1)%css, par(1)%rho, &
                  h0(kk), par(1)%thre, par(1)%the, par(1)%he, one/(par(1)%lambc*freezefac))

             ! there is a zero in between
             if ((tmp1*tmp2) < zero) then
                tmp1d3(kk) = rtbis_Tfrozen(tmp1d2(kk), dx(1), theta,par(1)%css, par(1)%rho, &
                     h0(kk), par(1)%thre, par(1)%the, par(1)%he, one/(par(1)%lambc*freezefac), &
                     Tsoil(1)-50._r_2, Tfreezing(kk))

                tmp1d4(kk) = thetalmax(tmp1d3(kk), S(1), par(1)%he, one/(par(1)%lambc*freezefac), &
                     par(1)%thre, par(1)%the) ! liquid content at new Tsoil
             else
                write(wlogn,*) "Found no solution for Tfrozen 2. ", kk, i
                write(wlogn,*) "Assume soil is totally frozen"
                tmp1d3(kk) = (tmp1d2(kk) + rhow*lambdaf*(theta*dx(1) +  h0(kk))) / &
                     (dx(1)*par(1)%css*par(1)%rho + rhow*csice*(theta*dx(1) + h0(kk)))
                tmp1d4(kk) = 0.0_r_2                
                write(wlogn,*) "frozen soil temperature: ", tmp1d3(kk)
             endif

             hice_tmp(kk) = hice(kk)
             hice(kk) = h0(kk)*var(1)%thetai/par(1)%thre
             vsnow(kk)%hsnow(1) = zero
             vsnow(kk)%hliq(1) = zero
             var(1)%thetal = max(tmp1d4(kk),zero)
             var(1)%thetai = theta - tmp1d4(kk)

             ! correct total energy stored in pond + soil
             Jcol_latent_S(kk) = Jcol_latent_S(kk)  - rhow*lambdaf*((hice(kk)-hice_tmp(kk)) + &
                  dx(1)*(var(1)%thetai-thetai(1)))


             Jcol_sensible(kk) = Jcol_sensible(kk) + &
                  JSoilLayer(tmp1d3(kk), &
                  dx(1), theta,par(1)%css, par(1)%rho, &
                  h0(kk), par(1)%thre, par(1)%the, &
                  par(1)%he, one/(par(1)%lambc*freezefac)) - &
                  tmp1d2(kk) - &
                  (-rhow*lambdaf*((hice(kk)-hice_tmp(kk)) + &
                  dx(1)*(var(1)%thetai-thetai(1))))
             Tsoil(1) = tmp1d3(kk)

             if (var(1)%thetai>zero) then
                var(1)%iice = 1
             endif
          endif

          vsnow(kk)%Jsensible(1) = (vsnow(kk)%hsnow(1)-vsnow(kk)%hliq(1))*rhow*csice*(vsnow(kk)%Tsn(1))
          vsnow(kk)%Jlatent(1) = (vsnow(kk)%hsnow(1)-vsnow(kk)%hliq(1))*rhow*(-lambdaf)
          vsnow(kk)%J = sum(vsnow(kk)%Jsensible(1:vsnow(kk)%nsnow)+vsnow(kk)%Jlatent(1:vsnow(kk)%nsnow))
          vsnow(kk)%deltaJ = vsnow(kk)%J - J0snow(kk)
          vsnow(kk)%FluxDivergence = vsnow(kk)%Qadv_rain + vsnow(kk)%Qadv_snow + vsnow(kk)%Qadv_vap + &
               vsnow(kk)%Qadv_melt + vsnow(kk)%Qcond_net  +  vsnow(kk)%Qadv_transfer

          vsnow(kk)%hsnow(1) = zero
          vsnow(kk)%hliq(1) = zero
          vsnow(kk)%depth(1) = zero


       endif ! termination of dedicated snow layer
    else

       if  (vsnow(kk)%nsnow_last==0) then ! snow layer initialisation (transfer pond/soil water to dedicated snow layer)
          vsnow(kk)%nsnow = 1
          vsnow(kk)%hsnow(1) = vsnow(kk)%wcol
          vsnow(kk)%depth(1) = vsnow(kk)%hsnow(1)*rhow/vsnow(kk)%dens(1)

          if (h0(kk)>vsnow(kk)%wcol) then ! extract new snow layer from pond
             ! total energy in new snow layer
             tmp1d1(kk) = (Tsoil(1))*rhow*vsnow(kk)%hsnow(1)*(cswat*(h0(kk)-hice(kk))/h0(kk) + &
                  csice*hice(kk)/h0(kk)) - rhow*lambdaf*vsnow(kk)%hsnow(1)*hice(kk)/h0(kk)
             ! correct total energy stored in pond + soil
             Jcol_latent_S(kk) = Jcol_latent_S(kk) + rhow*lambdaf*vsnow(kk)%hsnow(1)*hice(kk)/h0(kk)
             Jcol_sensible(kk) = Jcol_sensible(kk) - (Tsoil(1))*rhow*vsnow(kk)%hsnow(1) &
                  *(cswat*(h0(kk)-hice(kk))/h0(kk) +  csice*hice(kk)/h0(kk))

             ! total energy in snowpack totally frozen at 0degC
             tmp1d2(kk) = rhow*vsnow(kk)%hsnow(1)*(csice*zero - lambdaf)
             if (tmp1d1(kk)<=tmp1d2(kk)) then
                ! alll snow water frozen
                vsnow(kk)%hliq(1)=zero
                vsnow(kk)%tsn(1) = (tmp1d1(kk)+rhow*lambdaf*vsnow(kk)%hsnow(1))/(csice*rhow*vsnow(kk)%hsnow(1))
             else
                ! liquid snow water
                vsnow(kk)%hliq(1)=(tmp1d1(kk)-vsnow(kk)%hsnow(1)*rhow*(zero*csice-lambdaf))/ &
                     (rhow*(zero*cswat-zero*csice+lambdaf))
                vsnow(kk)%tsn(1) = zero
             endif

             h0(kk) = h0(kk) - vsnow(kk)%wcol
             hice(kk) = h0(kk)*var(1)%thetai/par(1)%thre

          else  ! extract new snow layer from soil ice + pond
             ! total energy in new snow layer (component extracted from soil ice)
             tmp1d1(kk) = (Tsoil(1))*rhow*(vsnow(kk)%hsnow(1)-h0(kk))*csice - &
                  rhow*lambdaf*(vsnow(kk)%hsnow(1)-h0(kk))
             h0_tmp(kk) = h0(kk)
             hice_tmp(kk) = hice(kk)
             if (h0(kk)>zero) then
                tmp1d1(kk) = tmp1d1(kk) + (Tsoil(1))*rhow*(cswat*(h0(kk)-hice(kk)) + &
                     csice*hice(kk)) - rhow*lambdaf*hice(kk)
             endif
             ! total energy in snowpack totally frozen at 0degC
             tmp1d2(kk) = rhow*vsnow(kk)%hsnow(1)*(csice*zero - lambdaf)
             if (tmp1d1(kk)<=tmp1d2(kk)) then
                ! all snow water frozen
                vsnow(kk)%hliq(1)=zero
                vsnow(kk)%tsn(1) = (tmp1d1(kk)+rhow*lambdaf*vsnow(kk)%hsnow(1))/(csice*rhow*vsnow(kk)%hsnow(1))
             else
                ! liquid snow water
                vsnow(kk)%hliq(1)=(tmp1d1(kk)-vsnow(kk)%hsnow(1)*rhow*(zero*csice-lambdaf))/ &
                     (rhow*(zero*cswat-zero*csice+lambdaf))
                vsnow(kk)%tsn(1) = zero
             endif

             ! correct soil moisture
             S(1) = S(1) - (vsnow(kk)%hsnow(1)-h0(kk))/dx(1)/par(1)%thre
             if (S(1).lt.one) then
                var(1)%isat= 0
             endif
             Tfreezing(kk) = Tfrz(S(1), par(1)%he, one/(par(1)%lambc*freezefac))
             if (S(1)<zero) then
                write(*,*) "error: over-extraction of soil water during snow pack init"
                write(*,*) "S(1), snow-col, deltaS"
                write(*,*) S(1) , vsnow(kk)%hsnow(1), - (vsnow(kk)%hsnow(1)-h0(kk))/dx(1)/par(1)%thre
             endif
             h0(kk) = zero
             hice(kk) = zero
             ! correct total energy stored in pond + soil
             ! correct soil temperature
             theta         = S(1)*(par(1)%thre) + (par(1)%the - par(1)%thre)
             ! total energy added to top soil layer
             tmp1d1(kk) = - tmp1d1(kk)
             ! total energy in old top soil layer
             tmp1d2(kk) = var(1)%csoil*dx(1)*(Tsoil(1)) -lambdaf*dx(1)*var(1)%thetai + &
                  (h0_tmp(kk)-hice_tmp(kk))*cswat*rhow*(Tsoil(1)) + &
                  hice_tmp(kk)*rhow*(csice*(Tsoil(1))-lambdaf)
             ! calculate energy in new top soil layer
             if (var(1)%iice==0) then
                var(1)%csoil = theta*rhow*cswat+par(1)%rho*par(1)%css
                tmp1d3(kk) = (tmp1d1(kk)+tmp1d2(kk))/((theta*rhow*cswat+par(1)%rho*par(1)%css)*dx(1)+ &
                     h0(kk)*rhow*cswat)
                Jcol_sensible(kk) = Jcol_sensible(kk) - &
                     var(1)%csoil*dx(1)*(Tsoil(1)) - &
                     (h0_tmp(kk))*cswat*rhow*(Tsoil(1)) + &
                     (theta*rhow*cswat+par(1)%rho*par(1)%css)*dx(1)*(tmp1d3(kk)) + &
                     (h0_tmp(kk))*cswat*rhow*(tmp1d3(kk))
                var(1)%csoil = theta*rhow*cswat+par(1)%rho*par(1)%css
                Tsoil(1) = tmp1d3(kk)
             else
                ! check if total energy in melt water and top soil (+pond) is enough for complete melting
                tmp1d3(kk) = var(1)%thetai*dx(1)*rhow*lambdaf + &
                     var(1)%thetai*h0_tmp(kk)/par(1)%thre*rhow*lambdaf

                tmp1d4(kk) = (Tfreezing(kk))*(rhow*cswat*(theta*dx(1)+h0(kk))+par(1)%rho*par(1)%css*dx(1))
                if ((tmp1d1(kk)+tmp1d2(kk)-tmp1d4(kk))>zero) then
                   !  complete melting
                   Jcol_latent_T(kk)    = Jcol_latent_T(kk) + tmp1d3(kk)
                   deltaJ_sensible_S(1) = zero
                   Tsoil(1) = var(1)%Tfrz + (tmp1d1(kk)+tmp1d2(kk)-tmp1d4(kk)) / &
                        (rhow*cswat*(theta*dx(1)+h0(kk))+par(1)%rho*par(1)%css*dx(1))
                   var(1)%iice   = 0
                   var(1)%thetai = zero
                   thetai(1)     = var(1)%thetai
                   var(1)%thetal = theta
                else

                   Jsoil = tmp1d1(kk)+tmp1d2(kk)! total energy in  soil layer
                   !check there is a zero
                   tmp1 = GTfrozen(Tsoil(1)-50._r_2, Jsoil, dx(1), theta,par(1)%css, par(1)%rho, &
                        h0(kk), par(1)%thre, par(1)%the, par(1)%he, one/(par(1)%lambc*freezefac))

                   tmp2 = GTFrozen(Tfreezing(kk), Jsoil, dx(1), theta,par(1)%css, par(1)%rho, &
                        h0(kk), par(1)%thre, par(1)%the, par(1)%he, one/(par(1)%lambc*freezefac))

                   ! there is a zero in between
                   if ((tmp1*tmp2) < zero) then
                      tmp1d3(kk) = rtbis_Tfrozen(tmp1d2(kk), dx(1), theta,par(1)%css, par(1)%rho, &
                           h0(kk), par(1)%thre, par(1)%the, par(1)%he, one/(par(1)%lambc*freezefac), &
                           Tsoil(1)-50._r_2, Tfreezing(kk))

                      tmp1d4(kk) = thetalmax(tmp1d3(kk), S(1), par(1)%he, one/(par(1)%lambc*freezefac), &
                           par(1)%thre, par(1)%the) ! liquid content at new Tsoil
                   else
                      write(wlogn,*) "Found no solution for Tfrozen 3. ", kk, i
                      write(wlogn,*) "Assume soil is totally frozen"
                      tmp1d3(kk) = (Jsoil + rhow*lambdaf*(theta*dx(1) +  h0(kk))) / &
                           (dx(1)*par(1)%css*par(1)%rho + rhow*csice*(theta*dx(1) + h0(kk)))
                      tmp1d4(kk) = 0.0_r_2                
                      write(wlogn,*) "frozen soil temperature: ", tmp1d3(kk)
                   endif

                   var(1)%thetal = max(tmp1d4(kk), zero)
                   var(1)%thetai = theta - tmp1d4(kk)
                   ! correct total energy stored in pond + soil
                   Jcol_latent_S(kk) = Jcol_latent_S(kk)- rhow*lambdaf*((hice(kk)-hice_tmp(kk)) + &
                        dx(1)*(var(1)%thetai-thetai(1)))

                   Jcol_sensible(kk) = Jcol_sensible(kk) - &
                        var(1)%csoil*dx(1)*(Tsoil(1)) - &
                        (h0_tmp(kk)-hice_tmp(kk))*cswat*rhow*(Tsoil(1)) - &
                        hice_tmp(kk)*rhow*(csice*(Tsoil(1))) + &
                        dx(1)*(tmp1d3(kk))*par(1)%rho*par(1)%css + &
                        (h0(kk)-hice(kk)+var(1)%thetal*dx(1))*cswat*rhow*(tmp1d3(kk)) + &
                        (hice_tmp(kk)+var(1)%thetai)*rhow*(csice*(tmp1d3(kk)))

                   thetai(1) = var(1)%thetai
                   Tsoil(1) = tmp1d3(kk)

                endif ! incomplete melting

             endif ! iice=1

          endif ! extract new snow layer from soil ice

          vsnow(kk)%deltaJsensible(1) = vsnow(kk)%deltaJsensible(1) + &
               (vsnow(kk)%hsnow(1)-vsnow(kk)%hliq(1))*rhow*csice*(vsnow(kk)%tsn(1)) + &
               vsnow(kk)%hliq(1)*rhow*cswat*(vsnow(kk)%tsn(1))

          vsnow(kk)%deltaJlatent(1) = vsnow(kk)%deltaJlatent(1) + &
               (vsnow(kk)%hsnow(1)-vsnow(kk)%hliq(1))*rhow*(-lambdaf)

          vsnow(kk)%Qadv_transfer = vsnow(kk)%Qadv_transfer + & ! soil to snow
               (vsnow(kk)%hsnow(1)-vsnow(kk)%hliq(1))*rhow*csice*(vsnow(kk)%tsn(1)) + &
               vsnow(kk)%hliq(1)*rhow*cswat*(vsnow(kk)%tsn(1)) + &
               (vsnow(kk)%hsnow(1)-vsnow(kk)%hliq(1))*rhow*(-lambdaf)

          vsnow(kk)%Qtransfer = vsnow(kk)%Qtransfer + vsnow(kk)%hsnow(1)
          qtransfer(kk) = vsnow(kk)%hsnow(1) ! transfer of water from snow to soil
          vsnow(kk)%Jsensible(1) = (vsnow(kk)%hsnow(1)-vsnow(kk)%hliq(1))*rhow*csice*(vsnow(kk)%Tsn(1))
          vsnow(kk)%Jlatent(1) = (vsnow(kk)%hsnow(1)-vsnow(kk)%hliq(1))*rhow*(-lambdaf)
          vsnow(kk)%J = sum(vsnow(kk)%Jsensible(1:vsnow(kk)%nsnow)+vsnow(kk)%Jlatent(1:vsnow(kk)%nsnow))
          vsnow(kk)%deltaJ = vsnow(kk)%J - J0snow(kk)
          vsnow(kk)%FluxDivergence = vsnow(kk)%Qadv_rain + vsnow(kk)%Qadv_snow + vsnow(kk)%Qadv_vap + &
               vsnow(kk)%Qadv_melt + vsnow(kk)%Qcond_net  +  vsnow(kk)%Qadv_transfer

       endif        ! snow layer initialisation

       do i=1, vsnow(kk)%nsnow
          vsnow(kk)%Jsensible(i) = (vsnow(kk)%hsnow(i)-vsnow(kk)%hliq(i))*rhow*csice*(vsnow(kk)%Tsn(i))
          vsnow(kk)%Jlatent(i) = (vsnow(kk)%hsnow(i)-vsnow(kk)%hliq(i))*rhow*(-lambdaf)
       enddo
       vsnow(kk)%J = sum(vsnow(kk)%Jsensible(1:vsnow(kk)%nsnow)+vsnow(kk)%Jlatent(1:vsnow(kk)%nsnow))
       vsnow(kk)%deltaJ = vsnow(kk)%J - J0snow(kk)
       vsnow(kk)%FluxDivergence = vsnow(kk)%Qadv_rain + vsnow(kk)%Qadv_snow + vsnow(kk)%Qadv_vap + &
            vsnow(kk)%Qadv_melt + vsnow(kk)%Qcond_net  +  vsnow(kk)%Qadv_transfer

       ! get snow melt from top layer
       if ((vsnow(kk)%hliq(1) - vsnow(kk)%hsnow(1)*vsnow(kk)%fsnowliq_max(1))>zero) then ! remove melt water
          qmelt(1) = max((vsnow(kk)%hliq(1) - vsnow(kk)%hsnow(1)*vsnow(kk)%fsnowliq_max(1)),zero)
          qmelt(1) = min(qmelt(1), max(0.9_r_2*(vsnow(kk)%hsnow(1)-snmin*(vsnow(kk)%dens(1)/rhow)),zero))
          vsnow(kk)%melt(1) = vsnow(kk)%melt(1) + qmelt(1)
          vsnow(kk)%hliq(1) = vsnow(kk)%hliq(1) - qmelt(1)
          vsnow(kk)%hsnow(1) = vsnow(kk)%hsnow(1) - qmelt(1)
          ! adjust depth and density for snow melt removal
          !tmp1d3(kk) = vsnow(kk)%dens(1)*(one-vsnow(kk)%hliq(1)/vsnow(kk)%hsnow(1))
          !vsnow(kk)%depth(1) = vsnow(kk)%depth(1) - qmelt(1)*rhow/tmp1d3(kk)
          ! vsnow(kk)%dens(1) = vsnow(kk)%hsnow(1)*rhow/vsnow(kk)%depth(1)
          !vsnow(kk)%dens(1) = max(vsnow(kk)%dens(1) - rhow*qmelt(1)/vsnow(kk)%depth(1), 50.0_r_2)
          do i=1, vsnow(kk)%nsnow
             vsnow(kk)%Jsensible(i) = (vsnow(kk)%hsnow(i)-vsnow(kk)%hliq(i))*rhow*csice*(vsnow(kk)%Tsn(i))
             vsnow(kk)%Jlatent(i) = (vsnow(kk)%hsnow(i)-vsnow(kk)%hliq(i))*rhow*(-lambdaf)
          enddo
       endif ! end remove snow melt from top layer

       ! add snow melt from above to layer below
       if (vsnow(kk)%nsnow>1) then

          ! simply add melt water from above if melt water already exists and total new amount doesn't exceed capacity
          if ((vsnow(kk)%hliq(vsnow(kk)%nsnow).gt.zero).and. &
               (vsnow(kk)%hliq(vsnow(kk)%nsnow) + qmelt(1)).le.(vsnow(kk)%hsnow(vsnow(kk)%nsnow) + &
               qmelt(1))*vsnow(kk)%fsnowliq_max(vsnow(kk)%nsnow)) then

             vsnow(kk)%hliq(vsnow(kk)%nsnow) = vsnow(kk)%hliq(vsnow(kk)%nsnow) + qmelt(1)
             vsnow(kk)%hsnow(vsnow(kk)%nsnow) = vsnow(kk)%hsnow(vsnow(kk)%nsnow) + qmelt(1)
             qmelt(vsnow(kk)%nsnow) = zero
             ! or convert excess to melt water
          elseif ((vsnow(kk)%hliq(vsnow(kk)%nsnow).gt.zero).and. &
               (vsnow(kk)%hliq(vsnow(kk)%nsnow) + qmelt(1)).gt.(vsnow(kk)%hsnow(vsnow(kk)%nsnow) + &
               qmelt(1))*vsnow(kk)%fsnowliq_max(vsnow(kk)%nsnow)) then

             tmp1d1(kk) = vsnow(kk)%hliq(vsnow(kk)%nsnow)
             tmp1d2(kk) = vsnow(kk)%hsnow(vsnow(kk)%nsnow)


             qmelt(vsnow(kk)%nsnow) =  tmp1d1(kk) + vsnow(kk)%melt(1) - &
                  vsnow(kk)%fsnowliq_max(vsnow(kk)%nsnow)*(tmp1d2(kk) - tmp1d1(kk))/ &
                  (1._r_2-vsnow(kk)%fsnowliq_max(vsnow(kk)%nsnow))

             vsnow(kk)%hliq(vsnow(kk)%nsnow) = vsnow(kk)%fsnowliq_max(vsnow(kk)%nsnow)*(tmp1d2(kk)  - tmp1d1(kk) )/ &
                  (1._r_2-vsnow(kk)%fsnowliq_max(vsnow(kk)%nsnow))


             vsnow(kk)%hsnow(vsnow(kk)%nsnow) = (tmp1d2(kk) - tmp1d1(kk)) + vsnow(kk)%hliq(vsnow(kk)%nsnow)

             ! or add melt water to completely frozen snowpack below
          else
             tmp1d1(kk) = (vsnow(kk)%hsnow(vsnow(kk)%nsnow)*(vsnow(kk)%tsn(vsnow(kk)%nsnow)*csice-lambdaf) &
                  /(vsnow(kk)%hsnow(vsnow(kk)%nsnow) + qmelt(1)) + lambdaf)/csice

             if (tmp1d1(kk).lt.zero) then ! snow pack remains completely frozen
                vsnow(kk)%tsn(vsnow(kk)%nsnow) = tmp1d1(kk)
                vsnow(kk)%hliq(vsnow(kk)%nsnow) = zero
                vsnow(kk)%hsnow(vsnow(kk)%nsnow) = vsnow(kk)%hsnow(vsnow(kk)%nsnow) + qmelt(1)
                qmelt(vsnow(kk)%nsnow) = zero

             else ! snowpack partially melted
                vsnow(kk)%hliq(vsnow(kk)%nsnow) = vsnow(kk)%hsnow(vsnow(kk)%nsnow)/ &
                     lambdaf*(csice*vsnow(kk)%tsn(vsnow(kk)%nsnow)-lambdaf) + &
                     (vsnow(kk)%hsnow(vsnow(kk)%nsnow) + qmelt(1))
                vsnow(kk)%tsn(vsnow(kk)%nsnow) = zero
                vsnow(kk)%hsnow(vsnow(kk)%nsnow) = vsnow(kk)%hsnow(vsnow(kk)%nsnow) + qmelt(1)
                qmelt(vsnow(kk)%nsnow) = zero


                if (vsnow(kk)%hliq(vsnow(kk)%nsnow).gt. &
                     vsnow(kk)%hsnow(vsnow(kk)%nsnow)*vsnow(kk)%fsnowliq_max(vsnow(kk)%nsnow)) then

                   tmp1d1(kk) = vsnow(kk)%hliq(vsnow(kk)%nsnow)
                   tmp1d2(kk) = vsnow(kk)%hsnow(vsnow(kk)%nsnow)

                   qmelt(vsnow(kk)%nsnow) =  vsnow(kk)%hliq(vsnow(kk)%nsnow) - &
                        vsnow(kk)%fsnowliq_max(vsnow(kk)%nsnow)*(vsnow(kk)%hsnow(vsnow(kk)%nsnow) - &
                        vsnow(kk)%hliq(vsnow(kk)%nsnow))/(1._r_2-vsnow(kk)%fsnowliq_max(vsnow(kk)%nsnow))

                   vsnow(kk)%hliq(vsnow(kk)%nsnow) = vsnow(kk)%fsnowliq_max(vsnow(kk)%nsnow)*(vsnow(kk)%hsnow(vsnow(kk)%nsnow) &
                        - vsnow(kk)%hliq(vsnow(kk)%nsnow))/ &
                        (1._r_2-vsnow(kk)%fsnowliq_max(vsnow(kk)%nsnow))

                   vsnow(kk)%hsnow(vsnow(kk)%nsnow) = vsnow(kk)%hsnow(vsnow(kk)%nsnow) -  &
                        tmp1d1(kk) + vsnow(kk)%hliq(vsnow(kk)%nsnow)
                   tmp1d3(kk) = vsnow(kk)%dens(vsnow(kk)%nsnow)*(one-vsnow(kk)%hliq(vsnow(kk)%nsnow)) &
                        /vsnow(kk)%hsnow(vsnow(kk)%nsnow)


                endif ! (vsnow(kk)%hliq(2).gt.vsnow(kk)%hsnow(2)*vsnow(kk)%fsnowliq_max(2))

             endif ! snowpack partially melted
             vsnow(kk)%melt(vsnow(kk)%nsnow) = vsnow(kk)%melt(vsnow(kk)%nsnow) + qmelt(vsnow(kk)%nsnow)
          endif ! (vsnow(kk)%hliq(2).gt.zero).and. &
          !(vsnow(kk)%hliq(2) + qmelt(1)).le.(vsnow(kk)%hsnow(2) + qmelt(1))*vsnow(kk)%fsnowliq_max(2))

       endif ! if (vsnow(kk)%nsnow>1)

       theta = S(1)*(par(1)%thre) + (par(1)%the - par(1)%thre)
       ! move snow melt from bottom snow layer to top soil + pond
       if (qmelt(vsnow(kk)%nsnow)>zero) then
          vsnow(kk)%Qmelt = vsnow(kk)%Qmelt + qmelt(vsnow(kk)%nsnow)
          if (melt_transfer) then
             h0_tmp(kk) = h0(kk)
             hice_tmp(kk) = hice(kk)

             ! total energy in melt-water is zero

             ! total energy in old top soil layer
             tmp1d2(kk) = JSoilLayer(Tsoil(1), &
                  dx(1), theta,par(1)%css, par(1)%rho, &
                  h0_tmp(kk), par(1)%thre, par(1)%the, &
                  par(1)%he, one/(par(1)%lambc*freezefac))

             !calculate new thetal, consistent with total energy and new pond height
             theta = S(1)*(par(1)%thre) + (par(1)%the - par(1)%thre)
             theta_tmp = theta
                      
             if (h0(kk)>zero) then
                h0(kk) = h0(kk) +  qmelt(vsnow(kk)%nsnow)
             elseif ((theta+qmelt(vsnow(kk)%nsnow)/dx(1))<par(1)%the) then
                h0(kk)=zero
                theta = theta+qmelt(vsnow(kk)%nsnow)/dx(1)
                S(1) = (theta - (par(1)%the - par(1)%thre) )/(par(1)%thre)
                if (S(1).lt.one) then
                   var(1)%isat= 0
                endif
             elseif ((theta+qmelt(vsnow(kk)%nsnow)/dx(1))>par(1)%the) then
                h0(kk) = ((theta+qmelt(vsnow(kk)%nsnow)/dx(1))-par(1)%the)*dx(1)

                S(1) = one
                var(1)%isat = 1
                ns(kk) = 0
                theta = par(1)%the
             endif
             Tfreezing(kk) = Tfrz(S(1), par(1)%he, one/(par(1)%lambc*freezefac))
             ! calculate energy in new top soil layer
             if (var(1)%iice==0) then
                var(1)%csoil = theta*rhow*cswat+par(1)%rho*par(1)%css
                tmp1d3(kk) = (tmp1d2(kk))/((theta*rhow*cswat+par(1)%rho*par(1)%css)*dx(1)+ &
                     h0(kk)*rhow*cswat)
                Jcol_sensible(kk) = Jcol_sensible(kk) - &
                     var(1)%csoil*dx(1)*(Tsoil(1)) - &
                     (h0_tmp(kk))*cswat*rhow*(Tsoil(1)) + &
                     (theta*rhow*cswat+par(1)%rho*par(1)%css)*dx(1)*(tmp1d3(kk)) + &
                     (h0_tmp(kk))*cswat*rhow*(tmp1d3(kk))

                var(1)%csoil = theta*rhow*cswat+par(1)%rho*par(1)%css
                Tsoil(1) = tmp1d3(kk)

             else
                ! check if total energy in melt water and top soil (+pond) is enough for complete melting
                tmp1d3(kk) = var(1)%thetai*dx(1)*rhow*lambdaf + &
                     var(1)%thetai*h0_tmp(kk)/par(1)%thre*rhow*lambdaf
                tmp1d4(kk) = (Tfreezing(kk))*(rhow*cswat*(theta*dx(1)+h0(kk))+par(1)%rho*par(1)%css*dx(1))
                if ((tmp1d2(kk)-tmp1d4(kk))>zero) then
                   !  complete melting
                   Jcol_latent_T(kk) = Jcol_latent_T(kk) + tmp1d3(kk)
                   deltaJ_sensible_S(1) = zero
                   ! Tsoil(1) = var(1)%Tfrz + (tmp1d1(kk)+tmp1d2(kk) -tmp1d4(kk))/ &
                   !      (rhow*cswat*(theta*dx(1)+h0(kk))+par(1)%rho*par(1)%css*dx(1))

                   Tsoil(1) = var(1)%Tfrz + (tmp1d2(kk) - tmp1d4(kk)) / &
                        (rhow*cswat*(theta*dx(1)+h0(kk))+par(1)%rho*par(1)%css*dx(1))

                   var(1)%iice = 0
                   var(1)%thetai = zero
                   var(1)%thetal = theta
                else
                 
                   ! frozen remaining frozen
                   Jsoil = tmp1d2(kk) ! total energy in  soil layer
                   !check there is a zero
                   tmp1 = GTfrozen(Tsoil(1)-50._r_2, Jsoil, dx(1), theta,par(1)%css, par(1)%rho, &
                        h0(kk), par(1)%thre, par(1)%the, par(1)%he, one/(par(1)%lambc*freezefac))

                   tmp2 = GTFrozen(Tfreezing(kk), Jsoil, dx(1), theta,par(1)%css, par(1)%rho, &
                        h0(kk), par(1)%thre, par(1)%the, par(1)%he, one/(par(1)%lambc*freezefac))

                   ! there is a zero in between
                   if ((tmp1*tmp2) < zero) then
                      tmp1d3(kk) = rtbis_Tfrozen(tmp1d2(kk), dx(1), theta,par(1)%css, par(1)%rho, &
                           h0(kk), par(1)%thre, par(1)%the, par(1)%he, one/(par(1)%lambc*freezefac), &
                           Tsoil(1)-50._r_2, Tfreezing(kk))

                      tmp1d4(kk) = thetalmax(tmp1d3(kk), S(1), par(1)%he, one/(par(1)%lambc*freezefac), &
                           par(1)%thre, par(1)%the) ! liquid content at new Tsoil
                   else
                      write(*,*) "Found no solution for Tfrozen 4.", irec, qmelt(1), h0(kk)
                      write(wlogn,*) "Assume soil is totally frozen"
                      tmp1d3(kk) = (Jsoil + rhow*lambdaf*(theta*dx(1) +  h0(kk))) / &
                           (dx(1)*par(1)%css*par(1)%rho + rhow*csice*(theta*dx(1) + h0(kk)))
                      tmp1d4(kk) = 0.0_r_2                
                      write(wlogn,*) "frozen soil temperature: ", tmp1d3(kk)
                   endif

                   var(1)%thetal = max(tmp1d4(kk), zero)
                   var(1)%thetai = theta - tmp1d4(kk)
                   hice_tmp(kk) = hice(kk)
                   hice = h0(kk)*var(1)%thetai/par(1)%thre

                   ! correct total energy stored in pond + soil
                   Jcol_latent_S(kk) = Jcol_latent_S(kk)- rhow*lambdaf*((hice(kk)-hice_tmp(kk)) + &
                        dx(1)*(var(1)%thetai-thetai(1)))

                   Jcol_sensible(kk) = Jcol_sensible(kk) - &
                        var(1)%csoil*dx(1)*(Tsoil(1)) - &
                        (h0_tmp(kk)-hice_tmp(kk))*cswat*rhow*(Tsoil(1)) - &
                        hice_tmp(kk)*rhow*(csice*(Tsoil(1))) + &
                        dx(1)*(tmp1d3(kk))*par(1)%rho*par(1)%css + &
                        (h0(kk)-hice(kk)+var(1)%thetal*dx(1))*cswat*rhow*(tmp1d3(kk)) + &
                        (hice_tmp(kk)+var(1)%thetai)*rhow*(csice*(tmp1d3(kk)))

                   thetai(1) = var(1)%thetai
                   Tsoil(1) = tmp1d3(kk)
                endif ! incomplete melting
             endif ! iice=1

             do i=1, vsnow(kk)%nsnow
                vsnow(kk)%Jsensible(i) = (vsnow(kk)%hsnow(i)-vsnow(kk)%hliq(i))*rhow*csice*(vsnow(kk)%Tsn(i))
                vsnow(kk)%Jlatent(i) = (vsnow(kk)%hsnow(i)-vsnow(kk)%hliq(i))*rhow*(-lambdaf)
             enddo
             vsnow(kk)%J = sum(vsnow(kk)%Jsensible(1:vsnow(kk)%nsnow)+vsnow(kk)%Jlatent(1:vsnow(kk)%nsnow))
             vsnow(kk)%deltaJ = vsnow(kk)%J - J0snow(kk)
             vsnow(kk)%FluxDivergence = vsnow(kk)%Qadv_rain + vsnow(kk)%Qadv_snow + vsnow(kk)%Qadv_vap + &
                  vsnow(kk)%Qadv_melt + vsnow(kk)%Qcond_net  +  vsnow(kk)%Qadv_transfer
          endif ! (melt_transfer==.true.)
       endif ! remove  melt water 

       do i=1, vsnow(kk)%nsnow
          vsnow(kk)%dens(i) = vsnow(kk)%hsnow(i)/vsnow(kk)%depth(i)*rhow
          if (vsnow(kk)%dens(i).lt.50._r_2) then
             vsnow(kk)%dens(i) = 50.0_r_2
             vsnow(kk)%depth(i) = vsnow(kk)%hsnow(i)*rhow/vsnow(kk)%dens(i)
          endif
          tmp1d3(kk) = (500._r_2*(vsnow(kk)%hsnow(i)-vsnow(kk)%hliq(i)) + rhow*vsnow(kk)%hliq(i))/vsnow(kk)%hsnow(i)
          if (vsnow(kk)%dens(i).gt.tmp1d3(kk)) then
             vsnow(kk)%dens(i) = tmp1d3(kk)
             vsnow(kk)%depth(i) = vsnow(kk)%hsnow(i)*rhow/vsnow(kk)%dens(i)
          endif

       enddo

       vsnow(kk)%totdepth = sum(vsnow(kk)%depth(1:vsnow(kk)%nsnow))

       ! adjust number of snow layers if required
       if (nsnow_max>1.and.vsnow(kk)%totdepth > 0.03_r_2) then
          tmp1d3(kk) = vsnow(kk)%depth(1)

          if (vsnow(kk)%totdepth.lt.0.06_r_2) then
             vsnow(kk)%depth(1) = vsnow(kk)%totdepth/2._r_2
          else
             vsnow(kk)%depth(1) = 0.03_r_2
          endif
          vsnow(kk)%depth(nsnow_max) = vsnow(kk)%totdepth -  vsnow(kk)%depth(1)

          if (vsnow(kk)%nsnow == 1) then
             !put excess into new layer below (initialise 2nd snow layer if necessary)
             vsnow(kk)%tsn(vsnow(kk)%nsnow+1) = vsnow(kk)%tsn(1)
             vsnow(kk)%dens(vsnow(kk)%nsnow+1) = vsnow(kk)%dens(1)
             vsnow(kk)%hsnow(1) = vsnow(kk)%depth(1)*(vsnow(kk)%dens(1)/rhow)
             vsnow(kk)%hsnow(vsnow(kk)%nsnow+1) = vsnow(kk)%depth(vsnow(kk)%nsnow+1)*(vsnow(kk)%dens(vsnow(kk)%nsnow+1)/rhow)
             vsnow(kk)%hliq(vsnow(kk)%nsnow+1)  = vsnow(kk)%hliq(1)*vsnow(kk)%hsnow(vsnow(kk)%nsnow+1)/(vsnow(kk)%hsnow(1)+ &
                  vsnow(kk)%hsnow(vsnow(kk)%nsnow+1))
             vsnow(kk)%hliq(1) = vsnow(kk)%hliq(1) - vsnow(kk)%hliq(vsnow(kk)%nsnow+1)
             vsnow(kk)%nsnow = 2
          else
             ! recalculate tsn, dens, hsnow, hliq according to new layer depths
             ! total energy of combined snow pack
             tmp1d1(kk) = rhow*(vsnow(kk)%hsnow(1)-vsnow(kk)%hliq(1))*(csice*vsnow(kk)%tsn(1) - lambdaf) + &
                  rhow*(vsnow(kk)%hsnow(vsnow(kk)%nsnow)-vsnow(kk)%hliq(vsnow(kk)%nsnow))* &
                  (csice*vsnow(kk)%tsn(vsnow(kk)%nsnow) - lambdaf)


             if (vsnow(kk)%depth(1).le.tmp1d3(kk)) then
                i= 1 ! if top layer decreases in depth, it retains old density
                j=2
             else
                i = 2 ! if bottom layer decreases in depth, it retains old density
                j=1
             endif
             vsnow(kk)%hliq(i)= vsnow(kk)%hliq(i)*vsnow(kk)%depth(i)*(vsnow(kk)%dens(i)/rhow)/vsnow(kk)%hsnow(i)
             vsnow(kk)%hsnow(j)= (vsnow(kk)%hsnow(i)+vsnow(kk)%hsnow(j)) - vsnow(kk)%depth(i)*(vsnow(kk)%dens(i)/rhow)
             vsnow(kk)%hsnow(i)= vsnow(kk)%depth(i)*(vsnow(kk)%dens(i)/rhow)

             ! energy content of augmented layer
             tmp1d2(kk) = tmp1d1(kk) - &
                  rhow*(vsnow(kk)%hsnow(i)-vsnow(kk)%hliq(i))*(csice*vsnow(kk)%tsn(i)- lambdaf)
             if (tmp1d2(kk).lt. -vsnow(kk)%hsnow(j)*rhow*lambdaf) then
                vsnow(kk)%hliq(j)= zero
                vsnow(kk)%tsn(j)= (tmp1d2(kk)/(rhow*vsnow(kk)%hsnow(j))+lambdaf)/csice
             else
                vsnow(kk)%tsn(j)= zero
                vsnow(kk)%hliq(j)= vsnow(kk)%hsnow(j)+ tmp1d2(kk)/(rhow*lambdaf)
             endif

             ! density of augmented layer
             vsnow(kk)%dens(j)= rhow * vsnow(kk)%hsnow(j)/vsnow(kk)%depth(j)

          endif

       elseif (nsnow_max>1.and.vsnow(kk)%nsnow == 2.and.(vsnow(kk)%totdepth).le.0.03_r_2) then
          ! combine top two snow layers into 1

          ! total energy of combined snow layer
          tmp1d1(kk) = rhow*(vsnow(kk)%hsnow(1)-vsnow(kk)%hliq(1))*(csice*vsnow(kk)%tsn(1) - lambdaf) + &
               rhow*(vsnow(kk)%hsnow(vsnow(kk)%nsnow)-vsnow(kk)%hliq(vsnow(kk)%nsnow))* &
               (csice*vsnow(kk)%tsn(vsnow(kk)%nsnow) - lambdaf)
          if (tmp1d1(kk).lt. -(vsnow(kk)%hsnow(1)+vsnow(kk)%hsnow(vsnow(kk)%nsnow))*rhow*lambdaf) then
             vsnow(kk)%hliq(1) = zero
             vsnow(kk)%tsn(1) = (tmp1d1(kk)/(rhow*(vsnow(kk)%hsnow(1)+vsnow(kk)%hsnow(vsnow(kk)%nsnow)))+lambdaf)/csice
          else
             vsnow(kk)%tsn(1) = zero
             vsnow(kk)%hliq(1) = (vsnow(kk)%hsnow(1)+vsnow(kk)%hsnow(vsnow(kk)%nsnow)) + tmp1d1(kk)/(rhow*lambdaf)
          endif

          vsnow(kk)%tsn(vsnow(kk)%nsnow) = zero
          vsnow(kk)%hliq(vsnow(kk)%nsnow) = zero
          vsnow(kk)%hsnow(1) = vsnow(kk)%hsnow(1) + vsnow(kk)%hsnow(vsnow(kk)%nsnow)
          vsnow(kk)%hsnow(vsnow(kk)%nsnow) = 0
          vsnow(kk)%depth(1) = vsnow(kk)%depth(1) + vsnow(kk)%depth(vsnow(kk)%nsnow)
          vsnow(kk)%depth(vsnow(kk)%nsnow) = 0
          vsnow(kk)%dens(1) = vsnow(kk)%hsnow(1)/vsnow(kk)%depth(1)*rhow
          vsnow(kk)%nsnow = 1
       endif

       if (vsnow(kk)%hsnow(1).lt.zero.or.vsnow(kk)%hsnow(nsnow_max).lt.zero) then
          write(wlogn,*) "hsnow<0. Set it to 0 (irec, kk, hsnow):", irec, kk, vsnow(kk)%hsnow(1)
          vsnow(kk)%hsnow(1) = zero
       endif

       do i=1, vsnow(kk)%nsnow
          vsnow(kk)%Dv(i) = Dva*((vsnow(kk)%tsn(i)+Tzero)/Tzero)**1.88_r_2 ! m2 s-1
          vsnow(kk)%sl(i) = slope_esat_ice(vsnow(kk)%tsn(i)) * Mw/thousand/Rgas/(vsnow(kk)%tsn(i)+Tzero)
          vsnow(kk)%kE(i)     = vsnow(kk)%Dv(i)*vsnow(kk)%sl(i)*thousand*lambdaf
          vsnow(kk)%kH(i) = 3.2217e-6_r_2 * vsnow(kk)%dens(i)**2
          vsnow(kk)%kth(i) = vsnow(kk)%kE(i) + vsnow(kk)%kH(i)
          vsnow(kk)%cv(i) = esat_ice(vsnow(kk)%tsn(i))*Mw/thousand/Rgas/(vsnow(kk)%tsn(i)+Tzero) ! m3 m-3
          vsnow(kk)%depth(i) = vsnow(kk)%hsnow(i)/(vsnow(kk)%dens(i)/rhow)
       enddo

    endif ! dedicated snow layer

  END SUBROUTINE snow_adjust

  !*********************************************************************************************************************

  SUBROUTINE isotope_vap(irec,isotopologue, n, nsnow, nsnow_last,  & ! scalar in
       ns, dx, deltaz, sig, dt, &     ! in soil
       Tsoil0, dTsoil, Sliqice, deltaSliqice, & ! in variables
       Sliq, deltaSliq,Sice, deltaSice, &
       Ts, Ta, &
       qsig, qlsig, qvsig,  & ! in fluxes
       qmelt, qtransfer, & ! melt water to soil and water transfer from soil to snow (+ve) or snow to soil (-ve)
       qprec, qprec_snow, qevap, qrunoff, qex, & ! in fluxes
       var_cv, var_Dv, &   ! in variables
       thetasat, thetar, tortuosity, deltacv, rbw, cva, civa, & ! in parameter
       cprec,cprec_snow, cali, & ! in iso
       ql0, qv0, & ! inout water
       ciso, cisoice, cisos, & ! inout iso
       qiso_in, qiso_out, qiso_evap, qiso_trans, qiso_liq_adv, qiso_vap_adv, qiso_liq_diff, qiso_vap_diff) ! out iso


    IMPLICIT NONE
    INTEGER(i_d),                    INTENT(IN)    :: irec
    INTEGER(i_d),                    INTENT(IN)    :: isotopologue ! which isotope
    INTEGER(i_d),                    INTENT(IN)    :: n            ! # of soil layers
    INTEGER(i_d),                    INTENT(IN)    :: nsnow        ! # of snow layers
    INTEGER(i_d),                    INTENT(IN)    :: nsnow_last   ! # of snow layers
    INTEGER(i_d),                    INTENT(IN)    :: ns           ! index of top of soil/snow column (-nsnow_max+1)
    REAL(r_2),    DIMENSION(ns:n),   INTENT(IN)    :: dx           ! soil/snow depths
    REAL(r_2),    DIMENSION(ns:n-1), INTENT(IN)    :: deltaz       ! soil/snow layer thickness
    REAL(r_2),                       INTENT(IN)    :: sig          ! implicit/explicit time steping constant
    REAL(r_2),                       INTENT(IN)    :: dt           ! time step
    REAL(r_2),    DIMENSION(ns:n),   INTENT(IN)    :: Tsoil0       ! soil/snow temperatures
    REAL(r_2),    DIMENSION(ns:n),   INTENT(IN)    :: dTsoil       ! soil/snow temperatures change
    REAL(r_2),    DIMENSION(ns:n),   INTENT(IN)    :: Sliqice      ! soil saturation
    REAL(r_2),    DIMENSION(ns:n),   INTENT(IN)    :: deltaSliqice ! soil sturation change
    REAL(r_2),    DIMENSION(ns:n),   INTENT(IN)    :: Sliq         ! soil saturation
    REAL(r_2),    DIMENSION(ns:n),   INTENT(IN)    :: deltaSliq    ! soil sturation change
    REAL(r_2),    DIMENSION(ns:n),   INTENT(IN)    :: Sice         ! soil saturation
    REAL(r_2),    DIMENSION(ns:n),   INTENT(IN)    :: deltaSice    ! soil sturation change
    REAL(r_2),                       INTENT(IN)    :: Ts           ! surface temperature
    REAL(r_2),                       INTENT(IN)    :: Ta           ! air temperature [degC]
    REAL(r_2),    DIMENSION(ns-1:n), INTENT(IN)    :: qsig         ! water flux
    REAL(r_2),    DIMENSION(ns-1:n), INTENT(INOUT) :: qlsig        ! liquid water flux
    REAL(r_2),    DIMENSION(ns-1:n), INTENT(INOUT) :: qvsig        ! vapour flux
    REAL(r_2),                       INTENT(IN)    :: qprec        ! liquid  precip
    REAL(r_2),                       INTENT(IN)    :: qprec_snow   ! solid precip
    REAL(r_2),                       INTENT(IN)    :: qevap        ! evaporation
    REAL(r_2),                       INTENT(IN)    :: qrunoff      ! runoff
    REAL(r_2),    DIMENSION(1:n),    INTENT(IN)    :: qex          ! root extraction
    REAL(r_2),                       INTENT(IN)    :: qmelt        !  melt water ! should eventually be dimensioned nsnow_max
    REAL(r_2),                       INTENT(IN)    :: qtransfer    ! water transfer from soil to snow (+ve) or snow to soil (-ve)
    ! TYPE(vars), DIMENSION(1:n),   INTENT(IN)    :: var           ! soil variables
    REAL(r_2), DIMENSION(ns:n),   INTENT(IN)    :: var_cv        ! soil/snow variables cv
    REAL(r_2), DIMENSION(ns:n),   INTENT(IN)    :: var_Dv        ! soil/snow variables Dv
    REAL(r_2), DIMENSION(ns:n),   INTENT(IN)    :: thetasat      ! saturation moisture
    REAL(r_2), DIMENSION(ns:n),   INTENT(IN)    :: thetar        ! residual moisture
    REAL(r_2), DIMENSION(ns:n),   INTENT(IN)    :: tortuosity    ! soil tortuosity
    REAL(r_2), DIMENSION(ns:n),   INTENT(INOUT) :: deltacv       !
    REAL(r_2),                    INTENT(IN)    :: rbw           ! boundary layer resistance
    REAL(r_2),                    INTENT(IN)    :: cva           !
    REAL(r_2),                    INTENT(IN)    :: civa          ! iso conc in air
    REAL(r_2),                    INTENT(IN)    :: cprec         ! iso conc in liquid precip
    REAL(r_2),                    INTENT(IN)    :: cprec_snow    ! iso conc in snow
    REAL(r_2),                    INTENT(IN)    :: cali          ! iso conc in alimenation water (from below)
    REAL(r_2),                    INTENT(INOUT) :: ql0           ! liquid flux into soil or snow surface
    REAL(r_2),                    INTENT(INOUT) :: qv0           ! vapour flux into soil or snow surface
    REAL(r_2), DIMENSION(ns:n),   INTENT(INOUT) :: ciso          ! iso conc in soil
    REAL(r_2), DIMENSION(ns:n),   INTENT(INOUT) :: cisoice       ! iso conc in soil
    REAL(r_2),                    INTENT(INOUT) :: cisos         ! iso conc on surface
    REAL(r_2),                    INTENT(OUT)   :: qiso_in       ! iso flux into soil
    REAL(r_2),                    INTENT(OUT)   :: qiso_out      ! iso flux out of soil
    REAL(r_2),                    INTENT(OUT)   :: qiso_evap     ! iso flux of evaporation
    REAL(r_2),                    INTENT(OUT)   :: qiso_trans    ! iso flux of transpiration
    REAL(r_2), DIMENSION(ns:n),   INTENT(OUT)   :: qiso_liq_adv  ! liquid iso flux in soil due to advection
    REAL(r_2), DIMENSION(ns:n),   INTENT(OUT)   :: qiso_vap_adv  ! vapour iso flux in soil due to advection
    REAL(r_2), DIMENSION(ns:n-1), INTENT(OUT)   :: qiso_liq_diff ! liquid iso flux in soil due to diffusion
    REAL(r_2), DIMENSION(ns:n-1), INTENT(OUT)   :: qiso_vap_diff ! vapour iso flux in soil due to diffusion

    ! Local variables

    REAL(r_2), DIMENSION(ns:n)   :: aa, bb, cc, dd, dc,  LHS, RHS
    REAL(r_2), DIMENSION(ns:n)   :: alphaplus, dalphaplusdT, deltaT, beta, deltabeta
    REAL(r_2), DIMENSION(ns:n)   :: alphaplus_liqice, dalphaplusdT_liqice
    ! diffusivities (liquid, liq-vap, coefft for D, surface H2O vapour, surface minor isotopologue vapour)
    REAL(r_2), DIMENSION(ns-1:n) :: Dl, Dv
    ! REAL(r_2)                    :: Dvs, Divs
    REAL(r_2)                    :: patm, nk, alphak, alphak_vdiff, alphak_ldiff
    REAL(r_2)                    :: cevapin, cevapout, qevapin, qevapout, dcevapoutdciso
    ! concentrations of advective fluxes and corresponding partial derivs wrt ciso
    REAL(r_2), DIMENSION(ns-1:n) :: cql, dcqldca, dcqldcb, cqv, dcqvdca, dcqvdcb
    REAL(r_2), DIMENSION(ns:n-1) :: wcql, wcqv
    REAL(r_2), DIMENSION(ns-1:n) :: betaqv, dbetaqv
    REAL(r_2), DIMENSION(ns-1:n) :: Dlmean, Dvmean, Dvbetamean, wl, wv
    REAL(r_2)                    :: coefA, coefB, coefC
    REAL(r_2)                    :: coefA_liqice, coefB_liqice, coefC_liqice
    REAL(r_2), DIMENSION(ns:n)   :: Seff, deltaSeff, S, Tsoil, cvsig, Sliqsig, Sicesig,qex_ss
    REAL(r_2), DIMENSION(ns:n)   :: thetaice, deltathetaice, dcice
    INTEGER(i_d)                 :: ns_ciso
    REAL(r_2)                    :: num, den, cv1
    REAL(r_2)                    :: alphaplus_s, alphaplus_a
    ! REAL(r_2)                    :: alphaplus_liqice
    REAL(r_2)                    :: cvs !, qevapL, qevapoutL, qevapinL
    ! REAL(r_2)                    :: cevapinL, cevapoutL, dcevapoutdcisoL, dcevapindcisoL
    REAL(r_2)                    :: w1, w2
    INTEGER(i_d), PARAMETER      :: formulation = 2  ! 1: betaql, betaqv of flux; 2: betaql, betaqv 0.5 of upper and lower
    REAL(r_2), DIMENSION(ns:n)   :: kfreeze, kfreeze2         ! combination of freezing variables
    REAL(r_2),   DIMENSION(ns:n) :: deltaS     !
    INTEGER(i_d), DIMENSION(1)   :: ii
    REAL(r_2)                    :: tmp
    REAL(r_2), DIMENSION(ns:n)   :: tmp1d1

    dcice = zero
    aa = zero
    bb = zero
    cc = zero
    dd = zero
    alphaplus = zero
    dalphaplusdT = zero
    alphaplus_liqice = zero
    dalphaplusdT_liqice = zero
    deltaT = zero
    Dl = zero
    Dv = zero
    cql = zero
    dcqldca = zero
    dcqldcb = zero
    cqv = zero
    dcqvdca = zero
    dcqvdcb = zero
    wcql = zero
    wcqv = zero
    beta = zero
    deltabeta = zero
    betaqv = zero
    dbetaqv = zero
    Dlmean = zero
    Dvmean = zero
    wl = zero
    wv = zero
    Seff = zero
    deltaSeff = zero
    S = zero
    Tsoil = zero
    cvsig = zero
    Sliqsig = zero
    Sicesig = zero
    thetaice = zero
    deltathetaice = zero
    dcice = zero
    kfreeze = zero
    kfreeze2 = zero
    qiso_liq_adv = zero
    qiso_vap_adv = zero
    qiso_liq_diff = zero
    qiso_vap_diff = zero
    dc = zero
    dcice = zero

    if (nsnow_last == 0) then
       ciso(0) = ciso(1)
       cisoice(0) = cisoice(1)
    endif

    ns_ciso = -nsnow+1 ! index of top layer (<=0 = snow)

    if ((nsnow_last.eq.0) .and. (ns_ciso.le.0)) then
       ciso(ns_ciso:0) = ciso(1)
       cisoice(ns_ciso:0) = cisoice(1)
    endif
    patm = one
    deltaT(ns_ciso:n) = dTsoil(ns_ciso:n)
    qex_ss(1:n) = qex(1:n)
    if (ns_ciso<1) then
       qex_ss(ns_ciso:0) = zero
    endif

    ! Tsoil to Tsoilsig
    Tsoil(ns_ciso:n) = Tsoil0(ns_ciso:n) + (sig-one)*deltaT(ns_ciso:n)

    ! equilibrium fractionation factors at the surface and in the soil
    coefA = zero
    coefB = zero
    coefC = zero
    if (isotopologue==1) then
       coefA = 24844.0_r_2
       coefB = -76.248_r_2
       coefC = 0.052612_r_2
       ! alphaplus_liqice = 1.0212_r_2
       coefA_liqice = 48888_r_2
       coefB_liqice = -203.10_r_2
       coefC_liqice = 0.2133_r_2
    endif
    if (isotopologue==2) then
       coefA = 1137.0_r_2
       coefB = -0.4156_r_2
       coefC = -0.0020667_r_2
       ! alphaplus_liqice = 1.00291_r_2
       coefA_liqice = 8312.5_r_2
       coefB_liqice = -49.192_r_2
       coefC_liqice = 0.0831_r_2
    endif
    !MC - Test
    ! alphaplus_a                    = one/exp(coefA/((Ta+Tzero)**2)+coefB/(Ta+Tzero)+coefC) ! at soil or litter surface
    ! alphaplus_s                    = one/exp(coefA/((Ts+Tzero)**2)+coefB/(Ts+Tzero)+coefC) ! at soil or litter surface
    ! alphaplus(ns_ciso:n)           = one/exp(coefA/((Tsoil(ns_ciso:n)+Tzero)**2)+coefB/(Tsoil(ns_ciso:n)+Tzero)+coefC)
    ! dalphaplusdT(ns_ciso:n)        = (two*coefA/(Tsoil(ns_ciso:n)+Tzero)**3 + coefB/(Tsoil(ns_ciso:n)+Tzero)**2) &
    !      / exp(coefA/((Tsoil(ns_ciso:n)+Tzero)**2)+coefB/(Tsoil(ns_ciso:n)+Tzero)+coefC)
    ! alphaplus_liqice(ns_ciso:n)    = one/exp(coefA_liqice/((Tsoil(ns_ciso:n)+Tzero)**2)+ &
    !      coefB_liqice/(Tsoil(ns_ciso:n)+Tzero)+coefC_liqice)
    ! dalphaplusdT_liqice(ns_ciso:n) = (two*coefA_liqice/(Tsoil(ns_ciso:n)+Tzero)**3 + &
    !      coefB_liqice/(Tsoil(ns_ciso:n)+Tzero)**2) &
    !      / exp(coefA_liqice/((Tsoil(ns_ciso:n)+Tzero)**2)+coefB_liqice/(Tsoil(ns_ciso:n)+Tzero)+coefC_liqice)
    ! alphaplus_liqice(ns_ciso:n) = one/alphaplus_liqice(ns_ciso:n)  ! needs to be > 1
    tmp = max(Ta, zero)
    alphaplus_a                    = one/exp(coefA/((tmp+Tzero)**2)+coefB/(tmp+Tzero)+coefC)
    tmp = max(Ts, zero)
    alphaplus_s                    = one/exp(coefA/((tmp+Tzero)**2)+coefB/(tmp+Tzero)+coefC)
    tmp1d1(ns_ciso:n) = max(Tsoil(ns_ciso:n), zero)
    alphaplus(ns_ciso:n)           = one/exp(coefA/((tmp1d1(ns_ciso:n)+Tzero)**2)+coefB/(tmp1d1(ns_ciso:n)+Tzero)+coefC)
    dalphaplusdT(ns_ciso:n)        = (two*coefA/(tmp1d1(ns_ciso:n)+Tzero)**3 + coefB/(tmp1d1(ns_ciso:n)+Tzero)**2) &
         / exp(coefA/((tmp1d1(ns_ciso:n)+Tzero)**2)+coefB/(tmp1d1(ns_ciso:n)+Tzero)+coefC)
    tmp1d1(ns_ciso:n) = min(Tsoil(ns_ciso:n), zero)
    alphaplus_liqice(ns_ciso:n)    = one/exp(coefA_liqice/((tmp1d1(ns_ciso:n)+Tzero)**2)+ &
         coefB_liqice/(tmp1d1(ns_ciso:n)+Tzero)+coefC_liqice)
    dalphaplusdT_liqice(ns_ciso:n) = (two*coefA_liqice/(tmp1d1(ns_ciso:n)+Tzero)**3 + &
         coefB_liqice/(tmp1d1(ns_ciso:n)+Tzero)**2) &
         / exp(coefA_liqice/((tmp1d1(ns_ciso:n)+Tzero)**2)+coefB_liqice/(tmp1d1(ns_ciso:n)+Tzero)+coefC_liqice)
    alphaplus_liqice(ns_ciso:n) = one/alphaplus_liqice(ns_ciso:n)  ! needs to be > 1

    ! if (experiment==1 .or. experiment==2.or. experiment==9.or. experiment==10 .or. experiment == 16) then
    if (experiment==1 .or. experiment==2.or. experiment==9.or. experiment==10) then
       alphaplus_s      = one
       alphaplus        = one
       dalphaplusdT     = zero
       alphaplus_liqice = one
    endif
    !MC - Test
    ! alphaplus_s      = one
    ! alphaplus        = one
    ! dalphaplusdT     = zero
    ! alphaplus_liqice = one

    ! beta
    beta(ns_ciso:n)      = alphaplus(ns_ciso:n)
    deltabeta(ns_ciso:n) = dalphaplusdT(ns_ciso:n)*deltaT(ns_ciso:n)
    beta(ns_ciso:n)      = beta(ns_ciso:n) + sig*deltabeta(ns_ciso:n) ! beta_sig

    ! adjust S and h and cv for sig of time-step
    S(ns_ciso:n)       = Sliqice(ns_ciso:n) + deltaSliqice(ns_ciso:n)*(sig-one) ! set S to Ssig. N.B. S = Sliqice
    Sliqsig(ns_ciso:n) = Sliq(ns_ciso:n) + deltaSliq(ns_ciso:n)*(sig-one)
    Sicesig(ns_ciso:n) = Sice(ns_ciso:n) + deltaSice(ns_ciso:n)*(sig-one)

    thetaice(1:n)      = Sice(1:n) * (thetasat(1:n)-thetar(1:n))
    deltathetaice(1:n) = deltaSice(1:n) * (thetasat(1:n)-thetar(1:n))

    kfreeze(ns_ciso:n)  = merge(alphaplus_liqice(ns_ciso:n)*deltaSice(ns_ciso:n)/Sice(ns_ciso:n), &
         zero, deltaSice(ns_ciso:n) > zero)
    kfreeze2(ns_ciso:n) = merge(deltaSice(ns_ciso:n)/Sice(ns_ciso:n)*(alphaplus_liqice(ns_ciso:n)*ciso(ns_ciso:n) - &
         cisoice(ns_ciso:n)), zero, deltaSice(ns_ciso:n) > zero)
    ! adapt kfreeze for totally frozen snow
    if (nsnow > 0) then
       where (.not. ((Sliq(ns_ciso:0) > zero) .or. (deltaSliq(ns_ciso:0) > zero))) ! solve for change in isotopes in ice
          kfreeze(ns_ciso:0)  = one
          kfreeze2(ns_ciso:0) = zero
          ciso(ns_ciso:0)     = cisoice(ns_ciso:0)
       endwhere
       ! composition of new ice is composition of snow
       if  (((Sliq(ns_ciso) > zero) .or. (deltaSliq(ns_ciso) > zero)) .and. (qprec_snow > zero)) then
          kfreeze(ns_ciso)  =  zero
          kfreeze2(ns_ciso) = deltaSice(ns_ciso)/Sice(ns_ciso)*(cprec_snow-cisoice(ns_ciso))
       endif
    endif
    ! composition of new ice in surface layer is composition of snow
    if  (((Sliq(ns_ciso) > zero) .or. (deltaSliq(ns_ciso) > zero)) .and. &
         (qprec_snow > zero) .and. (Sice(ns_ciso) > zero) .and. (deltaSice(ns_ciso) > zero)) then
       kfreeze(ns_ciso)  =  zero
       kfreeze2(ns_ciso) = deltaSice(ns_ciso)/Sice(ns_ciso)*(cprec_snow-cisoice(ns_ciso))
    endif
    if ((Sice(1) <= zero) .and. (deltaSice(1) > zero)) then
       write(*,*) "deltaSice inconistent with Sice"
    endif

    ! calculate deltacv which is consistent with deltaS and deltaSliqice
    deltaS(ns_ciso:n) = ( qlsig(ns_ciso-1:n-1) + qvsig(ns_ciso-1:n-1) &
         - qlsig(ns_ciso:n) - qvsig(ns_ciso:n) - qex_ss(ns_ciso:n) ) * &
         dt/(dx(ns_ciso:n)*thetasat(ns_ciso:n))

    ! modify for snow melt (including disappearance of snow pack) and/or transfer from top soil to new snow pack
    ! valid for 3 snow layers ???
    if ((qmelt.gt.zero) .or. ((qtransfer).gt.zero)) then
       deltaS(0) = (qlsig(-1) + qvsig(-1) - (qlsig(0)+qmelt-qtransfer) - qvsig(0))*dt/dx(0)/thetasat(0)
       deltaS(1) = (qlsig(0)+qmelt-qtransfer + qvsig(0) - qlsig(1)- qvsig(1)-qex_ss(1))*dt/dx(1)/thetasat(1)
    endif
    ! modify deltaS in top layer for precip
    deltaS(ns_ciso) = deltaS(ns_ciso) + (qprec+qprec_snow)*dt/dx(ns_ciso)/thetasat(ns_ciso)
    ! modify deltaS in top soil layer for runoff
    deltaS(1) = deltaS(1) - qrunoff*dt/thetasat(1)/dx(1)

    where (var_cv(ns_ciso:n).gt.zero)
       deltacv(ns_ciso:n) = (deltaS(ns_ciso:n) - deltaSliqice(ns_ciso:n) + var_cv(ns_ciso:n)*deltaSliqice(ns_ciso:n))/ &
            (1._r_2 - (S(ns_ciso:n) + deltaSliqice(ns_ciso:n)*(one-sig)))
    elsewhere
       deltacv(ns_ciso:n) = zero
    endwhere
    cvsig(ns_ciso:n) = var_cv(ns_ciso:n) + deltacv(ns_ciso:n)*(sig-one)

    where (S(ns_ciso:n)<one)
       Seff(ns_ciso:n) = Sliqsig(ns_ciso:n) + (one-S(ns_ciso:n))*cvsig(ns_ciso:n)*beta(ns_ciso:n) ! &
       Seff(ns_ciso:n) = Seff(ns_ciso:n)   + thetar(ns_ciso:n)/(thetasat(ns_ciso:n)-thetar(ns_ciso:n))
       deltaSeff(ns_ciso:n) = deltaSliq(ns_ciso:n) + beta(ns_ciso:n)*deltacv(ns_ciso:n) &
            + cvsig(ns_ciso:n)*deltabeta(ns_ciso:n) &
            - (S(ns_ciso:n)*beta(ns_ciso:n)*deltacv(ns_ciso:n) + cvsig(ns_ciso:n)*beta(ns_ciso:n)*deltaSliqice(ns_ciso:n) + &
            cvsig(ns_ciso:n)*S(ns_ciso:n)*deltabeta(ns_ciso:n))
    elsewhere
       Seff(ns_ciso:n)      = Sliqsig(ns_ciso:n)
       Seff(ns_ciso:n)      = Seff(ns_ciso:n) + thetar(ns_ciso:n)/thetasat(ns_ciso:n)
       deltaSeff(ns_ciso:n) = deltaSliq(ns_ciso:n)
    endwhere

    ! diffusional fractionation factor
    ! air
    if (isotopologue==1) alphak_vdiff = one / 1.0251_r_2   ! HDO diffusivity in air (Merlivat 1978)
    if (isotopologue==2) alphak_vdiff = one / 1.0285_r_2   ! H218O diffusivity in air (Merlivat 1978)
    if ((experiment >= 1 .and. experiment <= 5) .or. (experiment == 9 .or. experiment == 10)) alphak_vdiff = one
    ! if (experiment == 9 .or. experiment == 10 .or. experiment == 16) alphak_vdiff = one
    ! kinetic fractionation factor at the surface
    ! Dvs = Dva*1.e5_r_2/patm*((Ts+Tzero)/Tzero)**1.88_r_2 ! vapour diffusivity of water in air (m2s-1)
    ! if (isotopologue/=0) Divs = Dvs * alphak_vdiff  ! isotope diffusivity in air
    nk = min(S(1),one)*half + (one-min(S(1),one))*one
    if (experiment == 7 .or. experiment == 8) nk = one
    alphak = alphak_vdiff**nk ! = one/((Dvs/Divs)**nk)
    if ((experiment >= 1 .and. experiment <= 5) .or. (experiment == 9 .or. experiment == 10)) alphak = one
    ! if (experiment == 9 .or. experiment == 10 .or. experiment == 16) alphak = one

    ! liquid diffusivity in the pond and soil
    if (isotopologue==1) alphak_ldiff = one / 1.013_r_2
    if (isotopologue==2) alphak_ldiff = one / 1.026_r_2
    ! molecular diffusion of isotopes in normal liquid water (m2s-1)
    Dl(ns_ciso:n) = tortuosity(ns_ciso:n) * alphak_ldiff * 1.0e-7_r_2*exp(-577.0_r_2/((Tsoil(ns_ciso:n)+Tzero)-145._r_2))
    Dl(ns_ciso:n) = Dl(ns_ciso:n) * (min(Sliqsig(ns_ciso:n),one) * (thetasat(ns_ciso:n)-thetar(ns_ciso:n)) + thetar(ns_ciso:n))
    if ((experiment >= 1 .and. experiment <= 4) .or. (experiment == 9 .or. experiment == 10)) Dl = zero
    ! if (experiment == 9 .or. experiment == 10 .or. experiment == 16) Dl = zero

    ! vapour diffusivity in the soil
    Dv(ns_ciso:n) = var_Dv(ns_ciso:n) * alphak_vdiff  ! isotope diffusivity in soil air spaces Dvi
    !MC - Cannot remember if this was a test or not
    Dv(ns_ciso:n) = Dv(ns_ciso:n) * max(one-S(ns_ciso:n),zero)*(thetasat(ns_ciso:n)-thetar(ns_ciso:n))
    ! Dvi*cv always together on RHS
    Dv(ns_ciso:n) = Dv(ns_ciso:n) * cvsig(ns_ciso:n)

    ! upper boundary condition
    cvs = cva + qevap*rbw ! concentration of water vapour at soil/air interface (m3 (H2O liq)/ m3 (air))

    if (var_Dv(ns_ciso) /= zero) then
       cv1 = -qv0*(half*dx(ns_ciso))/var_Dv(ns_ciso) + cvs
    else
       cv1 = var_cv(ns_ciso)
    endif

    if (ql0 > zero) then
       w1 = zero
    else
       w1 = one
    endif
    if (qv0 > zero) then
       w2 = zero
    else
       w2 = one
    endif

    !MC - Test - alphak=1
    ! alphak = one
    if (nsnow .ge. 1) then
       !MC - Test
       ! cisos          = ciso(ns_ciso)
       ! dcevapoutdciso = one
       ! retain in incoming cisos
       num            = zero
       den            = zero
       dcevapoutdciso = alphak*alphaplus_s
    else
       num   = alphak_vdiff*var_Dv(1)/(half*dx(ns_ciso))*cv1*alphaplus(ns_ciso)*ciso(ns_ciso) &
            - ql0*ciso(ns_ciso)*w1 +civa*alphak/rbw + Dl(ns_ciso)*ciso(ns_ciso)/(half*dx(ns_ciso))
       den   = alphak*cvs*alphaplus_s/(rbw) + ql0*(one-w1) + &
            alphaplus(ns_ciso)*alphak_vdiff*cvs*var_Dv(1)/(half*dx(ns_ciso)) +Dl(ns_ciso)/(half*dx(ns_ciso))
       cisos = num/den
       dcevapoutdciso = alphak*alphaplus_s
       dcevapoutdciso = dcevapoutdciso*(alphak_vdiff*var_Dv(1)/(half*dx(ns_ciso))*var_cv(ns_ciso)*alphaplus(ns_ciso) - &
            ql0*w1  + Dl(ns_ciso)/(half*dx(ns_ciso)))/den
    endif

    qevapin  = cva/rbw
    qevapout = cvs/rbw
    !MC - Test
    ! cevapin  = alphak * civa/cva
    cevapin  = alphak * civa/cva / alphaplus_a
    cevapout = alphak * alphaplus_s * cisos

    ! concentrations_ciso of advective fluxes and corresponding partial derivs wrt ciso
    select case (formulation)
    case (1)
       cql(ns_ciso-1)         = ciso(ns_ciso)
       cql(ns_ciso:n-1)       = merge(ciso(ns_ciso:n-1), ciso(ns_ciso+1:n), qlsig(ns_ciso:n-1)>zero)
       dcqldca(ns_ciso-1:n-1) = merge(one,  zero, qlsig(ns_ciso-1:n-1)>zero)
       dcqldcb(ns_ciso-1:n-1) = merge(zero, one,  qlsig(ns_ciso-1:n-1)>zero)
       cql(n)                 = ciso(n)
       dcqldca(n)             = one
       dcqldcb(n)             = zero

       cqv(ns_ciso-1)         = ciso(ns_ciso) *beta(ns_ciso)
       cqv(ns_ciso:n-1)       = merge(ciso(ns_ciso:n-1), ciso(ns_ciso+1:n), qvsig(ns_ciso:n-1)>zero)
       dcqvdca(ns_ciso-1:n-1) = merge(one,  zero, qvsig(ns_ciso-1:n-1)>zero)
       dcqvdcb(ns_ciso-1:n-1) = merge(zero, one,  qvsig(ns_ciso-1:n-1)>zero)
       cqv(n)                 = ciso(n)
       dcqvdca(n)             = one
       dcqvdcb(n)             = zero

       betaqv(ns_ciso-1)      = beta(ns_ciso)  *alphak_vdiff
       betaqv(ns_ciso:n-1)    = merge(beta(ns_ciso:n-1), beta(ns_ciso:n), qvsig(ns_ciso:n-1)>0)   * alphak_vdiff
       betaqv(n)              = beta(n)  * alphak_vdiff

       dbetaqv(ns_ciso:n-1)   = (beta(ns_ciso:n-1) - beta(ns_ciso+1:n))/deltaz(ns_ciso:n-1)
       dbetaqv(ns_ciso-1)     = zero ! vh ????
       dbetaqv(n)             = zero

       ! modify for snow melt (including disappearance of snow pack) and/or transfer from top soil to new snow pack
       if ((qmelt.gt.zero).or.((qtransfer).gt.zero)) then
          cql(0)     = (qlsig(0)*cql(0)+qmelt*ciso(0)+ciso(1)*(-qtransfer))/(qlsig(0)+qmelt-qtransfer)
          dcqldcb(0) = merge(zero, one, qlsig(0)>zero)*qlsig(0)/(qlsig(0)+qmelt -qtransfer)
          dcqldca(0) = merge(one,  zero, qlsig(0)>zero)*qlsig(0)/(qlsig(0)+qmelt -qtransfer)
          qlsig(0)   = qlsig(0) + qmelt - qtransfer
       endif

       if (experiment == 16) then
          ! print*, 'Ha03 ', betaqv
          ! print*, 'Ha04 ', dbetaqv
          ! betaqv  = one
          ! dbetaqv = zero
       endif

    case (2)
       cql(ns_ciso-1)         = ciso(ns_ciso)
       cql(ns_ciso:n-1)       = half*(ciso(ns_ciso:n-1)+ciso(ns_ciso+1:n))
       dcqldca(ns_ciso-1:n-1) = half
       dcqldcb(ns_ciso-1:n-1) = half
       cql(n)                 = ciso(n)
       dcqldca(n)             = one
       dcqldcb(n)             = zero

       !cqv(ns_ciso-1)         = ciso(ns_ciso) *beta(ns_ciso)
       !cqv(ns_ciso:n-1)       = half*(ciso(ns_ciso:n-1)*beta(ns_ciso:n-1)+ciso(ns_ciso+1:n)*beta(ns_ciso+1:n))
       !dcqvdca(ns_ciso-1:n-1) = half * beta(ns_ciso-1:n-1)
       !dcqvdcb(ns_ciso-1:n-1) = half * beta(ns_ciso:n)
       !cqv(n)                 = ciso(n)*beta(n)
       !dcqvdca(n)             = one*beta(n)
       !dcqvdcb(n)             = zero

       cqv(ns_ciso-1)         = ciso(ns_ciso)
       cqv(ns_ciso:n-1)       = half*(ciso(ns_ciso:n-1)+ciso(ns_ciso+1:n))
       dcqvdca(ns_ciso-1:n-1) = half
       dcqvdcb(ns_ciso-1:n-1) = half
       cqv(n)                 = ciso(n)
       dcqvdca(n)             = one
       dcqvdcb(n)             = zero

       betaqv(ns_ciso-1)      = beta(ns_ciso)*alphak_vdiff
       betaqv(ns_ciso:n-1)    = half*(beta(ns_ciso:n-1)+beta(ns_ciso+1:n))*alphak_vdiff
       betaqv(n)              = beta(n)*alphak_vdiff
       ! print*, 'betaqv ', betaqv
       ! print*, 'beta ', beta
       ! print*, 'alpha ', alphak_vdiff, half

       dbetaqv(ns_ciso:n-1)   = (beta(ns_ciso:n-1) - beta(ns_ciso+1:n))/deltaz(ns_ciso:n-1)
       dbetaqv(ns_ciso-1)     = zero !! vh ???? !!
       dbetaqv(n)             = zero

       !MC - Test
       ! if (experiment == 16) then
       !    print*, 'Ha01 ', betaqv
       !    ! print*, 'Ha02 ', dbetaqv
       ! betaqv  = one
       !   dbetaqv = zero
       ! endif

       ! modify for snow melt (including disappearance of snow pack) and/or transfer from top soil to new snow pack
       if ((qmelt.gt.zero) .or. ((qtransfer).gt.zero)) then
          cql(0)      = (qlsig(0)*cql(0)+qmelt*ciso(0)+cisoice(1)*(-qtransfer))/(qlsig(0)+qmelt-qtransfer)
          dcqldcb(0)  =  half*qlsig(0)/(qlsig(0)+qmelt-qtransfer)
          dcqldca(0)  =  half*qlsig(0)/(qlsig(0)+qmelt-qtransfer)
          qlsig(0)    = qlsig(0) + qmelt - qtransfer
       endif

    case default
       write(*,*) "isotope_vap: illegal formulation [1-2]: ", formulation
       stop 2
    end select

    ! mean diffusivities
    wl(ns_ciso:n)         = dx(ns_ciso:n)
    Dlmean(ns_ciso-1:n-1) = (Dl(ns_ciso-1:n-1)*wl(ns_ciso-1:n-1) + &
         Dl(ns_ciso:n)*wl(ns_ciso:n))/(wl(ns_ciso-1:n-1)+wl(ns_ciso:n))
    DLmean(n)             = zero

    wv(ns_ciso:n)       = dx(ns_ciso:n)
    Dvmean(ns_ciso:n-1) = (Dv(ns_ciso:n-1)*wv(ns_ciso:n-1) + Dv(ns_ciso+1:n)*wv(ns_ciso+1:n))/(wv(ns_ciso:n-1)+wv(ns_ciso+1:n))
    Dvmean(n)           = zero

    Dvbetamean(ns_ciso:n-1) = (Dv(ns_ciso:n-1)*beta(ns_ciso:n-1)*wv(ns_ciso:n-1) &
         + Dv(ns_ciso+1:n)*beta(ns_ciso+1:n)*wv(ns_ciso+1:n))/(wv(ns_ciso:n-1)+wv(ns_ciso+1:n))
    Dvbetamean(n)           = zero

    ! coefficients of tridiagonal matrix
    aa(ns_ciso)     = zero

    aa(ns_ciso+1:n) = qlsig(ns_ciso:n-1)*dcqldca(ns_ciso:n-1) &
         + qvsig(ns_ciso:n-1)*betaqv(ns_ciso:n-1)*dcqvdca(ns_ciso:n-1) &
         + dbetaqv(ns_ciso:n-1)*dcqvdca(ns_ciso:n-1)*Dvmean(ns_ciso:n-1) &
         + Dlmean(ns_ciso:n-1)/deltaz(ns_ciso:n-1) &
         + Dvbetamean(ns_ciso:n-1)/deltaz(ns_ciso:n-1)
    !MC - thetasat-thetar?
    bb(ns_ciso)       = -(Seff(ns_ciso)+kfreeze(ns_ciso)*Sicesig(ns_ciso)) * &
         thetasat(ns_ciso)*dx(ns_ciso)/sig/dt & !!!vh!!! NB thetasat should be (thetasat-thetar)??
         - qevapout*dcevapoutdciso &
         - qlsig(ns_ciso)*dcqldca(ns_ciso) - qvsig(ns_ciso)*betaqv(ns_ciso)*dcqvdca(ns_ciso) &
         - dbetaqv(ns_ciso)*dcqvdca(ns_ciso)*Dvmean(ns_ciso) &
         - Dlmean(ns_ciso)/deltaz(ns_ciso) &
         - Dvbetamean(ns_ciso)/deltaz(ns_ciso) &
         - (qex_ss(ns_ciso))
    bb(ns_ciso+1:n-1) = -(Seff(ns_ciso+1:n-1)+kfreeze(ns_ciso+1:n-1)*Sicesig(ns_ciso+1:n-1)) * &
         thetasat(ns_ciso+1:n-1)*dx(ns_ciso+1:n-1)/sig/dt &
         + qlsig(ns_ciso:n-2)*dcqldcb(ns_ciso:n-2) &
         - qlsig(ns_ciso+1:n-1)*dcqldca(ns_ciso+1:n-1) &
         + qvsig(ns_ciso:n-2)*betaqv(ns_ciso:n-2)*dcqvdcb(ns_ciso:n-2) &
         - qvsig(ns_ciso+1:n-1)*betaqv(ns_ciso+1:n-1)*dcqvdca(ns_ciso+1:n-1) &
         - (Dlmean(ns_ciso:n-2)/deltaz(ns_ciso:n-2) + Dlmean(ns_ciso+1:n-1)/deltaz(ns_ciso+1:n-1)) &
         - (Dvbetamean(ns_ciso:n-2)/deltaz(ns_ciso:n-2) + Dvbetamean(ns_ciso+1:n-1)/deltaz(ns_ciso+1:n-1)) &
         - qex_ss(ns_ciso+1:n-1) &
         + dbetaqv(ns_ciso:n-2)*dcqvdcb(ns_ciso:n-2)*Dvmean(ns_ciso:n-2) &
         - dbetaqv(ns_ciso+1:n-1)*dcqvdca(ns_ciso+1:n-1)*Dvmean(ns_ciso+1:n-1)

    bb(1) = bb(1) - qrunoff

    bb(n) = -(Seff(n)+kfreeze(n)*Sicesig(n)) * &
         thetasat(n)*dx(n)/sig/dt &
         + qlsig(n-1)*dcqldcb(n-1)  +qvsig(n-1)*betaqv(n-1)*dcqvdcb(n-1) &
         + dbetaqv(n-1)*dcqvdcb(n-1) *Dvmean(n-1) &
         - qlsig(n)*dcqldca(n) &
         - Dlmean(n-1)/deltaz(n-1) &
         - Dvbetamean(n-1)/deltaz(n-1) &
         - qex_ss(n)

    cc(ns_ciso:n-1) = -qlsig(ns_ciso:n-1)*dcqldcb(ns_ciso:n-1)  - qvsig(ns_ciso:n-1)*betaqv(ns_ciso:n-1)*dcqvdcb(ns_ciso:n-1) &
         - dbetaqv(ns_ciso:n-1)*dcqvdcb(ns_ciso:n-1)*Dvmean(ns_ciso:n-1) &
         + Dlmean(ns_ciso:n-1)/deltaz(ns_ciso:n-1) &
         + Dvbetamean(ns_ciso:n-1)/deltaz(ns_ciso:n-1)

    cc(n)           = zero

    dd(ns_ciso) = thetasat(ns_ciso)*dx(ns_ciso)/sig/dt* &
         (ciso(ns_ciso)*deltaSeff(ns_ciso) + cisoice(ns_ciso)*deltaSice(ns_ciso) + &
         kfreeze2(ns_ciso)*Sicesig(ns_ciso)) &
         - qprec_snow*cprec_snow/sig - qprec*cprec/sig - qevapin*cevapin/sig + qevapout*cevapout/sig &
         + qlsig(ns_ciso)*cql(ns_ciso)/sig + qvsig(ns_ciso)*betaqv(ns_ciso)*cqv(ns_ciso)/sig &
         + dbetaqv(ns_ciso)*cqv(ns_ciso)*Dvmean(ns_ciso)/sig &
         + Dlmean(ns_ciso)/deltaz(ns_ciso)/sig*(ciso(ns_ciso) - ciso(ns_ciso+1)) &
         + Dvbetamean(ns_ciso)/deltaz(ns_ciso)/sig*(ciso(ns_ciso)-ciso(ns_ciso+1)) &
         + (qex_ss(ns_ciso))*ciso(ns_ciso)/sig

    dd(ns_ciso+1:n-1) = thetasat(ns_ciso+1:n-1)*dx(ns_ciso+1:n-1)/sig/dt* &
         (ciso(ns_ciso+1:n-1)*deltaSeff(ns_ciso+1:n-1) + cisoice(ns_ciso+1:n-1)*deltaSice(ns_ciso+1:n-1) + &
         kfreeze2(ns_ciso+1:n-1)*Sicesig(ns_ciso+1:n-1)) &
         - qlsig(ns_ciso:n-2)*cql(ns_ciso:n-2)/sig &
         - qvsig(ns_ciso:n-2)*betaqv(ns_ciso:n-2)*cqv(ns_ciso:n-2)/sig &
         - dbetaqv(ns_ciso:n-2)*cqv(ns_ciso:n-2)*Dvmean(ns_ciso:n-2)/sig &
         + qlsig(ns_ciso+1:n-1)*cql(ns_ciso+1:n-1)/sig &
         + qvsig(ns_ciso+1:n-1)*betaqv(ns_ciso+1:n-1)*cqv(ns_ciso+1:n-1)/sig &
         + dbetaqv(ns_ciso+1:n-1)*cqv(ns_ciso+1:n-1)*Dvmean(ns_ciso+1:n-1)/sig &
         - Dlmean(ns_ciso:n-2)/deltaz(ns_ciso:n-2)/sig*(ciso(ns_ciso:n-2) - ciso(ns_ciso+1:n-1)) &
         - Dvbetamean(ns_ciso:n-2)/deltaz(ns_ciso:n-2)/sig*(ciso(ns_ciso:n-2) - &
         ciso(ns_ciso+1:n-1)) &
         + Dlmean(ns_ciso+1:n-1)/deltaz(ns_ciso+1:n-1)/sig*(ciso(ns_ciso+1:n-1) - ciso(ns_ciso+2:n)) &
         + Dvbetamean(ns_ciso+1:n-1)/deltaz(ns_ciso+1:n-1)/sig*(ciso(ns_ciso+1:n-1) - &
         ciso(ns_ciso+2:n)) &
         + qex_ss(ns_ciso+1:n-1)*ciso(ns_ciso+1:n-1)/sig

    dd(1) = dd(1) + qrunoff*ciso(1)/sig

    dd(n) = thetasat(n)*dx(n)/sig/dt* &
         (ciso(n)*deltaSeff(n) + cisoice(n)*deltaSice(n) + &
         kfreeze2(n)*Sicesig(n)) &
         - qlsig(n-1)*cql(n-1)/sig &
         - qvsig(n-1)*betaqv(n-1)*cqv(n-1)/sig &
         - dbetaqv(n-1)*cqv(n-1)*Dvmean(n-1)/sig &
         + qlsig(n)*cql(n)/sig &
         - Dlmean(n-1)/deltaz(n-1)/sig*(ciso(n-1) - ciso(n)) &
         - Dvbetamean(n-1)/deltaz(n-1)/sig*(ciso(n-1) - ciso(n)) &
         + qex_ss(n)*ciso(n)/sig

    if (cali>zero .or. experiment==7 .or. experiment==8) then
       bb(n) = bb(n) + qlsig(n)*dcqldca(n)
       dd(n) = dd(n) + qlsig(n)*(cali-cql(n))/sig
    endif

    ! print*, 'aa ', aa(ns_ciso)
    ! print*, 'bb ', bb(ns_ciso)
    ! print*, 'cc ', cc(ns_ciso)
    ! print*, 'dd ', dd(ns_ciso)
    call tri(ns_ciso,n,aa,bb,cc,dd,dc)
    ! print*, 'After ', bb(ns_ciso) * dc(ns_ciso) + cc(ns_ciso) * dc(ns_ciso+1), dd(ns_ciso), &
    !      bb(ns_ciso) * dc(ns_ciso) + cc(ns_ciso) * dc(ns_ciso+1) - dd(ns_ciso)
    !MC - Why not sig*dc ... ?
    ! dcice(1:n) = sig*dc(1:n)*kfreeze(1:n) + kfreeze2(1:n)
    dcice(1:n) = dc(1:n)*kfreeze(1:n) + kfreeze2(1:n)
    if (nsnow.gt.0) then
       where ((Sliq(ns_ciso:0).gt.zero).or.(deltaSliq(ns_ciso:0).gt.zero))  ! solve for change in isotopes in liquid
          ! dcice(ns_ciso:0) = sig*dc(ns_ciso:0)*kfreeze(ns_ciso:0) + kfreeze2(ns_ciso:0)
          dcice(ns_ciso:0) = dc(ns_ciso:0)*kfreeze(ns_ciso:0) + kfreeze2(ns_ciso:0)
       elsewhere  ! solve for change in isotopes in ice
          dcice(ns_ciso:0) = dc(ns_ciso:0)
       endwhere
    endif

    ! check for mass balance
    !  if (ns_ciso == 1) then
    LHS (ns_ciso:n) = thetasat(ns_ciso:n)*dx(ns_ciso:n)/dt*(dc(ns_ciso:n)*Seff(ns_ciso:n) + &
         ciso(ns_ciso:n)*deltaSeff(ns_ciso:n)) + &
         thetasat(ns_ciso:n)*dx(ns_ciso:n)/dt*(dcice(ns_ciso:n)*Sicesig(ns_ciso:n) + &
         cisoice(ns_ciso:n)*deltaSice(ns_ciso:n))

    RHS(ns_ciso+1:n-1) = qlsig(ns_ciso:n-2)*(cql(ns_ciso:n-2) + sig*dc(ns_ciso:n-2)*dcqldca(ns_ciso:n-2) + &
         sig*dc(ns_ciso+1:n-1)*dcqldcb(ns_ciso:n-2)) &
         + qvsig(ns_ciso:n-2)*betaqv(ns_ciso:n-2)*(cqv(ns_ciso:n-2) + sig*dc(ns_ciso:n-2)*dcqvdca(ns_ciso:n-2) &
         + sig*dc(ns_ciso+1:n-1)*dcqvdcb(ns_ciso:n-2)) &
         -qlsig(ns_ciso+1:n-1)*(cql(ns_ciso+1:n-1) + sig*dc(ns_ciso+1:n-1)*dcqldca(ns_ciso+1:n-1) + &
         sig*dc(ns_ciso+2:n)*dcqldcb(ns_ciso+1:n-1) ) &
         - qvsig(ns_ciso+1:n-1)*betaqv(ns_ciso+1:n-1)*(cqv(ns_ciso+1:n-1) + sig*dc(ns_ciso+1:n-1)*dcqvdca(ns_ciso+1:n-1) + &
         sig*dc(ns_ciso+2:n)*dcqvdcb(ns_ciso+1:n-1)) &
         + Dlmean(ns_ciso:n-2)*((ciso(ns_ciso:n-2) + sig*dc(ns_ciso:n-2))- &
         (ciso(ns_ciso+1:n-1)+sig*dc(ns_ciso+1:n-1)))/deltaz(ns_ciso:n-2) &
         - Dlmean(ns_ciso+1:n-1)*((ciso(ns_ciso+1:n-1) + sig*dc(ns_ciso+1:n-1))- &
         (ciso(ns_ciso+2:n)+sig*dc(ns_ciso+2:n)))/deltaz(ns_ciso+1:n-1) &
         + Dvbetamean(ns_ciso:n-2)*((ciso(ns_ciso:n-2) + sig*dc(ns_ciso:n-2)) - &
         (ciso(ns_ciso+1:n-1) &
         + sig*dc(ns_ciso+1:n-1))*beta(ns_ciso+1:n-1))/deltaz(ns_ciso:n-2) &
         - Dvbetamean(ns_ciso+1:n-1)*((ciso(ns_ciso+1:n-1) + sig*dc(ns_ciso+1:n-1)) - &
         (ciso(ns_ciso+2:n) + sig*dc(ns_ciso+2:n)))/deltaz(ns_ciso+1:n-1) &
         - qex_ss(ns_ciso+1:n-1)*(ciso(ns_ciso+1:n-1) + sig*dc(ns_ciso+1:n-1)) &
         + dbetaqv(ns_ciso:n-2)*dcqvdcb(ns_ciso:n-2)*Dvmean(ns_ciso:n-2)*(ciso(ns_ciso:n-2) + sig*dc(ns_ciso:n-2)) &
         - dbetaqv(ns_ciso+1:n-1)*dcqvdcb(ns_ciso+1:n-1)*Dvmean(ns_ciso+1:n-1)*(ciso(ns_ciso+1:n-1) + sig*dc(ns_ciso+1:n-1))

    RHS(ns_ciso) = qprec_snow*cprec_snow + qprec*cprec - qevapout*(cevapout + sig*dc(ns_ciso)*dcevapoutdciso) &
         +qevapin*cevapin &
         -qlsig(ns_ciso)*(cql(ns_ciso) + sig*dc(ns_ciso)*dcqldca(ns_ciso) + sig*dc(ns_ciso+1)*dcqldcb(ns_ciso)) &
         - qvsig(ns_ciso)*betaqv(ns_ciso)*(cqv(ns_ciso) + sig*dc(ns_ciso)*dcqvdca(ns_ciso) + &
         sig*dc(ns_ciso+1)*dcqvdcb(ns_ciso)) &
         - Dlmean(ns_ciso)*((ciso(ns_ciso) + sig*dc(ns_ciso))-(ciso(ns_ciso+1)+sig*dc(ns_ciso+1)))/deltaz(ns_ciso) &
         - Dvbetamean(ns_ciso)*((ciso(ns_ciso) + sig*dc(ns_ciso)) - &
         (ciso(ns_ciso+1) + sig*dc(ns_ciso+1)))/deltaz(ns_ciso) &
         - (qex_ss(ns_ciso))*(ciso(ns_ciso) + sig*dc(ns_ciso)) &
         - dbetaqv(ns_ciso)*dcqvdcb(ns_ciso)*Dvmean(ns_ciso)*(ciso(ns_ciso) + sig*dc(ns_ciso))
    ! print*, 'Ha01 ', thetasat(ns_ciso)*dx(ns_ciso)/dt*(dc(ns_ciso)*Seff(ns_ciso) + &
    !      ciso(ns_ciso)*deltaSeff(ns_ciso)), &
    !      thetasat(ns_ciso)*dx(ns_ciso)/dt*(dcice(ns_ciso)*Sicesig(ns_ciso) + &
    !      (cisoice(ns_ciso:n)+sig*dcice(ns_ciso:n))*deltaSice(ns_ciso:n))
    ! print*, 'Ha02 ', qprec_snow*cprec_snow, &
    !      qprec*cprec, &
    !      qevapout*(cevapout + sig*dc(ns_ciso)*dcevapoutdciso), &
    !      qevapin*cevapin, &
    !      qlsig(ns_ciso)*(cql(ns_ciso) + sig*dc(ns_ciso)*dcqldca(ns_ciso) + sig*dc(ns_ciso+1)*dcqldcb(ns_ciso)), &
    !      qvsig(ns_ciso)*betaqv(ns_ciso)*(cqv(ns_ciso) + sig*dc(ns_ciso)*dcqvdca(ns_ciso) + &
    !      sig*dc(ns_ciso+1)*dcqvdcb(ns_ciso)), &
    !      Dlmean(ns_ciso)*((ciso(ns_ciso) + sig*dc(ns_ciso))-(ciso(ns_ciso+1)+sig*dc(ns_ciso+1)))/deltaz(ns_ciso), &
    !      Dvbetamean(ns_ciso)*((ciso(ns_ciso) + sig*dc(ns_ciso)) - &
    !      (ciso(ns_ciso+1) + sig*dc(ns_ciso+1)))/deltaz(ns_ciso), &
    !      (qex_ss(ns_ciso))*(ciso(ns_ciso) + sig*dc(ns_ciso)), &
    !      dbetaqv(ns_ciso)*dcqvdcb(ns_ciso)*Dvmean(ns_ciso)*(ciso(ns_ciso) + sig*dc(ns_ciso))
    ! print*, 'Ha03 ', qvsig(ns_ciso), betaqv(ns_ciso), cqv(ns_ciso), sig, dc(ns_ciso), dcqvdca(ns_ciso), &
    !      dc(ns_ciso+1), dcqvdcb(ns_ciso)
    ! print*, 'Ha04 ', cqv(ns_ciso), sig*dc(ns_ciso)*dcqvdca(ns_ciso), sig*dc(ns_ciso+1)*dcqvdcb(ns_ciso)


    RHS(1) = RHS(1) -qrunoff*(ciso(1) + sig*dc(1))

    RHS(n) =  qlsig(n-1)*(cql(n-1) + sig*dc(n-1)*dcqldca(n-1) + sig*dc(n)*dcqldcb(n-1)) &
         + qvsig(n-1)*betaqv(n-1)*(cqv(n-1) + sig*dc(n-1)*dcqvdca(n-1) + sig*dc(n)*dcqvdcb(n-1)) &
         -qsig(n)*(ciso(n) + sig*dc(n)) &
         + Dlmean(n-1)*((ciso(n-1) + sig*dc(n-1))-(ciso(n)+sig*dc(n)))/deltaz(n-1) &
         + Dvbetamean(n-1)*((ciso(n-1) + sig*dc(n-1)) - (ciso(n) + sig*dc(n)))/deltaz(n-1) &
         -qex_ss(n)*(ciso(n)+sig*dc(n)) &
         + dbetaqv(n)*dcqvdcb(n)*Dvmean(n)*(ciso(n) + sig*dc(n))

    if (cali>zero .or. experiment==7 .or. experiment==8) then
       RHS(n) = RHS(n) - qsig(n) * (cali - (ciso(n)+sig*dc(n)))
    endif

    ! if (any(abs(LHS(ns_ciso:n)-RHS(ns_ciso:n))>tiny(one))) then
    ! if (any(abs(LHS(ns_ciso:n)-RHS(ns_ciso:n))>epsilon(one))) then
    if (any(abs(LHS(ns_ciso:n)-RHS(ns_ciso:n))>1.e-6_r_2)) then
       write(*,*) 'abs(LHS-RHS) ', abs(LHS(ns_ciso:n)-RHS(ns_ciso:n))
       ii = maxloc(abs(LHS(ns_ciso:n)-RHS(ns_ciso:n)))
       write(*,*) 'Max of abs(LHS-RHS): ', maxval(abs(LHS(ns_ciso:n)-RHS(ns_ciso:n))), ii
       write(*,*) 'nsnow ', nsnow, irec, ns_ciso
       write(*,*) 'surface ', ciso(ns_ciso), dc(ns_ciso), dcice(ns_ciso)
       write(*,*) 'flux ', cprec_snow, cevapin, cevapout
       write(*,*) 'surface flux ', qv0, ql0, -qevap
       write(*,*) 'layer ', ciso(ii(1)), dc(ii(1)), dcice(ii(1))
       write(*,*) 'cqv', cqv
       write(*,*) 'cql', cql
       write(*,*) 'qprec', qprec, qprec_snow, kfreeze(ns_ciso), kfreeze2(ns_ciso)
       write(*,*) 'qmelt', qmelt, qtransfer, qrunoff, Seff(1)
       write(*,*) 'surface balance ', LHS(ns_ciso), RHS(ns_ciso), qprec_snow*cprec_snow
       write(*,*) 'layer balance ', LHS(ii(1)), RHS(ii(1))
       write(*,*) 'surface sigs ', nsnow_last, deltaSice(ns_ciso), deltaSliq(ns_ciso), cvsig(ns_ciso), deltacv(ns_ciso)
       stop 2
    endif

    !   endif ! 1==1

    ! isotopic fluxes
    qiso_in    = qprec_snow*cprec_snow + qprec*cprec +qevapin*cevapin
    qiso_out   = qevapout*(cevapout + sig*dc(ns_ciso)*dcevapoutdciso)
    qiso_evap  = qevapout*(cevapout + sig*dc(ns_ciso)*dcevapoutdciso) - qevapin*cevapin
    qiso_trans = sum(qex_ss(ns_ciso:n)*(ciso(ns_ciso:n)+sig*dc(ns_ciso:n)),1)

    qiso_liq_diff(ns_ciso:n-1) = Dlmean(ns_ciso:n-1)*((ciso(ns_ciso:n-1) + &
         sig*dc(ns_ciso:n-1))-(ciso(ns_ciso+1:n)+sig*dc(ns_ciso+1:n)))/deltaz(ns_ciso:n-1)
    qiso_vap_diff(ns_ciso:n-1) = Dvbetamean(ns_ciso:n-1)*((ciso(ns_ciso:n-1) + sig*dc(ns_ciso:n-1)) &
         - (ciso(ns_ciso+1:n) + sig*dc(ns_ciso+1:n)))/deltaz(ns_ciso :n-1)
    qiso_liq_adv(ns_ciso:n-1)  = qlsig(ns_ciso:n-1)*(cql(ns_ciso:n-1) + sig*dc(ns_ciso:n-1)*dcqldca(ns_ciso:n-1) + &
         sig*dc(ns_ciso+1:n)*dcqldcb(ns_ciso:n-1))
    qiso_liq_adv(n)      = qlsig(n)*(cql(n) + sig*dc(n))
    qiso_vap_adv(ns_ciso:n-1)  = qvsig(ns_ciso:n-1) * betaqv(ns_ciso:n-1) &
         * (cqv(ns_ciso:n-1) + sig*dc(ns_ciso:n-1)*dcqvdca(ns_ciso:n-1) + sig*dc(ns_ciso+1:n)*dcqvdcb(ns_ciso:n-1))
    qiso_vap_adv(n)      = zero

    ciso         = ciso + dc
    cisoice(1:n) = cisoice(1:n) + dcice(1:n)
    if (nsnow > 0) then
       cisoice(ns_ciso:0) = merge(cisoice(ns_ciso:0) + dcice(ns_ciso:0), & ! solve for change in isotopes in liquid
            cisoice(ns_ciso:n) + dc(ns_ciso:0), & ! solve for change in isotopes in ice
            (Sliq(ns_ciso:0) > zero) .or. (deltaSliq(ns_ciso:0) > zero))
    endif

    write(993,*) ciso(ns_ciso), cisos, civa/cva/alphaplus_a

  END SUBROUTINE isotope_vap

  !*********************************************************************************************************************

END MODULE sli_solve
