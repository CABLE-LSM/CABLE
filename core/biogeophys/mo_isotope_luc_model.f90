!> \file mo_isotope_luc_model.f90

!> \brief Generic isotopic land-use change model

!> \details Explicit solution of a generic land-use change model,
!> given all the fluxes, area changes between land use classes, source and sinks.
!> No fractionations included (yet).

!> \author Matthias Cuntz
!> \date Apr 2019

MODULE mo_isotope_luc_model

  implicit none

  private

  ! public routines
  public :: isotope_luc_model  ! Solves the next time step of a generic isotopic land-use change model

  ! ------------------------------------------------------------------

  !     NAME
  !         isotope_luc_model

  !     PURPOSE
  !>        \brief Generic isotopic land-use change model

  !>        \details Next explicit time step of the isotopic composition of a generic land-use change model,
  !>        given all fluxes, area changes between land use classes, sources and sinks.
  !>        The generic isotopic land-use model does not include any fractionations, yet.

  !>        The land-use change model for land-use class i is: \n
  !>            d(A(i)*C(i)) = -C(i)*sum(dA(i,:)) + sum(dA(:,i)*C(:)) + S(i) - T(i) \n
  !>        where \n
  !>        C [kg/m^2] are the carbon concentrations in the land-use classes, \n
  !>        A [m^2] are the areas of the land-use classes, \n
  !>        dA [m^2] are the area changes between land-use classes for the time step (A>0, A(i,i)=0), \n
  !>        S [kg] are sources other than contributions from other land use classes, and \n
  !>        T [kg] are sinks other than contributions to other land use classes. \n
  !>        S and T could be source and sinks between carbon pools, for example.

  !>        Note: \n
  !>             dA is change per time step in [m^2], NOT per second.\n
  !>             S and T are sources and sinks in total mass in [kg], NOT per area.

  !>        The land-use model is then solved from time step t to t+dt \n
  !>            A(i,t+dt)*C(i,t+dt) - A(i,t)*C(i,t) = -C(i,t)*sum(dA(i,:)) + sum(dA(:,i)*C(:,t)) + S(i,dt) - T(i,dt) \n
  !>            C(i,t+dt) = (C(i,t)*(A(i,t)-sum(dA(i,:))) + sum(dA(:,i)*C(:,t)) + S(i,dt) - T(i,dt)) / A(i,t+dt) \n
  !>        with S(i,dt), T(i,dt) the total sources and sinks during the time step.
  
  !>        The isotopic land-use model is then:
  !>            A(i,t+dt)*C'(i,t+dt) - A(i,t)*C'(i,t) =
  !>                -C'(i,t)*sum(dA(i,:)) + sum(dA(:,i)*C'(:,t)) + Rs(i,dt)*S(i,dt) - R(i,t)*T(i,dt) \n
  !>            C'(i,t+dt) = (C'(i,t)*(A(i,t)-sum(dA(i,:))) + sum(dA(:,i)*C'(:,t)) + Rs(i,dt)*S(i,dt) - R(i,t)*T(i,dt) ) /
  !>                        A(i,t+dt) \n
  !>        with C' the isotope carbon concentrations, R the carbon isotope ratio C'/C,
  !>        and Rs is the carbon isotope ratio of the source S during the time step.

  !>        This can also be written with carbon concentrations C at time t as \n
  !>            R(i,t) = C'(i,t) / C(i,t) \n
  !>            A(i,t+dt)*C'(i,t+dt) - A(i,t)*R(i,t)*C(i,t) = -R(i,t)*C(i,t)*sum(dA(i,:)) + sum(dA(:,i)*R(:,t)*C(:,t)) +
  !>                    Rs(i,dt)*S(i,dt) - R(i,t)*T(i,dt) \n
  !>            C'(i,t+dt) = (R(i,t)*C(i,t)*(A(i,t)-sum(dA(i,:))) + sum(dA(:,i)*R(:,t)*C(:,t)) +
  !>                          Rs(i,dt)*S(i,dt) - R(i,t)*T(i,dt)) / A(i,t+dt)

  !     CALLING SEQUENCE
  !         call isotope_luc_model(Ci, A, dA, C=C, S=S, Rs=Rs, T=T, Rt=Rt, At=At, trash=trash)

  !     INTENT
  !>        \param[inout] "real(dp) :: Ci(1:n)"
  !>                                   On entry: isotope concentrations in land-use classes at time step t-1 [kg/m^2]
  !>                                   On exit:  updated isotope concentrations of land-use classes at time step t [kg/m^2]
  !>        \param[in]    "real(dp) :: A(1:n)"                   Areas of land-use classes at time step t-1 [m^2]
  !>        \param[in]    "real(dp) :: dA(1:n,1:n)"              Change from each land-use class to the others
  !>                                   during the time step [m^2]
  !>        \param[in]    "real(dp), optional :: C(1:n)"         Non-isotope concentrations in land-use classes
  !>                                   at time step t-1 [kg/m^2] (default: not used)
  !>        \param[in]    "real(dp), optional :: S(1:n)"         Sources other than contributions from other land use classes,
  !>                                   for example from other carbon pool [kg/m^2] (default: 0.)
  !>        \param[in]    "real(dp), optional :: Rs(1:n)"        Isotope ratio of sources (default: 1.)
  !>        \param[in]    "real(dp), optional :: T(1:n)"         Sinks other than contributions from other land use classes,
  !>                                   for example to other carbon pool [kg/m^2] (default: 0.)
  !>        \param[in]    "real(dp), optional :: Rt(1:n)"        Isotope ratio of sinks
  !>                                   (default: isotope ratio of pool and land use class)
  !>        \param[in]    "real(dp), optional :: At(1:n)"         Areas of land-use classes at time step t [m^2]
  !>                                   (default: calculated from A and dA)
  !>        \param[inout] "real(dp), optional :: trash(1:n)"     Container to store possible inconsistencies,
  !>                                   might be numeric, between non-isotope and isotope model [kg/m^2].

  !     HISTORY
  !>        \author Written Matthias Cuntz
  !>        \date Apr 2019
  interface isotope_luc_model
     module procedure isotope_luc_model_1d
  end interface isotope_luc_model

  ! ------------------------------------------------------------------

contains

  ! ------------------------------------------------------------------

  subroutine isotope_luc_model_1d(Ci, A, dA, C, S, Rs, T, Rt, At, trash)

    use mo_kind,  only: dp, i4
    use mo_utils, only: eq !, ne

    implicit none

    real(dp), dimension(:),   intent(inout)           :: Ci     ! Iso land-use class
    real(dp), dimension(:),   intent(in)              :: A      ! Area land-use class
    real(dp), dimension(:,:), intent(in)              :: dA     ! Area changes between land-use classes
    real(dp), dimension(:),   intent(in),    optional :: C      ! Non-iso land-use class
    real(dp), dimension(:),   intent(in),    optional :: S      ! Source
    real(dp), dimension(:),   intent(in),    optional :: Rs     ! Isotope ration of source
    real(dp), dimension(:),   intent(in),    optional :: T      ! Sink
    real(dp), dimension(:),   intent(in),    optional :: Rt     ! Isotope ration of sink
    real(dp), dimension(:),   intent(in),    optional :: At     ! New area land-use class
    real(dp), dimension(:),   intent(inout), optional :: trash  ! garbage can for numerical inconsistencies

    ! Local variables
    integer(i4) :: i  ! counter
    integer(i4) :: nn ! number of pools
    real(dp), dimension(size(Ci,1)) :: R    ! Isotope ratio of pool
    ! defaults for optional inputs
    real(dp), dimension(size(Ci,1)) :: iC, iS, iRs, iT, iRt, itrash
    real(dp), dimension(size(Ci,1)) :: Anew ! A at t+dt

    ! Check sizes
    nn = size(Ci,1)
    if ( (size(A,1) /= nn) .or. (size(dA,1) /= nn) .or. (size(dA,2) /= nn) ) then
       write(*,*) 'Error isotope_luc_model_1d: non-fitting dimensions between isotopic concentrations,'
       write(*,*) '                            land areas and land-area changes.'
       write(*,*) '    size(Ci):   ', size(Ci,1)
       write(*,*) '    size(A):    ', size(A,1)
       write(*,*) '    size(dA,1): ', size(dA,1)
       write(*,*) '    size(dA,2): ', size(dA,2)
       stop 9
    endif

    ! Check dA >= 0
    if (any(dA < 0._dp)) then
       write(*,*) 'Error isotope_luc_model_1d: land area changes between land use classes must be >= 0.'
       write(*,*) '    dA: ', dA
       stop 9
    endif

    ! ! Check dA(i,:) == 0. if Ci(i) == 0.
    ! if (any(eq(Ci,0._dp))) then
    !    do i=1, nn
    !       if (eq(Ci(i),0._dp)) then
    !          if (any(ne(dA(i,:),0._dp))) then
    !             write(*,*) 'Error isotope_luc_model_1d: land area changes from land-use class i must be 0'
    !             write(*,*) '                            if isotope carbon concentration of land-use class is 0.'
    !             write(*,*) '    i, Ci(i): ', i, Ci(i)
    !             write(*,*) '       dA(i): ', dA(i,:)
    !             stop 9
    !          endif
    !       endif
    !    end do
    ! endif

    ! if (present(C)) then
    !    ! Check dA(i,:) == 0. if C(i) == 0.
    !    if (any(eq(C,0._dp))) then
    !       do i=1, nn
    !          if (eq(C(i),0._dp)) then
    !             if (any(ne(dA(i,:),0._dp))) then
    !                write(*,*) 'Error isotope_luc_model_1d: land area changes from land-use class i must be 0'
    !                write(*,*) '                            if carbon concentration of land-use class is 0.'
    !                write(*,*) '     i, C(i): ', i, C(i)
    !                write(*,*) '       dA(i): ', dA(i,:)
    !                stop 9
    !             endif
    !          endif
    !       end do
    !    endif
    ! endif

    !  Could be more checks for S, Rs, T, Rt, At
    
    ! Set optionals
    if (present(C)) then
       iC = C
    else
       iC = 1._dp
    endif
    if (present(S)) then
       iS = S
    else
       iS = 0._dp
    endif
    if (present(Rs)) then
       iRs = Rs
    else
       iRs = 1._dp
    endif
    if (present(T)) then
       iT = T
    else
       iT = 0._dp
    endif
    if (present(Rt)) then
       iRt = Rt
    else
       iRt = 1._dp
       where (iC > 0._dp) iRt = Ci / iC
    endif
    ! Land areas at t+dt
    if (present(At)) then
       Anew = At
    else
       Anew = A - sum(dA, dim=2) + sum(dA, dim=1)
    endif
    if (present(trash)) then
       itrash = trash
    else
       itrash = 0._dp
    endif

    ! Isotope ratio
    R(:) = 1._dp
    where (iC > 0._dp) R = Ci / iC

    ! isotopic LUC model
    Ci = R * iC * (A - sum(dA, dim=2)) + sum(dA * spread(R*iC, dim=2, ncopies=nn), dim=1) + iRs*iS - iRt*iT
    where (Anew > 0._dp)
       Ci = Ci / Anew
    elsewhere
       Ci     = 0._dp
       itrash = itrash + Ci
    endwhere

    ! Check final land-use classes
    ! Isotope land-use class became < 0.
    if (any(Ci < 0._dp)) then
       itrash = itrash + merge(abs(Ci), 0._dp, Ci < 0._dp)
       Ci = merge(0._dp, Ci, Ci < 0._dp)
    endif
    ! Non-isotope land-use class == 0. but isotope land-use class > 0.
    if (any(eq(iC,0._dp) .and. (Ci > 0._dp))) then
       itrash = itrash + merge(Ci, 0._dp, eq(iC,0._dp) .and. (Ci > 0._dp))
       Ci = merge(0._dp, Ci, eq(iC,0._dp) .and. (Ci > 0._dp))
    endif
    ! Non-isotope land-use class >0. but isotope land-use class == 0.
    ! ???

    if (present(trash)) trash = itrash

    return

  end subroutine isotope_luc_model_1d
  
  ! ------------------------------------------------------------------

END MODULE mo_isotope_luc_model
