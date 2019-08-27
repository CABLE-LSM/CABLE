!> \file mo_isotope_pool_model.f90

!> \brief Generic isotopic pool model

!> \details Explicit solution of a generic pool model,
!> given sink, sources and all the fluxes between all pools,
!> together with associated fractionation factors.

!> \author Matthias Cuntz
!> \date Apr 2019

MODULE mo_isotope_pool_model

  implicit none

  private

  public :: isotope_pool_model ! Solves the next time step of a generic isotopic pool model
  
  ! ------------------------------------------------------------------

  !     NAME
  !         isotope_pool_model

  !     PURPOSE
  !>        \brief Generic isotopic pool model

  !>        \details Next explicit time step of the isotopic composition of a generic pool model,
  !>        given sink, sources and all the fluxes between all pools,
  !>        together associated fractionation factors.

  !>        The pool model for pool i is: \n
  !>            dC(i)/dt = sum(F(:,i)) - sum(F(i,:)) + S(i) - T(i) \n
  !>        where C are the pools, F the fluxes between pools (F>0, F(i,i)=0),
  !>        S a pool-independent source term, and T a sink term, which can be pool-independent
  !>        or pool-dependent, i.e. T(i) = beta(i)*C(i)

  !>        The isotope model for pool i is then:\n
  !>            dC'(i)/dt = sum(alpha(:,i)*R(:)*F(:,i)) - R(i)*sum(alpha(i,:)*F(i,:))
  !>                       + Rs(i)*S(i) - alpha(i,i)*Rt(i)*T(i) \n
  !>        with C' the isotope pools. R are the isotope ratios of the pools at time t-dt,
  !>        and alpha are possible fractionation factors.
  !>        If Rt is not given, it is taken as the isotopic composition of the pool R(i).

  !>        All pools can have either dimension (n) or (n,m) with n pools and m land points (or fractions).
  !>        Fluxes and fractionation factors have dimensions (n,n) or (n,n,m).

  !     CALLING SEQUENCE
  !         call isotope_pool_model(dt, Ci, C, F, S=S, Rs=Rs, T=T, Rt=Rt, alpha=alpha, beta=beta, trash=trash)

  !     INTENT
  !>        \param[in]    "real(dp) :: dt"                       Time step
  !>        \param[inout] "real(dp) :: Ci(1:n[,1:m])"            Isotope concentrations in pools
  !>        \param[in]    "real(dp) :: C(1:n[,1:m])"             Non-isotope concentrations in pools at time step t-1
  !>        \param[in]    "real(dp) :: F(1:n,1:n[,1:m])"         Non-isotope fluxes between pools (F(i,i)=0) \n
  !>                                   F(i,j) positive flux from pool i to pool j \n
  !>                                   Isotope flux is alpha(i,j)*R(i)*F(i,j)
  !>        \param[in]    "real(dp), optional :: S(1:n[,1:m])"   Non-isotope source fluxes to pools (other than between pools) \n
  !>                                   Default: 0.
  !>        \param[in]    "real(dp), optional :: Rs(1:n[,1:m])"  Isotopic compositions of source fluxes \n
  !>                                   Any fractionations during source processes should be included in Rs. \n
  !>                                   Default: 1.
  !>        \param[in]    "real(dp), optional :: T(1:n[,1:m])"   Non-isotope sinks of pools (other than between pools) \n
  !>                                   Isotope flux is alpha(i,i)*Rt(i)*T(i) or alpha(i,i)*R(i)*T(i) \n
  !>                                   Default: 0.
  !>        \param[in]    "real(dp), optional :: Rt(1:n[,1:m])"  Isotopic compositions of sink fluxes. \n
  !>                                   If not given, the isotopic composition of the pool R(i) is taken. \n
  !>                                   Default: R(i)
  !>        \param[in]    "real(dp), optional :: alpha(1:n,1:n[,1:m])" Isotopic fractionation factors associated with
  !>                                   fluxes F between pools (alpha(i,j) i/=j) and of 
  !>                                   sinks T (other than between pools) (alpha(i,i)) \n
  !>                                   Default: 1.
  !>        \param[in]    "real(dp), optional :: beta(1:n[,1:m])"      Either T or beta can be given. \n
  !>                                   If beta is given then T(i) = beta(i)*C(i). \n
  !>                                   T supercedes beta, i.e. T will be taken if beta and T are given.
  !>        \param[inout] "real(dp), optional :: trash(1:n[,1:m])"     Container to store possible inconsistencies,
  !>                                   might be numeric, between non-isotope and isotope model. \n
  !>                                   Consistency checks are only performed if trash is given in call.
  !>        \param[inout] "logical,  optional :: trans"          Assumed order of pools and fluxes in 2D case \n
  !>                                   If .true.,  input order is [m,n] for pools and [m,n,n] for fluxes \n
  !>                                   If .false., input order is [n,m] for pools and [n,n,m] for fluxes \n
  !>                                   Default: .false.

  !     HISTORY
  !>        \author Written Matthias Cuntz
  !>        \date Apr 2019
  interface isotope_pool_model
     module procedure isotope_pool_model_1d, isotope_pool_model_2d
  end interface isotope_pool_model

  ! ------------------------------------------------------------------

  ! private routines
  private :: diag
  interface diag
     module procedure diag_2d, diag_3d
  end interface diag  

  ! ------------------------------------------------------------------

contains
  
  ! ------------------------------------------------------------------
  
  subroutine isotope_pool_model_1d(dt, Ci, C, F, S, Rs, T, Rt, alpha, beta, trash)

    use mo_kind,  only: dp, i4
    use mo_utils, only: eq, ne

    implicit none

    real(dp),                 intent(in)              :: dt     ! time step
    real(dp), dimension(:),   intent(inout)           :: Ci     ! Iso pool
    real(dp), dimension(:),   intent(in)              :: C      ! Non-iso pool
    real(dp), dimension(:,:), intent(in)              :: F      ! Fluxes between pools
    real(dp), dimension(:),   intent(in),    optional :: S      ! Sources not between pools
    real(dp), dimension(:),   intent(in),    optional :: Rs     ! Isotope ratio of S
    real(dp), dimension(:),   intent(in),    optional :: T      ! Sinks not between pools
    real(dp), dimension(:),   intent(in),    optional :: Rt     ! Isotope ratio of T
    real(dp), dimension(:,:), intent(in),    optional :: alpha  ! Fractionation factors sinks and fluxes between pools
    real(dp), dimension(:),   intent(in),    optional :: beta   ! Alternative to sinks = beta*Ct
    real(dp), dimension(:),   intent(inout), optional :: trash  ! garbage can for numerical inconsistencies

    ! Local variables
    integer(i4) :: i  ! counter
    integer(i4) :: nn ! number of pools
    real(dp), dimension(size(Ci,1)) :: R                 ! Isotope ratio of pool
    ! defaults for optional inputs
    real(dp), dimension(size(Ci,1))            :: iS, iRs, iT, iRt
    real(dp), dimension(size(Ci,1),size(Ci,1)) :: ialpha
    real(dp), dimension(size(Ci,1),size(Ci,1)) :: alphaF ! alpha*F
    real(dp), dimension(size(Ci,1)) :: sink              ! not-between pools sink
    real(dp), dimension(size(Ci,1)) :: source            ! not-between pools source
    real(dp), dimension(size(Ci,1)) :: isink             ! between pools sink
    real(dp), dimension(size(Ci,1)) :: isource           ! between pools source
    real(dp), dimension(size(Ci,1)) :: Cnew              ! New pool size for check

    ! Check sizes
    nn = size(Ci,1)
    if ( (size(C,1) /= nn) .or. (size(F,1) /= nn) .or. (size(F,2) /= nn) ) then
       write(*,*) 'Error isotope_pool_model_1d: non-fitting dimensions between isotopic pools,'
       write(*,*) '                             non-isotopic pools and fluxes.'
       write(*,*) '    size(Ci):  ', size(Ci,1)
       write(*,*) '    size(C):   ', size(C,1)
       write(*,*) '    size(F,1): ', size(F,1)
       write(*,*) '    size(F,2): ', size(F,2)
       stop 9
    endif

    ! Check F >= 0
    if (any(F < 0._dp)) then
       write(*,*) 'Error isotope_pool_model_1d: fluxes between pools must be >= 0.'
       write(*,*) '    F: ', F
       stop 9
    endif

    ! Check F(i,:) == 0. if C(i) == 0.
    if (any(eq(C,0._dp))) then
       do i=1, nn
          if (eq(C(i),0._dp)) then
             if (any(ne(F(i,:),0._dp))) then
                write(*,*) 'Error isotope_pool_model_1d: fluxes from pool i must be 0 if concentration in pool is 0.'
                write(*,*) '    i, C(i): ', i, C(i)
                write(*,*) '       F(i): ', F(i,:)
                stop 9
             endif
          endif
       end do
    endif

    ! Set optionals
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
       iRt = 1._dp ! could be zero as well to assure no isotopic sink if C=0
       where (C > 0._dp) iRt = Ci / C
    endif
    if (present(alpha)) then
       ialpha = alpha
    else
       ialpha = 1._dp
    endif
    if (present(alpha)) then
       ialpha = alpha
    else
       ialpha = 1._dp
    endif
    if (present(beta) .and. (.not. present(T))) then
       iT = beta * C
    endif

    ! Isotope ratio
    ! R(:) = 0._dp
    R(:) = 1._dp ! could be zero as well to assure no isotopic fluxes if initial C=0
    where (C > 0._dp) R = Ci / C

    ! alpha * F
    alphaF = ialpha * F

    ! between and not-between source and sinks
    sink    = diag(ialpha) * iRt * iT * dt
    source  = iRs * iS * dt
    isink   = sum(alphaF, dim=2) * R * dt
    isource = sum(alphaF * spread(R, dim=2, ncopies=nn), dim=1) * dt

    ! Explicit solution
    Ci = Ci - sink + source - isink + isource

    if (present(trash)) then
       ! Check final pools
       ! Isotope pool became < 0.
       if (any(Ci < 0._dp)) then
          trash = trash + merge(abs(Ci), 0._dp, Ci < 0._dp)
          Ci = merge(0._dp, Ci, Ci < 0._dp)
       endif
       ! Non-isotope pool == 0. but isotope pool > 0.
       Cnew = C - iT * dt + iS * dt - sum(F, dim=2)* dt + sum(F, dim=1) * dt
       if (any(eq(Cnew,0._dp) .and. (Ci > 0._dp))) then
          trash = trash + merge(Ci, 0._dp, eq(Cnew,0._dp) .and. (Ci > 0._dp))
          Ci = merge(0._dp, Ci, eq(Cnew,0._dp) .and. (Ci > 0._dp))
       endif
       ! Non-isotope pool >0. but isotope pool == 0.
       ! ???
    endif

    return

  end subroutine isotope_pool_model_1d

  
  subroutine isotope_pool_model_2d(dt, Ci, C, F, S, Rs, T, Rt, alpha, beta, trash, trans)

    use mo_kind,  only: dp, i4
    use mo_utils, only: eq, ne

    implicit none

    real(dp),                   intent(in)              :: dt     ! time step
    real(dp), dimension(:,:),   intent(inout)           :: Ci     ! Iso pool
    real(dp), dimension(:,:),   intent(in)              :: C      ! Non-iso pool
    real(dp), dimension(:,:,:), intent(in)              :: F      ! Fluxes between pools
    real(dp), dimension(:,:),   intent(in),    optional :: S      ! Sources not between pools
    real(dp), dimension(:,:),   intent(in),    optional :: Rs     ! Isotope ratio of S
    real(dp), dimension(:,:),   intent(in),    optional :: T      ! Sinks not between pools
    real(dp), dimension(:,:),   intent(in),    optional :: Rt     ! Isotope ratio of T
    real(dp), dimension(:,:,:), intent(in),    optional :: alpha  ! Fractionation factors sinks and fluxes between pools
    real(dp), dimension(:,:),   intent(in),    optional :: beta   ! Alternative to sinks = beta*Ct
    real(dp), dimension(:,:),   intent(inout), optional :: trash  ! garbage can for numerical inconsistencies
    logical,                    intent(in),    optional :: trans  ! transposed order pools and fluxes

   ! Local variables
    integer(i4) :: i, j  ! counter
    integer(i4) :: nland ! number of land points
    integer(i4) :: nn    ! number of pools
    real(dp), dimension(size(Ci,1),size(Ci,2)) :: R       ! Isotope ratio of pool
    ! defaults for optional inputs
    real(dp), dimension(size(Ci,1),size(Ci,2)) :: iS, iRs, iT, iRt
    real(dp), dimension(:,:,:), allocatable :: ialpha
    real(dp), dimension(:,:,:), allocatable :: alphaF     ! alpha*F
    real(dp), dimension(size(Ci,1),size(Ci,2)) :: sink    ! not-between pools sink
    real(dp), dimension(size(Ci,1),size(Ci,2)) :: source  ! not-between pools source
    real(dp), dimension(size(Ci,1),size(Ci,2)) :: isink   ! between pools sink
    real(dp), dimension(size(Ci,1),size(Ci,2)) :: isource ! between pools source
    real(dp), dimension(size(Ci,1),size(Ci,2)) :: Cnew    ! New pool size for check
    logical :: itrans                                     ! transposed order pools and fluxes

    ! Determine order of dimensions
    if (present(trans)) then
       itrans = trans
    else
       itrans = .false.
    endif

    ! Check sizes
    if (itrans) then
       nn    = size(Ci,2)
       nland = size(Ci,1)
       allocate(ialpha(nland,nn,nn))
       allocate(alphaF(nland,nn,nn))
    else
       nn    = size(Ci,1)
       nland = size(Ci,2)
       allocate(ialpha(nn,nn,nland))
       allocate(alphaF(nn,nn,nland))
    endif
    if (itrans) then
       if ( (size(C,2) /= nn) .or. (size(F,2) /= nn) .or. (size(F,3) /= nn) ) then
          write(*,*) 'Error isotope_pool_model_2d: non-fitting dimensions between isotopic pools,'
          write(*,*) '                             non-isotopic pools and fluxes.'
          write(*,*) '    size(Ci,2): ', size(Ci,2)
          write(*,*) '    size(C,2):  ', size(C,2)
          write(*,*) '    size(F,2):  ', size(F,2)
          write(*,*) '    size(F,3):  ', size(F,3)
          stop 9
       endif
       if ( (size(C,1) /= nland) .or. (size(F,1) /= nland) ) then
          write(*,*) 'Error isotope_pool_model_2d: non-fitting first dimensions between isotopic pools,'
          write(*,*) '                             non-isotopic pools and fluxes.'
          write(*,*) '    size(Ci,1): ', size(Ci,1)
          write(*,*) '    size(C,1):  ', size(C,1)
          write(*,*) '    size(F,1):  ', size(F,1)
          stop 9
       endif
    else
       if ( (size(C,1) /= nn) .or. (size(F,1) /= nn) .or. (size(F,2) /= nn) ) then
          write(*,*) 'Error isotope_pool_model_2d: non-fitting dimensions between isotopic pools,'
          write(*,*) '                             non-isotopic pools and fluxes.'
          write(*,*) '    size(Ci,1): ', size(Ci,1)
          write(*,*) '    size(C,1):  ', size(C,1)
          write(*,*) '    size(F,1):  ', size(F,1)
          write(*,*) '    size(F,2):  ', size(F,2)
          stop 9
       endif
       if ( (size(C,2) /= nland) .or. (size(F,3) /= nland) ) then
          write(*,*) 'Error isotope_pool_model_2d: non-fitting first dimensions between isotopic pools,'
          write(*,*) '                             non-isotopic pools and fluxes.'
          write(*,*) '    size(Ci,2): ', size(Ci,2)
          write(*,*) '    size(C,2):  ', size(C,2)
          write(*,*) '    size(F,3):  ', size(F,3)
          stop 9
       endif
    endif
    
    ! Check F >= 0
    if (any(F < 0._dp)) then
       write(*,*) 'Error isotope_pool_model_2d: fluxes between pools must be >= 0.'
       write(*,*) '    F: ', F
       stop 9
    endif

    ! Check F(i,:) == 0. if C(i) == 0.
    if (any(eq(C,0._dp))) then
       do j=1, nland
          do i=1, nn
             if (itrans) then
                if (eq(C(j,i),0._dp)) then
                   if (any(ne(F(j,i,:),0._dp))) then
                      write(*,*) 'Error isotope_pool_model_2d:'
                      write(*,*) '    fluxes from pool i at land point j must be 0 if concentration in pool is 0.'
                      write(*,*) '    i, j, C(j,i):   ', j, i, C(j,i)
                      write(*,*) '          F(j,i,:): ', F(j,i,:)
                      stop 9
                   endif
                endif
             else
                if (eq(C(i,j),0._dp)) then
                   if (any(ne(F(i,:,j),0._dp))) then
                      write(*,*) 'Error isotope_pool_model_2d:'
                      write(*,*) '    fluxes from pool i at land point j must be 0 if concentration in pool is 0.'
                      write(*,*) '    i, j, C(i,j):   ', i, j, C(i,j)
                      write(*,*) '          F(i,:,j): ', F(i,:,j)
                      stop 9
                   endif
                endif
             endif
          end do
       end do
    endif

    ! Set optionals
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
       ! iRt = 0._dp
       iRt = 1._dp ! could be zero as well to assure no isotopic sink if C=0
       where (C > 0._dp) iRt = Ci / C
    endif
    if (present(alpha)) then
       ialpha = alpha
    else
       ialpha = 1._dp
    endif
    if (present(alpha)) then
       ialpha = alpha
    else
       ialpha = 1._dp
    endif
    if (present(beta) .and. (.not. present(T))) then
       iT = beta * C
    endif

    ! Isotope ratio
    ! R(:,:) = 0._dp
    R(:,:) = 1._dp ! could be zero as well to assure no isotopic fluxes if initial C=0
    where (C > 0._dp) R = Ci / C

    ! alpha * F
    alphaF = ialpha * F

    ! between and not-between source and sinks
    sink    = diag(ialpha,trans=itrans) * iRt * iT * dt
    source  = iRs * iS * dt
    if (itrans) then
       isink   = sum(alphaF, dim=3) * R * dt
       isource = sum(alphaF * spread(R, dim=3, ncopies=nn), dim=2) * dt
    else
       isink   = sum(alphaF, dim=2) * R * dt
       isource = sum(alphaF * spread(R, dim=2, ncopies=nn), dim=1) * dt
    endif
    ! Explicit solution
    Ci = Ci - sink + source - isink + isource

    if (present(trash)) then
        ! Check final pools
        ! Isotope pool became < 0.
        if (any(Ci < 0._dp)) then
           trash = trash + merge(abs(Ci), 0._dp, Ci < 0._dp)
           Ci = merge(0._dp, Ci, Ci < 0._dp)
        endif
        ! Non-isotope pool == 0. but isotope pool > 0.
        if (itrans) then
           Cnew = C - iT * dt + iS * dt - sum(F, dim=3)* dt + sum(F, dim=2) * dt
        else
           Cnew = C - iT * dt + iS * dt - sum(F, dim=2)* dt + sum(F, dim=1) * dt
        endif
        if (any(eq(Cnew,0._dp) .and. (Ci > 0._dp))) then
           trash = trash + merge(Ci, 0._dp, eq(Cnew,0._dp) .and. (Ci > 0._dp))
           Ci = merge(0._dp, Ci, eq(Cnew,0._dp) .and. (Ci > 0._dp))
        endif
        ! Non-isotope pool >0. but isotope pool == 0.
        ! ???
     endif

    deallocate(ialpha)
    deallocate(alphaF)

    return

  end subroutine isotope_pool_model_2d
  
  ! ------------------------------------------------------------------

  ! Diagonal elements of a matrix
  function diag_2d(matrix)

    use mo_kind, only: dp, i4

    implicit none

    real(dp), dimension(:,:), intent(in) :: matrix
    real(dp), dimension(:), allocatable  :: diag_2d

    integer(i4) :: i

    if (size(matrix,1) /= size(matrix,2)) stop 'diag_2d: array must be squared matrix.'
    if (.not. allocated(diag_2d)) allocate(diag_2d(size(matrix,1)))

    forall(i=1:size(matrix,1)) diag_2d(i) = matrix(i,i)

  end function diag_2d

  ! Diagonal elements of the two first or two last dimensions of a matrix
  function diag_3d(matrix, trans)

    use mo_kind, only: dp, i4

    implicit none

    real(dp), dimension(:,:,:), intent(in)           :: matrix
    logical,                    intent(in), optional :: trans   ! two first or two last dimensions
    real(dp), dimension(:,:),   allocatable          :: diag_3d

    integer(i4) :: i, j
    logical :: itrans

    if (present(trans)) then
       itrans = trans
    else
       itrans = .false.
    endif

    if (itrans) then
       if (size(matrix,2) /= size(matrix,3)) &
            stop 'diag_3d: array must be squared matrix in the second and third dimensions if trans=.t.'
       if (.not. allocated(diag_3d)) allocate(diag_3d(size(matrix,1),size(matrix,2)))

       do j=1, size(matrix,1)
          forall(i=1:size(matrix,2)) diag_3d(j,i) = matrix(j,i,i)
       end do
    else
       if (size(matrix,1) /= size(matrix,2)) stop 'diag_3d: array must be squared matrix in the first and second dimensions.'
       if (.not. allocated(diag_3d)) allocate(diag_3d(size(matrix,1),size(matrix,3)))

       do j=1, size(matrix,3)
          forall(i=1:size(matrix,1)) diag_3d(i,j) = matrix(i,i,j)
       end do
    endif

  end function diag_3d
  
  ! ------------------------------------------------------------------

END MODULE mo_isotope_pool_model
