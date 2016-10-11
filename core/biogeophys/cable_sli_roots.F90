MODULE sli_roots

  USE cable_def_types_mod, ONLY: r_2, i_d
  USE sli_numbers,       ONLY: zero, half, one, two, e3

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: setroots, getrex ! routines

  ! b1, b2 - params used to get root distribution param b (Li et al., J Hydrolo 2001).
  REAL(r_2), PARAMETER :: b1     = 24.66
  REAL(r_2), PARAMETER :: b2     = 1.59
  REAL(r_2), PARAMETER :: lambda = 1.0

  INTERFACE setroots
     MODULE PROCEDURE setroots_1d
     MODULE PROCEDURE setroots_2d
  END INTERFACE setroots

  INTERFACE getrex
     MODULE PROCEDURE getrex_1d
     MODULE PROCEDURE getrex_2d
  END INTERFACE getrex

  ! Definitions:
  ! setroots - subroutine to set current root distribution based on Li, K.Y., De Jong, R. and J.B. Boisvert (2001).
  !            An exponential root-water-uptake model with water stress compensation. J. Hydrol. 252:189-204.
  ! getrex   - subroutine to get rate of water extraction from layers.

  !*********************************************************************************************************************

CONTAINS

  !*********************************************************************************************************************

  SUBROUTINE setroots_1d(x, F10, Zr, Fs)

    IMPLICIT NONE

    REAL(r_2), DIMENSION(:), INTENT(IN)  :: x
    REAL(r_2),               INTENT(IN)  :: F10
    REAL(r_2),               INTENT(IN)  :: Zr
    REAL(r_2), DIMENSION(:), INTENT(OUT) :: Fs
    ! Sets current weighted root length density distribn (Fs).
    !      Li et al. (J Hydrolo, 2001)
    ! Definitions of arguments:
    ! x(:) - depths to bottom of layers (cm).
    ! F10  - fraction of root length density in top 10% of the root zone
    ! Zr   - rooting depth (cm).
    INTEGER(i_d)                    :: ms
    REAL(r_2)                       :: b, extr
    REAL(r_2), DIMENSION(1:size(x)) :: ext0, ext1, xend, Fi
    REAL(r_2)                       :: tmp

    ms = size(x) ! soil layers

    xend(1)    = zero
    xend(2:ms) = x(1:ms-1)
    b          = b1 * F10**b2 / Zr ! root distrib param
    extr       = exp(-b*Zr)
    ext1(:)    = exp(-b*x(:))
    ext0(1)    = one
    ext0(2:ms) = ext1(1:ms-1)
    ! get fraction of root length density in layer i
    Fi(:) = (log((one+ext0(:))/(one+ext1(:))) + half*(ext0(:)-ext1(:))) / &
         (log(two/(one+extr)) + half*(one-extr))
    Fs(:) = exp(lambda*log(Fi(:))) ! weighted Fi
    where (xend(:) >= Zr) Fs(:) = zero

    ! ensure Fs sums to one
    tmp = sum(Fs(:))
    if (tmp > zero) then
       Fs(:) = Fs(:)/tmp
    else
       Fs(:) = zero
    endif

  END SUBROUTINE setroots_1d

  SUBROUTINE setroots_2d(x, F10, Zr, Fs)

    IMPLICIT NONE

    REAL(r_2), DIMENSION(:,:), INTENT(IN)  :: x
    REAL(r_2), DIMENSION(:),   INTENT(IN)  :: F10
    REAL(r_2), DIMENSION(:),   INTENT(IN)  :: Zr
    REAL(r_2), DIMENSION(:,:), INTENT(OUT) :: Fs

    INTEGER(i_d)                      :: i
    REAL(r_2), DIMENSION(1:size(x,2)) :: out

    do i=1, size(x,1) ! landpoints
       call setroots_1d(x(i,:), F10(i), Zr(i), out)
       Fs(i,:) = out
    end do

  END SUBROUTINE setroots_2d

  !*********************************************************************************************************************

  SUBROUTINE getrex_1d(S, rex, fws, Fs, thetaS, thetaw, Etrans, gamma, dx, dt)

    ! Lai and Katul formulation for root efficiency function vh 17/07/09
    ! changed to MCs maximum formulation

    IMPLICIT NONE

    REAL(r_2), DIMENSION(:), INTENT(IN)    :: S      ! relative saturation
    REAL(r_2), DIMENSION(:), INTENT(OUT)   :: rex    ! water extraction per layer
    REAL(r_2),               INTENT(INOUT) :: fws    ! stomatal limitation factor
    REAL(r_2), DIMENSION(:), INTENT(IN)    :: Fs     ! root length density
    REAL(r_2), DIMENSION(:), INTENT(IN)    :: thetaS ! saturation soil moisture
    REAL(r_2), DIMENSION(:), INTENT(IN)    :: thetaw ! soil moisture at wiliting point
    REAL(r_2),               INTENT(IN)    :: Etrans ! total transpiration
    REAL(r_2),               INTENT(IN)    :: gamma  ! skew of Li & Katul alpha2 function
    REAL(r_2), DIMENSION(:), INTENT(IN)    :: dx     ! layer thicknesses (m)
    REAL(r_2),               INTENT(IN)    :: dt

    ! Gets rate of water extraction compatible with CABLE stomatal conductance model
    ! theta(:) - soil moisture(m3 m-3)
    ! rex(:)   - rate of water extraction by roots from layers (cm/h).
    REAL(r_2), DIMENSION(1:size(S)) :: theta, lthetar, alpha_root, delta_root
    REAL(r_2)                       :: trex

    theta(:)   = S(:)*thetaS(:)
    lthetar(:) = log(max(theta(:)-thetaw(:),e3)/thetaS(:))

    where ((theta(:)-thetaw(:)) > e3)
       alpha_root(:) = exp( gamma/max(theta(:)-thetaw(:), e3) * lthetar(:) )
    elsewhere
       alpha_root(:) = zero
    endwhere

    where (Fs(:) > zero)
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
    !where ((rex*dt) > (theta(:)-0.01_r_2)*dx(:))
    where (((rex*dt) > (theta(:)-thetaw(:))*dx(:)) .and. ((rex*dt) > zero))
       alpha_root = alpha_root*(theta(:)-thetaw(:))*dx(:)/(rex*dt)
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
    if (any(((rex*dt) > (theta(:)-thetaw(:))*dx(:)) .and. ((rex*dt) > zero))) then
       fws = zero
       ! distribute extraction according to available water
       !rex(:) = (theta(:)-0.01_r_2)*dx(:)
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

  SUBROUTINE getrex_2d(S, rex, fws, Fs, thetaS, thetaw, Etrans, gamma, dx, dt)

    IMPLICIT NONE

    REAL(r_2), DIMENSION(:,:), INTENT(IN)    :: S      ! relative saturation
    REAL(r_2), DIMENSION(:,:), INTENT(OUT)   :: rex    ! water extraction per layer
    REAL(r_2), DIMENSION(:),   INTENT(INOUT) :: fws    ! stomatal limitation factor
    REAL(r_2), DIMENSION(:,:), INTENT(IN)    :: Fs     ! root length density
    REAL(r_2), DIMENSION(:,:), INTENT(IN)    :: thetaS ! saturation soil moisture
    REAL(r_2), DIMENSION(:,:), INTENT(IN)    :: thetaw ! soil moisture at wiliting point
    REAL(r_2), DIMENSION(:),   INTENT(IN)    :: Etrans ! total transpiration
    REAL(r_2), DIMENSION(:),   INTENT(IN)    :: gamma  ! skew of Li & Katul alpha2 function
    REAL(r_2), DIMENSION(:,:), INTENT(IN)    :: dx     ! layer thicknesses (m)
    REAL(r_2),                 INTENT(IN)    :: dt

    INTEGER(i_d)                      :: i
    REAL(r_2), DIMENSION(1:size(S,2)) :: rout
    REAL(r_2)                         :: fout

    do i=1, size(S,1) ! landpoints
       fout = fws(i)
       call getrex_1d(S(i,:), rout, fout, Fs(i,:), thetaS(i,:), thetaw(i,:), Etrans(i), gamma(i), dx(i,:), dt)
       rex(i,:) = rout
       fws(i)   = fout
    end do

  END SUBROUTINE getrex_2d

  !*********************************************************************************************************************

END MODULE sli_roots
