MODULE cbl_friction_vel_module

IMPLICIT NONE

PUBLIC comp_friction_vel
PUBLIC psim
PUBLIC psis

PRIVATE

CONTAINS

SUBROUTINE comp_friction_vel(friction_vel, iter, mp, CVONK, CUMIN, CPI_C,      &
                             zetar, zref_uv, zref_tq, z0m, ua )

IMPLICIT NONE

INTEGER, INTENT(IN) :: mp
REAL, INTENT(OUT) :: friction_vel(mp)   !canopy%us
INTEGER, INTENT(IN) :: iter
! physical constants
REAL, INTENT(IN) :: CVONK
REAL, INTENT(IN):: CUMIN   
! maths & other constants
REAL, INTENT(IN) :: CPI_C  

REAL, INTENT(IN) :: zetar(mp,iter)      !canopy%zetar
REAL, INTENT(IN) :: zref_uv(mp)         !rough%zref_uv
REAL, INTENT(IN) :: zref_tq(mp)         !rough%zref_tq
REAL, INTENT(IN) :: z0m(mp)             !rough%z0m
REAL, INTENT(IN) :: ua(mp)              !met%ua

!local vars
REAL :: lower_limit(mp), rescale(mp)
REAL :: psim_1(mp), psim_2(mp), psim_arg(mp)
REAL :: z_eff(mp)

!INH: Ticket #138 %us is defined based on U(rough%zref_uv)
! but zetar based on rough%zref_tq - changes to ensure consistency
!NB no RSL incorporated here

psim_1 = psim( zetar(:,iter) * zref_uv/zref_tq, mp, CPI_C   )

!rescale = CVONK * MAX( ua, SPREAD(CUMIN,1,mp ) ) 
rescale = CVONK * MAX( ua, CUMIN ) 
z_eff = zref_uv / z0m

psim_arg = zetar(:,iter) * z0m / zref_tq
psim_2 = psim( psim_arg, mp, CPI_C  )

lower_limit = rescale / ( LOG(z_eff) - psim_1 + psim_2 )

friction_vel= MIN( MAX(1.e-6, lower_limit), 10.0 )

RETURN
END SUBROUTINE comp_friction_vel



FUNCTION psim(zeta, mp, CPI_C ) RESULT(r)

! mrr, 16-sep-92 (from function psi: mrr, edinburgh 1977)
! computes integrated stability function psim(z/l) (z/l=zeta)
! for momentum, using the businger-dyer form for unstable cases
! and the Beljaars and Holtslag (1991) form for stable cases.

INTEGER, INTENT(IN) :: mp
REAL, INTENT(IN), DIMENSION(mp) ::  zeta       !
REAL, INTENT(IN) :: CPI_C  

! function result
REAL, DIMENSION(mp) :: r

REAL, PARAMETER ::                                                          &
  gu = 16.0,         & !
  gs = 5.0,          & ! 
  a  = 1.0,          & !
  b  = 0.667,        & !
  xc = 5.0,          & !
  d  = 0.35            !

REAL, DIMENSION(mp) ::                                                      &
  x,       & !
  z,       & !
  stable,  & !
  unstable   !
  
z = 0.5 + SIGN(0.5,zeta)    ! z=1 in stable, 0 in unstable

! Beljaars and Holtslag (1991) for stable
stable   = -a*zeta - b*(zeta - xc/d)*EXP( -d*zeta) - b*xc/d
x        = (1.0 + gu*ABS(zeta))**0.25
unstable = ALOG((1.0+x*x)*(1.0+x)**2/8) - 2.0*ATAN(x) + CPI_C*0.5

r        = z*stable + (1.0-z)*unstable

END FUNCTION psim

ELEMENTAL FUNCTION psis(zeta) RESULT(r)

! mrr, 16-sep-92 (from function psi: mrr, edinburgh 1977)
! computes integrated stability function psis(z/l) (z/l=zeta)
! for scalars, using the businger-dyer form for unstable cases
! and the webb form for stable cases. see paulson (1970).

REAL, INTENT(IN)     :: zeta

REAL, PARAMETER      ::                                                     &
     gu = 16.0,        & !
     gs = 5.0,         & !
     a = 1.0,          & !
     b = 0.667,        & !
     c = 5.0,          & !
     d = 0.35

REAL                 ::                                                     &
     r,                & !
     stzeta,           & !
     ustzeta,          & !
     z,                & !
     y,                & !
     stable,           & !
     unstable

z      = 0.5 + SIGN(0.5,zeta)    ! z=1 in stable, 0 in unstable

! Beljaars and Holtslag (1991) for stable
stzeta = MAX(0.,zeta)
stable = -(1.+2./3.*a*stzeta)**(3./2.) -  &
     b*(stzeta-c/d)*EXP(-d*stzeta) - b*c/d + 1.
y      = (1.0 + gu*ABS(zeta))**0.5
unstable = 2.0 * alog((1+y)*0.5)
r   = z*stable + (1.0-z)*unstable

END FUNCTION psis


END MODULE cbl_friction_vel_module
