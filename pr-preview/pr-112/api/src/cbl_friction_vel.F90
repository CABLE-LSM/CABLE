MODULE cbl_friction_vel_module
  !* This MODULE contains the SUBROUTINE [[comp_friction_vel]] and two
  ! FUNCTIONS ([[psim]] and [[psis]]) needed to
  ! evaluate the friction velocity over each land point/tile
  ! given the wind speed and the current estimate of the Monin-Obukhov
  ! stability parameter \(\xi\). Outputs `friction_vel`, `psim` and
  ! `psis` are used in [[define_canopy]].

IMPLICIT NONE

PUBLIC comp_friction_vel
PUBLIC psim
PUBLIC psis

PRIVATE

CONTAINS

SUBROUTINE comp_friction_vel(friction_vel, iter, mp, CVONK, CUMIN, CPI_C,      &
     zetar, zref_uv, zref_tq, z0m, ua )
  !*## Purpose
  !
  ! This SUBROUTINE evaluates the value of the friction velocity \(u_*\) as used
  ! during the iteration loop of the Monin-Obukhov (MO) similarity theory in
  ! [[define_canopy]]. \(u_*\) is used when quantifying the stability
  ! parameter `canopy%zetar` in [[update_zetar]], and when evaluating the
  ! resistance network components `rt0` and `rt1` in [[define_canopy]].
  ! The friction velocity quantifies the magnitude of turbulent mixing that
  ! is occurring as well as the flux of momentum from the atmosphere to
  ! the land.
  !
  ! The output `friction_vel` is known as `canopy%us` elsewhere in the code.
  !
  !
  !## Method
  !
  ! The basic formula is derived from the modified log-law for the
  ! vertical profile of wind speed near to the ground namely
  !
  ! \( u_{*} = U(z_{ref}) / ( \log [z_{ref}/z_{0m}] -
  !      \psi_m[\xi * z_{ref}/z_{refTq}]
  !    + \psi_m[\xi * z_{0m}/z_{refTq}] ) \)
  !
  ! with \(\psi_m\) the integrated similarity function given by [[psim]],
  ! \(U\) the wind speed at height \(z_{ref}\),
  ! \(\xi\) the current value of the stability parameter, and \(z_{0m}\) the
  ! roughness length.  
  ! A minimum value is applied to the input wind speed, `ua`= \(U(z_{ref})\),
  ! to assist with convergence of the MO iteration in light wind conditions.
  ! Small and large value limits are applied to the evaluated \(u_*\).
  !
  ! **NOTE** Most literature references assume that the reference levels
  ! for wind and temperature are equal inside the arguments to the
  ! stability functions \(\psi_m\) and \(\psi_s\) - in CABLE this need
  ! not be true.
  ! \(\xi\) is always defined relative to the reference level for temperature,
  ! in \(\psi_m\) the reference height for wind is needed.
  ! Consequently the calls to [[psim]] include conversion
  ! factors to account for the different reference levels.
  ! `zref_uv` and `zref_tq` passed to comp_friction_vel are the heights above
  ! the displacement height `rough%disp`.
  !
  !## References
  !
  ![Kowalczyk et al. (2006)](http://www.cmar.csiro.au/e-print/open/kowalczykea_2006a.pdf)
  ! - section 3.1, equations 1-9.

IMPLICIT NONE

INTEGER, INTENT(IN) :: mp               !! size of cable vector of land points (-)
REAL, INTENT(OUT) :: friction_vel(mp)   !! friction velocity (ms\(^{-1}\))
INTEGER, INTENT(IN) :: iter             !! MO iteration counter (-)
! physical constants
REAL, INTENT(IN) :: CVONK               !! von Karman constant (-)
REAL, INTENT(IN):: CUMIN                !! minimum value of wind speed (ms\(^{-1}\))
! maths & other constants
REAL, INTENT(IN) :: CPI_C               !! \(\pi\) (-)

REAL, INTENT(IN) :: zetar(mp,iter)      !! stability parameter - see [[update_zetar]] `canopy%zetar` (-)
REAL, INTENT(IN) :: zref_uv(mp)         !! reference height for wind `rough%zref_uv` (m)
REAL, INTENT(IN) :: zref_tq(mp)         !! reference height for temperature and humidity `rough%zref_tq` (m)
REAL, INTENT(IN) :: z0m(mp)             !! roughness length `rough%z0m` (m)
REAL, INTENT(IN) :: ua(mp)              !! wind speed `met%ua` (ms\(^{-1}\))

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
  !*## Purpose
  !
  ! Evaluates the integrated similarity function for momentum transfer,
  ! \(\psi_m(\xi)\).
  !
  !## Method
  !
  ! Uses the Businger-Dyer form for unstable conditions (\(\xi<0\)) and the
  ! Beljaars-Holtslag form for stable conditions (\(\xi>0\)).
  !
  ! This function is used in the evaluation of the friction velocity in
  ! [[comp_friction_vel]].
  !
  !## References
  ! - [Beljaars and Holtslag (1991)](https://doi.org/10.1175/1520-0450(1991)030<0327:FPOLSF>2.0.CO;2)
  ! - [Businger et al. (1971)](https://doi.org/10.1175/1520-0469(1971)028<0181:FPRITA>2.0.CO;2)
  ! - [Dyer (1974)](https://doi.org/10.1007/BF00240838)

! mrr, 16-sep-92 (from function psi: mrr, edinburgh 1977)
! computes integrated stability function psim(z/l) (z/l=zeta)
! for momentum, using the businger-dyer form for unstable cases
! and the Beljaars and Holtslag (1991) form for stable cases.

INTEGER, INTENT(IN) :: mp                !! size of cable vector of land points (-)
REAL, INTENT(IN), DIMENSION(mp) ::  zeta !! current value of stability parameter \(\xi\) (-)
REAL, INTENT(IN) :: CPI_C                !! \( \pi \) (-)

! function result
REAL, DIMENSION(mp) :: r                 ! result of FUNCTION

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
  !*## Purpose
  !
  ! Evaluates the integrated similarity function for turbulent transfer
  ! of scalars, \(\psi_s(\xi)\).
  !
  !## Method
  !
  ! Uses the Businger-Dyer form for unstable conditions (\(\xi<0\)) and the
  ! Beljaars-Holtslag form for stable conditions (\(\xi>0\)).
  !
  ! This function is used in the evaluation of the resistance network component
  ! `rt1usc` in [[define_canopy]].
  !
  !## References
  ! - [Beljaars and Holtslag (1991)](https://doi.org/10.1175/1520-0450(1991)030<0327:FPOLSF>2.0.CO;2)
  ! - [Businger et al. (1971)](https://doi.org/10.1175/1520-0469(1971)028<0181:FPRITA>2.0.CO;2)
  ! - [Dyer (1974)](https://doi.org/10.1007/BF00240838)

! mrr, 16-sep-92 (from function psi: mrr, edinburgh 1977)
! computes integrated stability function psis(z/l) (z/l=zeta)
! for scalars, using the businger-dyer form for unstable cases
! and the webb form for stable cases. see paulson (1970).

REAL, INTENT(IN)     :: zeta  !! current value of stability parameter \(\xi\) (-)

REAL, PARAMETER      ::                                                     &
     gu = 16.0,        & !
     gs = 5.0,         & !
     a = 1.0,          & !
     b = 0.667,        & !
     c = 5.0,          & !
     d = 0.35

REAL                 ::                                                     &
     r,                & !! OUT \(\psi_s (\xi) \)
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
