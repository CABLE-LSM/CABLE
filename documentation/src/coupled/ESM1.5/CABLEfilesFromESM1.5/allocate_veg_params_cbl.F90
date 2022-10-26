MODULE allocate_veg_params_mod

IMPLICIT NONE

PRIVATE

PUBLIC :: allocate_veg_parameter_type

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='allocate_veg_params_mod'

CONTAINS

SUBROUTINE allocate_veg_parameter_type(var, mp)

USE cable_params_mod,           ONLY: veg_parameter_type
USE cable_other_constants_mod,  ONLY: nrb !used as ambiguosly AND inconsistently as # rad. bands
USE cable_other_constants_mod,  ONLY: nscs ! number of soil carbon stores
USE cable_other_constants_mod,  ONLY: nvcs ! number of vegetation carbon stores

!-----------------------------------------------------------------------------
! Description:
!   Allocates variable arrays for vegetation parameters.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

TYPE(veg_parameter_type), INTENT(INOUT) :: var
INTEGER, INTENT(IN) :: mp

INTEGER ::                                                                    &
  error        = 0,                                                           &
                       ! Variable for trapping the error from each
                       ! individual call to allocate
  error_sum    = 0,                                                           &

                       ! Variable to track the sum of all errors
                       ! resulting from calls to allocate. Hence we
                       ! know that everything was successful if and
                       ! only if this is zero at the end
  errcode      = 101
                       ! Variable to use in error report
!jhan:FUDGE at least get this from right place
INTEGER, PARAMETER :: sm_levels = 6
INTEGER, PARAMETER :: zhook_in  = 0
INTEGER, PARAMETER :: zhook_out = 1
REAL :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ALLOCATE_VEG_PARAMETER_TYPE'

!End of header

ALLOCATE( var% iveg(mp), stat = error )
error_sum = error_sum + error
ALLOCATE( var% hc(mp), stat = error )
error_sum = error_sum + error
ALLOCATE( var% vlai(mp), stat = error )
error_sum = error_sum + error
ALLOCATE( var% xfang(mp), stat = error )
error_sum = error_sum + error
ALLOCATE( var%refl(mp,nrb), stat = error )
error_sum = error_sum + error
ALLOCATE( var%taul(mp,nrb), stat = error )
error_sum = error_sum + error

ALLOCATE( var% canst1(mp) )
ALLOCATE( var% dleaf(mp) )
ALLOCATE( var% ejmax(mp) )
!H!    ALLOCATE( var% iLU(mp) )
ALLOCATE( var% meth(mp) )
ALLOCATE( var% frac4(mp) )
ALLOCATE( var% xalbnir(mp) )
ALLOCATE( var% rp20(mp) )
ALLOCATE( var% rpcoef(mp) )
ALLOCATE( var% rs20(mp) )
ALLOCATE( var% shelrb(mp) )
ALLOCATE( var% vegcf(mp) )
ALLOCATE( var% tminvj(mp) )
!H!    ALLOCATE( var% toptvj(mp) )
ALLOCATE( var% tmaxvj(mp) )
ALLOCATE( var% vbeta(mp) )
ALLOCATE( var% vcmax(mp) )
ALLOCATE( var%extkn(mp) )
ALLOCATE( var%wai(mp) )
!H!    ALLOCATE( var%deciduous(mp) )
    !was nrb(=3), but never uses (:,3) in model
!H!    ALLOCATE( var%vlaimax(mp) )
ALLOCATE( var%a1gs(mp) )
ALLOCATE( var%d0gs(mp) )
ALLOCATE( var%alpha(mp) )
ALLOCATE( var%convex(mp) )
ALLOCATE( var%cfrd(mp) )
ALLOCATE( var%gswmin(mp) )
ALLOCATE( var%conkc0(mp) )
ALLOCATE( var%conko0(mp) )
ALLOCATE( var%ekc(mp) )
ALLOCATE( var%eko(mp) )
ALLOCATE( var% g0(mp) )   ! Ticket #56.
ALLOCATE( var% g1(mp) )   ! Ticket #56.

ALLOCATE( var%froot(mp,sm_levels) )

ALLOCATE ( var % rootbeta(mp) )
ALLOCATE ( var % GAMMA(mp) )
ALLOCATE ( var % f10(mp) )
ALLOCATE ( var % zr(mp) )
ALLOCATE ( var % clitt(mp) )

ALLOCATE(var % csoil(mp,nscs) )
ALLOCATE(var % cplant(mp,nvcs) )
ALLOCATE(var % ratecp(mp,nvcs) )
ALLOCATE(var % ratecs(mp,nscs) )
!H!    ALLOCATE ( var % disturbance_interval(mp,2) )
!H!    ALLOCATE ( var % disturbance_intensity(mp,2) )
!init to zero ~ first call or not. ALSO is redef as INTEGER. ALSO make SCALAR?
var%meth(:) = 0
IF (error_sum == 0) THEN
  var%iveg(:) = 0
  var%hc(:) = 0.0
  var%vlai(:) = 0.0
  var%refl(:,:) = 0.0
  var%taul(:,:) = 0.0
ELSE
  STOP
END IF

RETURN

END SUBROUTINE allocate_veg_parameter_type

END MODULE allocate_veg_params_mod
