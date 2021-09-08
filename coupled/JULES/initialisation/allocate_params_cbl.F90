MODULE allocate_veg_params_mod

!Common Non-science modules
USE yomhook,                  ONLY: lhook, dr_hook
USE ereport_mod,              ONLY: ereport
USE parkind1,                 ONLY: jprb, jpim

IMPLICIT NONE

PRIVATE

PUBLIC :: allocate_veg_parameter_type

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='allocate_veg_params_mod'

CONTAINS

SUBROUTINE allocate_veg_parameter_type(var, mp)

USE cable_params_mod,           ONLY: veg_parameter_type
!H! Clarify Components of impending radiation (direct and diffuse)
!H! AND radiation bands per component (NIR and VIS)
USE cable_other_constants_mod,  ONLY: nrb !used as ambiguosly AND inconsistently as # rad. bands
USE jules_soil_mod,             ONLY: sm_levels   ! number of soil levels
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

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ALLOCATE_VEG_PARAMETER_TYPE'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

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
  IF ( error_sum /= 0 )                                                       &
    CALL ereport("allocate_cable_arrays", errcode,                            &
                 "Error allocating CABLE model veg parameter arrays")
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE allocate_veg_parameter_type

END MODULE allocate_veg_params_mod


MODULE allocate_soil_params_mod

!Common Non-science modules
USE yomhook,                  ONLY: lhook, dr_hook
USE ereport_mod,              ONLY: ereport
USE parkind1,                 ONLY: jprb, jpim

IMPLICIT NONE

PRIVATE

PUBLIC :: allocate_soil_parameter_type

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='allocate_soil_params_mod'

CONTAINS

SUBROUTINE allocate_soil_parameter_type(var, mp)

USE cable_params_mod,           ONLY: soil_parameter_type
USE cable_other_constants_mod,  ONLY: nrb
USE jules_soil_mod,             ONLY: ms => sm_levels   ! number of soil levels

!-----------------------------------------------------------------------------
! Description:
!   Allocates soil parameter arrays.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

TYPE(soil_parameter_type), INTENT(INOUT) :: var
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

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ALLOCATE_SOIL_PARAMETER_TYPE'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ALLOCATE( var% isoilm(mp), stat = error )
error_sum = error_sum + error
ALLOCATE( var% albsoil(mp, nrb), stat = error )
error_sum = error_sum + error
ALLOCATE( var% albsoilf(mp), stat = error )
error_sum = error_sum + error
ALLOCATE( var% soilcol(mp), stat = error )
error_sum = error_sum + error
!H! from cable_define_types
ALLOCATE( var% bch(mp) )
ALLOCATE( var% c3(mp) )
ALLOCATE( var% clay(mp) )
ALLOCATE( var% css(mp) )
ALLOCATE( var% hsbh(mp) )
ALLOCATE( var% hyds(mp) )
ALLOCATE( var% i2bp3(mp) )
ALLOCATE( var% ibp2(mp) )
ALLOCATE( var% rhosoil(mp) )
ALLOCATE( var% sand(mp) )
ALLOCATE( var% sfc(mp) )
ALLOCATE( var% silt(mp) )
ALLOCATE( var% ssat(mp) )
ALLOCATE( var% sucs(mp) )
ALLOCATE( var% swilt(mp) )
ALLOCATE( var% zse(ms) )
ALLOCATE( var% zshh(ms+1) )
ALLOCATE( var% cnsd(mp) )
ALLOCATE( var% pwb_min(mp) )
!mrd561
!MD
!Aquifer properties
ALLOCATE( var%GWhyds_vec(mp) )
ALLOCATE( var%GWsucs_vec(mp) )
ALLOCATE( var%GWbch_vec(mp) )
ALLOCATE( var%GWssat_vec(mp) )
ALLOCATE( var%GWwatr(mp) )
var%GWwatr(:) = 0.05
ALLOCATE( var%GWz(mp) )
ALLOCATE( var%GWdz(mp) )
ALLOCATE( var%GWrhosoil_vec(mp) )
!soil properties (vary by layer)
ALLOCATE( var% zse_vec(mp,ms) )
ALLOCATE( var% heat_cap_lower_limit(mp,ms) )
ALLOCATE( var% css_vec(mp,ms) )
ALLOCATE( var% cnsd_vec(mp,ms) )
ALLOCATE( var%hyds_vec(mp,ms) )
ALLOCATE( var%sucs_vec(mp,ms) )
ALLOCATE( var%bch_vec(mp,ms) )
ALLOCATE( var%ssat_vec(mp,ms) )
ALLOCATE( var%watr(mp,ms) )
var%watr(:,:) = 0.05
ALLOCATE( var%sfc_vec(mp,ms) )
ALLOCATE( var%swilt_vec(mp,ms) )
ALLOCATE( var%sand_vec(mp,ms) )
ALLOCATE( var%clay_vec(mp,ms) )
ALLOCATE( var%silt_vec(mp,ms) )
ALLOCATE( var%org_vec(mp,ms) )
ALLOCATE( var%rhosoil_vec(mp,ms) )

ALLOCATE( var%drain_dens(mp) )
ALLOCATE( var%elev(mp) )
ALLOCATE( var%elev_std(mp) )
ALLOCATE( var%slope(mp) )
ALLOCATE( var%slope_std(mp) )

! Allocate variables for SLI soil model:
ALLOCATE ( var % nhorizons(mp) )
ALLOCATE ( var % ishorizon(mp,ms) )
ALLOCATE ( var % clitt(mp) )
ALLOCATE ( var % zeta(mp) )
ALLOCATE ( var % fsatmax(mp) )
!ALLOCATE ( var % swilt_vec(mp,ms) )
!ALLOCATE ( var % ssat_vec(mp,ms) )
!ALLOCATE ( var % sfc_vec(mp,ms) )
IF ( .NOT. (ASSOCIATED(var % swilt_vec))) ALLOCATE ( var % swilt_vec(mp,ms) )
IF ( .NOT. (ASSOCIATED(var % ssat_vec))) ALLOCATE ( var % ssat_vec(mp,ms) )
IF ( .NOT. (ASSOCIATED(var % sfc_vec))) ALLOCATE ( var % sfc_vec(mp,ms) )



IF (error_sum == 0) THEN
  var%isoilm(:) = 0
  var%albsoil(:,:) = 0.0
  var%albsoilf(:) = 0.0
  var%soilcol(:) = 0.0
ELSE
  IF ( error_sum /= 0 )                                                       &
    CALL ereport("allocate_cable_arrays", errcode,                            &
                 "Error allocating CABLE model soil parameter arrays")
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE allocate_soil_parameter_type

END MODULE allocate_soil_params_mod
