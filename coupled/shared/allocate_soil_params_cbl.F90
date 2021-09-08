MODULE allocate_soil_params_mod

IMPLICIT NONE

PRIVATE

PUBLIC :: allocate_soil_parameter_type

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='allocate_soil_params_mod'

CONTAINS

SUBROUTINE allocate_soil_parameter_type(var, mp)

USE cable_params_mod,           ONLY: soil_parameter_type
USE cable_other_constants_mod,  ONLY: nrb

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

!jhan:FUDGE at least get this from right place
INTEGER, PARAMETER :: ms = 6
INTEGER, PARAMETER :: zhook_in  = 0
INTEGER, PARAMETER :: zhook_out = 1
REAL :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ALLOCATE_SOIL_PARAMETER_TYPE'

!End of header

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
  STOP
END IF

RETURN

END SUBROUTINE allocate_soil_parameter_type

END MODULE allocate_soil_params_mod
