MODULE progs_cbl_vars_mod 

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PROGS_CBL_VARS_MOD'
! Tiled soil prognostics to be initialized from IO
TYPE progs_cbl_vars_data_type

  REAL, ALLOCATABLE, PUBLIC ::                                                &
    SoilTemp_CABLE(:,:,:),                                                    &
    SoilMoisture_CABLE(:,:,:),                                                &
    FrozenSoilFrac_CABLE(:,:,:),                                              &
    SnowDepth_CABLE(:,:,:),                                                   &
    SnowMass_CABLE(:,:,:),                                                    &
    SnowTemp_CABLE(:,:,:),                                                    &
    SnowDensity_CABLE(:,:,:),                                                 &
    ThreeLayerSnowFlag_CABLE(:,:),                                            &
    OneLyrSnowDensity_CABLE(:,:),                                             &
    SnowAge_CABLE(:,:),                                                       &
    snowOsurft(:,:)

END TYPE progs_cbl_vars_data_type

TYPE progs_cbl_vars_type

  REAL, POINTER, PUBLIC ::                                                    &
    SoilTemp_CABLE(:,:,:),                                                    &
    SoilMoisture_CABLE(:,:,:),                                                &
    FrozenSoilFrac_CABLE(:,:,:),                                              &
    SnowDepth_CABLE(:,:,:),                                                   &
    SnowMass_CABLE(:,:,:),                                                    &
    SnowTemp_CABLE(:,:,:),                                                    &
    SnowDensity_CABLE(:,:,:),                                                 &
    ThreeLayerSnowFlag_CABLE(:,:),                                            &
    OneLyrSnowDensity_CABLE(:,:),                                             &
    SnowAge_CABLE(:,:),                                                       &
    snowOsurft(:,:)
END TYPE progs_cbl_vars_type

CONTAINS

!===============================================================================
SUBROUTINE progs_cbl_vars_alloc(land_pts, nsurft, sm_levels, lsm_id, cable, progs_cbl_vars_data )

!Replacements for the argument list
USE grid_constants_cbl_mod,   ONLY: tsl

!Common Non-science modules
USE parkind1,                 ONLY: jprb, jpim
USE yomhook,                  ONLY: lhook, dr_hook
USE jules_print_mgr,          ONLY: jules_message, jules_print, PrNorm
USE ereport_mod,              ONLY: ereport

IMPLICIT NONE

!Arguments
INTEGER, INTENT(IN) :: land_pts, nsurft,sm_levels
INTEGER, INTENT(IN) :: lsm_id, cable

TYPE(progs_cbl_vars_data_type), INTENT(IN OUT) :: progs_cbl_vars_data

!-----------------------------------------------------------------------
! Local variables for error trapping
!-----------------------------------------------------------------------
INTEGER ::                                                                    &
  error        = 0,                                                           &
                       ! Variable for trapping the error from each
                       ! individual call to allocate
  error_sum    = 0,                                                           &

                       ! Variable to track the sum of all errors
                       ! resulting from calls to allocate. Hence we
                       ! know that everything was successful if and
                       ! only if this is zero at the end
  errcode 
                       ! Variable to use in error report

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PROGS_CBL_VARS_ALLOC'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( lsm_id == cable ) THEN

  ! CABLE vars to be initialized via JULES i/o
  ALLOCATE( progs_cbl_vars_data%SoilTemp_CABLE(land_pts, nsurft, sm_levels), stat = error )
  error_sum = error_sum + error

  ALLOCATE( progs_cbl_vars_data%SoilMoisture_CABLE(land_pts, nsurft, sm_levels), stat = error )
  error_sum = error_sum + error

  ALLOCATE( progs_cbl_vars_data%FrozenSoilFrac_CABLE(land_pts, nsurft, sm_levels), stat = error )
  error_sum = error_sum + error

  ALLOCATE( progs_cbl_vars_data%SnowDepth_CABLE(land_pts, nsurft, tsl), stat = error )
  error_sum = error_sum + error

  ALLOCATE( progs_cbl_vars_data%SnowMass_CABLE(land_pts,nsurft, tsl),stat = error )
  error_sum = error_sum + error

  ALLOCATE( progs_cbl_vars_data%SnowTemp_CABLE(land_pts, nsurft, tsl), stat = error )
  error_sum = error_sum + error

  ALLOCATE( progs_cbl_vars_data%SnowDensity_CABLE(land_pts,nsurft,tsl),stat = error)
  error_sum = error_sum + error

  ALLOCATE( progs_cbl_vars_data%ThreeLayerSnowFlag_CABLE(land_pts,nsurft),stat = error)
  error_sum = error_sum + error

  ALLOCATE( progs_cbl_vars_data%OneLyrSnowDensity_CABLE(land_pts,nsurft),stat = error)
  error_sum = error_sum + error

  ALLOCATE( progs_cbl_vars_data%SnowAge_CABLE(land_pts, nsurft), stat = error )
  error_sum = error_sum + error

  ALLOCATE( progs_cbl_vars_data%snowOsurft(land_pts, nsurft), stat = error )
  error_sum = error_sum + error

  !-----------------------------------------------------------------------
  ! Write out an error if there was one
  !-----------------------------------------------------------------------
  ! Check for error.
  IF ( error_sum /= 0 )                                                       &
  CALL ereport(RoutineName, errcode,                                          &
              "Error allocating CABLE prognostic array")

  progs_cbl_vars_data%SoilTemp_CABLE        = 0.0
  progs_cbl_vars_data%SoilMoisture_CABLE    = 0.0
  progs_cbl_vars_data%FrozenSoilFrac_CABLE  = 0.0
  progs_cbl_vars_data%SnowDepth_CABLE       = 0.0
  progs_cbl_vars_data%SnowMass_CABLE        = 0.0
  progs_cbl_vars_data%SnowTemp_CABLE        = 0.0
  progs_cbl_vars_data%SnowDensity_CABLE     = 0.0
  progs_cbl_vars_data%SnowAge_CABLE         = 0.0
  progs_cbl_vars_data%snowOsurft            = 0.0
  progs_cbl_vars_data%ThreeLayerSnowFlag_CABLE  = 0.0
  progs_cbl_vars_data%OneLyrSnowDensity_CABLE   = 0.0

ELSE 
  ! CABLE vars to be initialized via JULES i/o
  ALLOCATE( progs_cbl_vars_data%SoilTemp_CABLE(1, 1, 1), stat = error )
  error_sum = error_sum + error

  ALLOCATE( progs_cbl_vars_data%SoilMoisture_CABLE(1, 1, 1), stat = error )
  error_sum = error_sum + error

  ALLOCATE( progs_cbl_vars_data%FrozenSoilFrac_CABLE(1, 1, 1), stat = error )
  error_sum = error_sum + error

  ALLOCATE( progs_cbl_vars_data%SnowDepth_CABLE(1, 1, 1), stat = error )
  error_sum = error_sum + error

  ALLOCATE( progs_cbl_vars_data%SnowMass_CABLE(1,1, 1),stat = error )
  error_sum = error_sum + error

  ALLOCATE( progs_cbl_vars_data%SnowTemp_CABLE(1, 1, 1), stat = error )
  error_sum = error_sum + error

  ALLOCATE( progs_cbl_vars_data%SnowDensity_CABLE(1,1,1),stat = error)
  error_sum = error_sum + error

  ALLOCATE( progs_cbl_vars_data%ThreeLayerSnowFlag_CABLE(1,1),stat = error)
  error_sum = error_sum + error

  ALLOCATE( progs_cbl_vars_data%OneLyrSnowDensity_CABLE(1,1),stat = error)
  error_sum = error_sum + error

  ALLOCATE( progs_cbl_vars_data%SnowAge_CABLE(1, 1), stat = error )
  error_sum = error_sum + error

  ALLOCATE( progs_cbl_vars_data%snowOsurft(1, 1), stat = error )
  error_sum = error_sum + error

  !-----------------------------------------------------------------------
  ! Write out an error if there was one
  !-----------------------------------------------------------------------
  ! Check for error.
  IF ( error_sum /= 0 )                                                       &
  CALL ereport(RoutineName, errcode,                                          &
              "Error allocating CABLE prognostic array")

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE progs_cbl_vars_alloc

!===============================================================================
SUBROUTINE progs_cbl_vars_dealloc(progs_cbl_vars_data )

  !Common Non-science modules
USE parkind1,                 ONLY: jprb, jpim
USE yomhook,                  ONLY: lhook, dr_hook
  
IMPLICIT NONE
  
!Arguments
TYPE(progs_cbl_vars_data_type), INTENT(IN OUT) :: progs_cbl_vars_data
  
!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
  
CHARACTER(LEN=*), PARAMETER :: RoutineName='PROGS_CBL_VARS_DEALLOC'
  
!End of header
  
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
  
  
! CABLE vars to be initialized via JULES i/o
DEALLOCATE( progs_cbl_vars_data%snowOsurft)
DEALLOCATE( progs_cbl_vars_data%SnowAge_CABLE)
DEALLOCATE( progs_cbl_vars_data%OneLyrSnowDensity_CABLE)
DEALLOCATE( progs_cbl_vars_data%ThreeLayerSnowFlag_CABLE)
DEALLOCATE( progs_cbl_vars_data%SnowDensity_CABLE)
DEALLOCATE( progs_cbl_vars_data%SnowTemp_CABLE)
DEALLOCATE( progs_cbl_vars_data%SnowMass_CABLE)
DEALLOCATE( progs_cbl_vars_data%SnowDepth_CABLE)
DEALLOCATE( progs_cbl_vars_data%FrozenSoilFrac_CABLE)
DEALLOCATE( progs_cbl_vars_data%SoilMoisture_CABLE)
DEALLOCATE( progs_cbl_vars_data%SoilTemp_CABLE)
  
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE progs_cbl_vars_dealloc
  
!===============================================================================
SUBROUTINE progs_cbl_vars_assoc(progs_cbl_vars, progs_cbl_vars_data )

!No USE statements other than Dr Hook
USE parkind1,    ONLY: jprb, jpim
USE yomhook,     ONLY: lhook, dr_hook

IMPLICIT NONE

!Arguments
TYPE(progs_cbl_vars_type), INTENT(IN OUT) :: progs_cbl_vars
TYPE(progs_cbl_vars_data_type), INTENT(IN OUT), TARGET :: progs_cbl_vars_data

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PROGS_CBL_VARS_ASSOC'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL progs_cbl_vars_nullify(progs_cbl_vars)

progs_cbl_vars%SoilTemp_CABLE            => progs_cbl_vars_data%SoilTemp_CABLE
progs_cbl_vars%SoilMoisture_CABLE        =>                                   &
        progs_cbl_vars_data%SoilMoisture_CABLE    
progs_cbl_vars%FrozenSoilFrac_CABLE      =>                                   &
        progs_cbl_vars_data%FrozenSoilFrac_CABLE  
progs_cbl_vars%SnowDepth_CABLE           => progs_cbl_vars_data%SnowDepth_CABLE      
progs_cbl_vars%SnowMass_CABLE            => progs_cbl_vars_data%SnowMass_CABLE      
progs_cbl_vars%SnowTemp_CABLE            => progs_cbl_vars_data%SnowTemp_CABLE     
progs_cbl_vars%SnowDensity_CABLE         =>                                   &
        progs_cbl_vars_data%SnowDensity_CABLE 
progs_cbl_vars%SnowAge_CABLE             => progs_cbl_vars_data%SnowAge_CABLE    
progs_cbl_vars%snowOsurft                => progs_cbl_vars_data%snowOsurft    
progs_cbl_vars%ThreeLayerSnowFlag_CABLE  =>                                   &
        progs_cbl_vars_data%ThreeLayerSnowFlag_CABLE
progs_cbl_vars%OneLyrSnowDensity_CABLE   =>                                   &
        progs_cbl_vars_data%OneLyrSnowDensity_CABLE

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE progs_cbl_vars_assoc

!===============================================================================
SUBROUTINE progs_cbl_vars_nullify(progs_cbl_vars)

  !No USE statements other than Dr Hook
USE parkind1,    ONLY: jprb, jpim
USE yomhook,     ONLY: lhook, dr_hook
  
IMPLICIT NONE
  
!Arguments
TYPE(progs_cbl_vars_type), INTENT(IN OUT) :: progs_cbl_vars
  
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
  
CHARACTER(LEN=*), PARAMETER :: RoutineName='PROGS_CBL_VARS_NULLIFY'
  
!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

NULLIFY(progs_cbl_vars%SoilTemp_CABLE)        
NULLIFY(progs_cbl_vars%SoilMoisture_CABLE)        
NULLIFY(progs_cbl_vars%FrozenSoilFrac_CABLE)    
NULLIFY(progs_cbl_vars%SnowDepth_CABLE)             
NULLIFY(progs_cbl_vars%SnowMass_CABLE)              
NULLIFY(progs_cbl_vars%SnowTemp_CABLE)             
NULLIFY(progs_cbl_vars%SnowDensity_CABLE)      
NULLIFY(progs_cbl_vars%SnowAge_CABLE)             
NULLIFY(progs_cbl_vars%snowOsurft)                
NULLIFY(progs_cbl_vars%ThreeLayerSnowFlag_CABLE) 
NULLIFY(progs_cbl_vars%OneLyrSnowDensity_CABLE)   
  
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE progs_cbl_vars_nullify
  
END MODULE progs_cbl_vars_mod 

