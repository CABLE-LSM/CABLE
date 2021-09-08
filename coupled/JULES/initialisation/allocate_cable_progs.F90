MODULE allocate_cable_progs_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ALLOCATE_CABLE_ARRAYS_MOD'

CONTAINS
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine ALLOCATE_JULES_ARRAYS
!
! Description: Routine that allocates memory to the JULES arrays
! This assume that the values in the jules_surface_types module have been set
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.
!
!   Code Owner: Please refer to ModuleLeaders.txt and UM file CodeOwners.txt
!   This file belongs in section: Land

SUBROUTINE allocate_cable_progs()

!Replacements for the argument list
USE ancil_info,               ONLY: land_pts, nsurft
USE jules_soil_mod,           ONLY: sm_levels
USE cable_prognostic_info_mod
USE prognostics,              ONLY: snow_surft

!Common Non-science modules
USE parkind1,                 ONLY: jprb, jpim
USE yomhook,                  ONLY: lhook, dr_hook
USE jules_print_mgr,          ONLY: jules_message, jules_print, PrNorm
USE ereport_mod,              ONLY: ereport
USE cable_surface_types_mod,           ONLY: msn_cable

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Allocates the model arrays using sizes determined during initialisation
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!
! Modules are ordered alphabetically. Please respect it if you add something
! new.
!
! The allocation statements are also grouped into 'JULES', 'UM' and 'Common'
! according to whether they are needed for just JULES, the UM or both.
!
!-----------------------------------------------------------------------------

#if defined(UM_JULES)
! Input variables for dimensioning
INTEGER, INTENT(IN) ::                                                        &
  land_pts,                                                                   &
                       ! Number of land points
  nsurft,                                                                     &
                       ! Number of surface tiles
  sm_levels
                       ! Number of soil layers
#endif

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

  temp_size,                                                                  &
  temp_tiles,                                                                 &
  temp_layers,                                                                &
                       ! For storing the size of array to allocate for variables
                       ! that are sometimes set to size 1.
                       ! Removes some duplicate allocate statements
  errcode 
                       ! Variable to use in error report

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ALLOCATE_JULES_ARRAYS'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
errcode = 101

! CABLE vars to be initialized via JULES i/o
ALLOCATE( SoilTemp_CABLE(land_pts, nsurft, sm_levels), stat = error )
error_sum = error_sum + error

ALLOCATE( SoilMoisture_CABLE(land_pts, nsurft, sm_levels), stat = error )
error_sum = error_sum + error

ALLOCATE( FrozenSoilFrac_CABLE(land_pts, nsurft, sm_levels), stat = error )
error_sum = error_sum + error

ALLOCATE( SnowDepth_CABLE(land_pts, nsurft, msn_cable), stat = error )
error_sum = error_sum + error

ALLOCATE( SnowMass_CABLE(land_pts,nsurft, msn_cable),stat = error )
error_sum = error_sum + error

ALLOCATE( SnowTemp_CABLE(land_pts, nsurft, msn_cable), stat = error )
error_sum = error_sum + error

ALLOCATE( SnowDensity_CABLE(land_pts,nsurft,msn_cable),stat = error)
error_sum = error_sum + error

ALLOCATE( ThreeLayerSnowFlag_CABLE(land_pts,nsurft),stat = error)
error_sum = error_sum + error

ALLOCATE( OneLyrSnowDensity_CABLE(land_pts,nsurft),stat = error)
error_sum = error_sum + error

ALLOCATE( SnowAge_CABLE(land_pts, nsurft), stat = error )
error_sum = error_sum + error

ALLOCATE( snowOsurft(land_pts, nsurft), stat = error )
error_sum = error_sum + error

snowOsurft = snow_surft

!-----------------------------------------------------------------------
! Write out an error if there was one
!-----------------------------------------------------------------------
! Check for error.
IF ( error_sum /= 0 )                                                         &
  CALL ereport("allocate_cable_arrays", errcode,                              &
               "Error allocating CABLE model arrays")

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE allocate_cable_progs

END MODULE allocate_cable_progs_mod
