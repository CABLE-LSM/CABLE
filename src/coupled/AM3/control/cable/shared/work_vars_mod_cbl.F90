!******************************COPYRIGHT********************************************
! (c) CSIRO 2022.
! All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms and
! conditions set out therein.
!
! [Met Office Ref SC0237]
!******************************COPYRIGHT********************************************
MODULE work_vars_mod_cbl

!------------------------------------------------------------------------------
! Description:
!   Declares/(de)allocates/assigns CABLE "working" variables
!   These are vars reqested to be kept across the surf_couple* pathways
!   and/or timesteps. Some will be elevated to prognostics, others will be
!   removed via rewriting of the algorithm where requested
!
! This MODULE is USEd by:
!      cable_fields_mod.F90,
!      surf_couple_explicit_mod.F90,
!      surf_couple_extra_mod.F90,
!      surf_couple_implicit_mod.F90,
!      control.F90,
!      init_cable_working_vars.F90,
!      init.F90
!
! This MODULE contains 4 public Subroutines:
!      alloc_work_vars_cbl,
!      dealloc_work_vars_cbl,
!      assoc_work_vars_cbl,
!      nullify_assoc_work_vars_cbl
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!------------------------------------------------------------------------------

USE cable_canopy_type_mod,    ONLY: canopy_type
USE cable_canopy_type_mod,    ONLY: canopy_data_type
USE cable_air_type_mod,       ONLY: air_type
USE cable_air_type_mod,       ONLY: air_data_type
USE cable_met_type_mod,       ONLY: met_type
USE cable_met_type_mod,       ONLY: met_data_type
USE cable_balances_type_mod,  ONLY: balances_type
USE cable_balances_type_mod,  ONLY: balances_data_type
USE cable_soil_type_mod,      ONLY: soil_type
USE cable_soil_type_mod,      ONLY: soil_data_type
USE cable_veg_type_mod,       ONLY: veg_type
USE cable_veg_type_mod,       ONLY: veg_data_type
USE cable_soil_snow_type_mod, ONLY: soil_snow_type
USE cable_soil_snow_type_mod, ONLY: soil_snow_data_type
USE cable_radiation_type_mod, ONLY: radiation_type
USE cable_radiation_type_mod, ONLY: radiation_data_type
USE cable_roughness_type_mod, ONLY: roughness_type
USE cable_roughness_type_mod, ONLY: roughness_data_type
USE cable_climate_type_mod,   ONLY: climate_type
USE cable_climate_type_mod,   ONLY: climate_data_type
USE cable_bgc_pool_type_mod,  ONLY: bgc_pool_type
USE cable_bgc_pool_type_mod,  ONLY: bgc_pool_data_type
USE cable_sum_flux_type_mod,  ONLY: sum_flux_type
USE cable_sum_flux_type_mod,  ONLY: sum_flux_data_type

IMPLICIT NONE

PUBLIC :: alloc_work_vars_cbl
PUBLIC :: assoc_work_vars_cbl
PUBLIC :: work_vars_data_type
PUBLIC :: work_vars_type
PRIVATE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='WORK_VARS_MOD_CBL'

TYPE :: work_vars_data_type

  !fields returned @ hydrology (surf_couple_extra) level of JULES
  !computed in CABLE @ implicit (surf_couple_implicit) level
  REAL, ALLOCATABLE, PUBLIC :: snow_surft(:,:)
  REAL, ALLOCATABLE, PUBLIC :: lying_snow(:)
  REAL, ALLOCATABLE, PUBLIC :: surf_roff(:)
  REAL, ALLOCATABLE, PUBLIC :: sub_surf_roff(:)
  REAL, ALLOCATABLE, PUBLIC :: tot_tfall(:)
  REAL, ALLOCATABLE, PUBLIC :: sw_down_ij(:,:,:) 
  REAL, ALLOCATABLE, PUBLIC :: reducedLAIdue2snow(:)
!check dims
REAL, ALLOCATABLE, PUBLIC :: tot_wb_lake(:)

END TYPE work_vars_data_type

TYPE :: work_vars_type

  !fields returned @ hydrology (surf_couple_extra) level of JULES
  !computed in CABLE @ implicit (surf_couple_implicit) level
  REAL, POINTER, PUBLIC :: snow_surft(:,:)
  REAL, POINTER, PUBLIC :: lying_snow(:)
  REAL, POINTER, PUBLIC :: surf_roff(:)
  REAL, POINTER, PUBLIC :: sub_surf_roff(:)
  REAL, POINTER, PUBLIC :: tot_tfall(:)
  REAL, POINTER, PUBLIC :: sw_down_ij(:,:,:)  
  REAL, POINTER, PUBLIC :: reducedLAIdue2snow(:)
  
  TYPE(canopy_type)     :: canopy
  TYPE(air_type)        :: air
  TYPE(met_type)        :: met
  TYPE(balances_type)   :: bal
  TYPE(soil_type)       :: soil
  TYPE(veg_type)        :: veg
  TYPE(soil_snow_type)  :: ssnow
  TYPE(radiation_type)  :: rad
  TYPE(roughness_type)  :: rough
  TYPE(climate_type)    :: climate 
  TYPE(bgc_pool_type)   :: bgc
  TYPE(sum_flux_type)   :: sum_flux

END TYPE work_vars_type

! data arrays need to be declared outside of the work% TYPE
TYPE(canopy_data_type),    PUBLIC, TARGET :: canopy_data
TYPE(air_data_type),       PUBLIC, TARGET :: air_data
TYPE(met_data_type),       PUBLIC, TARGET :: met_data
TYPE(balances_data_type),  PUBLIC, TARGET :: bal_data
TYPE(veg_data_type),       PUBLIC, TARGET :: veg_data
TYPE(soil_data_type),      PUBLIC, TARGET :: soil_data
TYPE(soil_snow_data_type), PUBLIC, TARGET :: soil_snow_data
TYPE(radiation_data_type), PUBLIC, TARGET :: rad_data
TYPE(roughness_data_type), PUBLIC, TARGET :: rough_data
TYPE(climate_data_type),   PUBLIC, TARGET :: climate_data 
TYPE(bgc_pool_data_type),  PUBLIC, TARGET :: bgc_data
TYPE(sum_flux_data_type),  PUBLIC, TARGET :: sum_flux_data

CONTAINS

!===============================================================================
SUBROUTINE alloc_work_vars_cbl( row_length, rows, land_pts, nsurft, sm_levels, &
                                lsm_id, cable, work_data_cbl )

! Description:
!   Allocate the CABLE work data variables in the derived type structure

!Replacements for the argument list
USE grid_constants_mod_cbl,       ONLY: nsnl

!Common Non-science modules
USE parkind1,                 ONLY: jprb, jpim
USE yomhook,                  ONLY: lhook, dr_hook
USE jules_print_mgr,          ONLY: jules_message, jules_print, PrNorm
USE ereport_mod,              ONLY: ereport

IMPLICIT NONE

!Arguments
INTEGER, INTENT(IN) :: row_length, rows 
INTEGER, INTENT(IN) :: land_pts, nsurft,sm_levels
INTEGER, INTENT(IN) :: lsm_id, cable

TYPE(work_vars_data_type), INTENT(IN OUT) :: work_data_cbl

!-----------------------------------------------------------------------
! Local variables for error trapping
!-----------------------------------------------------------------------
INTEGER ::                                                                     &
  ERROR        = 0,                                                            &
                       ! Variable for trapping the error from each
                       ! individual call to allocate
  error_sum    = 0,                                                            &

                       ! Variable to track the sum of all errors
                       ! resulting from calls to allocate. Hence we
                       ! know that everything was successful if and
                       ! only if this is zero at the end
  errcode
                       ! Variable to use in error report

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ALLOC_WORK_VARS_CBL'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! CABLE vars to be initialized
IF ( lsm_id == cable ) THEN

  ALLOCATE( work_data_cbl%snow_surft(land_pts, nsurft), STAT = ERROR )
  ALLOCATE( work_data_cbl%lying_snow(land_pts),        STAT = ERROR )
  ALLOCATE( work_data_cbl%surf_roff(land_pts),         STAT = ERROR )
  ALLOCATE( work_data_cbl%tot_tfall(land_pts),         STAT = ERROR )
  ALLOCATE( work_data_cbl%sub_surf_roff(land_pts),     STAT = ERROR )
  !CM3 
  ALLOCATE( work_data_cbl%sw_down_ij(row_length, rows,4), STAT = ERROR )
  ALLOCATE( work_data_cbl%reducedLAIdue2snow(land_pts * nsurft), STAT = ERROR )
  !mp=sum(surft_pts) would do it
  error_sum = error_sum + ERROR

  !-----------------------------------------------------------------------
  ! Write out an error if there was one
  !-----------------------------------------------------------------------
  IF ( error_sum /= 0 )                                                        &
  CALL ereport(RoutineName, errcode,                                           &
              "Error allocating CABLE work type")

ELSE

  ALLOCATE( work_data_cbl%snow_surft(1,1), STAT = ERROR )
  ALLOCATE( work_data_cbl%lying_snow(1),  STAT = ERROR )
  ALLOCATE( work_data_cbl%surf_roff(1),   STAT = ERROR )
  ALLOCATE( work_data_cbl%tot_tfall(1),   STAT = ERROR )
  ALLOCATE( work_data_cbl%sub_surf_roff(1), STAT = ERROR )
  !CM3 
  ALLOCATE( work_data_cbl%sw_down_ij(1,1,1), STAT = ERROR )
  error_sum = error_sum + ERROR

  !-----------------------------------------------------------------------
  ! Write out an error if there was one
  !-----------------------------------------------------------------------
  IF ( error_sum /= 0 )                                                        &
  CALL ereport(RoutineName, errcode,                                           &
              "Error allocating CABLE work type")

END IF

work_data_cbl% snow_surft     (:,:)       = 0.0
work_data_cbl% lying_snow     (:)         = 0.0
work_data_cbl% surf_roff      (:)         = 0.0
work_data_cbl% sub_surf_roff  (:)         = 0.0
work_data_cbl% tot_tfall      (:)         = 0.0
work_data_cbl% sw_down_ij(:,:,:)          = 0.0
work_data_cbl%reducedLAIdue2snow(:)       = 0.0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE alloc_work_vars_cbl

!===============================================================================
SUBROUTINE dealloc_work_vars_cbl(work_data_cbl )

! Description:
!   Deallocate the CABLE work data variables in the derived type structure

!Common Non-science modules
USE parkind1,                 ONLY: jprb, jpim
USE yomhook,                  ONLY: lhook, dr_hook

IMPLICIT NONE

!Arguments
TYPE(work_vars_data_type), INTENT(IN OUT) :: work_data_cbl

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DEALLOC_WORK_VARS_CBL'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! CABLE vars to be initialized via JULES i/o
DEALLOCATE( work_data_cbl %snow_surft     )
DEALLOCATE( work_data_cbl %lying_snow    )
DEALLOCATE( work_data_cbl %surf_roff     )
DEALLOCATE( work_data_cbl %sub_surf_roff )
DEALLOCATE( work_data_cbl %tot_tfall     )
DEALLOCATE( work_data_cbl %sw_down_ij    )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE dealloc_work_vars_cbl

!===============================================================================
SUBROUTINE assoc_work_vars_cbl(work_cbl, work_data_cbl )

! Description:
!   Associate the CABLE work pointers in the derived type structure

!No USE statements other than Dr Hook
USE parkind1,    ONLY: jprb, jpim
USE yomhook,     ONLY: lhook, dr_hook

IMPLICIT NONE

!Arguments
TYPE(work_vars_type),      INTENT(IN OUT)         :: work_cbl
TYPE(work_vars_data_type), INTENT(IN OUT), TARGET :: work_data_cbl

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASSOC_WORK_VARS_CBL'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL nullify_assoc_work_vars_cbl(work_cbl)

work_cbl% snow_surft         => work_data_cbl% snow_surft
work_cbl% lying_snow         => work_data_cbl% lying_snow
work_cbl% surf_roff          => work_data_cbl% surf_roff
work_cbl% sub_surf_roff      => work_data_cbl% sub_surf_roff
work_cbl% tot_tfall          => work_data_cbl% tot_tfall
work_cbl% sw_down_ij         => work_data_cbl% sw_down_ij
work_cbl% reducedLAIdue2snow => work_data_cbl% reducedLAIdue2snow


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE assoc_work_vars_cbl

!===============================================================================
SUBROUTINE nullify_assoc_work_vars_cbl(work_cbl)

! Description:
!   Nullify the CABLE work pointers in the derived type structure

  !No USE statements other than Dr Hook
USE parkind1,    ONLY: jprb, jpim
USE yomhook,     ONLY: lhook, dr_hook

IMPLICIT NONE

!Arguments
TYPE(work_vars_type), INTENT(IN OUT) :: work_cbl

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='NULLIFY_ASSOC_WORK_VARS_CBL'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

NULLIFY( work_cbl% snow_surft     )
NULLIFY( work_cbl% lying_snow     )
NULLIFY( work_cbl% surf_roff      )
NULLIFY( work_cbl% sub_surf_roff  )
NULLIFY( work_cbl% tot_tfall      )
NULLIFY( work_cbl% sw_down_ij     )
NULLIFY( work_cbl% reducedLAIdue2snow )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE nullify_assoc_work_vars_cbl

END MODULE work_vars_mod_cbl

