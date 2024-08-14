#if defined(UM_JULES)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Wrapper module containing subroutines for reading cable namelists
!
MODULE read_cable_namelists_mod

! Description:
!  Contains read_cable_<namelist> and read_<urban_namelist> subroutines
!  for reading namelists into cable during a UM-JULES job.
!
! Method:
!  The unit number holding the namelist is passed as the sole argument
!  to each file.
!
! Code Owner: Please refer to ModuleLeaders.txt and UM file CodeOwners.txt
! This file belongs in section: top_level
!
! Code Description:
!  Language: FORTRAN 95.
!  This code is written to UMDP3 v8.5 programming standards.


USE umPrintMgr     ,  ONLY:                                                    &
  PrintStatus, PrStatus_Oper

USE UM_ParCore,       ONLY:                                                    &
  mype

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='READ_CABLE_NAMELISTS_MOD'

CONTAINS

! *********************************************************************

SUBROUTINE read_cable_surface_types (unitnumber)

! Description:
!  Read the cable_SURFACE_TYPES namelist

USE cable_surface_types_mod,  ONLY:                                            &
  print_nlist_cable_surface_types,                                             &
  check_cable_surface_types, read_nml_cable_surface_types,                     &
  set_derived_variables_cable_surface_types


USE land_tile_ids_mod,     ONLY: surface_type_ids_jls => surface_type_ids
USE land_tile_ids_mod_cbl, ONLY: set_surface_type_ids_cbl
USE land_tile_ids_mod_cbl, ONLY: surface_type_ids_cbl => surface_type_ids

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: unitnumber

REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_cable_SURFACE_TYPES'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL read_nml_cable_surface_types(unitnumber)
CALL set_derived_variables_cable_surface_types()
IF (PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
  CALL print_nlist_cable_surface_types()
END IF
! Set the surface_type_ids array and carry out additional checks
CALL set_surface_type_ids_cbl()
surface_type_ids_jls = surface_type_ids_cbl
CALL check_cable_surface_types()


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE read_cable_surface_types

SUBROUTINE read_cable_model_environment(unitnumber)

! Description:
! Read the cable namelist

USE cable_model_env_mod, ONLY: read_nml_cable_model_env
USE cable_model_env_mod, ONLY: set_derived_variables_cable_model_env

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: unitnumber

REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_CABLE_MODEL_ENVIRONMENT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL read_nml_cable_model_env(unitnumber)

CALL set_derived_variables_cable_model_env()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE read_cable_model_environment


END MODULE read_cable_namelists_mod
#endif
