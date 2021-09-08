MODULE init_cable_progs_mod

USE input_mod, ONLY: fill_variables_from_file

USE logging_mod, ONLY: log_info, log_warn, log_fatal

IMPLICIT NONE

PRIVATE
PUBLIC  init_cable_progs

CONTAINS

SUBROUTINE init_cable_progs()

USE io_constants, ONLY: max_sdf_name_len, max_file_name_len, namelist_unit

USE string_utils_mod, ONLY: to_string

USE templating_mod, ONLY: tpl_has_var_name, tpl_substitute_var

USE model_interface_mod, ONLY: identifier_len, populate_var, get_var_id

USE ancil_info, ONLY: land_pts!, soil_pts, soil_index, sm_levels


IMPLICIT NONE

!------------------------------------------------------------------------------
! Description:
!   
!   Reads in information about CABLE prognostic variables for their
!   initialisation
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!------------------------------------------------------------------------------
! Work variables
INTEGER, PARAMETER :: nCABLE_VARS  = 10

INTEGER :: nvars_required      ! The number of prognostic variables that are
                               ! required in this configuration
  
CHARACTER(LEN=identifier_len) :: required_vars(nCABLE_VARS)
                               ! The variable identifiers of the required
                               ! variables

INTEGER :: nvars_file       ! The number of variables that will be set
                            ! from the given file (template?)

! Variables passed to fill_variables_from_file
CHARACTER(LEN=identifier_len) :: file_var(nCABLE_VARS)
                      ! The variable identifiers of the variables to set
                      ! from file
CHARACTER(LEN=max_sdf_name_len) :: file_var_name(nCABLE_VARS)
                      ! The name of each variable in the file

CHARACTER(LEN=max_sdf_name_len) :: file_tpl_name(nCABLE_VARS)
                      ! The name to substitute in a template for each
                      ! variable

INTEGER :: i,l  ! Index variables

INTEGER :: error  ! Error indicator

!-----------------------------------------------------------------------------
! Definition of the cable_progs namelist
!-----------------------------------------------------------------------------
CHARACTER(LEN=max_file_name_len) :: FILE
                      ! The name of the file (or variable name template) to
                      ! use for variables that need to be filled from file

INTEGER :: nvars      ! The number of variables in this section

CHARACTER(LEN=identifier_len) :: var(nCABLE_VARS)
                      ! The variable identifiers of the variables
LOGICAL :: use_file(nCABLE_VARS)
                      !   T - the variable uses the file
                      !   F - the variable is set using a constant value
CHARACTER(LEN=max_sdf_name_len) :: var_name(nCABLE_VARS)
                      ! The name of each variable in the file
CHARACTER(LEN=max_sdf_name_len) :: tpl_name(nCABLE_VARS)
                      ! The name to substitute in a template for each
                      ! variable
REAL :: const_val(nCABLE_VARS)
INTEGER:: iconst_val(nCABLE_VARS)
                      ! The constant value to use for each variable if
                      ! use_file = F for that variable

NAMELIST  / cable_progs/ FILE, nvars, use_file, var, var_name

!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
nvars_required = 0
nvars_file     = 0
nvars          = 0
use_file(:)    = .TRUE.  ! Default is for every variable to be read from file

!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
CALL log_info("init_cable_progs", "Reading CABLE_PROGS namelist...")

! First, we read the cable_progs namelist
READ(namelist_unit, NML = cable_progs, IOSTAT = error)

IF ( error /= 0 )                                                             &
   CALL log_fatal("init_cable_progs",                                         &
                  "Error reading namelist CABLE_PROGS" //                     &
                  "(IOSTAT=" // TRIM(to_string(error)) // ")")

!-----------------------------------------------------------------------------
! Set up CABLE prognostics using namelist values
!-----------------------------------------------------------------------------
! Set up the required variables
! All the CABLE variables are always required for CABLE runs
nvars_required = nCABLE_VARS
required_vars(:) = (/                                                         &
                        'ThreeLayerSnowFlag_CABLE',                           &
                        'OneLyrSnowDensity_CABLE ',                           &
                        'SnowAge_CABLE           ',                           &
                        'SnowDensity_CABLE       ',                           &
                        'SnowMass_CABLE          ',                           &
                        'SnowDepth_CABLE         ',                           &
                        'SnowTemp_CABLE          ',                           &
                        'FrozenSoilFrac_CABLE    ',                           &
                        'SoilMoisture_CABLE      ',                           &
                        'SoilTemp_CABLE          '                            &
                        /)

!-----------------------------------------------------------------------------
! Check that all the required variables are there
!-----------------------------------------------------------------------------

DO i = 1,nvars_required
  IF ( .NOT. ANY(var(1:nvars) == TRIM(required_vars(i))) )                    &
  !IF ( trim( var(1) ) /= TRIM( required_vars(1) ) )                         &
    CALL log_fatal("init_cable_progs",                                        &
                   "No value given for required variable '" //                &
                   TRIM(required_vars(i)) // "'")
END DO


!-----------------------------------------------------------------------------
! Check which variables we will be using and partition them into variables
! set to constant values and variables set from file
!-----------------------------------------------------------------------------
DO i = 1,nvars
  !-----------------------------------------------------------------------------
  ! If the variable is one of the required vars, then we will be using it
  !-----------------------------------------------------------------------------
  IF ( ANY(required_vars(1:nvars_required) == var(i)) ) THEN
    IF ( use_file(i) ) THEN
      CALL log_info("init_cable_progs",                                       &
                    "'" // TRIM(var(i)) // "' will be read from file")

      ! If the variable will be filled from file, register it here
      nvars_file = nvars_file + 1
      file_var(nvars_file) = var(i)
      file_var_name(nvars_file) = var_name(i)
      file_tpl_name(nvars_file) = tpl_name(i)
    ELSE
      ! If the variable is being set as a constant, just populate it here
      CALL log_info("init_cable_progs",                                       &
                    "'" // TRIM(var(i)) // "' will be set to a constant")

      CALL populate_var(get_var_id(var(i)), const_val = const_val(i))
    END IF
  ELSE
    ! If the variable is not a required variable, warn about not using it
    CALL log_warn("init_cable_progs",                                         &
                  "Provided variable '" // TRIM(var(i)) //                    &
                  "' is not required, so will be ignored")
  END IF
END DO

!-----------------------------------------------------------------------------
! Set variables from file
!-----------------------------------------------------------------------------
IF ( nvars_file > 0 ) THEN
  IF ( tpl_has_var_name(FILE) ) THEN
    ! We are using a file name template, so loop through the variables setting
    ! one from each file
    DO i = 1,nvars_file
        
      CALL fill_variables_from_file(                                          &
        tpl_substitute_var(FILE, file_tpl_name(i)),                           &
        (/ file_var(i) /), (/ file_var_name(i) /)                             &
      )
    END DO
  ELSE
    ! We are not using a file name template, so set all variables from the same
    ! file
       
    CALL fill_variables_from_file(                                            &
      FILE,file_var(1:nvars_file), file_var_name(1:nvars_file)                &
    )
  END IF
END IF

RETURN

END SUBROUTINE init_cable_progs

END MODULE init_cable_progs_mod
