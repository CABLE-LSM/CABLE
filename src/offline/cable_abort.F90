!==============================================================================
! This source code is part of the
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
!
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) is located
! in each directory containing CABLE code.
!
! ==============================================================================
! Purpose: Error management for CABLE offline
!
! Contact: Bernard.Pak@csiro.au
!
! History: Developed by Gab Abramowitz and Harvey Davies
!
!
! ==============================================================================

MODULE cable_abort_module

  USE iso_fortran_env, ONLY: error_unit
  USE cable_mpi_mod, ONLY: mpi_grp_t

  IMPLICIT NONE

  TYPE(mpi_grp_t), PRIVATE :: mpi_grp_global

CONTAINS

  SUBROUTINE cable_abort_module_init(mpi_grp)
    !! Initialise abort module
    TYPE(mpi_grp_t), intent(in) :: mpi_grp
    mpi_grp_global = mpi_grp
  END SUBROUTINE

  SUBROUTINE cable_abort(message, file, line)
    !! Print the error message and stop the code
    CHARACTER(LEN=*), INTENT(IN) :: message !! Error message
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: file
    INTEGER, INTENT(IN), OPTIONAL :: line
    CHARACTER(5) :: line_string

    IF (present(file) .AND. present(line)) THEN
      WRITE (line_string, "(I5)") line
      WRITE (error_unit, *) file // ":" // trim(adjustl(line_string)) // ": " // message
    ELSE
      WRITE (error_unit, *) message
    END IF
    call mpi_grp_global%abort()

  END SUBROUTINE

  !==============================================================================
  !
  ! Name: nc_abort
  !
  ! Purpose: For NETCDF errors. Prints an error message then stops the code
  !
  ! CALLed from: get_restart_data
  !              extrarestart
  !              get_default_lai
  !              open_met_file
  !              get_met_data
  !              close_met_file
  !              load_parameters
  !              open_output_file
  !              write_output
  !              close_output_file
  !              create_restart
  !              read_gridinfo
  !              readpar_i
  !              readpar_r
  !              readpar_rd
  !              readpar_r2
  !              readpar_r2d
  !              define_output_variable_r1
  !              define_output_variable_r2
  !              define_output_parameter_r1
  !              define_output_parameter_r2
  !              write_output_variable_r1
  !              write_output_variable_r2
  !              write_output_parameter_r1
  !              write_output_parameter_r1d
  !              write_output_parameter_r2
  !              write_output_parameter_r2d
  !
  ! MODULEs used: netcdf
  !
  !==============================================================================

  SUBROUTINE nc_abort(ok, message)

    USE netcdf

    ! Input arguments
    CHARACTER(LEN=*), INTENT(IN) :: message
    INTEGER, INTENT(IN) :: ok

    WRITE (*, *) message ! error from subroutine
    WRITE (*, *) NF90_STRERROR(ok) ! netcdf error details

    STOP

  END SUBROUTINE nc_abort

END MODULE cable_abort_module
