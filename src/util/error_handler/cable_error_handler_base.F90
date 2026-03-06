! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

module cable_error_handler_base_mod
  !* This module defines the base error handler type for CABLE.
  ! It provides a default implementation of the error handling behaviour, which
  ! can be extended to provide custom error handling.
  use iso_fortran_env, only: error_unit
  implicit none
  private

  public :: cable_error_handler_base_t

  integer, parameter, public :: DEFAULT_ERROR_CODE = 999 !! Default error code to use when none is provided

  type cable_error_handler_base_t
    !* Base error handler type for CABLE.
    ! This type provides a default implementation of the error handling behaviour.
  contains
    procedure :: build_error_message => cable_error_handler_base_build_error_message
    procedure :: abort => cable_error_handler_base_abort
  end type cable_error_handler_base_t

contains

  function cable_error_handler_base_build_error_message(this, message, file, line, error_code) result(error_message)
    !! Build an error message string.
    class(cable_error_handler_base_t), intent(inout) :: this
    character(len=*), intent(in) :: message !! Error message to display
    character(len=*), intent(in) :: file !! Source file where the error occurred
    integer, intent(in) :: line !! Line number where the error occurred
    integer, intent(in), optional :: error_code !! Optional error code
    character(len=:), allocatable :: error_message
    character(5) :: line_string

    write (line_string, "(I5)") line
    error_message = "Error: " // file // ":" // "L" // trim(adjustl(line_string)) // ": " // message

  end function cable_error_handler_base_build_error_message

  subroutine cable_error_handler_base_abort(this, message, file, line, error_code)
    !! Default implementation of the abort procedure for the base error handler.
    class(cable_error_handler_base_t), intent(inout) :: this
    character(len=*), intent(in) :: message !! Error message to display
    character(len=*), intent(in) :: file !! Source file where the error occurred
    integer, intent(in) :: line !! Line number where the error occurred
    integer, intent(in), optional :: error_code !! Optional error code

    write(unit=error_unit, fmt="(A)") this%build_error_message(message, file, line, DEFAULT_ERROR_CODE)
    error stop DEFAULT_ERROR_CODE

  end subroutine

end module
