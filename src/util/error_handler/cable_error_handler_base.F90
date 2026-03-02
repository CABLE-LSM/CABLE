module cable_error_handler_base_mod
  use iso_fortran_env, only: error_unit
  implicit none
  private
  public :: cable_error_handler_base_t

  integer, parameter, public :: DEFAULT_ERROR_CODE = 999

  type cable_error_handler_base_t
  contains
    procedure :: error_message => cable_error_handler_base_build_error_message
    procedure :: abort => cable_error_handler_base_abort
  end type cable_error_handler_base_t

contains

  function cable_error_handler_base_build_error_message(this, message, file, line, error_code) result(error_message)
    class(cable_error_handler_base_t), intent(inout) :: this
    character(len=*), intent(in) :: message
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer, intent(in), optional :: error_code
    character(len=:), allocatable :: error_message
    character(5) :: line_string

    write (line_string, "(I5)") line
    error_message = "Error: " // file // ":" // "L" // trim(adjustl(line_string)) // ": " // message

  end function cable_error_handler_base_build_error_message

  subroutine cable_error_handler_base_abort(this, message, file, line, error_code)
    class(cable_error_handler_base_t), intent(inout) :: this
    character(len=*), intent(in) :: message
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer, intent(in), optional :: error_code

    write(error_unit, *) this%error_message(message, file, line, DEFAULT_ERROR_CODE)
    error stop DEFAULT_ERROR_CODE

  end subroutine

end module
