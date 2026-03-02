module cable_error_handler_mpi_mod
  use iso_fortran_env, only: error_unit
  use cable_error_handler_base_mod, only: cable_error_handler_base_t
  use cable_error_handler_base_mod, only: DEFAULT_ERROR_CODE
  use cable_mpi_mod, only: mpi_grp_t
  implicit none
  private

  public :: cable_error_handler_mpi_t

  type, extends(cable_error_handler_base_t) :: cable_error_handler_mpi_t
    type(mpi_grp_t) :: mpi_grp
  contains
    procedure :: abort => cable_error_handler_mpi_abort
  end type cable_error_handler_mpi_t

contains

  subroutine cable_error_handler_mpi_abort(this, message, file, line, error_code)
    class(cable_error_handler_mpi_t), intent(inout) :: this
    character(len=*), intent(in) :: message
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer, intent(in), optional :: error_code

    integer :: err_code

    if (present(error_code)) then
      err_code = error_code
    else
      err_code = DEFAULT_ERROR_CODE
    end if

    write(error_unit, *) this%error_message(message, file, line, err_code)
    call this%mpi_grp%abort(err_code)

  end subroutine

end module
