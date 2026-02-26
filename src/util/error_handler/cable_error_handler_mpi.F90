! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

module cable_error_handler_mpi_mod
  !* This module defines an MPI-aware error handler for CABLE.
  ! It extends the base error handler to provide functionality for aborting an MPI program.
  use iso_fortran_env, only: error_unit
  use cable_error_handler_base_mod, only: cable_error_handler_base_t
  use cable_error_handler_base_mod, only: DEFAULT_ERROR_CODE
  use cable_mpi_mod, only: mpi_grp_t
  implicit none
  private

  public :: cable_error_handler_mpi_t

  type, extends(cable_error_handler_base_t) :: cable_error_handler_mpi_t
    !* MPI-aware error handler type for CABLE.
    ! This type extends the base error handler to provide functionality for aborting an MPI program.
    type(mpi_grp_t) :: mpi_grp
  contains
    procedure :: abort => cable_error_handler_mpi_abort
  end type cable_error_handler_mpi_t

contains

  subroutine cable_error_handler_mpi_abort(this, message, file, line, error_code)
    !! Implementation of the abort procedure for the MPI-aware error handler.
    class(cable_error_handler_mpi_t), intent(inout) :: this
    character(len=*), intent(in) :: message !! Error message to display
    character(len=*), intent(in) :: file !! Source file where the error occurred
    integer, intent(in) :: line !! Line number where the error occurred
    integer, intent(in), optional :: error_code !! Optional error code

    integer :: err_code

    if (present(error_code)) then
      err_code = error_code
    else
      err_code = DEFAULT_ERROR_CODE
    end if

    write(unit=error_unit, fmt="(A)") this%build_error_message(message, file, line, err_code)
    call this%mpi_grp%abort(err_code)

  end subroutine

end module
