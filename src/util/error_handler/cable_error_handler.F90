! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

module cable_error_handler_mod
  !* This module provides error handling functionality that can be used
  ! throughout the CABLE codebase.
  !
  ! Error handling behaviour is controlled internally by either the global error
  ! handler instance (`error_handler_global`) or the fallback error handler
  ! instance (`error_handler_fallback`). The global error handler instance can be
  ! set via `cable_error_handler_set` to provide custom error handling behaviour,
  ! while the fallback error handler instance provides a default implementation of
  ! the error handling behaviour. The global error handler is polymorphic which
  ! allows for customising the error handling behaviour dynamically at runtime. To
  ! do this we can create a new type that extends `cable_error_handler_base_t` with
  ! the new error handling behaviour and set the global error handler to an
  ! instance of the extended type.

  use cable_error_handler_base_mod, only: cable_error_handler_base_t

  implicit none
  private

  public :: cable_error_handler_base_t
  public :: cable_error_handler_set
  public :: cable_error_handler_free
  public :: cable_abort

  type(cable_error_handler_base_t), target :: error_handler_fallback = cable_error_handler_base_t()

  class(cable_error_handler_base_t), allocatable, target :: error_handler_global

contains

  subroutine cable_error_handler_set(new_error_handler)
    !! Set the global error handler instance.
    class(cable_error_handler_base_t), intent(in) :: new_error_handler
      !! New error handler instance to set as the global error handler.
    error_handler_global = new_error_handler
  end subroutine

  subroutine cable_error_handler_free()
    !! Free the global error handler instance.
    if (allocated(error_handler_global)) deallocate(error_handler_global)
  end subroutine

  subroutine cable_abort(message, file, line, error_code)
    !! Abort CABLE with an error message.
    character(len=*), intent(in) :: message !! Error message to display
    character(len=*), intent(in) :: file !! Source file where the error occurred
    integer, intent(in) :: line !! Line number where the error occurred
    integer, intent(in), optional :: error_code !! Optional error code
    class(cable_error_handler_base_t), pointer :: error_handler

    error_handler => error_handler_fallback
    if (allocated(error_handler_global)) error_handler => error_handler_global

    call error_handler%abort(message, file, line, error_code)

  end subroutine

end module
