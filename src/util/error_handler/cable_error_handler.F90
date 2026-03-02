module cable_error_handler_mod
  use cable_error_handler_base_mod, only: cable_error_handler_base_t
  implicit none
  private

  public :: cable_error_handler_base_t
  public :: cable_error_handler_mod_init
  public :: cable_error_handler_mod_end
  public :: cable_abort

  type(cable_error_handler_base_t), target :: error_handler_fallback = cable_error_handler_base_t()

  !> Error handler instance
  ! The class keyword allows for the error handler to be called polymorphically.
  ! This can be useful in scenarios where we want to customise the error handling
  ! behaviour dynamically at runtime. In those cases, we can create a new type
  ! that extends cable_error_handler_base_t and initialise the global error
  ! handler with the extended type.
  class(cable_error_handler_base_t), allocatable, target :: error_handler

contains

  subroutine cable_error_handler_mod_init(new_error_handler)
    class(cable_error_handler_base_t), intent(in), optional :: new_error_handler
    if (present(new_error_handler)) error_handler = new_error_handler
  end subroutine

  subroutine cable_error_handler_mod_end()
    if (allocated(error_handler)) deallocate(error_handler)
  end subroutine

  subroutine cable_abort(message, file, line, error_code)
    character(len=*), intent(in) :: message
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer, intent(in), optional :: error_code
    class(cable_error_handler_base_t), pointer :: err_handler

    err_handler => error_handler_fallback
    if (allocated(error_handler)) err_handler => error_handler

    call err_handler%abort(message, file, line, error_code)

  end subroutine

end module
