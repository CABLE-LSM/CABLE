module file_utils
  use fortuno_interface_m, only: global_comm
  use cable_mpi_mod, only: mpi_grp_t
  implicit none

  private

  public :: &
    file_exists, &
    file_delete, &
    file_delete_collective

contains

  function file_exists(file_name)
    character(len=*), intent(in) :: file_name
    logical :: file_exists
    inquire(file=trim(file_name), exist=file_exists)
  end function

  subroutine file_delete(file_name)
    character(len=*), intent(in) :: file_name
    integer :: file_unit
    if (.not. file_exists(file_name)) return
    open(file=file_name, newunit=file_unit)
    close(file_unit, status="delete")
  end subroutine

  subroutine file_delete_collective(file_name)
    character(len=*), intent(in) :: file_name
    type(mpi_grp_t) :: mpi_grp
    mpi_grp = mpi_grp_t(global_comm())
    if (mpi_grp%rank == 0) call file_delete(file_name)
  end subroutine

end module file_utils
