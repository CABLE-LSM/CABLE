submodule (cable_netcdf_mod) cable_netcdf_internal
  use cable_netcdf_nf90_mod
  use cable_netcdf_pio_mod
  implicit none

contains

  module subroutine cable_netcdf_mod_init(mpi_grp)
    type(mpi_grp_t), intent(in) :: mpi_grp
    if (mpi_grp%size > 1) then
      cable_netcdf_io_handler = cable_netcdf_pio_io_t(mpi_grp)
    else
      cable_netcdf_io_handler = cable_netcdf_nf90_io_t()
    end if
    call cable_netcdf_io_handler%init()
  end subroutine

  module subroutine cable_netcdf_mod_end()
    call cable_netcdf_io_handler%finalise()
  end subroutine

  module function cable_netcdf_create_file(path, iotype, mode) result(file)
    character(len=*), intent(in) :: path
    integer, intent(in) :: iotype
    integer, intent(in), optional :: mode
    class(cable_netcdf_file_t), allocatable :: file
    file = cable_netcdf_io_handler%create_file(path, iotype, mode)
  end function

  module function cable_netcdf_open_file(path, iotype, mode) result(file)
    character(len=*), intent(in) :: path
    integer, intent(in) :: iotype
    integer, intent(in), optional :: mode
    class(cable_netcdf_file_t), allocatable :: file
    file = cable_netcdf_io_handler%open_file(path, iotype, mode)
  end function

  module function cable_netcdf_create_decomp(compmap, dims, type) result(decomp)
    integer, intent(in) :: compmap(:), dims(:), type
    class(cable_netcdf_decomp_t), allocatable :: decomp
    decomp = cable_netcdf_io_handler%create_decomp(compmap, dims, type)
  end function

end submodule cable_netcdf_internal
