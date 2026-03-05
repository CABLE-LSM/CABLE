! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

module cable_netcdf_pio_mod
  !! A stub implementation of the netCDF file handling interface for when PIO is unavailable.
  use cable_netcdf_mod
  use cable_mpi_mod, only: mpi_grp_t
  use cable_abort_module, only: cable_abort
  use cable_netcdf_stub_types_mod, only: cable_netcdf_stub_io_t
  implicit none

  type, extends(cable_netcdf_stub_io_t) :: cable_netcdf_pio_io_t
  end type

  interface cable_netcdf_pio_io_t
    procedure cable_netcdf_pio_io_constructor
  end interface

contains

  function cable_netcdf_pio_io_constructor(mpi_grp) result(this)
    type(cable_netcdf_pio_io_t) :: this
    type(mpi_grp_t), intent(in) :: mpi_grp
    call cable_abort("Error instantiating cable_netcdf_pio_io_t: PIO support not available", file=__FILE__, line=__LINE__)
    this = cable_netcdf_pio_io_t()
  end function

end module cable_netcdf_pio_mod
