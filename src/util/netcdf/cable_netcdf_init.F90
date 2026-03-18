! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

submodule (cable_netcdf_mod) cable_netcdf_init_smod
  !! Submodule for initialising the I/O handler in cable_netcdf_mod.

  use cable_netcdf_nf90_mod
  use cable_netcdf_pio_mod

  implicit none

contains

  !> Module initialization procedure for the cable_netcdf_mod module.
  !! This procedure should be called before any other procedures in this module.
  module subroutine cable_netcdf_mod_init(mpi_grp)
    type(mpi_grp_t), intent(in) :: mpi_grp
      !* The MPI group for the set of processes that will be performing netCDF I/O operations.
      ! All procedures in this module should be called collectively by all
      ! processes in this group.
    if (mpi_grp%size > 1) then
      cable_netcdf_io_handler = cable_netcdf_pio_io_t(mpi_grp)
    else
      cable_netcdf_io_handler = cable_netcdf_nf90_io_t()
    end if
    call cable_netcdf_io_handler%init()
  end subroutine

end submodule cable_netcdf_init_smod
