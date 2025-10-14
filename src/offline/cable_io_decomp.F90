module cable_io_decomp_mod
  use cable_def_types_mod, only: mp, mp_global
  use cable_def_types_mod, only: mland, mland_global

  use cable_io_vars_module, only: io_decomp_t
  use cable_io_vars_module, only: xdimsize, ydimsize
  use cable_io_vars_module, only: land_x, land_y
  use cable_io_vars_module, only: landpt
  use cable_io_vars_module, only: land_decomp_start
  use cable_io_vars_module, only: patch_decomp_start

  use cable_netcdf_decomp_util_mod, only: dim_spec_t
  use cable_netcdf_decomp_util_mod, only: io_decomp_grid
  use cable_netcdf_decomp_util_mod, only: io_decomp_grid_patch
  use cable_netcdf_decomp_util_mod, only: io_decomp_land
  use cable_netcdf_decomp_util_mod, only: io_decomp_land_patch
  use cable_netcdf_decomp_util_mod, only: io_decomp_patch

  use cable_netcdf_mod, only: CABLE_NETCDF_FLOAT

  implicit none
  private

  public :: cable_io_decomp_init

contains

  subroutine cable_io_decomp_init(io_decomp)
    type(io_decomp_t), intent(out) :: io_decomp

    call io_decomp_grid( &
      land_x, &
      land_y, &
      mem_shape_spec=[dim_spec_t('land', mland)], &
      grid_shape_spec=[dim_spec_t('x', xdimsize), dim_spec_t('y', ydimsize)], &
      type=CABLE_NETCDF_FLOAT, &
      decomp=io_decomp%grid_real32_1d &
    )

  end subroutine

end module
