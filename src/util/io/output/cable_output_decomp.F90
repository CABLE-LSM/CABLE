module cable_output_decomp_mod

  use cable_error_handler_mod, only: cable_abort

  use cable_array_utils_mod, only: array_eq

  use cable_def_types_mod, only: mp
  use cable_def_types_mod, only: mp_global
  use cable_def_types_mod, only: mland
  use cable_def_types_mod, only: mland_global
  use cable_def_types_mod, only: ms
  use cable_def_types_mod, only: msn
  use cable_def_types_mod, only: nrb
  use cable_def_types_mod, only: ncs
  use cable_def_types_mod, only: ncp

  use cable_io_vars_module, only: xdimsize
  use cable_io_vars_module, only: ydimsize
  use cable_io_vars_module, only: max_vegpatches
  use cable_io_vars_module, only: land_x, land_y
  use cable_io_vars_module, only: landpt
  use cable_io_vars_module, only: land_decomp_start
  use cable_io_vars_module, only: patch_decomp_start

  use cable_netcdf_mod, only: cable_netcdf_decomp_t
  use cable_netcdf_mod, only: CABLE_NETCDF_INT
  use cable_netcdf_mod, only: CABLE_NETCDF_FLOAT
  use cable_netcdf_mod, only: CABLE_NETCDF_DOUBLE

  use cable_netcdf_decomp_util_mod, only: io_decomp_land_to_x_y
  use cable_netcdf_decomp_util_mod, only: io_decomp_patch_to_x_y_patch
  use cable_netcdf_decomp_util_mod, only: io_decomp_land_to_land
  use cable_netcdf_decomp_util_mod, only: io_decomp_patch_to_land_patch
  use cable_netcdf_decomp_util_mod, only: io_decomp_patch_to_patch

  use cable_output_types_mod, only: cable_output_variable_t
  use cable_output_types_mod, only: cable_output_profile_t

  implicit none
  private

  public :: allocate_decompositions
  public :: deallocate_decompositions
  public :: associate_decomp_int32
  public :: associate_decomp_real32
  public :: associate_decomp_real64

  type :: cable_output_decomp_t
    class(cable_netcdf_decomp_t), allocatable :: land_int32
    class(cable_netcdf_decomp_t), allocatable :: land_real32
    class(cable_netcdf_decomp_t), allocatable :: land_real64
    class(cable_netcdf_decomp_t), allocatable :: land_soil_int32
    class(cable_netcdf_decomp_t), allocatable :: land_soil_real32
    class(cable_netcdf_decomp_t), allocatable :: land_soil_real64
    class(cable_netcdf_decomp_t), allocatable :: land_snow_int32
    class(cable_netcdf_decomp_t), allocatable :: land_snow_real32
    class(cable_netcdf_decomp_t), allocatable :: land_snow_real64
    class(cable_netcdf_decomp_t), allocatable :: land_rad_int32
    class(cable_netcdf_decomp_t), allocatable :: land_rad_real32
    class(cable_netcdf_decomp_t), allocatable :: land_rad_real64
    class(cable_netcdf_decomp_t), allocatable :: land_plantcarbon_int32
    class(cable_netcdf_decomp_t), allocatable :: land_plantcarbon_real32
    class(cable_netcdf_decomp_t), allocatable :: land_plantcarbon_real64
    class(cable_netcdf_decomp_t), allocatable :: land_soilcarbon_int32
    class(cable_netcdf_decomp_t), allocatable :: land_soilcarbon_real32
    class(cable_netcdf_decomp_t), allocatable :: land_soilcarbon_real64
    class(cable_netcdf_decomp_t), allocatable :: patch_int32
    class(cable_netcdf_decomp_t), allocatable :: patch_real32
    class(cable_netcdf_decomp_t), allocatable :: patch_real64
    class(cable_netcdf_decomp_t), allocatable :: patch_soil_int32
    class(cable_netcdf_decomp_t), allocatable :: patch_soil_real32
    class(cable_netcdf_decomp_t), allocatable :: patch_soil_real64
    class(cable_netcdf_decomp_t), allocatable :: patch_snow_int32
    class(cable_netcdf_decomp_t), allocatable :: patch_snow_real32
    class(cable_netcdf_decomp_t), allocatable :: patch_snow_real64
    class(cable_netcdf_decomp_t), allocatable :: patch_rad_int32
    class(cable_netcdf_decomp_t), allocatable :: patch_rad_real32
    class(cable_netcdf_decomp_t), allocatable :: patch_rad_real64
    class(cable_netcdf_decomp_t), allocatable :: patch_plantcarbon_int32
    class(cable_netcdf_decomp_t), allocatable :: patch_plantcarbon_real32
    class(cable_netcdf_decomp_t), allocatable :: patch_plantcarbon_real64
    class(cable_netcdf_decomp_t), allocatable :: patch_soilcarbon_int32
    class(cable_netcdf_decomp_t), allocatable :: patch_soilcarbon_real32
    class(cable_netcdf_decomp_t), allocatable :: patch_soilcarbon_real64
  end type

  type(cable_output_decomp_t), target :: decomp_grid_x_y
  type(cable_output_decomp_t), target :: decomp_grid_land
  type(cable_output_decomp_t), target :: decomp_grid_restart

contains

  subroutine allocate_decompositions()

    decomp_grid_x_y%land_int32               = io_decomp_land_to_x_y(land_x, land_y, [mland],      [xdimsize, ydimsize],      CABLE_NETCDF_INT)
    decomp_grid_x_y%land_real32              = io_decomp_land_to_x_y(land_x, land_y, [mland],      [xdimsize, ydimsize],      CABLE_NETCDF_FLOAT)
    decomp_grid_x_y%land_real64              = io_decomp_land_to_x_y(land_x, land_y, [mland],      [xdimsize, ydimsize],      CABLE_NETCDF_DOUBLE)
    decomp_grid_x_y%land_soil_int32          = io_decomp_land_to_x_y(land_x, land_y, [mland, ms],  [xdimsize, ydimsize, ms],  CABLE_NETCDF_INT)
    decomp_grid_x_y%land_soil_real32         = io_decomp_land_to_x_y(land_x, land_y, [mland, ms],  [xdimsize, ydimsize, ms],  CABLE_NETCDF_FLOAT)
    decomp_grid_x_y%land_soil_real64         = io_decomp_land_to_x_y(land_x, land_y, [mland, ms],  [xdimsize, ydimsize, ms],  CABLE_NETCDF_DOUBLE)
    decomp_grid_x_y%land_snow_int32          = io_decomp_land_to_x_y(land_x, land_y, [mland, msn], [xdimsize, ydimsize, msn], CABLE_NETCDF_INT)
    decomp_grid_x_y%land_snow_real32         = io_decomp_land_to_x_y(land_x, land_y, [mland, msn], [xdimsize, ydimsize, msn], CABLE_NETCDF_FLOAT)
    decomp_grid_x_y%land_snow_real64         = io_decomp_land_to_x_y(land_x, land_y, [mland, msn], [xdimsize, ydimsize, msn], CABLE_NETCDF_DOUBLE)
    decomp_grid_x_y%land_rad_int32           = io_decomp_land_to_x_y(land_x, land_y, [mland, nrb], [xdimsize, ydimsize, nrb], CABLE_NETCDF_INT)
    decomp_grid_x_y%land_rad_real32          = io_decomp_land_to_x_y(land_x, land_y, [mland, nrb], [xdimsize, ydimsize, nrb], CABLE_NETCDF_FLOAT)
    decomp_grid_x_y%land_rad_real64          = io_decomp_land_to_x_y(land_x, land_y, [mland, nrb], [xdimsize, ydimsize, nrb], CABLE_NETCDF_DOUBLE)
    decomp_grid_x_y%land_plantcarbon_int32   = io_decomp_land_to_x_y(land_x, land_y, [mland, ncp], [xdimsize, ydimsize, ncp], CABLE_NETCDF_INT)
    decomp_grid_x_y%land_plantcarbon_real32  = io_decomp_land_to_x_y(land_x, land_y, [mland, ncp], [xdimsize, ydimsize, ncp], CABLE_NETCDF_FLOAT)
    decomp_grid_x_y%land_plantcarbon_real64  = io_decomp_land_to_x_y(land_x, land_y, [mland, ncp], [xdimsize, ydimsize, ncp], CABLE_NETCDF_DOUBLE)
    decomp_grid_x_y%land_soilcarbon_int32    = io_decomp_land_to_x_y(land_x, land_y, [mland, ncs], [xdimsize, ydimsize, ncs], CABLE_NETCDF_INT)
    decomp_grid_x_y%land_soilcarbon_real32   = io_decomp_land_to_x_y(land_x, land_y, [mland, ncs], [xdimsize, ydimsize, ncs], CABLE_NETCDF_FLOAT)
    decomp_grid_x_y%land_soilcarbon_real64   = io_decomp_land_to_x_y(land_x, land_y, [mland, ncs], [xdimsize, ydimsize, ncs], CABLE_NETCDF_DOUBLE)
    decomp_grid_x_y%patch_int32              = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, [mp],      [xdimsize, ydimsize, max_vegpatches],      CABLE_NETCDF_INT)
    decomp_grid_x_y%patch_real32             = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, [mp],      [xdimsize, ydimsize, max_vegpatches],      CABLE_NETCDF_FLOAT)
    decomp_grid_x_y%patch_real64             = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, [mp],      [xdimsize, ydimsize, max_vegpatches],      CABLE_NETCDF_DOUBLE)
    decomp_grid_x_y%patch_soil_int32         = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, [mp, ms],  [xdimsize, ydimsize, max_vegpatches, ms],  CABLE_NETCDF_INT)
    decomp_grid_x_y%patch_soil_real32        = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, [mp, ms],  [xdimsize, ydimsize, max_vegpatches, ms],  CABLE_NETCDF_FLOAT)
    decomp_grid_x_y%patch_soil_real64        = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, [mp, ms],  [xdimsize, ydimsize, max_vegpatches, ms],  CABLE_NETCDF_DOUBLE)
    decomp_grid_x_y%patch_snow_int32         = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, [mp, msn], [xdimsize, ydimsize, max_vegpatches, msn], CABLE_NETCDF_INT)
    decomp_grid_x_y%patch_snow_real32        = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, [mp, msn], [xdimsize, ydimsize, max_vegpatches, msn], CABLE_NETCDF_FLOAT)
    decomp_grid_x_y%patch_snow_real64        = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, [mp, msn], [xdimsize, ydimsize, max_vegpatches, msn], CABLE_NETCDF_DOUBLE)
    decomp_grid_x_y%patch_rad_int32          = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, [mp, nrb], [xdimsize, ydimsize, max_vegpatches, nrb], CABLE_NETCDF_INT)
    decomp_grid_x_y%patch_rad_real32         = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, [mp, nrb], [xdimsize, ydimsize, max_vegpatches, nrb], CABLE_NETCDF_FLOAT)
    decomp_grid_x_y%patch_rad_real64         = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, [mp, nrb], [xdimsize, ydimsize, max_vegpatches, nrb], CABLE_NETCDF_DOUBLE)
    decomp_grid_x_y%patch_plantcarbon_int32  = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, [mp, ncp], [xdimsize, ydimsize, max_vegpatches, ncp], CABLE_NETCDF_INT)
    decomp_grid_x_y%patch_plantcarbon_real32 = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, [mp, ncp], [xdimsize, ydimsize, max_vegpatches, ncp], CABLE_NETCDF_FLOAT)
    decomp_grid_x_y%patch_plantcarbon_real64 = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, [mp, ncp], [xdimsize, ydimsize, max_vegpatches, ncp], CABLE_NETCDF_DOUBLE)
    decomp_grid_x_y%patch_soilcarbon_int32   = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, [mp, ncs], [xdimsize, ydimsize, max_vegpatches, ncs], CABLE_NETCDF_INT)
    decomp_grid_x_y%patch_soilcarbon_real32  = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, [mp, ncs], [xdimsize, ydimsize, max_vegpatches, ncs], CABLE_NETCDF_FLOAT)
    decomp_grid_x_y%patch_soilcarbon_real64  = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, [mp, ncs], [xdimsize, ydimsize, max_vegpatches, ncs], CABLE_NETCDF_DOUBLE)

    decomp_grid_land%land_int32               = io_decomp_land_to_land(land_decomp_start, [mland],      [mland_global],      CABLE_NETCDF_INT)
    decomp_grid_land%land_real32              = io_decomp_land_to_land(land_decomp_start, [mland],      [mland_global],      CABLE_NETCDF_FLOAT)
    decomp_grid_land%land_real64              = io_decomp_land_to_land(land_decomp_start, [mland],      [mland_global],      CABLE_NETCDF_DOUBLE)
    decomp_grid_land%land_soil_int32          = io_decomp_land_to_land(land_decomp_start, [mland, ms],  [mland_global, ms],  CABLE_NETCDF_INT)
    decomp_grid_land%land_soil_real32         = io_decomp_land_to_land(land_decomp_start, [mland, ms],  [mland_global, ms],  CABLE_NETCDF_FLOAT)
    decomp_grid_land%land_soil_real64         = io_decomp_land_to_land(land_decomp_start, [mland, ms],  [mland_global, ms],  CABLE_NETCDF_DOUBLE)
    decomp_grid_land%land_snow_int32          = io_decomp_land_to_land(land_decomp_start, [mland, msn], [mland_global, msn], CABLE_NETCDF_INT)
    decomp_grid_land%land_snow_real32         = io_decomp_land_to_land(land_decomp_start, [mland, msn], [mland_global, msn], CABLE_NETCDF_FLOAT)
    decomp_grid_land%land_snow_real64         = io_decomp_land_to_land(land_decomp_start, [mland, msn], [mland_global, msn], CABLE_NETCDF_DOUBLE)
    decomp_grid_land%land_rad_int32           = io_decomp_land_to_land(land_decomp_start, [mland, nrb], [mland_global, nrb], CABLE_NETCDF_INT)
    decomp_grid_land%land_rad_real32          = io_decomp_land_to_land(land_decomp_start, [mland, nrb], [mland_global, nrb], CABLE_NETCDF_FLOAT)
    decomp_grid_land%land_rad_real64          = io_decomp_land_to_land(land_decomp_start, [mland, nrb], [mland_global, nrb], CABLE_NETCDF_DOUBLE)
    decomp_grid_land%land_plantcarbon_int32   = io_decomp_land_to_land(land_decomp_start, [mland, ncp], [mland_global, ncp], CABLE_NETCDF_INT)
    decomp_grid_land%land_plantcarbon_real32  = io_decomp_land_to_land(land_decomp_start, [mland, ncp], [mland_global, ncp], CABLE_NETCDF_FLOAT)
    decomp_grid_land%land_plantcarbon_real64  = io_decomp_land_to_land(land_decomp_start, [mland, ncp], [mland_global, ncp], CABLE_NETCDF_DOUBLE)
    decomp_grid_land%land_soilcarbon_int32    = io_decomp_land_to_land(land_decomp_start, [mland, ncs], [mland_global, ncs], CABLE_NETCDF_INT)
    decomp_grid_land%land_soilcarbon_real32   = io_decomp_land_to_land(land_decomp_start, [mland, ncs], [mland_global, ncs], CABLE_NETCDF_FLOAT)
    decomp_grid_land%land_soilcarbon_real64   = io_decomp_land_to_land(land_decomp_start, [mland, ncs], [mland_global, ncs], CABLE_NETCDF_DOUBLE)
    decomp_grid_land%patch_int32              = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, [mp],      [mland_global, max_vegpatches],      CABLE_NETCDF_INT)
    decomp_grid_land%patch_real32             = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, [mp],      [mland_global, max_vegpatches],      CABLE_NETCDF_FLOAT)
    decomp_grid_land%patch_real64             = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, [mp],      [mland_global, max_vegpatches],      CABLE_NETCDF_DOUBLE)
    decomp_grid_land%patch_soil_int32         = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, [mp, ms],  [mland_global, max_vegpatches, ms],  CABLE_NETCDF_INT)
    decomp_grid_land%patch_soil_real32        = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, [mp, ms],  [mland_global, max_vegpatches, ms],  CABLE_NETCDF_FLOAT)
    decomp_grid_land%patch_soil_real64        = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, [mp, ms],  [mland_global, max_vegpatches, ms],  CABLE_NETCDF_DOUBLE)
    decomp_grid_land%patch_snow_int32         = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, [mp, msn], [mland_global, max_vegpatches, msn], CABLE_NETCDF_INT)
    decomp_grid_land%patch_snow_real32        = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, [mp, msn], [mland_global, max_vegpatches, msn], CABLE_NETCDF_FLOAT)
    decomp_grid_land%patch_snow_real64        = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, [mp, msn], [mland_global, max_vegpatches, msn], CABLE_NETCDF_DOUBLE)
    decomp_grid_land%patch_rad_int32          = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, [mp, nrb], [mland_global, max_vegpatches, nrb], CABLE_NETCDF_INT)
    decomp_grid_land%patch_rad_real32         = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, [mp, nrb], [mland_global, max_vegpatches, nrb], CABLE_NETCDF_FLOAT)
    decomp_grid_land%patch_rad_real64         = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, [mp, nrb], [mland_global, max_vegpatches, nrb], CABLE_NETCDF_DOUBLE)
    decomp_grid_land%patch_plantcarbon_int32  = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, [mp, ncp], [mland_global, max_vegpatches, ncp], CABLE_NETCDF_INT)
    decomp_grid_land%patch_plantcarbon_real32 = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, [mp, ncp], [mland_global, max_vegpatches, ncp], CABLE_NETCDF_FLOAT)
    decomp_grid_land%patch_plantcarbon_real64 = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, [mp, ncp], [mland_global, max_vegpatches, ncp], CABLE_NETCDF_DOUBLE)
    decomp_grid_land%patch_soilcarbon_int32   = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, [mp, ncs], [mland_global, max_vegpatches, ncs], CABLE_NETCDF_INT)
    decomp_grid_land%patch_soilcarbon_real32  = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, [mp, ncs], [mland_global, max_vegpatches, ncs], CABLE_NETCDF_FLOAT)
    decomp_grid_land%patch_soilcarbon_real64  = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, [mp, ncs], [mland_global, max_vegpatches, ncs], CABLE_NETCDF_DOUBLE)

    decomp_grid_restart%patch_int32              = io_decomp_patch_to_patch(patch_decomp_start, [mp],      [mp_global],      CABLE_NETCDF_INT)
    decomp_grid_restart%patch_real32             = io_decomp_patch_to_patch(patch_decomp_start, [mp],      [mp_global],      CABLE_NETCDF_FLOAT)
    decomp_grid_restart%patch_real64             = io_decomp_patch_to_patch(patch_decomp_start, [mp],      [mp_global],      CABLE_NETCDF_DOUBLE)
    decomp_grid_restart%patch_soil_int32         = io_decomp_patch_to_patch(patch_decomp_start, [mp, ms],  [mp_global, ms],  CABLE_NETCDF_INT)
    decomp_grid_restart%patch_soil_real32        = io_decomp_patch_to_patch(patch_decomp_start, [mp, ms],  [mp_global, ms],  CABLE_NETCDF_FLOAT)
    decomp_grid_restart%patch_soil_real64        = io_decomp_patch_to_patch(patch_decomp_start, [mp, ms],  [mp_global, ms],  CABLE_NETCDF_DOUBLE)
    decomp_grid_restart%patch_snow_int32         = io_decomp_patch_to_patch(patch_decomp_start, [mp, msn], [mp_global, msn], CABLE_NETCDF_INT)
    decomp_grid_restart%patch_snow_real32        = io_decomp_patch_to_patch(patch_decomp_start, [mp, msn], [mp_global, msn], CABLE_NETCDF_FLOAT)
    decomp_grid_restart%patch_snow_real64        = io_decomp_patch_to_patch(patch_decomp_start, [mp, msn], [mp_global, msn], CABLE_NETCDF_DOUBLE)
    decomp_grid_restart%patch_rad_int32          = io_decomp_patch_to_patch(patch_decomp_start, [mp, nrb], [mp_global, nrb], CABLE_NETCDF_INT)
    decomp_grid_restart%patch_rad_real32         = io_decomp_patch_to_patch(patch_decomp_start, [mp, nrb], [mp_global, nrb], CABLE_NETCDF_FLOAT)
    decomp_grid_restart%patch_rad_real64         = io_decomp_patch_to_patch(patch_decomp_start, [mp, nrb], [mp_global, nrb], CABLE_NETCDF_DOUBLE)
    decomp_grid_restart%patch_plantcarbon_int32  = io_decomp_patch_to_patch(patch_decomp_start, [mp, ncp], [mp_global, ncp], CABLE_NETCDF_INT)
    decomp_grid_restart%patch_plantcarbon_real32 = io_decomp_patch_to_patch(patch_decomp_start, [mp, ncp], [mp_global, ncp], CABLE_NETCDF_FLOAT)
    decomp_grid_restart%patch_plantcarbon_real64 = io_decomp_patch_to_patch(patch_decomp_start, [mp, ncp], [mp_global, ncp], CABLE_NETCDF_DOUBLE)
    decomp_grid_restart%patch_soilcarbon_int32   = io_decomp_patch_to_patch(patch_decomp_start, [mp, ncs], [mp_global, ncs], CABLE_NETCDF_INT)
    decomp_grid_restart%patch_soilcarbon_real32  = io_decomp_patch_to_patch(patch_decomp_start, [mp, ncs], [mp_global, ncs], CABLE_NETCDF_FLOAT)
    decomp_grid_restart%patch_soilcarbon_real64  = io_decomp_patch_to_patch(patch_decomp_start, [mp, ncs], [mp_global, ncs], CABLE_NETCDF_DOUBLE)

  end subroutine

  subroutine deallocate_decompositions()

    deallocate(decomp_grid_x_y%land_int32)
    deallocate(decomp_grid_x_y%land_real32)
    deallocate(decomp_grid_x_y%land_real64)
    deallocate(decomp_grid_x_y%land_soil_int32)
    deallocate(decomp_grid_x_y%land_soil_real32)
    deallocate(decomp_grid_x_y%land_soil_real64)
    deallocate(decomp_grid_x_y%land_snow_int32)
    deallocate(decomp_grid_x_y%land_snow_real32)
    deallocate(decomp_grid_x_y%land_snow_real64)
    deallocate(decomp_grid_x_y%land_rad_int32)
    deallocate(decomp_grid_x_y%land_rad_real32)
    deallocate(decomp_grid_x_y%land_rad_real64)
    deallocate(decomp_grid_x_y%land_plantcarbon_int32)
    deallocate(decomp_grid_x_y%land_plantcarbon_real32)
    deallocate(decomp_grid_x_y%land_plantcarbon_real64)
    deallocate(decomp_grid_x_y%land_soilcarbon_int32)
    deallocate(decomp_grid_x_y%land_soilcarbon_real32)
    deallocate(decomp_grid_x_y%land_soilcarbon_real64)
    deallocate(decomp_grid_x_y%patch_int32)
    deallocate(decomp_grid_x_y%patch_real32)
    deallocate(decomp_grid_x_y%patch_real64)
    deallocate(decomp_grid_x_y%patch_soil_int32)
    deallocate(decomp_grid_x_y%patch_soil_real32)
    deallocate(decomp_grid_x_y%patch_soil_real64)
    deallocate(decomp_grid_x_y%patch_snow_int32)
    deallocate(decomp_grid_x_y%patch_snow_real32)
    deallocate(decomp_grid_x_y%patch_snow_real64)
    deallocate(decomp_grid_x_y%patch_rad_int32)
    deallocate(decomp_grid_x_y%patch_rad_real32)
    deallocate(decomp_grid_x_y%patch_rad_real64)
    deallocate(decomp_grid_x_y%patch_plantcarbon_int32)
    deallocate(decomp_grid_x_y%patch_plantcarbon_real32)
    deallocate(decomp_grid_x_y%patch_plantcarbon_real64)
    deallocate(decomp_grid_x_y%patch_soilcarbon_int32)
    deallocate(decomp_grid_x_y%patch_soilcarbon_real32)
    deallocate(decomp_grid_x_y%patch_soilcarbon_real64)

    deallocate(decomp_grid_land%land_int32)
    deallocate(decomp_grid_land%land_real32)
    deallocate(decomp_grid_land%land_real64)
    deallocate(decomp_grid_land%land_soil_int32)
    deallocate(decomp_grid_land%land_soil_real32)
    deallocate(decomp_grid_land%land_soil_real64)
    deallocate(decomp_grid_land%land_snow_int32)
    deallocate(decomp_grid_land%land_snow_real32)
    deallocate(decomp_grid_land%land_snow_real64)
    deallocate(decomp_grid_land%land_rad_int32)
    deallocate(decomp_grid_land%land_rad_real32)
    deallocate(decomp_grid_land%land_rad_real64)
    deallocate(decomp_grid_land%land_plantcarbon_int32)
    deallocate(decomp_grid_land%land_plantcarbon_real32)
    deallocate(decomp_grid_land%land_plantcarbon_real64)
    deallocate(decomp_grid_land%land_soilcarbon_int32)
    deallocate(decomp_grid_land%land_soilcarbon_real32)
    deallocate(decomp_grid_land%land_soilcarbon_real64)
    deallocate(decomp_grid_land%patch_int32)
    deallocate(decomp_grid_land%patch_real32)
    deallocate(decomp_grid_land%patch_real64)
    deallocate(decomp_grid_land%patch_soil_int32)
    deallocate(decomp_grid_land%patch_soil_real32)
    deallocate(decomp_grid_land%patch_soil_real64)
    deallocate(decomp_grid_land%patch_snow_int32)
    deallocate(decomp_grid_land%patch_snow_real32)
    deallocate(decomp_grid_land%patch_snow_real64)
    deallocate(decomp_grid_land%patch_rad_int32)
    deallocate(decomp_grid_land%patch_rad_real32)
    deallocate(decomp_grid_land%patch_rad_real64)
    deallocate(decomp_grid_land%patch_plantcarbon_int32)
    deallocate(decomp_grid_land%patch_plantcarbon_real32)
    deallocate(decomp_grid_land%patch_plantcarbon_real64)
    deallocate(decomp_grid_land%patch_soilcarbon_int32)
    deallocate(decomp_grid_land%patch_soilcarbon_real32)
    deallocate(decomp_grid_land%patch_soilcarbon_real64)

    deallocate(decomp_grid_restart%patch_int32)
    deallocate(decomp_grid_restart%patch_real32)
    deallocate(decomp_grid_restart%patch_real64)
    deallocate(decomp_grid_restart%patch_soil_int32)
    deallocate(decomp_grid_restart%patch_soil_real32)
    deallocate(decomp_grid_restart%patch_soil_real64)
    deallocate(decomp_grid_restart%patch_snow_int32)
    deallocate(decomp_grid_restart%patch_snow_real32)
    deallocate(decomp_grid_restart%patch_snow_real64)
    deallocate(decomp_grid_restart%patch_rad_int32)
    deallocate(decomp_grid_restart%patch_rad_real32)
    deallocate(decomp_grid_restart%patch_rad_real64)
    deallocate(decomp_grid_restart%patch_plantcarbon_int32)
    deallocate(decomp_grid_restart%patch_plantcarbon_real32)
    deallocate(decomp_grid_restart%patch_plantcarbon_real64)
    deallocate(decomp_grid_restart%patch_soilcarbon_int32)
    deallocate(decomp_grid_restart%patch_soilcarbon_real32)
    deallocate(decomp_grid_restart%patch_soilcarbon_real64)

  end subroutine

  subroutine associate_decomp_int32(output_profile, output_var, decomp)
    type(cable_output_profile_t), intent(in) :: output_profile
    type(cable_output_variable_t), intent(in) :: output_var
    class(cable_netcdf_decomp_t), pointer, intent(inout) :: decomp
    type(cable_output_decomp_t), pointer :: output_decomp

    select case (output_profile%grid_type)
    case ("restart")
      call associate_decomp_restart_int32(output_var, decomp)
      return
    case ("mask")
      output_decomp => decomp_grid_x_y
    case ("land")
      output_decomp => decomp_grid_land
    case default
      call cable_abort("Unsupported grid type for output profile " // output_profile%file_name, __FILE__, __LINE__)
    end select

    if (array_eq(output_var%data_shape(:)%size, [mp])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_int32
      else
        decomp => output_decomp%land_int32
      end if
    else if (array_eq(output_var%data_shape(:)%size, [mp, ms])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_soil_int32
      else
        decomp => output_decomp%land_soil_int32
      end if
    else if (array_eq(output_var%data_shape(:)%size, [mp, msn])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_snow_int32
      else
        decomp => output_decomp%land_snow_int32
      end if
    else if (array_eq(output_var%data_shape(:)%size, [mp, nrb])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_rad_int32
      else
        decomp => output_decomp%land_rad_int32
      end if
    else if (array_eq(output_var%data_shape(:)%size, [mp, ncp])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_plantcarbon_int32
      else
        decomp => output_decomp%land_plantcarbon_int32
      end if
    else if (array_eq(output_var%data_shape(:)%size, [mp, ncs])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_soilcarbon_int32
      else
        decomp => output_decomp%land_soilcarbon_int32
      end if
    else
      call cable_abort("Unsupported data shape for output variable " // output_var%field_name, __FILE__, __LINE__)
    end if

  end subroutine associate_decomp_int32

  subroutine associate_decomp_real32(output_profile, output_var, decomp)
    type(cable_output_profile_t), intent(in) :: output_profile
    type(cable_output_variable_t), intent(in) :: output_var
    class(cable_netcdf_decomp_t), pointer, intent(inout) :: decomp
    type(cable_output_decomp_t), pointer :: output_decomp

    select case (output_profile%grid_type)
    case ("restart")
      call associate_decomp_restart_real32(output_var, decomp)
      return
    case ("mask")
      output_decomp => decomp_grid_x_y
    case ("land")
      output_decomp => decomp_grid_land
    case default
      call cable_abort("Unsupported grid type for output profile " // output_profile%file_name, __FILE__, __LINE__)
    end select

    if (array_eq(output_var%data_shape(:)%size, [mp])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_real32
      else
        decomp => output_decomp%land_real32
      end if
    else if (array_eq(output_var%data_shape(:)%size, [mp, ms])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_soil_real32
      else
        decomp => output_decomp%land_soil_real32
      end if
    else if (array_eq(output_var%data_shape(:)%size, [mp, msn])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_snow_real32
      else
        decomp => output_decomp%land_snow_real32
      end if
    else if (array_eq(output_var%data_shape(:)%size, [mp, nrb])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_rad_real32
      else
        decomp => output_decomp%land_rad_real32
      end if
    else if (array_eq(output_var%data_shape(:)%size, [mp, ncp])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_plantcarbon_real32
      else
        decomp => output_decomp%land_plantcarbon_real32
      end if
    else if (array_eq(output_var%data_shape(:)%size, [mp, ncs])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_soilcarbon_real32
      else
        decomp => output_decomp%land_soilcarbon_real32
      end if
    else
      call cable_abort("Unsupported data shape for output variable " // output_var%field_name, __FILE__, __LINE__)
    end if

  end subroutine associate_decomp_real32

  subroutine associate_decomp_real64(output_profile, output_var, decomp)
    type(cable_output_profile_t), intent(in) :: output_profile
    type(cable_output_variable_t), intent(in) :: output_var
    class(cable_netcdf_decomp_t), pointer, intent(inout) :: decomp
    type(cable_output_decomp_t), pointer :: output_decomp

    select case (output_profile%grid_type)
    case ("restart")
      call associate_decomp_restart_real64(output_var, decomp)
      return
    case ("mask")
      output_decomp => decomp_grid_x_y
    case ("land")
      output_decomp => decomp_grid_land
    case default
      call cable_abort("Unsupported grid type for output profile " // output_profile%file_name, __FILE__, __LINE__)
    end select

    if (array_eq(output_var%data_shape(:)%size, [mp])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_real64
      else
        decomp => output_decomp%land_real64
      end if
    else if (array_eq(output_var%data_shape(:)%size, [mp, ms])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_soil_real64
      else
        decomp => output_decomp%land_soil_real64
      end if
    else if (array_eq(output_var%data_shape(:)%size, [mp, msn])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_snow_real64
      else
        decomp => output_decomp%land_snow_real64
      end if
    else if (array_eq(output_var%data_shape(:)%size, [mp, nrb])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_rad_real64
      else
        decomp => output_decomp%land_rad_real64
      end if
    else if (array_eq(output_var%data_shape(:)%size, [mp, ncp])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_plantcarbon_real64
      else
        decomp => output_decomp%land_plantcarbon_real64
      end if
    else if (array_eq(output_var%data_shape(:)%size, [mp, ncs])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_soilcarbon_real64
      else
        decomp => output_decomp%land_soilcarbon_real64
      end if
    else
      call cable_abort("Unsupported data shape for output variable " // output_var%field_name, __FILE__, __LINE__)
    end if

  end subroutine associate_decomp_real64

  subroutine associate_decomp_restart_int32(output_var, decomp)
    type(cable_output_variable_t), intent(in) :: output_var
    class(cable_netcdf_decomp_t), pointer, intent(inout) :: decomp
    type(cable_output_decomp_t), pointer :: output_decomp

    if (array_eq(output_var%data_shape(:)%size, [mp])) then
      decomp => decomp_grid_restart%patch_int32
    else if (array_eq(output_var%data_shape(:)%size, [mp, ms])) then
      decomp => decomp_grid_restart%patch_soil_int32
    else if (array_eq(output_var%data_shape(:)%size, [mp, msn])) then
      decomp => decomp_grid_restart%patch_snow_int32
    else if (array_eq(output_var%data_shape(:)%size, [mp, nrb])) then
      decomp => decomp_grid_restart%patch_rad_int32
    else if (array_eq(output_var%data_shape(:)%size, [mp, ncp])) then
      decomp => decomp_grid_restart%patch_plantcarbon_int32
    else if (array_eq(output_var%data_shape(:)%size, [mp, ncs])) then
      decomp => decomp_grid_restart%patch_soilcarbon_int32
    else
      call cable_abort("Unsupported data shape for output variable " // output_var%field_name, __FILE__, __LINE__)
    end if

  end subroutine associate_decomp_restart_int32

  subroutine associate_decomp_restart_real32(output_var, decomp)
    type(cable_output_variable_t), intent(in) :: output_var
    class(cable_netcdf_decomp_t), pointer, intent(inout) :: decomp

    if (array_eq(output_var%data_shape(:)%size, [mp])) then
      decomp => decomp_grid_restart%patch_real32
    else if (array_eq(output_var%data_shape(:)%size, [mp, ms])) then
      decomp => decomp_grid_restart%patch_soil_real32
    else if (array_eq(output_var%data_shape(:)%size, [mp, msn])) then
      decomp => decomp_grid_restart%patch_snow_real32
    else if (array_eq(output_var%data_shape(:)%size, [mp, nrb])) then
      decomp => decomp_grid_restart%patch_rad_real32
    else if (array_eq(output_var%data_shape(:)%size, [mp, ncp])) then
      decomp => decomp_grid_restart%patch_plantcarbon_real32
    else if (array_eq(output_var%data_shape(:)%size, [mp, ncs])) then
      decomp => decomp_grid_restart%patch_soilcarbon_real32
    else
      call cable_abort("Unsupported data shape for output variable " // output_var%field_name, __FILE__, __LINE__)
    end if

  end subroutine associate_decomp_restart_real32

  subroutine associate_decomp_restart_real64(output_var, decomp)
    type(cable_output_variable_t), intent(in) :: output_var
    class(cable_netcdf_decomp_t), pointer, intent(inout) :: decomp

    if (array_eq(output_var%data_shape(:)%size, [mp])) then
      decomp => decomp_grid_restart%patch_real64
    else if (array_eq(output_var%data_shape(:)%size, [mp, ms])) then
      decomp => decomp_grid_restart%patch_soil_real64
    else if (array_eq(output_var%data_shape(:)%size, [mp, msn])) then
      decomp => decomp_grid_restart%patch_snow_real64
    else if (array_eq(output_var%data_shape(:)%size, [mp, nrb])) then
      decomp => decomp_grid_restart%patch_rad_real64
    else if (array_eq(output_var%data_shape(:)%size, [mp, ncp])) then
      decomp => decomp_grid_restart%patch_plantcarbon_real64
    else if (array_eq(output_var%data_shape(:)%size, [mp, ncs])) then
      decomp => decomp_grid_restart%patch_soilcarbon_real64
    else
      call cable_abort("Unsupported data shape for output variable " // output_var%field_name, __FILE__, __LINE__)
    end if

  end subroutine associate_decomp_restart_real64

end module
