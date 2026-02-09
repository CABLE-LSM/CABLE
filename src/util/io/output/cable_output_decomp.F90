module cable_output_decomp_mod

  use cable_abort_module, only: cable_abort

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

  use cable_netcdf_decomp_util_mod, only: dim_spec_t
  use cable_netcdf_decomp_util_mod, only: io_decomp_land_to_x_y
  use cable_netcdf_decomp_util_mod, only: io_decomp_patch_to_x_y_patch
  use cable_netcdf_decomp_util_mod, only: io_decomp_land_to_land
  use cable_netcdf_decomp_util_mod, only: io_decomp_patch_to_land_patch
  use cable_netcdf_decomp_util_mod, only: io_decomp_patch_to_patch

  use cable_output_types_mod, only: cable_output_variable_t
  use cable_output_types_mod, only: cable_output_profile_t
  use cable_output_types_mod, only: CABLE_OUTPUT_DIM_PATCH
  use cable_output_types_mod, only: CABLE_OUTPUT_DIM_SOIL
  use cable_output_types_mod, only: CABLE_OUTPUT_DIM_SNOW
  use cable_output_types_mod, only: CABLE_OUTPUT_DIM_RAD
  use cable_output_types_mod, only: CABLE_OUTPUT_DIM_PLANTCARBON
  use cable_output_types_mod, only: CABLE_OUTPUT_DIM_SOILCARBON

  use cable_output_utils_mod, only: data_shape_eq

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

    type(dim_spec_t), allocatable :: mem_shape_land(:)
    type(dim_spec_t), allocatable :: mem_shape_land_soil(:)
    type(dim_spec_t), allocatable :: mem_shape_land_snow(:)
    type(dim_spec_t), allocatable :: mem_shape_land_rad(:)
    type(dim_spec_t), allocatable :: mem_shape_land_plantcarbon(:)
    type(dim_spec_t), allocatable :: mem_shape_land_soilcarbon(:)
    type(dim_spec_t), allocatable :: mem_shape_patch(:)
    type(dim_spec_t), allocatable :: mem_shape_patch_soil(:)
    type(dim_spec_t), allocatable :: mem_shape_patch_snow(:)
    type(dim_spec_t), allocatable :: mem_shape_patch_rad(:)
    type(dim_spec_t), allocatable :: mem_shape_patch_plantcarbon(:)
    type(dim_spec_t), allocatable :: mem_shape_patch_soilcarbon(:)

    type(dim_spec_t), allocatable :: var_shape_x_y(:)
    type(dim_spec_t), allocatable :: var_shape_x_y_soil(:)
    type(dim_spec_t), allocatable :: var_shape_x_y_snow(:)
    type(dim_spec_t), allocatable :: var_shape_x_y_rad(:)
    type(dim_spec_t), allocatable :: var_shape_x_y_plantcarbon(:)
    type(dim_spec_t), allocatable :: var_shape_x_y_soilcarbon(:)
    type(dim_spec_t), allocatable :: var_shape_x_y_patch(:)
    type(dim_spec_t), allocatable :: var_shape_x_y_patch_soil(:)
    type(dim_spec_t), allocatable :: var_shape_x_y_patch_snow(:)
    type(dim_spec_t), allocatable :: var_shape_x_y_patch_rad(:)
    type(dim_spec_t), allocatable :: var_shape_x_y_patch_plantcarbon(:)
    type(dim_spec_t), allocatable :: var_shape_x_y_patch_soilcarbon(:)

    type(dim_spec_t), allocatable :: var_shape_land(:)
    type(dim_spec_t), allocatable :: var_shape_land_soil(:)
    type(dim_spec_t), allocatable :: var_shape_land_snow(:)
    type(dim_spec_t), allocatable :: var_shape_land_rad(:)
    type(dim_spec_t), allocatable :: var_shape_land_plantcarbon(:)
    type(dim_spec_t), allocatable :: var_shape_land_soilcarbon(:)
    type(dim_spec_t), allocatable :: var_shape_land_patch(:)
    type(dim_spec_t), allocatable :: var_shape_land_patch_soil(:)
    type(dim_spec_t), allocatable :: var_shape_land_patch_snow(:)
    type(dim_spec_t), allocatable :: var_shape_land_patch_rad(:)
    type(dim_spec_t), allocatable :: var_shape_land_patch_plantcarbon(:)
    type(dim_spec_t), allocatable :: var_shape_land_patch_soilcarbon(:)

    type(dim_spec_t), allocatable :: var_shape_patch(:)
    type(dim_spec_t), allocatable :: var_shape_patch_soil(:)
    type(dim_spec_t), allocatable :: var_shape_patch_snow(:)
    type(dim_spec_t), allocatable :: var_shape_patch_rad(:)
    type(dim_spec_t), allocatable :: var_shape_patch_plantcarbon(:)
    type(dim_spec_t), allocatable :: var_shape_patch_soilcarbon(:)

    mem_shape_land              = [dim_spec_t('land', mland)]
    mem_shape_land_soil         = [dim_spec_t('land', mland), dim_spec_t('soil', ms)]
    mem_shape_land_snow         = [dim_spec_t('land', mland), dim_spec_t('snow', msn)]
    mem_shape_land_rad          = [dim_spec_t('land', mland), dim_spec_t('rad', nrb)]
    mem_shape_land_plantcarbon  = [dim_spec_t('land', mland), dim_spec_t('plantcarbon', ncp)]
    mem_shape_land_soilcarbon   = [dim_spec_t('land', mland), dim_spec_t('soilcarbon', ncs)]
    mem_shape_patch             = [dim_spec_t('patch', mp)]
    mem_shape_patch_soil        = [dim_spec_t('patch', mp), dim_spec_t('soil', ms)]
    mem_shape_patch_snow        = [dim_spec_t('patch', mp), dim_spec_t('snow', msn)]
    mem_shape_patch_rad         = [dim_spec_t('patch', mp), dim_spec_t('rad', nrb)]
    mem_shape_patch_plantcarbon = [dim_spec_t('patch', mp), dim_spec_t('plantcarbon', ncp)]
    mem_shape_patch_soilcarbon  = [dim_spec_t('patch', mp), dim_spec_t('soilcarbon', ncs)]

    var_shape_x_y                    = [dim_spec_t('x', xdimsize), dim_spec_t('y', ydimsize)]
    var_shape_x_y_soil               = [dim_spec_t('x', xdimsize), dim_spec_t('y', ydimsize), dim_spec_t('soil', ms)]
    var_shape_x_y_snow               = [dim_spec_t('x', xdimsize), dim_spec_t('y', ydimsize), dim_spec_t('snow', msn)]
    var_shape_x_y_rad                = [dim_spec_t('x', xdimsize), dim_spec_t('y', ydimsize), dim_spec_t('rad', nrb)]
    var_shape_x_y_plantcarbon        = [dim_spec_t('x', xdimsize), dim_spec_t('y', ydimsize), dim_spec_t('plantcarbon', ncp)]
    var_shape_x_y_soilcarbon         = [dim_spec_t('x', xdimsize), dim_spec_t('y', ydimsize), dim_spec_t('soilcarbon', ncs)]
    var_shape_x_y_patch              = [dim_spec_t('x', xdimsize), dim_spec_t('y', ydimsize), dim_spec_t('patch', max_vegpatches)]
    var_shape_x_y_patch_soil         = [dim_spec_t('x', xdimsize), dim_spec_t('y', ydimsize), dim_spec_t('patch', max_vegpatches), dim_spec_t('soil', ms)]
    var_shape_x_y_patch_snow         = [dim_spec_t('x', xdimsize), dim_spec_t('y', ydimsize), dim_spec_t('patch', max_vegpatches), dim_spec_t('snow', msn)]
    var_shape_x_y_patch_rad          = [dim_spec_t('x', xdimsize), dim_spec_t('y', ydimsize), dim_spec_t('patch', max_vegpatches), dim_spec_t('rad', nrb)]
    var_shape_x_y_patch_plantcarbon  = [dim_spec_t('x', xdimsize), dim_spec_t('y', ydimsize), dim_spec_t('patch', max_vegpatches), dim_spec_t('plantcarbon', ncp)]
    var_shape_x_y_patch_soilcarbon   = [dim_spec_t('x', xdimsize), dim_spec_t('y', ydimsize), dim_spec_t('patch', max_vegpatches), dim_spec_t('soilcarbon', ncs)]

    var_shape_land                   = [dim_spec_t('land', mland_global)]
    var_shape_land_soil              = [dim_spec_t('land', mland_global), dim_spec_t('soil', ms)]
    var_shape_land_snow              = [dim_spec_t('land', mland_global), dim_spec_t('snow', msn)]
    var_shape_land_rad               = [dim_spec_t('land', mland_global), dim_spec_t('rad', nrb)]
    var_shape_land_plantcarbon       = [dim_spec_t('land', mland_global), dim_spec_t('plantcarbon', ncp)]
    var_shape_land_soilcarbon        = [dim_spec_t('land', mland_global), dim_spec_t('soilcarbon', ncs)]
    var_shape_land_patch             = [dim_spec_t('land', mland_global)]
    var_shape_land_patch_soil        = [dim_spec_t('land', mland_global), dim_spec_t('patch', max_vegpatches), dim_spec_t('soil', ms)]
    var_shape_land_patch_snow        = [dim_spec_t('land', mland_global), dim_spec_t('patch', max_vegpatches), dim_spec_t('snow', msn)]
    var_shape_land_patch_rad         = [dim_spec_t('land', mland_global), dim_spec_t('patch', max_vegpatches), dim_spec_t('rad', nrb)]
    var_shape_land_patch_plantcarbon = [dim_spec_t('land', mland_global), dim_spec_t('patch', max_vegpatches), dim_spec_t('plantcarbon', ncp)]
    var_shape_land_patch_soilcarbon  = [dim_spec_t('land', mland_global), dim_spec_t('patch', max_vegpatches), dim_spec_t('soilcarbon', ncs)]

    var_shape_patch                  = [dim_spec_t('patch', mp_global)]
    var_shape_patch_soil             = [dim_spec_t('patch', mp_global), dim_spec_t('soil', ms)]
    var_shape_patch_snow             = [dim_spec_t('patch', mp_global), dim_spec_t('snow', msn)]
    var_shape_patch_rad              = [dim_spec_t('patch', mp_global), dim_spec_t('rad', nrb)]
    var_shape_patch_plantcarbon      = [dim_spec_t('patch', mp_global), dim_spec_t('plantcarbon', ncp)]
    var_shape_patch_soilcarbon       = [dim_spec_t('patch', mp_global), dim_spec_t('soilcarbon', ncs)]

    decomp_grid_x_y%land_int32               = io_decomp_land_to_x_y(land_x, land_y, mem_shape_land, var_shape_x_y, CABLE_NETCDF_INT)
    decomp_grid_x_y%land_real32              = io_decomp_land_to_x_y(land_x, land_y, mem_shape_land, var_shape_x_y, CABLE_NETCDF_FLOAT)
    decomp_grid_x_y%land_real64              = io_decomp_land_to_x_y(land_x, land_y, mem_shape_land, var_shape_x_y, CABLE_NETCDF_DOUBLE)
    decomp_grid_x_y%land_soil_int32          = io_decomp_land_to_x_y(land_x, land_y, mem_shape_land_soil, var_shape_x_y_soil, CABLE_NETCDF_INT)
    decomp_grid_x_y%land_soil_real32         = io_decomp_land_to_x_y(land_x, land_y, mem_shape_land_soil, var_shape_x_y_soil, CABLE_NETCDF_FLOAT)
    decomp_grid_x_y%land_soil_real64         = io_decomp_land_to_x_y(land_x, land_y, mem_shape_land_soil, var_shape_x_y_soil, CABLE_NETCDF_DOUBLE)
    decomp_grid_x_y%land_snow_int32          = io_decomp_land_to_x_y(land_x, land_y, mem_shape_land_snow, var_shape_x_y_snow, CABLE_NETCDF_INT)
    decomp_grid_x_y%land_snow_real32         = io_decomp_land_to_x_y(land_x, land_y, mem_shape_land_snow, var_shape_x_y_snow, CABLE_NETCDF_FLOAT)
    decomp_grid_x_y%land_snow_real64         = io_decomp_land_to_x_y(land_x, land_y, mem_shape_land_snow, var_shape_x_y_snow, CABLE_NETCDF_DOUBLE)
    decomp_grid_x_y%land_rad_int32           = io_decomp_land_to_x_y(land_x, land_y, mem_shape_land_rad, var_shape_x_y_rad, CABLE_NETCDF_INT)
    decomp_grid_x_y%land_rad_real32          = io_decomp_land_to_x_y(land_x, land_y, mem_shape_land_rad, var_shape_x_y_rad, CABLE_NETCDF_FLOAT)
    decomp_grid_x_y%land_rad_real64          = io_decomp_land_to_x_y(land_x, land_y, mem_shape_land_rad, var_shape_x_y_rad, CABLE_NETCDF_DOUBLE)
    decomp_grid_x_y%land_plantcarbon_int32   = io_decomp_land_to_x_y(land_x, land_y, mem_shape_land_plantcarbon, var_shape_x_y_plantcarbon, CABLE_NETCDF_INT)
    decomp_grid_x_y%land_plantcarbon_real32  = io_decomp_land_to_x_y(land_x, land_y, mem_shape_land_plantcarbon, var_shape_x_y_plantcarbon, CABLE_NETCDF_FLOAT)
    decomp_grid_x_y%land_plantcarbon_real64  = io_decomp_land_to_x_y(land_x, land_y, mem_shape_land_plantcarbon, var_shape_x_y_plantcarbon, CABLE_NETCDF_DOUBLE)
    decomp_grid_x_y%land_soilcarbon_int32    = io_decomp_land_to_x_y(land_x, land_y, mem_shape_land_soilcarbon, var_shape_x_y_soilcarbon, CABLE_NETCDF_INT)
    decomp_grid_x_y%land_soilcarbon_real32   = io_decomp_land_to_x_y(land_x, land_y, mem_shape_land_soilcarbon, var_shape_x_y_soilcarbon, CABLE_NETCDF_FLOAT)
    decomp_grid_x_y%land_soilcarbon_real64   = io_decomp_land_to_x_y(land_x, land_y, mem_shape_land_soilcarbon, var_shape_x_y_soilcarbon, CABLE_NETCDF_DOUBLE)
    decomp_grid_x_y%patch_int32              = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch, var_shape_x_y_patch, CABLE_NETCDF_INT)
    decomp_grid_x_y%patch_real32             = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch, var_shape_x_y_patch, CABLE_NETCDF_FLOAT)
    decomp_grid_x_y%patch_real64             = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch, var_shape_x_y_patch, CABLE_NETCDF_DOUBLE)
    decomp_grid_x_y%patch_soil_int32         = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_soil, var_shape_x_y_patch_soil, CABLE_NETCDF_INT)
    decomp_grid_x_y%patch_soil_real32        = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_soil, var_shape_x_y_patch_soil, CABLE_NETCDF_FLOAT)
    decomp_grid_x_y%patch_soil_real64        = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_soil, var_shape_x_y_patch_soil, CABLE_NETCDF_DOUBLE)
    decomp_grid_x_y%patch_snow_int32         = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_snow, var_shape_x_y_patch_snow, CABLE_NETCDF_INT)
    decomp_grid_x_y%patch_snow_real32        = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_snow, var_shape_x_y_patch_snow, CABLE_NETCDF_FLOAT)
    decomp_grid_x_y%patch_snow_real64        = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_snow, var_shape_x_y_patch_snow, CABLE_NETCDF_DOUBLE)
    decomp_grid_x_y%patch_rad_int32          = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_rad, var_shape_x_y_patch_rad, CABLE_NETCDF_INT)
    decomp_grid_x_y%patch_rad_real32         = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_rad, var_shape_x_y_patch_rad, CABLE_NETCDF_FLOAT)
    decomp_grid_x_y%patch_rad_real64         = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_rad, var_shape_x_y_patch_rad, CABLE_NETCDF_DOUBLE)
    decomp_grid_x_y%patch_plantcarbon_int32  = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_plantcarbon, var_shape_x_y_patch_plantcarbon, CABLE_NETCDF_INT)
    decomp_grid_x_y%patch_plantcarbon_real32 = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_plantcarbon, var_shape_x_y_patch_plantcarbon, CABLE_NETCDF_FLOAT)
    decomp_grid_x_y%patch_plantcarbon_real64 = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_plantcarbon, var_shape_x_y_patch_plantcarbon, CABLE_NETCDF_DOUBLE)
    decomp_grid_x_y%patch_soilcarbon_int32   = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_soilcarbon, var_shape_x_y_patch_soilcarbon, CABLE_NETCDF_INT)
    decomp_grid_x_y%patch_soilcarbon_real32  = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_soilcarbon, var_shape_x_y_patch_soilcarbon, CABLE_NETCDF_FLOAT)
    decomp_grid_x_y%patch_soilcarbon_real64  = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_soilcarbon, var_shape_x_y_patch_soilcarbon, CABLE_NETCDF_DOUBLE)

    decomp_grid_land%land_int32               = io_decomp_land_to_land(land_decomp_start, mem_shape_land, var_shape_land, CABLE_NETCDF_INT)
    decomp_grid_land%land_real32              = io_decomp_land_to_land(land_decomp_start, mem_shape_land, var_shape_land, CABLE_NETCDF_FLOAT)
    decomp_grid_land%land_real64              = io_decomp_land_to_land(land_decomp_start, mem_shape_land, var_shape_land, CABLE_NETCDF_DOUBLE)
    decomp_grid_land%land_soil_int32          = io_decomp_land_to_land(land_decomp_start, mem_shape_land_soil, var_shape_land_soil, CABLE_NETCDF_INT)
    decomp_grid_land%land_soil_real32         = io_decomp_land_to_land(land_decomp_start, mem_shape_land_soil, var_shape_land_soil, CABLE_NETCDF_FLOAT)
    decomp_grid_land%land_soil_real64         = io_decomp_land_to_land(land_decomp_start, mem_shape_land_soil, var_shape_land_soil, CABLE_NETCDF_DOUBLE)
    decomp_grid_land%land_snow_int32          = io_decomp_land_to_land(land_decomp_start, mem_shape_land_snow, var_shape_land_snow, CABLE_NETCDF_INT)
    decomp_grid_land%land_snow_real32         = io_decomp_land_to_land(land_decomp_start, mem_shape_land_snow, var_shape_land_snow, CABLE_NETCDF_FLOAT)
    decomp_grid_land%land_snow_real64         = io_decomp_land_to_land(land_decomp_start, mem_shape_land_snow, var_shape_land_snow, CABLE_NETCDF_DOUBLE)
    decomp_grid_land%land_rad_int32           = io_decomp_land_to_land(land_decomp_start, mem_shape_land_rad, var_shape_land_rad, CABLE_NETCDF_INT)
    decomp_grid_land%land_rad_real32          = io_decomp_land_to_land(land_decomp_start, mem_shape_land_rad, var_shape_land_rad, CABLE_NETCDF_FLOAT)
    decomp_grid_land%land_rad_real64          = io_decomp_land_to_land(land_decomp_start, mem_shape_land_rad, var_shape_land_rad, CABLE_NETCDF_DOUBLE)
    decomp_grid_land%land_plantcarbon_int32   = io_decomp_land_to_land(land_decomp_start, mem_shape_land_plantcarbon, var_shape_land_plantcarbon, CABLE_NETCDF_INT)
    decomp_grid_land%land_plantcarbon_real32  = io_decomp_land_to_land(land_decomp_start, mem_shape_land_plantcarbon, var_shape_land_plantcarbon, CABLE_NETCDF_FLOAT)
    decomp_grid_land%land_plantcarbon_real64  = io_decomp_land_to_land(land_decomp_start, mem_shape_land_plantcarbon, var_shape_land_plantcarbon, CABLE_NETCDF_DOUBLE)
    decomp_grid_land%land_soilcarbon_int32    = io_decomp_land_to_land(land_decomp_start, mem_shape_land_soilcarbon, var_shape_land_soilcarbon, CABLE_NETCDF_INT)
    decomp_grid_land%land_soilcarbon_real32   = io_decomp_land_to_land(land_decomp_start, mem_shape_land_soilcarbon, var_shape_land_soilcarbon, CABLE_NETCDF_FLOAT)
    decomp_grid_land%land_soilcarbon_real64   = io_decomp_land_to_land(land_decomp_start, mem_shape_land_soilcarbon, var_shape_land_soilcarbon, CABLE_NETCDF_DOUBLE)
    decomp_grid_land%patch_int32              = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch, var_shape_land_patch, CABLE_NETCDF_INT)
    decomp_grid_land%patch_real32             = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch, var_shape_land_patch, CABLE_NETCDF_FLOAT)
    decomp_grid_land%patch_real64             = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch, var_shape_land_patch, CABLE_NETCDF_DOUBLE)
    decomp_grid_land%patch_soil_int32         = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_soil, var_shape_land_patch_soil, CABLE_NETCDF_INT)
    decomp_grid_land%patch_soil_real32        = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_soil, var_shape_land_patch_soil, CABLE_NETCDF_FLOAT)
    decomp_grid_land%patch_soil_real64        = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_soil, var_shape_land_patch_soil, CABLE_NETCDF_DOUBLE)
    decomp_grid_land%patch_snow_int32         = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_snow, var_shape_land_patch_snow, CABLE_NETCDF_INT)
    decomp_grid_land%patch_snow_real32        = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_snow, var_shape_land_patch_snow, CABLE_NETCDF_FLOAT)
    decomp_grid_land%patch_snow_real64        = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_snow, var_shape_land_patch_snow, CABLE_NETCDF_DOUBLE)
    decomp_grid_land%patch_rad_int32          = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_rad, var_shape_land_patch_rad, CABLE_NETCDF_INT)
    decomp_grid_land%patch_rad_real32         = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_rad, var_shape_land_patch_rad, CABLE_NETCDF_FLOAT)
    decomp_grid_land%patch_rad_real64         = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_rad, var_shape_land_patch_rad, CABLE_NETCDF_DOUBLE)
    decomp_grid_land%patch_plantcarbon_int32  = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_plantcarbon, var_shape_land_patch_plantcarbon, CABLE_NETCDF_INT)
    decomp_grid_land%patch_plantcarbon_real32 = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_plantcarbon, var_shape_land_patch_plantcarbon, CABLE_NETCDF_FLOAT)
    decomp_grid_land%patch_plantcarbon_real64 = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_plantcarbon, var_shape_land_patch_plantcarbon, CABLE_NETCDF_DOUBLE)
    decomp_grid_land%patch_soilcarbon_int32   = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_soilcarbon, var_shape_land_patch_soilcarbon, CABLE_NETCDF_INT)
    decomp_grid_land%patch_soilcarbon_real32  = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_soilcarbon, var_shape_land_patch_soilcarbon, CABLE_NETCDF_FLOAT)
    decomp_grid_land%patch_soilcarbon_real64  = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_soilcarbon, var_shape_land_patch_soilcarbon, CABLE_NETCDF_DOUBLE)

    decomp_grid_restart%patch_int32              = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch, var_shape_patch, CABLE_NETCDF_INT)
    decomp_grid_restart%patch_real32             = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch, var_shape_patch, CABLE_NETCDF_FLOAT)
    decomp_grid_restart%patch_real64             = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch, var_shape_patch, CABLE_NETCDF_DOUBLE)
    decomp_grid_restart%patch_soil_int32         = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch_soil, var_shape_patch_soil, CABLE_NETCDF_INT)
    decomp_grid_restart%patch_soil_real32        = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch_soil, var_shape_patch_soil, CABLE_NETCDF_FLOAT)
    decomp_grid_restart%patch_soil_real64        = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch_soil, var_shape_patch_soil, CABLE_NETCDF_DOUBLE)
    decomp_grid_restart%patch_snow_int32         = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch_snow, var_shape_patch_snow, CABLE_NETCDF_INT)
    decomp_grid_restart%patch_snow_real32        = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch_snow, var_shape_patch_snow, CABLE_NETCDF_FLOAT)
    decomp_grid_restart%patch_snow_real64        = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch_snow, var_shape_patch_snow, CABLE_NETCDF_DOUBLE)
    decomp_grid_restart%patch_rad_int32          = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch_rad, var_shape_patch_rad, CABLE_NETCDF_INT)
    decomp_grid_restart%patch_rad_real32         = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch_rad, var_shape_patch_rad, CABLE_NETCDF_FLOAT)
    decomp_grid_restart%patch_rad_real64         = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch_rad, var_shape_patch_rad, CABLE_NETCDF_DOUBLE)
    decomp_grid_restart%patch_plantcarbon_int32  = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch_plantcarbon, var_shape_patch_plantcarbon, CABLE_NETCDF_INT)
    decomp_grid_restart%patch_plantcarbon_real32 = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch_plantcarbon, var_shape_patch_plantcarbon, CABLE_NETCDF_FLOAT)
    decomp_grid_restart%patch_plantcarbon_real64 = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch_plantcarbon, var_shape_patch_plantcarbon, CABLE_NETCDF_DOUBLE)
    decomp_grid_restart%patch_soilcarbon_int32   = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch_soilcarbon, var_shape_patch_soilcarbon, CABLE_NETCDF_INT)
    decomp_grid_restart%patch_soilcarbon_real32  = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch_soilcarbon, var_shape_patch_soilcarbon, CABLE_NETCDF_FLOAT)
    decomp_grid_restart%patch_soilcarbon_real64  = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch_soilcarbon, var_shape_patch_soilcarbon, CABLE_NETCDF_DOUBLE)

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

    if (data_shape_eq(output_var%data_shape, [CABLE_OUTPUT_DIM_PATCH])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_int32
      else
        decomp => output_decomp%land_int32
      end if
    else if (data_shape_eq(output_var%data_shape, [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_SOIL])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_soil_int32
      else
        decomp => output_decomp%land_soil_int32
      end if
    else if (data_shape_eq(output_var%data_shape, [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_SNOW])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_snow_int32
      else
        decomp => output_decomp%land_snow_int32
      end if
    else if (data_shape_eq(output_var%data_shape, [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_RAD])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_rad_int32
      else
        decomp => output_decomp%land_rad_int32
      end if
    else if (data_shape_eq(output_var%data_shape, [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_PLANTCARBON])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_plantcarbon_int32
      else
        decomp => output_decomp%land_plantcarbon_int32
      end if
    else if (data_shape_eq(output_var%data_shape, [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_SOILCARBON])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_soilcarbon_int32
      else
        decomp => output_decomp%land_soilcarbon_int32
      end if
    else
      call cable_abort("Unsupported data shape for output variable " // output_var%name, __FILE__, __LINE__)
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

    if (data_shape_eq(output_var%data_shape, [CABLE_OUTPUT_DIM_PATCH])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_real32
      else
        decomp => output_decomp%land_real32
      end if
    else if (data_shape_eq(output_var%data_shape, [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_SOIL])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_soil_real32
      else
        decomp => output_decomp%land_soil_real32
      end if
    else if (data_shape_eq(output_var%data_shape, [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_SNOW])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_snow_real32
      else
        decomp => output_decomp%land_snow_real32
      end if
    else if (data_shape_eq(output_var%data_shape, [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_RAD])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_rad_real32
      else
        decomp => output_decomp%land_rad_real32
      end if
    else if (data_shape_eq(output_var%data_shape, [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_PLANTCARBON])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_plantcarbon_real32
      else
        decomp => output_decomp%land_plantcarbon_real32
      end if
    else if (data_shape_eq(output_var%data_shape, [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_SOILCARBON])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_soilcarbon_real32
      else
        decomp => output_decomp%land_soilcarbon_real32
      end if
    else
      call cable_abort("Unsupported data shape for output variable " // output_var%name, __FILE__, __LINE__)
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

    if (data_shape_eq(output_var%data_shape, [CABLE_OUTPUT_DIM_PATCH])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_real64
      else
        decomp => output_decomp%land_real64
      end if
    else if (data_shape_eq(output_var%data_shape, [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_SOIL])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_soil_real64
      else
        decomp => output_decomp%land_soil_real64
      end if
    else if (data_shape_eq(output_var%data_shape, [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_SNOW])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_snow_real64
      else
        decomp => output_decomp%land_snow_real64
      end if
    else if (data_shape_eq(output_var%data_shape, [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_RAD])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_rad_real64
      else
        decomp => output_decomp%land_rad_real64
      end if
    else if (data_shape_eq(output_var%data_shape, [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_PLANTCARBON])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_plantcarbon_real64
      else
        decomp => output_decomp%land_plantcarbon_real64
      end if
    else if (data_shape_eq(output_var%data_shape, [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_SOILCARBON])) then
      if (output_var%reduction_method == "none") then
        decomp => output_decomp%patch_soilcarbon_real64
      else
        decomp => output_decomp%land_soilcarbon_real64
      end if
    else
      call cable_abort("Unsupported data shape for output variable " // output_var%name, __FILE__, __LINE__)
    end if

  end subroutine associate_decomp_real64

  subroutine associate_decomp_restart_int32(output_var, decomp)
    type(cable_output_variable_t), intent(in) :: output_var
    class(cable_netcdf_decomp_t), pointer, intent(inout) :: decomp
    type(cable_output_decomp_t), pointer :: output_decomp

    if (data_shape_eq(output_var%data_shape, [CABLE_OUTPUT_DIM_PATCH])) then
      decomp => decomp_grid_restart%patch_int32
    else if (data_shape_eq(output_var%data_shape, [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_SOIL])) then
      decomp => decomp_grid_restart%patch_soil_int32
    else if (data_shape_eq(output_var%data_shape, [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_SNOW])) then
      decomp => decomp_grid_restart%patch_snow_int32
    else if (data_shape_eq(output_var%data_shape, [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_RAD])) then
      decomp => decomp_grid_restart%patch_rad_int32
    else if (data_shape_eq(output_var%data_shape, [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_PLANTCARBON])) then
      decomp => decomp_grid_restart%patch_plantcarbon_int32
    else if (data_shape_eq(output_var%data_shape, [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_SOILCARBON])) then
      decomp => decomp_grid_restart%patch_soilcarbon_int32
    else
      call cable_abort("Unsupported data shape for output variable " // output_var%name, __FILE__, __LINE__)
    end if

  end subroutine associate_decomp_restart_int32

  subroutine associate_decomp_restart_real32(output_var, decomp)
    type(cable_output_variable_t), intent(in) :: output_var
    class(cable_netcdf_decomp_t), pointer, intent(inout) :: decomp

    if (data_shape_eq(output_var%data_shape, [CABLE_OUTPUT_DIM_PATCH])) then
      decomp => decomp_grid_restart%patch_real32
    else if (data_shape_eq(output_var%data_shape, [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_SOIL])) then
      decomp => decomp_grid_restart%patch_soil_real32
    else if (data_shape_eq(output_var%data_shape, [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_SNOW])) then
      decomp => decomp_grid_restart%patch_snow_real32
    else if (data_shape_eq(output_var%data_shape, [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_RAD])) then
      decomp => decomp_grid_restart%patch_rad_real32
    else if (data_shape_eq(output_var%data_shape, [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_PLANTCARBON])) then
      decomp => decomp_grid_restart%patch_plantcarbon_real32
    else if (data_shape_eq(output_var%data_shape, [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_SOILCARBON])) then
      decomp => decomp_grid_restart%patch_soilcarbon_real32
    else
      call cable_abort("Unsupported data shape for output variable " // output_var%name, __FILE__, __LINE__)
    end if

  end subroutine associate_decomp_restart_real32

  subroutine associate_decomp_restart_real64(output_var, decomp)
    type(cable_output_variable_t), intent(in) :: output_var
    class(cable_netcdf_decomp_t), pointer, intent(inout) :: decomp

    if (data_shape_eq(output_var%data_shape, [CABLE_OUTPUT_DIM_PATCH])) then
      decomp => decomp_grid_restart%patch_real64
    else if (data_shape_eq(output_var%data_shape, [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_SOIL])) then
      decomp => decomp_grid_restart%patch_soil_real64
    else if (data_shape_eq(output_var%data_shape, [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_SNOW])) then
      decomp => decomp_grid_restart%patch_snow_real64
    else if (data_shape_eq(output_var%data_shape, [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_RAD])) then
      decomp => decomp_grid_restart%patch_rad_real64
    else if (data_shape_eq(output_var%data_shape, [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_PLANTCARBON])) then
      decomp => decomp_grid_restart%patch_plantcarbon_real64
    else if (data_shape_eq(output_var%data_shape, [CABLE_OUTPUT_DIM_PATCH, CABLE_OUTPUT_DIM_SOILCARBON])) then
      decomp => decomp_grid_restart%patch_soilcarbon_real64
    else
      call cable_abort("Unsupported data shape for output variable " // output_var%name, __FILE__, __LINE__)
    end if

  end subroutine associate_decomp_restart_real64

end module
