module cable_io_decomp_mod
  use cable_def_types_mod, only: mp, mp_global
  use cable_def_types_mod, only: mland, mland_global
  use cable_def_types_mod, only: ms
  use cable_def_types_mod, only: msn
  use cable_def_types_mod, only: nrb
  use cable_def_types_mod, only: ncp
  use cable_def_types_mod, only: ncs

  use cable_io_vars_module, only: xdimsize, ydimsize
  use cable_io_vars_module, only: land_x, land_y
  use cable_io_vars_module, only: landpt
  use cable_io_vars_module, only: max_vegpatches
  use cable_io_vars_module, only: land_decomp_start
  use cable_io_vars_module, only: patch_decomp_start
  use cable_io_vars_module, only: output
  use cable_io_vars_module, only: metGrid

  use cable_netcdf_decomp_util_mod, only: dim_spec_t
  use cable_netcdf_decomp_util_mod, only: io_decomp_land_to_x_y
  use cable_netcdf_decomp_util_mod, only: io_decomp_patch_to_x_y_patch
  use cable_netcdf_decomp_util_mod, only: io_decomp_land_to_land
  use cable_netcdf_decomp_util_mod, only: io_decomp_patch_to_land_patch
  use cable_netcdf_decomp_util_mod, only: io_decomp_patch_to_patch

  use cable_netcdf_mod, only: cable_netcdf_decomp_t, CABLE_NETCDF_FLOAT

  implicit none
  private

  public :: &
    io_decomp_t, &
    cable_io_decomp_init

  type io_decomp_t
    class(cable_netcdf_decomp_t), allocatable :: patch_to_x_y_patch_real32
    class(cable_netcdf_decomp_t), allocatable :: patch_soil_to_x_y_patch_soil_real32
    class(cable_netcdf_decomp_t), allocatable :: patch_snow_to_x_y_patch_snow_real32
    class(cable_netcdf_decomp_t), allocatable :: patch_rad_to_x_y_patch_rad_real32
    class(cable_netcdf_decomp_t), allocatable :: patch_plantcarbon_to_x_y_patch_plantcarbon_real32
    class(cable_netcdf_decomp_t), allocatable :: patch_soilcarbon_to_x_y_patch_soilcarbon_real32

    class(cable_netcdf_decomp_t), allocatable :: patch_to_land_patch_real32
    class(cable_netcdf_decomp_t), allocatable :: patch_soil_to_land_patch_soil_real32
    class(cable_netcdf_decomp_t), allocatable :: patch_snow_to_land_patch_snow_real32
    class(cable_netcdf_decomp_t), allocatable :: patch_rad_to_land_patch_rad_real32
    class(cable_netcdf_decomp_t), allocatable :: patch_plantcarbon_to_land_patch_plantcarbon_real32
    class(cable_netcdf_decomp_t), allocatable :: patch_soilcarbon_to_land_patch_soilcarbon_real32

    class(cable_netcdf_decomp_t), allocatable :: patch_to_patch_real32
    class(cable_netcdf_decomp_t), allocatable :: patch_soil_to_patch_soil_real32
    class(cable_netcdf_decomp_t), allocatable :: patch_snow_to_patch_snow_real32
    class(cable_netcdf_decomp_t), allocatable :: patch_rad_to_patch_rad_real32
    class(cable_netcdf_decomp_t), allocatable :: patch_plantcarbon_to_patch_plantcarbon_real32
    class(cable_netcdf_decomp_t), allocatable :: patch_soilcarbon_to_patch_soilcarbon_real32

    class(cable_netcdf_decomp_t), allocatable :: land_to_x_y_real32
    class(cable_netcdf_decomp_t), allocatable :: land_soil_to_x_y_soil_real32
    class(cable_netcdf_decomp_t), allocatable :: land_snow_to_x_y_snow_real32
    class(cable_netcdf_decomp_t), allocatable :: land_rad_to_x_y_rad_real32
    class(cable_netcdf_decomp_t), allocatable :: land_plantcarbon_to_x_y_plantcarbon_real32
    class(cable_netcdf_decomp_t), allocatable :: land_soilcarbon_to_x_y_soilcarbon_real32

    class(cable_netcdf_decomp_t), allocatable :: land_to_land_real32
    class(cable_netcdf_decomp_t), allocatable :: land_soil_to_land_soil_real32
    class(cable_netcdf_decomp_t), allocatable :: land_snow_to_land_snow_real32
    class(cable_netcdf_decomp_t), allocatable :: land_rad_to_land_rad_real32
    class(cable_netcdf_decomp_t), allocatable :: land_plantcarbon_to_land_plantcarbon_real32
    class(cable_netcdf_decomp_t), allocatable :: land_soilcarbon_to_land_soilcarbon_real32

    class(cable_netcdf_decomp_t), pointer :: output_base_real32
    class(cable_netcdf_decomp_t), pointer :: output_base_soil_real32
    class(cable_netcdf_decomp_t), pointer :: output_base_snow_real32
    class(cable_netcdf_decomp_t), pointer :: output_base_rad_real32
    class(cable_netcdf_decomp_t), pointer :: output_base_plantcarbon_real32
    class(cable_netcdf_decomp_t), pointer :: output_base_soilcarbon_real32
    class(cable_netcdf_decomp_t), pointer :: output_base_patch_real32
    class(cable_netcdf_decomp_t), pointer :: output_base_patch_soil_real32
    class(cable_netcdf_decomp_t), pointer :: output_base_patch_snow_real32
    class(cable_netcdf_decomp_t), pointer :: output_base_patch_rad_real32
    class(cable_netcdf_decomp_t), pointer :: output_base_patch_plantcarbon_real32
    class(cable_netcdf_decomp_t), pointer :: output_base_patch_soilcarbon_real32
  end type io_decomp_t

contains

  subroutine cable_io_decomp_init(io_decomp)
    type(io_decomp_t), intent(out), target :: io_decomp

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

    logical :: requires_land_output_grid, requires_x_y_output_grid

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

    io_decomp%land_to_x_y_real32                                 = io_decomp_land_to_x_y(land_x, land_y, mem_shape_land, var_shape_x_y, CABLE_NETCDF_FLOAT)
    io_decomp%land_soil_to_x_y_soil_real32                       = io_decomp_land_to_x_y(land_x, land_y, mem_shape_land_soil, var_shape_x_y_soil, CABLE_NETCDF_FLOAT)
    io_decomp%land_snow_to_x_y_snow_real32                       = io_decomp_land_to_x_y(land_x, land_y, mem_shape_land_snow, var_shape_x_y_snow, CABLE_NETCDF_FLOAT)
    io_decomp%land_rad_to_x_y_rad_real32                         = io_decomp_land_to_x_y(land_x, land_y, mem_shape_land_rad, var_shape_x_y_rad, CABLE_NETCDF_FLOAT)
    io_decomp%land_plantcarbon_to_x_y_plantcarbon_real32         = io_decomp_land_to_x_y(land_x, land_y, mem_shape_land_plantcarbon, var_shape_x_y_plantcarbon, CABLE_NETCDF_FLOAT)
    io_decomp%land_soilcarbon_to_x_y_soilcarbon_real32           = io_decomp_land_to_x_y(land_x, land_y, mem_shape_land_soilcarbon, var_shape_x_y_soilcarbon, CABLE_NETCDF_FLOAT)

    io_decomp%land_to_land_real32                                = io_decomp_land_to_land(land_decomp_start, mem_shape_land, var_shape_land, CABLE_NETCDF_FLOAT)
    io_decomp%land_soil_to_land_soil_real32                      = io_decomp_land_to_land(land_decomp_start, mem_shape_land_soil, var_shape_land_soil, CABLE_NETCDF_FLOAT)
    io_decomp%land_snow_to_land_snow_real32                      = io_decomp_land_to_land(land_decomp_start, mem_shape_land_snow, var_shape_land_snow, CABLE_NETCDF_FLOAT)
    io_decomp%land_rad_to_land_rad_real32                        = io_decomp_land_to_land(land_decomp_start, mem_shape_land_rad, var_shape_land_rad, CABLE_NETCDF_FLOAT)
    io_decomp%land_plantcarbon_to_land_plantcarbon_real32        = io_decomp_land_to_land(land_decomp_start, mem_shape_land_plantcarbon, var_shape_land_plantcarbon, CABLE_NETCDF_FLOAT)
    io_decomp%land_soilcarbon_to_land_soilcarbon_real32          = io_decomp_land_to_land(land_decomp_start, mem_shape_land_soilcarbon, var_shape_land_soilcarbon, CABLE_NETCDF_FLOAT)

    io_decomp%patch_to_x_y_patch_real32                          = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch, var_shape_x_y_patch, CABLE_NETCDF_FLOAT)
    io_decomp%patch_soil_to_x_y_patch_soil_real32                = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_soil, var_shape_x_y_soil, CABLE_NETCDF_FLOAT)
    io_decomp%patch_snow_to_x_y_patch_snow_real32                = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_snow, var_shape_x_y_snow, CABLE_NETCDF_FLOAT)
    io_decomp%patch_rad_to_x_y_patch_rad_real32                  = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_rad, var_shape_x_y_rad, CABLE_NETCDF_FLOAT)
    io_decomp%patch_plantcarbon_to_x_y_patch_plantcarbon_real32  = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_plantcarbon, var_shape_x_y_plantcarbon, CABLE_NETCDF_FLOAT)
    io_decomp%patch_soilcarbon_to_x_y_patch_soilcarbon_real32    = io_decomp_patch_to_x_y_patch(land_x, land_y, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_soilcarbon, var_shape_x_y_soilcarbon, CABLE_NETCDF_FLOAT)

    io_decomp%patch_to_land_patch_real32                         = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch, var_shape_land_patch, CABLE_NETCDF_FLOAT)
    io_decomp%patch_soil_to_land_patch_soil_real32               = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_soil, var_shape_land_soil, CABLE_NETCDF_FLOAT)
    io_decomp%patch_snow_to_land_patch_snow_real32               = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_snow, var_shape_land_snow, CABLE_NETCDF_FLOAT)
    io_decomp%patch_rad_to_land_patch_rad_real32                 = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_rad, var_shape_land_rad, CABLE_NETCDF_FLOAT)
    io_decomp%patch_plantcarbon_to_land_patch_plantcarbon_real32 = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_plantcarbon, var_shape_land_plantcarbon, CABLE_NETCDF_FLOAT)
    io_decomp%patch_soilcarbon_to_land_patch_soilcarbon_real32   = io_decomp_patch_to_land_patch(land_decomp_start, landpt(:)%cstart, landpt(:)%nap, mem_shape_patch_soilcarbon, var_shape_land_soilcarbon, CABLE_NETCDF_FLOAT)

    io_decomp%patch_to_patch_real32                              = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch, var_shape_patch, CABLE_NETCDF_FLOAT)
    io_decomp%patch_soil_to_patch_soil_real32                    = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch_soil, var_shape_patch_soil, CABLE_NETCDF_FLOAT)
    io_decomp%patch_snow_to_patch_snow_real32                    = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch_snow, var_shape_patch_snow, CABLE_NETCDF_FLOAT)
    io_decomp%patch_rad_to_patch_rad_real32                      = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch_rad, var_shape_patch_rad, CABLE_NETCDF_FLOAT)
    io_decomp%patch_plantcarbon_to_patch_plantcarbon_real32      = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch_plantcarbon, var_shape_patch_plantcarbon, CABLE_NETCDF_FLOAT)
    io_decomp%patch_soilcarbon_to_patch_soilcarbon_real32        = io_decomp_patch_to_patch(patch_decomp_start, mem_shape_patch_soilcarbon, var_shape_patch_soilcarbon, CABLE_NETCDF_FLOAT)

    requires_x_y_output_grid = ( &
      ( &
        output%grid == 'default' .AND. metGrid == 'mask' &
      ) .OR. ( &
        output%grid == 'mask' .OR. output%grid == 'ALMA' &
      ) &
    )
    if (requires_x_y_output_grid) then
      io_decomp%output_base_real32                   => io_decomp%land_to_x_y_real32
      io_decomp%output_base_soil_real32              => io_decomp%land_soil_to_x_y_soil_real32
      io_decomp%output_base_snow_real32              => io_decomp%land_snow_to_x_y_snow_real32
      io_decomp%output_base_rad_real32               => io_decomp%land_rad_to_x_y_rad_real32
      io_decomp%output_base_plantcarbon_real32       => io_decomp%land_plantcarbon_to_x_y_plantcarbon_real32
      io_decomp%output_base_soilcarbon_real32        => io_decomp%land_soilcarbon_to_x_y_soilcarbon_real32
      io_decomp%output_base_patch_real32             => io_decomp%patch_to_x_y_patch_real32
      io_decomp%output_base_patch_soil_real32        => io_decomp%patch_soil_to_x_y_patch_soil_real32
      io_decomp%output_base_patch_snow_real32        => io_decomp%patch_snow_to_x_y_patch_snow_real32
      io_decomp%output_base_patch_rad_real32         => io_decomp%patch_rad_to_x_y_patch_rad_real32
      io_decomp%output_base_patch_plantcarbon_real32 => io_decomp%patch_plantcarbon_to_x_y_patch_plantcarbon_real32
      io_decomp%output_base_patch_soilcarbon_real32  => io_decomp%patch_soilcarbon_to_x_y_patch_soilcarbon_real32
    end if

    requires_land_output_grid = ( &
      output%grid == 'land' .OR. (output%grid == 'default' .AND. metGrid == 'land') &
    )
    if (requires_land_output_grid) then
      io_decomp%output_base_real32                   => io_decomp%land_to_land_real32
      io_decomp%output_base_soil_real32              => io_decomp%land_soil_to_land_soil_real32
      io_decomp%output_base_snow_real32              => io_decomp%land_snow_to_land_snow_real32
      io_decomp%output_base_rad_real32               => io_decomp%land_rad_to_land_rad_real32
      io_decomp%output_base_plantcarbon_real32       => io_decomp%land_plantcarbon_to_land_plantcarbon_real32
      io_decomp%output_base_soilcarbon_real32        => io_decomp%land_soilcarbon_to_land_soilcarbon_real32
      io_decomp%output_base_patch_real32             => io_decomp%patch_to_land_patch_real32
      io_decomp%output_base_patch_soil_real32        => io_decomp%patch_soil_to_land_patch_soil_real32
      io_decomp%output_base_patch_snow_real32        => io_decomp%patch_snow_to_land_patch_snow_real32
      io_decomp%output_base_patch_rad_real32         => io_decomp%patch_rad_to_land_patch_rad_real32
      io_decomp%output_base_patch_plantcarbon_real32 => io_decomp%patch_plantcarbon_to_land_patch_plantcarbon_real32
      io_decomp%output_base_patch_soilcarbon_real32  => io_decomp%patch_soilcarbon_to_land_patch_soilcarbon_real32
    end if

  end subroutine

end module
