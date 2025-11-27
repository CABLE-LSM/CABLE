module cable_output_definitions_mod
  use cable_abort_module, only: cable_abort

  use cable_def_types_mod, only: canopy_type

  use cable_io_vars_module, only: metGrid

  use cable_netcdf_mod, only: CABLE_NETCDF_FLOAT

  use aggregator_mod, only: new_aggregator

  use cable_netcdf_mod, only: cable_netcdf_decomp_t
  use cable_netcdf_mod, only: MAX_LEN_DIM => CABLE_NETCDF_MAX_STR_LEN_DIM

  use cable_io_decomp_mod, only: io_decomp_t

  use cable_output_prototype_v2_mod, only: requires_x_y_output_grid
  use cable_output_prototype_v2_mod, only: requires_land_output_grid
  use cable_output_prototype_v2_mod, only: output, patchout
  use cable_output_prototype_v2_mod, only: cable_output_add_variable

  use cable_checks_module, only: ranges

  implicit none
  private

  public :: cable_output_definitions_set

contains

  subroutine cable_output_definitions_set(io_decomp, canopy)
    class(io_decomp_t), intent(in), target :: io_decomp
    type(canopy_type), intent(inout) :: canopy

    character(len=MAX_LEN_DIM), allocatable :: base_dims(:)

    class(cable_netcdf_decomp_t), pointer :: output_decomp_base_int32
    class(cable_netcdf_decomp_t), pointer :: output_decomp_base_real32
    class(cable_netcdf_decomp_t), pointer :: output_decomp_base_real64
    class(cable_netcdf_decomp_t), pointer :: output_decomp_base_soil_int32
    class(cable_netcdf_decomp_t), pointer :: output_decomp_base_soil_real32
    class(cable_netcdf_decomp_t), pointer :: output_decomp_base_soil_real64
    class(cable_netcdf_decomp_t), pointer :: output_decomp_base_snow_int32
    class(cable_netcdf_decomp_t), pointer :: output_decomp_base_snow_real32
    class(cable_netcdf_decomp_t), pointer :: output_decomp_base_snow_real64
    class(cable_netcdf_decomp_t), pointer :: output_decomp_base_rad_int32
    class(cable_netcdf_decomp_t), pointer :: output_decomp_base_rad_real32
    class(cable_netcdf_decomp_t), pointer :: output_decomp_base_rad_real64
    class(cable_netcdf_decomp_t), pointer :: output_decomp_base_plantcarbon_int32
    class(cable_netcdf_decomp_t), pointer :: output_decomp_base_plantcarbon_real32
    class(cable_netcdf_decomp_t), pointer :: output_decomp_base_plantcarbon_real64
    class(cable_netcdf_decomp_t), pointer :: output_decomp_base_soilcarbon_int32
    class(cable_netcdf_decomp_t), pointer :: output_decomp_base_soilcarbon_real32
    class(cable_netcdf_decomp_t), pointer :: output_decomp_base_soilcarbon_real64
    class(cable_netcdf_decomp_t), pointer :: output_decomp_base_patch_int32
    class(cable_netcdf_decomp_t), pointer :: output_decomp_base_patch_real32
    class(cable_netcdf_decomp_t), pointer :: output_decomp_base_patch_real64
    class(cable_netcdf_decomp_t), pointer :: output_decomp_base_patch_soil_int32
    class(cable_netcdf_decomp_t), pointer :: output_decomp_base_patch_soil_real32
    class(cable_netcdf_decomp_t), pointer :: output_decomp_base_patch_soil_real64
    class(cable_netcdf_decomp_t), pointer :: output_decomp_base_patch_snow_int32
    class(cable_netcdf_decomp_t), pointer :: output_decomp_base_patch_snow_real32
    class(cable_netcdf_decomp_t), pointer :: output_decomp_base_patch_snow_real64
    class(cable_netcdf_decomp_t), pointer :: output_decomp_base_patch_rad_int32
    class(cable_netcdf_decomp_t), pointer :: output_decomp_base_patch_rad_real32
    class(cable_netcdf_decomp_t), pointer :: output_decomp_base_patch_rad_real64
    class(cable_netcdf_decomp_t), pointer :: output_decomp_base_patch_plantcarbon_int32
    class(cable_netcdf_decomp_t), pointer :: output_decomp_base_patch_plantcarbon_real32
    class(cable_netcdf_decomp_t), pointer :: output_decomp_base_patch_plantcarbon_real64
    class(cable_netcdf_decomp_t), pointer :: output_decomp_base_patch_soilcarbon_int32
    class(cable_netcdf_decomp_t), pointer :: output_decomp_base_patch_soilcarbon_real32
    class(cable_netcdf_decomp_t), pointer :: output_decomp_base_patch_soilcarbon_real64

    if (requires_x_y_output_grid(output%grid, metGrid)) then
      base_dims = ["x", "y"]
      output_decomp_base_int32                    => io_decomp%land_to_x_y_int32
      output_decomp_base_real32                   => io_decomp%land_to_x_y_real32
      output_decomp_base_real64                   => io_decomp%land_to_x_y_real64
      output_decomp_base_soil_int32               => io_decomp%land_soil_to_x_y_soil_int32
      output_decomp_base_soil_real32              => io_decomp%land_soil_to_x_y_soil_real32
      output_decomp_base_soil_real64              => io_decomp%land_soil_to_x_y_soil_real64
      output_decomp_base_snow_int32               => io_decomp%land_snow_to_x_y_snow_int32
      output_decomp_base_snow_real32              => io_decomp%land_snow_to_x_y_snow_real32
      output_decomp_base_snow_real64              => io_decomp%land_snow_to_x_y_snow_real64
      output_decomp_base_rad_int32                => io_decomp%land_rad_to_x_y_rad_int32
      output_decomp_base_rad_real32               => io_decomp%land_rad_to_x_y_rad_real32
      output_decomp_base_rad_real64               => io_decomp%land_rad_to_x_y_rad_real64
      output_decomp_base_plantcarbon_int32        => io_decomp%land_plantcarbon_to_x_y_plantcarbon_int32
      output_decomp_base_plantcarbon_real32       => io_decomp%land_plantcarbon_to_x_y_plantcarbon_real32
      output_decomp_base_plantcarbon_real64       => io_decomp%land_plantcarbon_to_x_y_plantcarbon_real64
      output_decomp_base_soilcarbon_int32         => io_decomp%land_soilcarbon_to_x_y_soilcarbon_int32
      output_decomp_base_soilcarbon_real32        => io_decomp%land_soilcarbon_to_x_y_soilcarbon_real32
      output_decomp_base_soilcarbon_real64        => io_decomp%land_soilcarbon_to_x_y_soilcarbon_real64
      output_decomp_base_patch_int32              => io_decomp%patch_to_x_y_patch_int32
      output_decomp_base_patch_real32             => io_decomp%patch_to_x_y_patch_real32
      output_decomp_base_patch_real64             => io_decomp%patch_to_x_y_patch_real64
      output_decomp_base_patch_soil_int32         => io_decomp%patch_soil_to_x_y_patch_soil_int32
      output_decomp_base_patch_soil_real32        => io_decomp%patch_soil_to_x_y_patch_soil_real32
      output_decomp_base_patch_soil_real64        => io_decomp%patch_soil_to_x_y_patch_soil_real64
      output_decomp_base_patch_snow_int32         => io_decomp%patch_snow_to_x_y_patch_snow_int32
      output_decomp_base_patch_snow_real32        => io_decomp%patch_snow_to_x_y_patch_snow_real32
      output_decomp_base_patch_snow_real64        => io_decomp%patch_snow_to_x_y_patch_snow_real64
      output_decomp_base_patch_rad_int32          => io_decomp%patch_rad_to_x_y_patch_rad_int32
      output_decomp_base_patch_rad_real32         => io_decomp%patch_rad_to_x_y_patch_rad_real32
      output_decomp_base_patch_rad_real64         => io_decomp%patch_rad_to_x_y_patch_rad_real64
      output_decomp_base_patch_plantcarbon_int32  => io_decomp%patch_plantcarbon_to_x_y_patch_plantcarbon_int32
      output_decomp_base_patch_plantcarbon_real32 => io_decomp%patch_plantcarbon_to_x_y_patch_plantcarbon_real32
      output_decomp_base_patch_plantcarbon_real64 => io_decomp%patch_plantcarbon_to_x_y_patch_plantcarbon_real64
      output_decomp_base_patch_soilcarbon_int32   => io_decomp%patch_soilcarbon_to_x_y_patch_soilcarbon_int32
      output_decomp_base_patch_soilcarbon_real32  => io_decomp%patch_soilcarbon_to_x_y_patch_soilcarbon_real32
      output_decomp_base_patch_soilcarbon_real64  => io_decomp%patch_soilcarbon_to_x_y_patch_soilcarbon_real64
    else if (requires_land_output_grid(output%grid, metGrid)) then
      base_dims = ["land"]
      output_decomp_base_int32                    => io_decomp%land_to_land_int32
      output_decomp_base_real32                   => io_decomp%land_to_land_real32
      output_decomp_base_real64                   => io_decomp%land_to_land_real64
      output_decomp_base_soil_int32               => io_decomp%land_soil_to_land_soil_int32
      output_decomp_base_soil_real32              => io_decomp%land_soil_to_land_soil_real32
      output_decomp_base_soil_real64              => io_decomp%land_soil_to_land_soil_real64
      output_decomp_base_snow_int32               => io_decomp%land_snow_to_land_snow_int32
      output_decomp_base_snow_real32              => io_decomp%land_snow_to_land_snow_real32
      output_decomp_base_snow_real64              => io_decomp%land_snow_to_land_snow_real64
      output_decomp_base_rad_int32                => io_decomp%land_rad_to_land_rad_int32
      output_decomp_base_rad_real32               => io_decomp%land_rad_to_land_rad_real32
      output_decomp_base_rad_real64               => io_decomp%land_rad_to_land_rad_real64
      output_decomp_base_plantcarbon_int32        => io_decomp%land_plantcarbon_to_land_plantcarbon_int32
      output_decomp_base_plantcarbon_real32       => io_decomp%land_plantcarbon_to_land_plantcarbon_real32
      output_decomp_base_plantcarbon_real64       => io_decomp%land_plantcarbon_to_land_plantcarbon_real64
      output_decomp_base_soilcarbon_int32         => io_decomp%land_soilcarbon_to_land_soilcarbon_int32
      output_decomp_base_soilcarbon_real32        => io_decomp%land_soilcarbon_to_land_soilcarbon_real32
      output_decomp_base_soilcarbon_real64        => io_decomp%land_soilcarbon_to_land_soilcarbon_real64
      output_decomp_base_patch_int32              => io_decomp%patch_to_land_patch_int32
      output_decomp_base_patch_real32             => io_decomp%patch_to_land_patch_real32
      output_decomp_base_patch_real64             => io_decomp%patch_to_land_patch_real64
      output_decomp_base_patch_soil_int32         => io_decomp%patch_soil_to_land_patch_soil_int32
      output_decomp_base_patch_soil_real32        => io_decomp%patch_soil_to_land_patch_soil_real32
      output_decomp_base_patch_soil_real64        => io_decomp%patch_soil_to_land_patch_soil_real64
      output_decomp_base_patch_snow_int32         => io_decomp%patch_snow_to_land_patch_snow_int32
      output_decomp_base_patch_snow_real32        => io_decomp%patch_snow_to_land_patch_snow_real32
      output_decomp_base_patch_snow_real64        => io_decomp%patch_snow_to_land_patch_snow_real64
      output_decomp_base_patch_rad_int32          => io_decomp%patch_rad_to_land_patch_rad_int32
      output_decomp_base_patch_rad_real32         => io_decomp%patch_rad_to_land_patch_rad_real32
      output_decomp_base_patch_rad_real64         => io_decomp%patch_rad_to_land_patch_rad_real64
      output_decomp_base_patch_plantcarbon_int32  => io_decomp%patch_plantcarbon_to_land_patch_plantcarbon_int32
      output_decomp_base_patch_plantcarbon_real32 => io_decomp%patch_plantcarbon_to_land_patch_plantcarbon_real32
      output_decomp_base_patch_plantcarbon_real64 => io_decomp%patch_plantcarbon_to_land_patch_plantcarbon_real64
      output_decomp_base_patch_soilcarbon_int32   => io_decomp%patch_soilcarbon_to_land_patch_soilcarbon_int32
      output_decomp_base_patch_soilcarbon_real32  => io_decomp%patch_soilcarbon_to_land_patch_soilcarbon_real32
      output_decomp_base_patch_soilcarbon_real64  => io_decomp%patch_soilcarbon_to_land_patch_soilcarbon_real64
    else
      call cable_abort("Error: Unable to determine output grid type", __FILE__, __LINE__)
    end if

    call cable_output_add_variable( &
      name="Qh", &
      dims=[base_dims, "patch", "time"], &
      var_type=CABLE_NETCDF_FLOAT, &
      units="W/m^2", &
      long_name="Surface sensible heat flux", &
      range=ranges%Qh, &
      active=output%Qh .and. (output%patch .OR. patchout%Qh), &
      decomp=output_decomp_base_patch_real32, &
      aggregator=new_aggregator( &
        source_data=canopy%fh, &
        method="mean" &
      ) &
    )

    call cable_output_add_variable( &
      name="Qh", &
      dims=[base_dims, "time"], &
      var_type=CABLE_NETCDF_FLOAT, &
      units="W/m^2", &
      long_name="Surface sensible heat flux", &
      range=ranges%Qh, &
      active=output%Qh .and. .not. (output%patch .OR. patchout%Qh), &
      reduction_method="grid_cell_average", &
      decomp=output_decomp_base_real32, &
      aggregator=new_aggregator( &
        source_data=canopy%fh, &
        method="mean" &
      ) &
    )

    call cable_output_add_variable( &
      name="Tmx", &
      dims=[base_dims, "patch", "time"], &
      var_type=CABLE_NETCDF_FLOAT, &
      units="oC", &
      long_name="averaged daily maximum screen-level T", &
      active=( &
        output%Tex .and. &
        output%averaging == "monthly" .and. &
        (output%patch .OR. patchout%Tex) &
      ), &
      decomp=output_decomp_base_patch_real32, &
      range=ranges%Tscrn, &
      accumulation_frequency="daily", &
      aggregator=new_aggregator( &
        source_data=canopy%tscrn_max_daily%storage, &
        method="mean" &
      ) &
    )

    call cable_output_add_variable( &
      name="Tmx", &
      dims=[base_dims, "time"], &
      var_type=CABLE_NETCDF_FLOAT, &
      units="oC", &
      long_name="averaged daily maximum screen-level T", &
      active=( &
        output%Tex .and. &
        output%averaging == "monthly" .and. &
        .not. (output%patch .OR. patchout%Tex) &
      ), &
      reduction_method="grid_cell_average", &
      decomp=output_decomp_base_real32, &
      range=ranges%Tscrn, &
      accumulation_frequency="daily", &
      aggregator=new_aggregator( &
        source_data=canopy%tscrn_max_daily%storage, &
        method="mean" &
      ) &
    )

  end subroutine cable_output_definitions_set

end module
