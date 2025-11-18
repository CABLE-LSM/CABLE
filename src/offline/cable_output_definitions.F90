module cable_output_definitions_mod
  use iso_fortran_env, only: real32

  use cable_abort_module, only: cable_abort

  use cable_def_types_mod, only: canopy_type

  use cable_io_vars_module, only: metGrid

  use cable_netcdf_mod, only: CABLE_NETCDF_FLOAT

  use aggregator_mod, only: new_aggregator
  use aggregator_mod, only: aggregator_real32_1d_t

  use cable_netcdf_mod, only: MAX_LEN_DIM => CABLE_NETCDF_MAX_STR_LEN_DIM

  use cable_output_prototype_v2_mod, only: requires_x_y_output_grid
  use cable_output_prototype_v2_mod, only: requires_land_output_grid
  use cable_output_prototype_v2_mod, only: cable_output_add_variable
  use cable_output_prototype_v2_mod, only: cable_output_aggregator_t
  use cable_output_prototype_v2_mod, only: cable_output_add_aggregator
  use cable_output_prototype_v2_mod, only: output_options, patchout_options
  use cable_output_prototype_v2_mod, only: CABLE_OUTPUT_SHAPE_TYPE_BASE
  use cable_output_prototype_v2_mod, only: CABLE_OUTPUT_SHAPE_TYPE_BASE_SOIL
  use cable_output_prototype_v2_mod, only: CABLE_OUTPUT_SHAPE_TYPE_BASE_SNOW
  use cable_output_prototype_v2_mod, only: CABLE_OUTPUT_SHAPE_TYPE_BASE_RAD
  use cable_output_prototype_v2_mod, only: CABLE_OUTPUT_SHAPE_TYPE_BASE_PLANTCARBON
  use cable_output_prototype_v2_mod, only: CABLE_OUTPUT_SHAPE_TYPE_BASE_SOILCARBON
  use cable_output_prototype_v2_mod, only: CABLE_OUTPUT_SHAPE_TYPE_BASE_PATCH
  use cable_output_prototype_v2_mod, only: CABLE_OUTPUT_SHAPE_TYPE_BASE_PATCH_SOIL
  use cable_output_prototype_v2_mod, only: CABLE_OUTPUT_SHAPE_TYPE_BASE_PATCH_SNOW
  use cable_output_prototype_v2_mod, only: CABLE_OUTPUT_SHAPE_TYPE_BASE_PATCH_RAD
  use cable_output_prototype_v2_mod, only: CABLE_OUTPUT_SHAPE_TYPE_BASE_PATCH_PLANTCARBON
  use cable_output_prototype_v2_mod, only: CABLE_OUTPUT_SHAPE_TYPE_BASE_PATCH_SOILCARBON

  use cable_checks_module, only: ranges ! TODO(Sean): pass ranges via an argument rather than use module

  implicit none
  private

  public :: cable_output_definitions_set

contains

  subroutine cable_output_definitions_set(canopy)
    type(canopy_type), intent(inout) :: canopy

    character(len=MAX_LEN_DIM), allocatable :: base_dims(:)

    if (requires_x_y_output_grid(output_options%grid, metGrid)) then
      base_dims = ["x", "y"]
    else if (requires_land_output_grid(output_options%grid, metGrid)) then
      base_dims = ["land"]
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
      active=output_options%Qh .and. (output_options%patch .OR. patchout_options%Qh), &
      grid_cell_averaging=.false., &
      shape_type=CABLE_OUTPUT_SHAPE_TYPE_BASE_PATCH, &
      accumulation_frequency="all", &
      aggregation_frequency=output_options%averaging, &
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
      active=output_options%Qh .and. .not. (output_options%patch .OR. patchout_options%Qh), &
      grid_cell_averaging=.true., &
      shape_type=CABLE_OUTPUT_SHAPE_TYPE_BASE, &
      accumulation_frequency="all", &
      aggregation_frequency=output_options%averaging, &
      aggregator=new_aggregator( &
        source_data=canopy%fh, &
        method="mean" &
      ) &
    )

    add_variable_Tmx: block
      real(kind=real32), pointer :: tdaymx(:)
      type(cable_output_aggregator_t), target :: tdaymx_intermediate_aggregator

      if (output_options%Tex .and. output_options%averaging == "monthly") then
        ! Create an intermmediate aggregator to compute daily maximum T
        call cable_output_add_aggregator( &
          aggregator=new_aggregator( &
            source_data=canopy%tscrn, &
            method="max" &
          ), &
          accumulation_frequency="all", &
          aggregation_frequency="daily", &
          output_aggregator=tdaymx_intermediate_aggregator &
        )
        select type(aggregator => tdaymx_intermediate_aggregator%aggregator_handle%aggregator)
        type is (aggregator_real32_1d_t)
          ! This is required to ensure that the storage for tdaymx is allocated.
          call aggregator%init()
          tdaymx => aggregator%storage
        end select
      else
        tdaymx => canopy%tscrn ! dummy assignment when Tmx is not needed
      end if

      call cable_output_add_variable( &
        name="Tmx", &
        dims=[base_dims, "time"], &
        var_type=CABLE_NETCDF_FLOAT, &
        units="oC", &
        long_name="averaged daily maximum screen-level T", &
        active=( &
          output_options%Tex .and. &
          output_options%averaging == "monthly" .and. &
          .not. (output_options%patch .OR. patchout_options%Tex) &
        ), &
        grid_cell_averaging=.true., &
        shape_type=CABLE_OUTPUT_SHAPE_TYPE_BASE, &
        range=ranges%Tscrn, &
        accumulation_frequency="daily", &
        aggregation_frequency=output_options%averaging, &
        aggregator=new_aggregator( &
          source_data=tdaymx, &
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
          output_options%Tex .and. &
          output_options%averaging == "monthly" .and. &
          (output_options%patch .OR. patchout_options%Tex) &
        ), &
        grid_cell_averaging=.false., &
        shape_type=CABLE_OUTPUT_SHAPE_TYPE_BASE_PATCH, &
        range=ranges%Tscrn, &
        accumulation_frequency="daily", &
        aggregation_frequency=output_options%averaging, &
        aggregator=new_aggregator( &
          source_data=tdaymx, &
          method="mean" &
        ) &
      )

    end block add_variable_Tmx

  end subroutine cable_output_definitions_set

end module
