module cable_output_definitions_mod
  use cable_def_types_mod, only: canopy_type
  use cable_def_types_mod, only: soil_parameter_type
  use cable_def_types_mod, only: mvtype, mstype

  use cable_netcdf_mod, only: CABLE_NETCDF_INT
  use cable_netcdf_mod, only: CABLE_NETCDF_FLOAT

  use aggregator_mod, only: new_aggregator

  use cable_io_vars_module, only: output, patchout
  use cable_io_vars_module, only: landpt_global
  use cable_io_vars_module, only: patch
  use cable_io_vars_module, only: lat_all, lon_all
  use cable_io_vars_module, only: latitude, longitude

  use cable_output_types_mod, only: cable_output_variable_t
  use cable_output_types_mod, only: attribute => cable_output_attribute_t
  use cable_output_types_mod, only: DIM_PATCH => CABLE_OUTPUT_DIM_PATCH
  use cable_output_types_mod, only: DIM_SOIL => CABLE_OUTPUT_DIM_SOIL
  use cable_output_types_mod, only: DIM_RAD => CABLE_OUTPUT_DIM_RAD
  use cable_output_types_mod, only: DIM_LAND => CABLE_OUTPUT_DIM_LAND
  use cable_output_types_mod, only: DIM_LAND_GLOBAL => CABLE_OUTPUT_DIM_LAND_GLOBAL
  use cable_output_types_mod, only: DIM_X => CABLE_OUTPUT_DIM_X
  use cable_output_types_mod, only: DIM_Y => CABLE_OUTPUT_DIM_Y

  use cable_checks_module, only: ranges

  use cable_abort_module, only: cable_abort

  implicit none
  private

  public :: core_outputs
  public :: coordinate_variables

contains

  function core_outputs(canopy, soil) result(output_variables)
    type(canopy_type), intent(inout) :: canopy
    type(soil_parameter_type), intent(in) :: soil

    type(cable_output_variable_t), allocatable :: output_variables(:)

    output_variables = [ &
      cable_output_variable_t( &
        name="isoil", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_INT, &
        range=ranges%isoil, &
        active=output%isoil, &
        patchout=output%patch .or. patchout%isoil, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        restart=.true., &
        aggregator=new_aggregator(soil%isoilm), &
        metadata=[ &
          attribute("units", "-"), &
          attribute("long_name", "Soil type") &
        ] &
      ), &
      cable_output_variable_t( &
        name="swilt", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%swilt, &
        active=output%swilt, &
        patchout=output%patch .or. patchout%swilt, &
        reduction_method="grid_cell_average", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(soil%swilt), &
        metadata=[ &
          attribute("units", "1"), &
          attribute("long_name", "") &
        ] &
      ), &
      cable_output_variable_t( &
        name="albsoil", &
        data_shape=[DIM_PATCH, DIM_RAD], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%albsoil, &
        active=output%albsoil, &
        patchout=output%patch .or. patchout%albsoil, &
        reduction_method="first_patch_in_grid_cell", &
        aggregation_method="point", &
        parameter=.true., &
        aggregator=new_aggregator(soil%albsoil), &
        metadata=[ &
          attribute("units", "1"), &
          attribute("long_name", "") &
        ] &
      ), &
      cable_output_variable_t( &
        name="Qh", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=ranges%Qh, &
        active=output%Qh, &
        patchout=output%patch .or. patchout%Qh, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        aggregator=new_aggregator(canopy%fh), &
        metadata=[ &
           attribute("units", "W/m^2"), &
           attribute("long_name", "Surface sensible heat flux") &
        ] &
      ), &
      cable_output_variable_t( &
        name="Tmx", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        active=output%Tex .and. output%averaging == "monthly", &
        patchout=output%patch .or. patchout%Tex, &
        reduction_method="grid_cell_average", &
        aggregation_method="mean", &
        range=ranges%Tscrn, &
        accumulation_frequency="daily", &
        aggregator=new_aggregator(canopy%tscrn_max_daily%aggregated_data), &
        metadata=[ &
          attribute("units", "oC"), &
          attribute("long_name", "averaged daily maximum screen-level T") &
        ] &
      ), &
      cable_output_variable_t( &
        name="nap", &
        data_shape=[DIM_LAND_GLOBAL], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=[-huge(0.0), huge(0.0)], &
        active=.false., &
        restart=.true., &
        distributed=.false., &
        aggregation_method="point", &
        aggregator=new_aggregator(landpt_global(:)%nap), &
        metadata=[ &
          attribute("units", ""), &
          attribute("long_name", "") &
        ] &
      ), &
      cable_output_variable_t( &
        name="patchfrac", &
        data_shape=[DIM_PATCH], &
        var_type=CABLE_NETCDF_FLOAT, &
        range=[0.0, 1.0], &
        active=.false., &
        restart=.true., &
        distributed=.true., &
        aggregation_method="point", &
        aggregator=new_aggregator(patch(:)%frac), &
        metadata=[ &
          attribute("units", ""), &
          attribute("long_name", "Fraction of vegetated grid cell area occupied by a vegetation/soil patch") &
        ] &
      ), &
      cable_output_variable_t( &
        name="mvtype", &
        var_type=CABLE_NETCDF_FLOAT, &
        range=[-huge(0.0), huge(0.0)], &
        active=.false., &
        restart=.true., &
        distributed=.false., &
        aggregation_method="point", &
        aggregator=new_aggregator(mvtype), &
        metadata=[ &
          attribute("units", ""), &
          attribute("long_name", "Number of vegetation types") &
        ] &
      ), &
      cable_output_variable_t( &
        name="mstype", &
        var_type=CABLE_NETCDF_FLOAT, &
        range=[-huge(0.0), huge(0.0)], &
        active=.false., &
        restart=.true., &
        distributed=.false., &
        aggregation_method="point", &
        aggregator=new_aggregator(mstype), &
        metadata=[ &
          attribute("units", ""), &
          attribute("long_name", "Number of soil types") &
        ] &
      ) &
    ]

  end function core_outputs

  function coordinate_variables(grid_type) result(output_variables)
    character(len=*), intent(in) :: grid_type

    type(cable_output_variable_t), allocatable :: output_variables(:)

    select case (grid_type)
    case ("restart")
      output_variables = [ &
        cable_output_variable_t( &
          name="latitude", &
          data_shape=[DIM_LAND_GLOBAL], &
          var_type=CABLE_NETCDF_FLOAT, &
          range=[-90.0, 90.0], &
          active=.false., &
          restart=.true., &
          distributed=.false., &
          aggregation_method="point", &
          aggregator=new_aggregator(latitude), &
          metadata=[ &
            attribute("units", "degrees_north"), &
            attribute("long_name", "") &
          ] &
        ), &
        cable_output_variable_t( &
          name="longitude", &
          data_shape=[DIM_LAND_GLOBAL], &
          var_type=CABLE_NETCDF_FLOAT, &
          range=[-huge(0.0), huge(0.0)], & ! TODO(Sean): this depends on the met forcing input?
          active=.false., &
          restart=.true., &
          distributed=.false., &
          aggregation_method="point", &
          aggregator=new_aggregator(longitude), &
          metadata=[ &
            attribute("units", "degrees_east"), &
            attribute("long_name", "") &
          ] &
        ) &
      ]
    case ("mask")
      output_variables = [ &
        ! Define latitude and longitude variable (ALMA):
        cable_output_variable_t( &
          name="latitude", &
          data_shape=[DIM_X, DIM_Y], &
          var_type=CABLE_NETCDF_FLOAT, &
          range=[-90.0, 90.0], &
          active=.true., &
          parameter=.true., &
          distributed=.false., &
          aggregation_method="point", &
          aggregator=new_aggregator(lat_all), &
          metadata=[ &
            attribute("units", "degrees_north"), &
            attribute("long_name", "") &
          ] &
        ), &
        cable_output_variable_t( &
          name="longitude", &
          data_shape=[DIM_X, DIM_Y], &
          var_type=CABLE_NETCDF_FLOAT, &
          range=[-huge(0.0), huge(0.0)], & ! TODO(Sean): this depends on the met forcing input?
          active=.true., &
          parameter=.true., &
          distributed=.false., &
          aggregation_method="point", &
          aggregator=new_aggregator(lon_all), &
          metadata=[ &
           attribute("units", "degrees_east"), &
           attribute("long_name", "") &
          ] &
        ), &
        ! Write "cordinate variables" to enable reading by GrADS:
        cable_output_variable_t( &
          name="x", &
          data_shape=[DIM_X], &
          var_type=CABLE_NETCDF_FLOAT, &
          range=[-huge(0.0), huge(0.0)], & ! TODO(Sean): this depends on the met forcing input?
          active=.true., &
          parameter=.true., &
          distributed=.false., &
          aggregation_method="point", &
          aggregator=new_aggregator(lon_all(:, 1)), &
          metadata=[ &
            attribute("units", "degrees_east"), &
            attribute("long_name", ""), &
            attribute("comment", "x coordinate variable for GrADS compatibility") &
          ] &
        ), &
        cable_output_variable_t( &
          name="y", &
          data_shape=[DIM_Y], &
          var_type=CABLE_NETCDF_FLOAT, &
          range=[-90.0, 90.0], & ! TODO(Sean): this depends on the met forcing input?
          active=.true., &
          parameter=.true., &
          distributed=.false., &
          aggregation_method="point", &
          aggregator=new_aggregator(lat_all(1, :)), &
          metadata=[ &
            attribute("units", "degrees_north"), &
            attribute("long_name", ""), &
            attribute("comment", "y coordinate variable for GrADS compatibility") &
          ] &
        ) &
      ]
    case ("land")
      output_variables = [ &
        cable_output_variable_t( &
          name="local_lat", &
          data_shape=[DIM_LAND_GLOBAL], &
          var_type=CABLE_NETCDF_FLOAT, &
          range=[-90.0, 90.0], &
          active=.true., &
          parameter=.true., &
          distributed=.false., &
          aggregation_method="point", &
          aggregator=new_aggregator(latitude), &
          metadata=[ &
            attribute("units", "degrees_north"), &
            attribute("long_name", "") &
          ] &
        ), &
        cable_output_variable_t( &
          name="local_lon", &
          data_shape=[DIM_LAND_GLOBAL], &
          var_type=CABLE_NETCDF_FLOAT, &
          range=[-huge(0.0), huge(0.0)], & ! TODO(Sean): this depends on the met forcing input?
          active=.true., &
          parameter=.true., &
          distributed=.false., &
          aggregation_method="point", &
          aggregator=new_aggregator(longitude), &
          metadata=[ &
            attribute("units", "degrees_east"), &
            attribute("long_name", "") &
          ] &
        ) &
      ]
    case default
      call cable_abort("Unexpected grid type '" // grid_type // "'", __FILE__, __LINE__)
    end select

  end function

end module
