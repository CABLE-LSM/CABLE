! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

submodule (cable_output_mod) cable_output_internal

  use cable_error_handler_mod, only: cable_abort
  use cable_netcdf_mod, only: cable_netcdf_decomp_t
  use cable_netcdf_mod, only: CABLE_NETCDF_FLOAT

  use aggregator_mod, only: new_aggregator
  use aggregator_mod, only: aggregator_int32_0d_t
  use aggregator_mod, only: aggregator_int32_1d_t
  use aggregator_mod, only: aggregator_int32_2d_t
  use aggregator_mod, only: aggregator_int32_3d_t
  use aggregator_mod, only: aggregator_real32_0d_t
  use aggregator_mod, only: aggregator_real32_1d_t
  use aggregator_mod, only: aggregator_real32_2d_t
  use aggregator_mod, only: aggregator_real32_3d_t
  use aggregator_mod, only: aggregator_real64_0d_t
  use aggregator_mod, only: aggregator_real64_1d_t
  use aggregator_mod, only: aggregator_real64_2d_t
  use aggregator_mod, only: aggregator_real64_3d_t

  use cable_checks_module, only: check_range

  use cable_io_vars_module, only: lat_all, lon_all
  use cable_io_vars_module, only: latitude, longitude

  implicit none

  interface
    module subroutine cable_output_define_variables(output_stream, restart)
      type(cable_output_stream_t), intent(inout) :: output_stream
      logical, intent(in), optional :: restart
    end subroutine
  end interface

  interface
    module subroutine cable_output_write_variable(output_stream, output_variable, patch, landpt, frame, restart)
      type(cable_output_stream_t), intent(inout) :: output_stream
      type(cable_output_variable_t), intent(inout), target :: output_variable
      type(patch_type), intent(in), optional :: patch(:)
      type(land_type), intent(in), optional :: landpt(:)
      integer, intent(in), optional :: frame
      logical, intent(in), optional :: restart
    end subroutine
  end interface

  interface
    module subroutine cable_output_decomp_init()
    end subroutine
    module subroutine cable_output_decomp_free()
    end subroutine
    module subroutine cable_output_decomp_associate(output_stream, output_var, decomp)
      type(cable_output_stream_t), intent(in) :: output_stream
      type(cable_output_variable_t), intent(in) :: output_var
      class(cable_netcdf_decomp_t), pointer, intent(inout) :: decomp
    end subroutine
    module subroutine cable_output_reduction_buffers_init()
    end subroutine
    module subroutine cable_output_reduction_buffers_free()
    end subroutine
  end interface

  interface cable_output_reduction_buffers_associate
    module subroutine cable_output_reduction_buffers_associate_1d_int32(output_var, temp_buffer)
      type(cable_output_variable_t), intent(inout) :: output_var
      integer(kind=int32), pointer, intent(inout) :: temp_buffer(:)
    end subroutine
    module subroutine cable_output_reduction_buffers_associate_1d_real32(output_var, temp_buffer)
      type(cable_output_variable_t), intent(inout) :: output_var
      real(kind=real32), pointer, intent(inout) :: temp_buffer(:)
    end subroutine
    module subroutine cable_output_reduction_buffers_associate_1d_real64(output_var, temp_buffer)
      type(cable_output_variable_t), intent(inout) :: output_var
      real(kind=real64), pointer, intent(inout) :: temp_buffer(:)
    end subroutine
    module subroutine cable_output_reduction_buffers_associate_2d_int32(output_var, temp_buffer)
      type(cable_output_variable_t), intent(inout) :: output_var
      integer(kind=int32), pointer, intent(inout) :: temp_buffer(:,:)
    end subroutine
    module subroutine cable_output_reduction_buffers_associate_2d_real32(output_var, temp_buffer)
      type(cable_output_variable_t), intent(inout) :: output_var
      real(kind=real32), pointer, intent(inout) :: temp_buffer(:,:)
    end subroutine
    module subroutine cable_output_reduction_buffers_associate_2d_real64(output_var, temp_buffer)
      type(cable_output_variable_t), intent(inout) :: output_var
      real(kind=real64), pointer, intent(inout) :: temp_buffer(:,:)
    end subroutine
    module subroutine cable_output_reduction_buffers_associate_3d_int32(output_var, temp_buffer)
      type(cable_output_variable_t), intent(inout) :: output_var
      integer(kind=int32), pointer, intent(inout) :: temp_buffer(:,:,:)
    end subroutine
    module subroutine cable_output_reduction_buffers_associate_3d_real32(output_var, temp_buffer)
      type(cable_output_variable_t), intent(inout) :: output_var
      real(kind=real32), pointer, intent(inout) :: temp_buffer(:,:,:)
    end subroutine
    module subroutine cable_output_reduction_buffers_associate_3d_real64(output_var, temp_buffer)
      type(cable_output_variable_t), intent(inout) :: output_var
      real(kind=real64), pointer, intent(inout) :: temp_buffer(:,:,:)
    end subroutine
  end interface

contains

  logical function frequency_is_greater_than(freq_a, freq_b) result(freq_a_greater_than_b)
    character(len=*), intent(in) :: freq_a
    character(len=*), intent(in) :: freq_b

    integer :: period_in_hours_a, period_in_hours_b

    select case (freq_a)
    case ("all")
      freq_a_greater_than_b = freq_b /= "all"
    case ("user")
      read(freq_a(5:7), *) period_in_hours_a
      if (freq_b == "user") then
        read(freq_b(5:7), *) period_in_hours_b
        freq_a_greater_than_b = period_in_hours_a < period_in_hours_b
      else
        freq_a_greater_than_b = freq_b /= "all"
      end if
    case ("daily")
      freq_a_greater_than_b = freq_b == "monthly"
    case ("monthly")
      freq_a_greater_than_b = .false.
    case default
      call cable_abort("Unexpected sampling frequency '" // freq_a, __FILE__, __LINE__)
    end select

  end function frequency_is_greater_than

  subroutine check_variable_range(output_variable, time_index, met)
    type(cable_output_variable_t), intent(in) :: output_variable
    integer, intent(in) :: time_index
    type(met_type), intent(in), optional :: met

    select type (aggregator => output_variable%aggregator)
    type is (aggregator_int32_0d_t)
      call check_range(output_variable%field_name, [real(aggregator%source_data)], output_variable%range, time_index, met)
    type is (aggregator_int32_1d_t)
      call check_range(output_variable%field_name, real(aggregator%source_data), output_variable%range, time_index, met)
    type is (aggregator_int32_2d_t)
      call check_range(output_variable%field_name, real(aggregator%source_data), output_variable%range, time_index, met)
    type is (aggregator_int32_3d_t)
      call check_range(output_variable%field_name, real(aggregator%source_data), output_variable%range, time_index, met)
    type is (aggregator_real32_0d_t)
      call check_range(output_variable%field_name, [aggregator%source_data], output_variable%range, time_index, met)
    type is (aggregator_real32_1d_t)
      call check_range(output_variable%field_name, aggregator%source_data, output_variable%range, time_index, met)
    type is (aggregator_real32_2d_t)
      call check_range(output_variable%field_name, aggregator%source_data, output_variable%range, time_index, met)
    type is (aggregator_real32_3d_t)
      call check_range(output_variable%field_name, aggregator%source_data, output_variable%range, time_index, met)
    type is (aggregator_real64_0d_t)
      call check_range(output_variable%field_name, [real(aggregator%source_data)], output_variable%range, time_index, met)
    type is (aggregator_real64_1d_t)
      call check_range(output_variable%field_name, real(aggregator%source_data), output_variable%range, time_index, met)
    type is (aggregator_real64_2d_t)
      call check_range(output_variable%field_name, real(aggregator%source_data), output_variable%range, time_index, met)
    type is (aggregator_real64_3d_t)
      call check_range(output_variable%field_name, real(aggregator%source_data), output_variable%range, time_index, met)
    class default
      call cable_abort("Unexpected aggregator type", __FILE__, __LINE__)
    end select

  end subroutine check_variable_range

  function coordinate_variables_list(grid_type) result(coord_variables)
    character(len=*), intent(in) :: grid_type

    type(cable_output_variable_t), allocatable :: coord_variables(:)
    type(cable_output_variable_t), allocatable :: mask_coord_variables(:)

    type(cable_output_dim_t) :: dim_x, dim_y, dim_land_global

    dim_x = cable_output_get_dimension("x")
    dim_y = cable_output_get_dimension("y")
    dim_land_global = cable_output_get_dimension("land_global")

    mask_coord_variables = [ &
      cable_output_variable_t( &
        field_name="lat_all", &
        netcdf_name="latitude", &
        data_shape=[dim_x, dim_y], &
        var_type=CABLE_NETCDF_FLOAT, &
        parameter=.true., &
        distributed=.false., &
        aggregator=new_aggregator(lat_all), &
        metadata=[cable_output_attribute_t("units", "degrees_north")] &
      ), &
      cable_output_variable_t( &
        field_name="lon_all", &
        netcdf_name="longitude", &
        data_shape=[dim_x, dim_y], &
        parameter=.true., &
        distributed=.false., &
        aggregator=new_aggregator(lon_all), &
        metadata=[cable_output_attribute_t("units", "degrees_east")] &
      ), &
      cable_output_variable_t( &
        field_name="x", &
        data_shape=[dim_x], &
        parameter=.true., &
        distributed=.false., &
        aggregator=new_aggregator(lon_all(:, 1)), &
        metadata=[ &
          cable_output_attribute_t("units", "degrees_east"), &
          cable_output_attribute_t("comment", "x coordinate variable for GrADS compatibility") &
        ] &
      ), &
      cable_output_variable_t( &
        field_name="y", &
        data_shape=[dim_y], &
        parameter=.true., &
        distributed=.false., &
        aggregator=new_aggregator(lat_all(1, :)), &
        metadata=[ &
          cable_output_attribute_t("units", "degrees_north"), &
          cable_output_attribute_t("comment", "y coordinate variable for GrADS compatibility") &
        ] &
      ) &
    ]

    select case (grid_type)
    case ("restart")
      coord_variables = [ &
        cable_output_variable_t( &
          field_name="latitude", &
          data_shape=[dim_land_global], &
          distributed=.false., &
          aggregator=new_aggregator(latitude), &
          metadata=[cable_output_attribute_t("units", "degrees_north")] &
        ), &
        cable_output_variable_t( &
          field_name="longitude", &
          data_shape=[dim_land_global], &
          distributed=.false., &
          aggregator=new_aggregator(longitude), &
          metadata=[cable_output_attribute_t("units", "degrees_east")] &
        ) &
      ]
    case ("mask")
      coord_variables = mask_coord_variables
    case ("land")
      coord_variables = [ &
        mask_coord_variables, &
        cable_output_variable_t( &
          field_name="local_lat", &
          data_shape=[dim_land_global], &
          parameter=.true., &
          distributed=.false., &
          aggregator=new_aggregator(latitude), &
          metadata=[cable_output_attribute_t("units", "degrees_north")] &
        ), &
        cable_output_variable_t( &
          field_name="local_lon", &
          data_shape=[dim_land_global], &
          parameter=.true., &
          distributed=.false., &
          aggregator=new_aggregator(longitude), &
          metadata=[cable_output_attribute_t("units", "degrees_east")] &
        ) &
      ]
    case default
      call cable_abort("Unexpected grid type '" // grid_type // "'", __FILE__, __LINE__)
    end select

  end function coordinate_variables_list

end submodule cable_output_internal
