! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

submodule (cable_output_mod) cable_output_internal_smod
  !* Internal interfaces and procedures for [[cable_output_mod]].
  !
  ! This module declares interfaces for the procedures that are used by
  ! [[cable_output_impl]], as well as various utilities used in other parts of
  ! the output system.

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
    !! Interfaces for procedures used by [[cable_output_impl]].

    module subroutine cable_output_decomp_init()
      !! Intialises I/O decompositions used in the output system.
    end subroutine

    module subroutine cable_output_decomp_free()
      !! Deallocates I/O decompositions used in the output system.
    end subroutine

    module subroutine cable_output_decomp_associate(output_stream, output_var, decomp)
      !* Associates an I/O decomposition pointer with the appropriate I/O
      ! decomposition, taking into account the output variable shape and type, and
      ! the output stream grid type.
      type(cable_output_stream_t), intent(in) :: output_stream
        !! The output stream for which to associate the decomposition.
      type(cable_output_variable_t), intent(in) :: output_var
        !! The output variable for which to associate the decomposition.
      class(cable_netcdf_decomp_t), pointer, intent(inout) :: decomp
        !! The decomposition pointer to associate.
    end subroutine

    module subroutine cable_output_define_stream(output_stream, restart)
      !* Defines all variables, dimensions and attributes for a given output
      ! stream.
      type(cable_output_stream_t), intent(inout) :: output_stream
        !! The output stream to define.
      logical, intent(in), optional :: restart
        !* Whether this is a restart stream definition. Set to `.false.` by
        ! default.
    end subroutine

    module subroutine cable_output_reduction_buffers_init()
      !! Initialises the buffers used for performing grid reductions in the output system.
    end subroutine

    module subroutine cable_output_reduction_buffers_free()
      !! Deallocates the buffers used for performing grid reductions in the output system.
    end subroutine

    module subroutine cable_output_write_variable(output_stream, output_variable, patch, landpt, frame, restart)
      !! Writes a variable to the output stream.
      type(cable_output_stream_t), intent(inout) :: output_stream !! The output stream to write to.
      type(cable_output_variable_t), intent(inout), target :: output_variable !! The variable to write.
      type(patch_type), intent(in), optional :: patch(:)
        !! The patch type instance for performing grid reductions over the patch dimension if required.
      type(land_type), intent(in), optional :: landpt(:)
        !! The land type instance for performing grid reductions over the patch dimension if required.
      integer, intent(in), optional :: frame !! The frame or unlimited dimension index to write at.
      logical, intent(in), optional :: restart !! Whether this is a restart stream write.
    end subroutine

  end interface

  interface cable_output_reduction_buffers_associate
    !* Interface for associating a pointer array with the the appropriate
    ! reduction buffer, taking into account the output variable shape, type and
    ! reduction method.
    module subroutine cable_output_reduction_buffers_associate_1d_int32(output_var, temp_buffer)
      !! The reduction buffer association subroutine for 1D 32-bit integer variables.
      type(cable_output_variable_t), intent(inout) :: output_var
        !! The output variable for which to associate the reduction buffer.
      integer(kind=int32), pointer, intent(inout) :: temp_buffer(:)
        !! The pointer array to associate with the appropriate reduction buffer.
    end subroutine
    module subroutine cable_output_reduction_buffers_associate_1d_real32(output_var, temp_buffer)
      !! The reduction buffer association subroutine for 1D 32-bit real variables.
      type(cable_output_variable_t), intent(inout) :: output_var
        !! The output variable for which to associate the reduction buffer.
      real(kind=real32), pointer, intent(inout) :: temp_buffer(:)
        !! The pointer array to associate with the appropriate reduction buffer.
    end subroutine
    module subroutine cable_output_reduction_buffers_associate_1d_real64(output_var, temp_buffer)
      !! The reduction buffer association subroutine for 1D 64-bit real variables.
      type(cable_output_variable_t), intent(inout) :: output_var
        !! The output variable for which to associate the reduction buffer.
      real(kind=real64), pointer, intent(inout) :: temp_buffer(:)
        !! The pointer array to associate with the appropriate reduction buffer.
    end subroutine
    module subroutine cable_output_reduction_buffers_associate_2d_int32(output_var, temp_buffer)
      !! The reduction buffer association subroutine for 2D 32-bit integer variables.
      type(cable_output_variable_t), intent(inout) :: output_var
        !! The output variable for which to associate the reduction buffer.
      integer(kind=int32), pointer, intent(inout) :: temp_buffer(:,:)
        !! The pointer array to associate with the appropriate reduction buffer.
    end subroutine
    module subroutine cable_output_reduction_buffers_associate_2d_real32(output_var, temp_buffer)
      !! The reduction buffer association subroutine for 2D 32-bit real variables.
      type(cable_output_variable_t), intent(inout) :: output_var
        !! The output variable for which to associate the reduction buffer.
      real(kind=real32), pointer, intent(inout) :: temp_buffer(:,:)
        !! The pointer array to associate with the appropriate reduction buffer.
    end subroutine
    module subroutine cable_output_reduction_buffers_associate_2d_real64(output_var, temp_buffer)
      !! The reduction buffer association subroutine for 2D 64-bit real variables.
      type(cable_output_variable_t), intent(inout) :: output_var
        !! The output variable for which to associate the reduction buffer.
      real(kind=real64), pointer, intent(inout) :: temp_buffer(:,:)
        !! The pointer array to associate with the appropriate reduction buffer.
    end subroutine
    module subroutine cable_output_reduction_buffers_associate_3d_int32(output_var, temp_buffer)
      !! The reduction buffer association subroutine for 3D 32-bit integer variables.
      type(cable_output_variable_t), intent(inout) :: output_var
        !! The output variable for which to associate the reduction buffer.
      integer(kind=int32), pointer, intent(inout) :: temp_buffer(:,:,:)
        !! The pointer array to associate with the appropriate reduction buffer.
    end subroutine
    module subroutine cable_output_reduction_buffers_associate_3d_real32(output_var, temp_buffer)
      !! The reduction buffer association subroutine for 3D 32-bit real variables.
      type(cable_output_variable_t), intent(inout) :: output_var
        !! The output variable for which to associate the reduction buffer.
      real(kind=real32), pointer, intent(inout) :: temp_buffer(:,:,:)
        !! The pointer array to associate with the appropriate reduction buffer.
    end subroutine
    module subroutine cable_output_reduction_buffers_associate_3d_real64(output_var, temp_buffer)
      !! The reduction buffer association subroutine for 3D 64-bit real variables.
      type(cable_output_variable_t), intent(inout) :: output_var
        !! The output variable for which to associate the reduction buffer.
      real(kind=real64), pointer, intent(inout) :: temp_buffer(:,:,:)
        !! The pointer array to associate with the appropriate reduction buffer.
    end subroutine
  end interface

contains

  function native_to_netcdf_dimensions(native_dimension, grid_type, reduction_method) result(netcdf_dimensions)
    !* Returns the netCDF dimension(s) corresponding to a given output
    ! variable dimension, taking into account the output grid type and reduction
    ! method. This function is used to determine the dimensions of netCDF
    ! variables based on the in-memory data shapes of CABLE variables as
    ! described by `cable_output_dim_t` instances.
    type(cable_output_dim_t), intent(in) :: native_dimension
      !! The in-memory dimension.
    character(len=*), intent(in) :: grid_type
      !* The output grid type. See [[allowed_grid_types]] for the available
      ! grid types.
    character(len=*), intent(in) :: reduction_method
      !* The reduction method applied to the variable. See
      ! [[allowed_reduction_methods]] for the available reduction methods.
    type(cable_output_dim_t), allocatable :: netcdf_dimensions(:)

    select case (native_dimension%name())
    case (NATIVE_DIM_NAME_PATCH)
      select case (grid_type)
      case ("restart")
        netcdf_dimensions = [cable_output_dim_t("mp", mp_global)]
      case ("land")
        if (reduction_method == "none") then
          netcdf_dimensions = [ &
            cable_output_dim_t("land", mland_global), &
            cable_output_dim_t("patch", max_vegpatches) &
          ]
        else
          netcdf_dimensions = [cable_output_dim_t("land", mland_global)]
        end if
      case ("mask")
        if (reduction_method == "none") then
          netcdf_dimensions = [ &
            cable_output_dim_t("x", xdimsize), &
            cable_output_dim_t("y", ydimsize), &
            cable_output_dim_t("patch", max_vegpatches) &
          ]
        else
          netcdf_dimensions = [ &
            cable_output_dim_t("x", xdimsize), &
            cable_output_dim_t("y", ydimsize) &
          ]
        end if
      case default
        call cable_abort("Unable to determine output grid type.", __FILE__, __LINE__)
      end select
    case (NATIVE_DIM_NAME_PATCH_GLOBAL)
      netcdf_dimensions = [cable_output_dim_t("mp", mp_global)]
    case (NATIVE_DIM_NAME_PATCH_GRID_CELL)
      netcdf_dimensions = [cable_output_dim_t("patch", max_vegpatches)]
    case (NATIVE_DIM_NAME_LAND)
      select case (grid_type)
      case ("restart")
        netcdf_dimensions = [cable_output_dim_t("mland", mland_global)]
      case ("land")
        netcdf_dimensions = [cable_output_dim_t("land", mland_global)]
      case ("mask")
        netcdf_dimensions = [ &
          cable_output_dim_t("x", xdimsize), &
          cable_output_dim_t("y", ydimsize) &
        ]
      case default
        call cable_abort("Unable to determine output grid type.", __FILE__, __LINE__)
      end select
    case (NATIVE_DIM_NAME_LAND_GLOBAL)
      if (grid_type == "restart") then
        netcdf_dimensions = [cable_output_dim_t("mland", mland_global)]
      else
        netcdf_dimensions = [cable_output_dim_t("land", mland_global)]
      end if
    case default
      netcdf_dimensions = [native_dimension]
    end select

  end function native_to_netcdf_dimensions

  subroutine check_variable_range(output_variable, time_index, met)
    !* Checks whether the value(s) of an output variable are within their
    ! specified range of physical values.
    !
    ! Note that range checks are done on the native diagnostic (not to be
    ! confused with the netCDF variable which may have different units).
    type(cable_output_variable_t), intent(in) :: output_variable
      !! The output variable for which to check the range.
    integer, intent(in) :: time_index
      !! The current time step index, used for error messages.
    type(met_type), intent(in), optional :: met
      !! The met_type instance containing the current meteorological conditions, used for error messages.

    select type (aggregator => output_variable%aggregator)
    type is (aggregator_int32_0d_t)
      call check_range(output_variable%field_name, [real(aggregator%source_data)], output_variable%range_native, time_index, met)
    type is (aggregator_int32_1d_t)
      call check_range(output_variable%field_name, real(aggregator%source_data), output_variable%range_native, time_index, met)
    type is (aggregator_int32_2d_t)
      call check_range(output_variable%field_name, real(aggregator%source_data), output_variable%range_native, time_index, met)
    type is (aggregator_int32_3d_t)
      call check_range(output_variable%field_name, real(aggregator%source_data), output_variable%range_native, time_index, met)
    type is (aggregator_real32_0d_t)
      call check_range(output_variable%field_name, [aggregator%source_data], output_variable%range_native, time_index, met)
    type is (aggregator_real32_1d_t)
      call check_range(output_variable%field_name, aggregator%source_data, output_variable%range_native, time_index, met)
    type is (aggregator_real32_2d_t)
      call check_range(output_variable%field_name, aggregator%source_data, output_variable%range_native, time_index, met)
    type is (aggregator_real32_3d_t)
      call check_range(output_variable%field_name, aggregator%source_data, output_variable%range_native, time_index, met)
    type is (aggregator_real64_0d_t)
      call check_range(output_variable%field_name, [real(aggregator%source_data)], output_variable%range_native, time_index, met)
    type is (aggregator_real64_1d_t)
      call check_range(output_variable%field_name, real(aggregator%source_data), output_variable%range_native, time_index, met)
    type is (aggregator_real64_2d_t)
      call check_range(output_variable%field_name, real(aggregator%source_data), output_variable%range_native, time_index, met)
    type is (aggregator_real64_3d_t)
      call check_range(output_variable%field_name, real(aggregator%source_data), output_variable%range_native, time_index, met)
    class default
      call cable_abort("Unexpected aggregator type", __FILE__, __LINE__)
    end select

  end subroutine check_variable_range

  function coordinate_variables_list(grid_type) result(coord_variables)
    !* Returns a list of coordinate variables to be included in an output stream
    ! based on the output grid type.
    character(len=*), intent(in) :: grid_type
      !! The output grid type. See [[allowed_grid_types]] for the available grid types.

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

end submodule cable_output_internal_smod
