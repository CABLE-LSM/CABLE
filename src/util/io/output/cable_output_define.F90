submodule (cable_output_mod:cable_output_internal) cable_output_define

  use cable_common_module, only: filename

  use cable_def_types_mod, only: mp_global
  use cable_def_types_mod, only: mland_global

  use cable_io_vars_module, only: xdimsize
  use cable_io_vars_module, only: ydimsize
  use cable_io_vars_module, only: max_vegpatches
  use cable_io_vars_module, only: timeunits
  use cable_io_vars_module, only: time_coord
  use cable_io_vars_module, only: calendar

  use cable_netcdf_mod, only: CABLE_NETCDF_UNLIMITED
  use cable_netcdf_mod, only: CABLE_NETCDF_INT
  use cable_netcdf_mod, only: CABLE_NETCDF_FLOAT
  use cable_netcdf_mod, only: CABLE_NETCDF_DOUBLE

  implicit none

contains

  integer function netcdf_var_type(output_variable, use_native_type)
    type(cable_output_variable_t), intent(in) :: output_variable
    logical, intent(in), optional :: use_native_type
    logical :: native_type

    native_type = .false.
    if (present(use_native_type)) native_type = use_native_type

    if (.not. native_type .and. output_variable%var_type /= CABLE_OUTPUT_VAR_TYPE_UNDEFINED) then
      netcdf_var_type = output_variable%var_type
      return
    end if

    select case (output_variable%aggregator%type())
    case ("int32")
      netcdf_var_type = CABLE_NETCDF_INT
    case ("real32")
      netcdf_var_type = CABLE_NETCDF_FLOAT
    case ("real64")
      netcdf_var_type = CABLE_NETCDF_DOUBLE
    case default
      call cable_abort("Unable to infer variable type for variable " // trim(output_variable%field_name), __FILE__, __LINE__)
    end select

  end function

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


  function infer_netcdf_dimensions(output_stream, output_variable, time_axis) result(netcdf_dimensions)
    type(cable_output_stream_t), intent(in) :: output_stream
    type(cable_output_variable_t), intent(in) :: output_variable
    logical, intent(in), optional :: time_axis

    type(cable_output_dim_t), allocatable :: netcdf_dimensions(:)
    integer :: i

    allocate(netcdf_dimensions(0))
    if (allocated(output_variable%data_shape)) then
      netcdf_dimensions = [( &
        native_to_netcdf_dimensions( &
          native_dimension=output_variable%data_shape(i), &
          grid_type=output_stream%grid_type, &
          reduction_method=output_variable%reduction_method &
        ), &
        i = 1, size(output_variable%data_shape) &
      )]
    end if

    if (present(time_axis)) then; if (time_axis) then
        netcdf_dimensions = [netcdf_dimensions, cable_output_dim_t("time", CABLE_NETCDF_UNLIMITED)]
    end if; end if

  end function infer_netcdf_dimensions

  module subroutine cable_output_define_variables(output_stream, restart)
    type(cable_output_stream_t), intent(inout) :: output_stream
    logical, intent(in), optional :: restart

    logical :: restart_local
    integer :: i, j

    type(cable_output_variable_t), allocatable :: all_output_variables(:)
    type(cable_output_dim_t), allocatable :: required_dimensions(:), netcdf_dimensions(:)
    character(len=128) :: variable_name

    restart_local = .false.
    if (present(restart)) restart_local = restart

    all_output_variables = [ &
      output_stream%coordinate_variables, &
      output_stream%output_variables &
    ]

    do i = 1, size(all_output_variables)
      associate(output_var => all_output_variables(i))
        if (.not. allocated(output_var%data_shape)) cycle
        netcdf_dimensions = infer_netcdf_dimensions( &
          output_stream, &
          output_var, &
          time_axis=(.not. (restart_local .or. output_var%parameter)) &
        )
        if (.not. allocated(required_dimensions)) then
          required_dimensions = netcdf_dimensions
        else
          required_dimensions = [ &
            required_dimensions, &
            pack(netcdf_dimensions, [( &
              all(netcdf_dimensions(j)%name() /= required_dimensions(:)%name()), &
              j = 1, &
              size(netcdf_dimensions) &
            )]) &
          ]
        end if
      end associate
    end do

    do i = 1, size(required_dimensions)
      if (required_dimensions(i)%name() == "time") cycle
      call output_stream%output_file%def_dims([required_dimensions(i)%name()], [required_dimensions(i)%size()])
    end do

    if (output_stream%grid_type == "restart") then
      call output_stream%output_file%def_dims(["time"], [1])
    else
      call output_stream%output_file%def_dims(["time"], [CABLE_NETCDF_UNLIMITED])
    end if

    call output_stream%output_file%def_var("time", CABLE_NETCDF_DOUBLE, ["time"])
    call output_stream%output_file%put_att("time", "units", timeunits)
    call output_stream%output_file%put_att("time", "coordinate", time_coord)
    call output_stream%output_file%put_att("time", "calendar", calendar)

    do i = 1, size(output_stream%coordinate_variables)
      associate(coord_var => output_stream%coordinate_variables(i))
        variable_name = coord_var%get_netcdf_name()
        netcdf_dimensions = infer_netcdf_dimensions(output_stream, coord_var)
        call output_stream%output_file%def_var( &
          var_name=variable_name, &
          dim_names=netcdf_dimensions(:)%name(), &
          type=netcdf_var_type(coord_var) &
        )
        if (allocated(coord_var%metadata)) then
          do j = 1, size(coord_var%metadata)
            call output_stream%output_file%put_att(variable_name, coord_var%metadata(j)%name, coord_var%metadata(j)%value)
          end do
        end if
      end associate
    end do

    do i = 1, size(output_stream%output_variables)
      associate(output_var => output_stream%output_variables(i))
        variable_name = output_var%get_netcdf_name()
        if (restart_local) variable_name = output_var%field_name
        netcdf_dimensions = infer_netcdf_dimensions( &
          output_stream, &
          output_var, &
          time_axis=(.not. (restart_local .or. output_var%parameter)) &
        )
        call output_stream%output_file%def_var( &
          var_name=variable_name, &
          dim_names=netcdf_dimensions(:)%name(), &
          type=netcdf_var_type(output_var, use_native_type=restart_local) &
        )
        if (allocated(output_var%metadata)) then
          do j = 1, size(output_var%metadata)
            call output_stream%output_file%put_att(variable_name, output_var%metadata(j)%name, output_var%metadata(j)%value)
          end do
        end if
        select case (netcdf_var_type(output_var, use_native_type=restart_local))
        case (CABLE_NETCDF_INT)
          call output_stream%output_file%put_att(variable_name, "_FillValue", CABLE_OUTPUT_FILL_VALUE_INT32)
          call output_stream%output_file%put_att(variable_name, "missing_value", CABLE_OUTPUT_FILL_VALUE_INT32)
        case (CABLE_NETCDF_FLOAT)
          call output_stream%output_file%put_att(variable_name, "_FillValue", CABLE_OUTPUT_FILL_VALUE_REAL32)
          call output_stream%output_file%put_att(variable_name, "missing_value", CABLE_OUTPUT_FILL_VALUE_REAL32)
        case (CABLE_NETCDF_DOUBLE)
          call output_stream%output_file%put_att(variable_name, "_FillValue", CABLE_OUTPUT_FILL_VALUE_REAL64)
          call output_stream%output_file%put_att(variable_name, "missing_value", CABLE_OUTPUT_FILL_VALUE_REAL64)
        end select
      end associate
    end do

    if (.not. restart_local) call set_global_attributes(output_stream)

  end subroutine cable_output_define_variables

  subroutine set_global_attributes(output_stream)
    type(cable_output_stream_t), intent(inout) :: output_stream

    character(32) :: todaydate, nowtime
    integer :: i

    if (allocated(output_stream%metadata)) then
      do i = 1, size(output_stream%metadata)
        call output_stream%output_file%put_att( &
          att_name=output_stream%metadata(i)%name, &
          att_value=output_stream%metadata(i)%value &
        )
      end do
    end if

    call date_and_time(todaydate, nowtime)
    todaydate = todaydate(1:4) // "/" // todaydate(5:6) // "/" // todaydate(7:8)
    nowtime = nowtime(1:2) // ":" // nowtime(3:4) // ":" // nowtime(5:6)
    call output_stream%output_file%put_att("Production", trim(todaydate) // " at " // trim(nowtime))
    call output_stream%output_file%put_att("Source", "CABLE LSM output file")
    call output_stream%output_file%put_att("CABLE_input_file", trim(filename%met))

    select case (output_stream%sampling_frequency)
    case ("user")
       call output_stream%output_file%put_att("Output_averaging", TRIM(output_stream%sampling_frequency(5:7)) // "-hourly output")
    case ("all")
       call output_stream%output_file%put_att("Output_averaging", "all timesteps recorded")
    case ("daily")
       call output_stream%output_file%put_att("Output_averaging", "daily")
    case ("monthly")
       call output_stream%output_file%put_att("Output_averaging", "monthly")
    case default
       call cable_abort("Invalid sampling frequency '" // output_stream%sampling_frequency // "'", __FILE__, __LINE__)
    end select

  end subroutine set_global_attributes

end submodule cable_output_define
