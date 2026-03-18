module cable_output_utils_mod

  use cable_common_module, only: filename

  use cable_def_types_mod, only: mp_global
  use cable_def_types_mod, only: mland_global

  use cable_io_vars_module, only: xdimsize
  use cable_io_vars_module, only: ydimsize
  use cable_io_vars_module, only: max_vegpatches
  use cable_io_vars_module, only: timeunits
  use cable_io_vars_module, only: time_coord
  use cable_io_vars_module, only: calendar

  use cable_error_handler_mod, only: cable_abort

  use cable_netcdf_mod, only: CABLE_NETCDF_UNLIMITED
  use cable_netcdf_mod, only: CABLE_NETCDF_INT
  use cable_netcdf_mod, only: CABLE_NETCDF_FLOAT
  use cable_netcdf_mod, only: CABLE_NETCDF_DOUBLE

  use cable_output_types_mod, only: cable_output_dim_t
  use cable_output_types_mod, only: cable_output_variable_t
  use cable_output_types_mod, only: cable_output_profile_t
  use cable_output_types_mod, only: FILL_VALUE_INT32
  use cable_output_types_mod, only: FILL_VALUE_REAL32
  use cable_output_types_mod, only: FILL_VALUE_REAL64
  use cable_output_types_mod, only: CABLE_OUTPUT_VAR_TYPE_UNDEFINED

  use aggregator_mod, only: new_aggregator

  implicit none
  private

  public :: check_duplicate_variable_names
  public :: check_sampling_frequency
  public :: var_name
  public :: define_variables
  public :: set_global_attributes

  integer, parameter :: MAX_LEN_DIM = 32

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

  elemental function var_name(output_variable)
    type(cable_output_variable_t), intent(in) :: output_variable
    character(len=128) :: var_name

    if (len_trim(output_variable%netcdf_name) > 0) then
      var_name = output_variable%netcdf_name
    else
      var_name = output_variable%field_name
    end if

  end function var_name

  subroutine check_duplicate_variable_names(output_profile)
    type(cable_output_profile_t), intent(in) :: output_profile
    integer :: i

    do i = 1, size(output_profile%output_variables)
      associate(output_var => output_profile%output_variables(i))
        if (count(var_name(output_var) == var_name(output_profile%output_variables(:))) > 1) then
          call cable_abort("Duplicate variable name found: " // trim(var_name(output_var)), __FILE__, __LINE__)
        end if
      end associate
    end do

  end subroutine

  subroutine check_sampling_frequency(output_profile)
    type(cable_output_profile_t), intent(in) :: output_profile

    integer :: i, sampling_period_in_hours, accumulation_period_in_hours

    character(len=256) :: err_message

    do i = 1, size(output_profile%output_variables)
      associate(output_var => output_profile%output_variables(i))
        err_message = ( &
          "Invalid combination of sampling frequency '" // output_profile%sampling_frequency // &
          "' with accumulation frequency '" // output_var%accumulation_frequency // "' for variable '" // &
          output_var%field_name // "' in file '" // output_profile%file_name // "'" &
        )
        select case (output_profile%sampling_frequency)
        case ("all")
          if (output_var%accumulation_frequency /= "all") call cable_abort(err_message, __FILE__, __LINE__)
        case ("user")
          read(output_profile%sampling_frequency(5:7), *) sampling_period_in_hours
          if (output_var%accumulation_frequency == "user") then
            read(output_var%accumulation_frequency(5:7), *) accumulation_period_in_hours
            if (sampling_period_in_hours < accumulation_period_in_hours) then
              call cable_abort(err_message, __FILE__, __LINE__)
            end if
          else if (output_var%accumulation_frequency /= "all") then
            call cable_abort(err_message, __FILE__, __LINE__)
          end if
        case ("daily")
          if (.not. any(output_var%accumulation_frequency == ["all  ", "daily", "user "])) then
            call cable_abort(err_message, __FILE__, __LINE__)
          end if
        case ("monthly")
          if (.not. any(output_var%accumulation_frequency == ["all    ", "daily  ", "user   ", "monthly"])) then
            call cable_abort(err_message, __FILE__, __LINE__)
          end if
        case default
          call cable_abort("Invalid sampling frequency '" // output_profile%sampling_frequency // &
            "' for variable '" // output_var%field_name // "' in file '" // output_profile%file_name // "'", __FILE__, __LINE__)
        end select
      end associate
    end do

  end subroutine check_sampling_frequency

  function infer_netcdf_dimensions(output_profile, output_variable, time_axis) result(netcdf_dimensions)
    type(cable_output_profile_t), intent(in) :: output_profile
    type(cable_output_variable_t), intent(in) :: output_variable
    logical, intent(in), optional :: time_axis

    type(cable_output_dim_t), allocatable :: netcdf_dimensions(:)
    integer :: i

    allocate(netcdf_dimensions(0))
    if (allocated(output_variable%data_shape)) then
      do i = 1, size(output_variable%data_shape)
        if (output_variable%data_shape(i)%name == "patch") then
          select case (output_profile%grid_type)
          case ("restart")
            netcdf_dimensions = [netcdf_dimensions, cable_output_dim_t("mp", mp_global)]
          case ("land")
            if (output_variable%reduction_method == "none") then
              netcdf_dimensions = [ &
                netcdf_dimensions, &
                cable_output_dim_t("land", mland_global), &
                cable_output_dim_t("patch", max_vegpatches) &
              ]
            else
              netcdf_dimensions = [ &
                netcdf_dimensions, &
                cable_output_dim_t("land", mland_global) &
              ]
            end if
          case ("mask")
            if (output_variable%reduction_method == "none") then
              netcdf_dimensions = [ &
                netcdf_dimensions, &
                cable_output_dim_t("x", xdimsize), &
                cable_output_dim_t("y", ydimsize), &
                cable_output_dim_t("patch", max_vegpatches) &
              ]
            else
              netcdf_dimensions = [ &
                netcdf_dimensions, &
                cable_output_dim_t("x", xdimsize), &
                cable_output_dim_t("y", ydimsize) &
              ]
            end if
          case default
            call cable_abort("Unable to determine output grid type.", __FILE__, __LINE__)
          end select
        else if (output_variable%data_shape(i)%name == "land") then
          select case (output_profile%grid_type)
          case ("restart")
            netcdf_dimensions = [netcdf_dimensions, cable_output_dim_t("mland", mland_global)]
          case ("land")
            netcdf_dimensions = [ &
              netcdf_dimensions, &
              cable_output_dim_t("land", mland_global) &
            ]
          case ("mask")
            netcdf_dimensions = [ &
              netcdf_dimensions, &
              cable_output_dim_t("x", xdimsize), &
              cable_output_dim_t("y", ydimsize) &
            ]
          case default
            call cable_abort("Unable to determine output grid type.", __FILE__, __LINE__)
          end select
        else if (output_variable%data_shape(i)%name == "land_global") then
          if (output_profile%grid_type == "restart") then
            netcdf_dimensions = [netcdf_dimensions, cable_output_dim_t("mland", mland_global)]
          else
            netcdf_dimensions = [netcdf_dimensions, cable_output_dim_t("land", mland_global)]
          end if
        else
          netcdf_dimensions = [netcdf_dimensions, output_variable%data_shape(i)]
        end if
      end do
    end if

    if (present(time_axis)) then; if (time_axis) then
        netcdf_dimensions = [netcdf_dimensions, cable_output_dim_t("time", CABLE_NETCDF_UNLIMITED)]
    end if; end if

  end function infer_netcdf_dimensions

  subroutine define_variables(output_profile, restart)
    type(cable_output_profile_t), intent(inout) :: output_profile
    logical, intent(in), optional :: restart

    logical :: restart_local
    integer :: i, j

    type(cable_output_variable_t), allocatable :: all_output_variables(:)
    type(cable_output_dim_t), allocatable :: required_dimensions(:), netcdf_dimensions(:)
    character(len=128) :: variable_name

    restart_local = .false.
    if (present(restart)) restart_local = restart

    all_output_variables = [ &
      output_profile%coordinate_variables, &
      output_profile%output_variables &
    ]

    do i = 1, size(all_output_variables)
      associate(output_var => all_output_variables(i))
        if (.not. allocated(output_var%data_shape)) cycle
        netcdf_dimensions = infer_netcdf_dimensions( &
          output_profile, &
          output_var, &
          time_axis=(.not. (restart_local .or. output_var%parameter)) &
        )
        if (.not. allocated(required_dimensions)) then
          required_dimensions = netcdf_dimensions
        else
          required_dimensions = [ &
            required_dimensions, &
            pack(netcdf_dimensions, [( &
              all(netcdf_dimensions(j)%name /= required_dimensions(:)%name), &
              j = 1, &
              size(netcdf_dimensions) &
            )]) &
          ]
        end if
      end associate
    end do

    do i = 1, size(required_dimensions)
      if (required_dimensions(i)%name == "time") cycle
      call output_profile%output_file%def_dims([required_dimensions(i)%name], [required_dimensions(i)%size])
    end do

    if (output_profile%grid_type == "restart") then
      call output_profile%output_file%def_dims(["time"], [1])
    else
      call output_profile%output_file%def_dims(["time"], [CABLE_NETCDF_UNLIMITED])
    end if

    call output_profile%output_file%def_var("time", CABLE_NETCDF_DOUBLE, ["time"])
    call output_profile%output_file%put_att("time", "units", timeunits)
    call output_profile%output_file%put_att("time", "coordinate", time_coord)
    call output_profile%output_file%put_att("time", "calendar", calendar)

    do i = 1, size(output_profile%coordinate_variables)
      associate(coord_var => output_profile%coordinate_variables(i))
        variable_name = var_name(coord_var)
        netcdf_dimensions = infer_netcdf_dimensions(output_profile, coord_var)
        call output_profile%output_file%def_var( &
          var_name=variable_name, &
          dim_names=netcdf_dimensions(:)%name, &
          type=netcdf_var_type(coord_var) &
        )
        if (allocated(coord_var%metadata)) then
          do j = 1, size(coord_var%metadata)
            call output_profile%output_file%put_att(variable_name, coord_var%metadata(j)%name, coord_var%metadata(j)%value)
          end do
        end if
      end associate
    end do

    do i = 1, size(output_profile%output_variables)
      associate(output_var => output_profile%output_variables(i))
        variable_name = var_name(output_var)
        if (restart_local) variable_name = output_var%field_name
        netcdf_dimensions = infer_netcdf_dimensions( &
          output_profile, &
          output_var, &
          time_axis=(.not. (restart_local .or. output_var%parameter)) &
        )
        call output_profile%output_file%def_var( &
          var_name=variable_name, &
          dim_names=netcdf_dimensions(:)%name, &
          type=netcdf_var_type(output_var, use_native_type=restart_local) &
        )
        if (allocated(output_var%metadata)) then
          do j = 1, size(output_var%metadata)
            call output_profile%output_file%put_att(variable_name, output_var%metadata(j)%name, output_var%metadata(j)%value)
          end do
        end if
        select case (netcdf_var_type(output_var, use_native_type=restart_local))
        case (CABLE_NETCDF_INT)
          call output_profile%output_file%put_att(variable_name, "_FillValue", FILL_VALUE_INT32)
          call output_profile%output_file%put_att(variable_name, "missing_value", FILL_VALUE_INT32)
        case (CABLE_NETCDF_FLOAT)
          call output_profile%output_file%put_att(variable_name, "_FillValue", FILL_VALUE_REAL32)
          call output_profile%output_file%put_att(variable_name, "missing_value", FILL_VALUE_REAL32)
        case (CABLE_NETCDF_DOUBLE)
          call output_profile%output_file%put_att(variable_name, "_FillValue", FILL_VALUE_REAL64)
          call output_profile%output_file%put_att(variable_name, "missing_value", FILL_VALUE_REAL64)
        end select
      end associate
    end do

  end subroutine define_variables

  subroutine set_global_attributes(output_profile)
    type(cable_output_profile_t), intent(inout) :: output_profile

    character(32) :: todaydate, nowtime
    integer :: i

    if (allocated(output_profile%metadata)) then
      do i = 1, size(output_profile%metadata)
        call output_profile%output_file%put_att( &
          att_name=output_profile%metadata(i)%name, &
          att_value=output_profile%metadata(i)%value &
        )
      end do
    end if

    call date_and_time(todaydate, nowtime)
    todaydate = todaydate(1:4) // "/" // todaydate(5:6) // "/" // todaydate(7:8)
    nowtime = nowtime(1:2) // ":" // nowtime(3:4) // ":" // nowtime(5:6)
    call output_profile%output_file%put_att("Production", trim(todaydate) // " at " // trim(nowtime))
    call output_profile%output_file%put_att("Source", "CABLE LSM output file")
    call output_profile%output_file%put_att("CABLE_input_file", trim(filename%met))

    select case (output_profile%sampling_frequency)
    case ("user")
       call output_profile%output_file%put_att("Output_averaging", TRIM(output_profile%sampling_frequency(5:7)) // "-hourly output")
    case ("all")
       call output_profile%output_file%put_att("Output_averaging", "all timesteps recorded")
    case ("daily")
       call output_profile%output_file%put_att("Output_averaging", "daily")
    case ("monthly")
       call output_profile%output_file%put_att("Output_averaging", "monthly")
    case default
       call cable_abort("Invalid sampling frequency '" // output_profile%sampling_frequency // "'", __FILE__, __LINE__)
    end select

  end subroutine set_global_attributes

end module
