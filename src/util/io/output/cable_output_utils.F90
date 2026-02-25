module cable_output_utils_mod

  use cable_common_module, only: filename

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
  use cable_io_vars_module, only: timeunits
  use cable_io_vars_module, only: time_coord
  use cable_io_vars_module, only: calendar
  use cable_io_vars_module, only: lat_all, lon_all
  use cable_io_vars_module, only: latitude, longitude

  use cable_abort_module, only: cable_abort

  use cable_netcdf_mod, only: MAX_LEN_DIM => CABLE_NETCDF_MAX_STR_LEN_DIM
  use cable_netcdf_mod, only: CABLE_NETCDF_UNLIMITED
  use cable_netcdf_mod, only: CABLE_NETCDF_INT
  use cable_netcdf_mod, only: CABLE_NETCDF_FLOAT
  use cable_netcdf_mod, only: CABLE_NETCDF_DOUBLE

  use cable_output_types_mod, only: cable_output_dim_t
  use cable_output_types_mod, only: cable_output_variable_t
  use cable_output_types_mod, only: cable_output_profile_t
  use cable_output_types_mod, only: attribute => cable_output_attribute_t
  use cable_output_types_mod, only: CABLE_OUTPUT_DIM_PATCH
  use cable_output_types_mod, only: CABLE_OUTPUT_DIM_SOIL
  use cable_output_types_mod, only: CABLE_OUTPUT_DIM_SNOW
  use cable_output_types_mod, only: CABLE_OUTPUT_DIM_RAD
  use cable_output_types_mod, only: CABLE_OUTPUT_DIM_PLANTCARBON
  use cable_output_types_mod, only: CABLE_OUTPUT_DIM_SOILCARBON
  use cable_output_types_mod, only: CABLE_OUTPUT_DIM_LAND
  use cable_output_types_mod, only: CABLE_OUTPUT_DIM_LAND_GLOBAL
  use cable_output_types_mod, only: CABLE_OUTPUT_DIM_X
  use cable_output_types_mod, only: CABLE_OUTPUT_DIM_Y
  use cable_output_types_mod, only: FILL_VALUE_INT32
  use cable_output_types_mod, only: FILL_VALUE_REAL32
  use cable_output_types_mod, only: FILL_VALUE_REAL64
  use cable_output_types_mod, only: CABLE_OUTPUT_VAR_TYPE_UNDEFINED

  use aggregator_mod, only: new_aggregator

  implicit none
  private

  public :: coordinate_variables
  public :: check_duplicate_variable_names
  public :: check_sampling_frequency
  public :: dim_size
  public :: var_name
  public :: infer_dim_names
  public :: define_variables
  public :: set_global_attributes
  public :: data_shape_eq

contains

  function coordinate_variables(grid_type) result(coord_variables)
    character(len=*), intent(in) :: grid_type

    type(cable_output_variable_t), allocatable :: coord_variables(:)

    select case (grid_type)
    case ("restart")
      coord_variables = [ &
        cable_output_variable_t( &
          field_name="latitude", &
          data_shape=[CABLE_OUTPUT_DIM_LAND_GLOBAL], &
          distributed=.false., &
          aggregator=new_aggregator(latitude), &
          metadata=[attribute("units", "degrees_north")] &
        ), &
        cable_output_variable_t( &
          field_name="longitude", &
          data_shape=[CABLE_OUTPUT_DIM_LAND_GLOBAL], &
          distributed=.false., &
          aggregator=new_aggregator(longitude), &
          metadata=[attribute("units", "degrees_east")] &
        ) &
      ]
    case ("mask")
      coord_variables = [ &
        cable_output_variable_t( &
          field_name="lat_all", &
          netcdf_name="latitude", &
          data_shape=[CABLE_OUTPUT_DIM_X, CABLE_OUTPUT_DIM_Y], &
          var_type=CABLE_NETCDF_FLOAT, &
          parameter=.true., &
          distributed=.false., &
          aggregator=new_aggregator(lat_all), &
          metadata=[attribute("units", "degrees_north")] &
        ), &
        cable_output_variable_t( &
          field_name="lon_all", &
          netcdf_name="longitude", &
          data_shape=[CABLE_OUTPUT_DIM_X, CABLE_OUTPUT_DIM_Y], &
          parameter=.true., &
          distributed=.false., &
          aggregator=new_aggregator(lon_all), &
          metadata=[attribute("units", "degrees_east")] &
        ), &
        cable_output_variable_t( &
          field_name="x", &
          data_shape=[CABLE_OUTPUT_DIM_X], &
          parameter=.true., &
          distributed=.false., &
          aggregator=new_aggregator(lon_all(:, 1)), &
          metadata=[ &
            attribute("units", "degrees_east"), &
            attribute("comment", "x coordinate variable for GrADS compatibility") &
          ] &
        ), &
        cable_output_variable_t( &
          field_name="y", &
          data_shape=[CABLE_OUTPUT_DIM_Y], &
          parameter=.true., &
          distributed=.false., &
          aggregator=new_aggregator(lat_all(1, :)), &
          metadata=[ &
            attribute("units", "degrees_north"), &
            attribute("comment", "y coordinate variable for GrADS compatibility") &
          ] &
        ) &
      ]
    case ("land")
      coord_variables = [ &
        cable_output_variable_t( &
          field_name="local_lat", &
          data_shape=[CABLE_OUTPUT_DIM_LAND_GLOBAL], &
          parameter=.true., &
          distributed=.false., &
          aggregator=new_aggregator(latitude), &
          metadata=[attribute("units", "degrees_north")] &
        ), &
        cable_output_variable_t( &
          field_name="local_lon", &
          data_shape=[CABLE_OUTPUT_DIM_LAND_GLOBAL], &
          parameter=.true., &
          distributed=.false., &
          aggregator=new_aggregator(longitude), &
          metadata=[attribute("units", "degrees_east")] &
        ) &
      ]
    case default
      call cable_abort("Unexpected grid type '" // grid_type // "'", __FILE__, __LINE__)
    end select

  end function coordinate_variables

  logical function data_shape_eq(shape1, shape2)
    type(cable_output_dim_t), dimension(:), intent(in) :: shape1, shape2
    data_shape_eq = size(shape1) == size(shape2) .and. all(shape1 == shape2)
  end function

  function dim_size(dims)
    type(cable_output_dim_t), intent(in) :: dims(:)
    integer, allocatable :: dim_size(:)
    integer :: i

    allocate(dim_size(size(dims)))
    do i = 1, size(dims)
      select case (dims(i)%value)
      case (CABLE_OUTPUT_DIM_PATCH%value)
        dim_size(i) = mp
      case (CABLE_OUTPUT_DIM_SOIL%value)
        dim_size(i) = ms
      case (CABLE_OUTPUT_DIM_SNOW%value)
        dim_size(i) = msn
      case (CABLE_OUTPUT_DIM_RAD%value)
        dim_size(i) = nrb
      case (CABLE_OUTPUT_DIM_PLANTCARBON%value)
        dim_size(i) = ncp
      case (CABLE_OUTPUT_DIM_SOILCARBON%value)
        dim_size(i) = ncs
      case (CABLE_OUTPUT_DIM_LAND%value)
        dim_size(i) = mland
      case (CABLE_OUTPUT_DIM_LAND_GLOBAL%value)
        dim_size(i) = mland_global
      case (CABLE_OUTPUT_DIM_X%value)
        dim_size(i) = xdimsize
      case (CABLE_OUTPUT_DIM_Y%value)
        dim_size(i) = ydimsize
      case default
        call cable_abort("Unexpected dimension type", __FILE__, __LINE__)
      end select
    end do

  end function

  integer function netcdf_var_type(output_variable)
    type(cable_output_variable_t), intent(in) :: output_variable

    if (output_variable%var_type /= CABLE_OUTPUT_VAR_TYPE_UNDEFINED) then
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

  function infer_dim_names(output_profile, output_variable, time_axis) result(dim_names)
    type(cable_output_profile_t), intent(in) :: output_profile
    type(cable_output_variable_t), intent(in) :: output_variable
    logical, intent(in), optional :: time_axis

    character(MAX_LEN_DIM), allocatable :: dim_names(:)
    integer :: j

    allocate(dim_names(0))
    if (allocated(output_variable%data_shape)) then
      do j = 1, size(output_variable%data_shape)
        select case (output_variable%data_shape(j)%value)
        case (CABLE_OUTPUT_DIM_PATCH%value)
          select case (output_profile%grid_type)
          case ("restart")
            dim_names = [dim_names, "mp"]
          case ("land")
            if (output_variable%reduction_method == "none") then
              dim_names = [dim_names, "land", "patch"]
            else
              dim_names = [dim_names, "land"]
            end if
          case ("mask")
            if (output_variable%reduction_method == "none") then
              dim_names = [dim_names, "x", "y", "patch"]
            else
              dim_names = [dim_names, "x", "y"]
            end if
          case default
            call cable_abort("Unexpected grid type '" // output_profile%grid_type // &
              "' for variable '" // output_variable%field_name // "'", __FILE__, __LINE__)
          end select
        case (CABLE_OUTPUT_DIM_LAND%value)
          select case (output_profile%grid_type)
          case ("restart")
            dim_names = [dim_names, "mland"]
          case ("land")
            dim_names = [dim_names, "land"]
          case ("mask")
            dim_names = [dim_names, "x", "y"]
          case default
            call cable_abort("Unexpected grid type '" // output_profile%grid_type // &
              "' for variable '" // output_variable%field_name // "'", __FILE__, __LINE__)
          end select
        case (CABLE_OUTPUT_DIM_LAND_GLOBAL%value)
          if (output_profile%grid_type == "restart") then
            dim_names = [dim_names, "mland"]
          else
            dim_names = [dim_names, "land"]
          end if
        case (CABLE_OUTPUT_DIM_SOIL%value)
          dim_names = [dim_names, "soil"]
        case (CABLE_OUTPUT_DIM_SNOW%value)
          dim_names = [dim_names, "snow"]
        case (CABLE_OUTPUT_DIM_RAD%value)
          dim_names = [dim_names, "rad"]
        case (CABLE_OUTPUT_DIM_PLANTCARBON%value)
          if (output_profile%grid_type == "restart") then
            dim_names = [dim_names, "plant_carbon_pools"]
          else
            dim_names = [dim_names, "plantcarbon"]
          end if
        case (CABLE_OUTPUT_DIM_SOILCARBON%value)
          if (output_profile%grid_type == "restart") then
            dim_names = [dim_names, "soil_carbon_pools"]
          else
            dim_names = [dim_names, "soilcarbon"]
          end if
        case (CABLE_OUTPUT_DIM_X%value)
          dim_names = [dim_names, "x"]
        case (CABLE_OUTPUT_DIM_Y%value)
          dim_names = [dim_names, "y"]
        case default
          call cable_abort("Unexpected data shape for variable " // output_variable%field_name, __FILE__, __LINE__)
        end select
      end do
    end if

    if (present(time_axis)) then; if (time_axis) then
        dim_names = [dim_names, "time"]
    end if; end if

  end function

  function infer_cell_methods(output_variable) result(cell_methods)
    type(cable_output_variable_t), intent(in) :: output_variable
    character(len=256) :: cell_methods
    character(len=256) :: cell_methods_time, cell_methods_area

    if (.not. output_variable%parameter) then
      select case (output_variable%aggregation_method)
      case ("point")
        cell_methods_time = "time: point"
      case ("mean")
        cell_methods_time = "time: mean"
      case ("sum")
        cell_methods_time = "time: sum"
      case ("min")
        cell_methods_time = "time: minimum"
      case ("max")
        cell_methods_time = "time: maximum"
      case default
        call cable_abort("Unexpected aggregation method '" // output_variable%aggregation_method // &
          "' for variable '" // output_variable%field_name // "'", __FILE__, __LINE__)
      end select
    else
      cell_methods_time = ""
    end if

    select case (output_variable%reduction_method)
    case ("none")
      ! TODO(Sean): the cell_method for this case should be `area: point where
      ! pft` where `pft` is a string-valued auxiliary coordinate variable
      ! describing the labels for all patches:
      cell_methods_area = ""
    case ("first_patch_in_grid_cell")
      ! TODO(Sean): the cell_method for this case should be `area: point where
      ! <area_type>` where `<area_type>` is the area type of the first patch in
      ! the grid cell:
      cell_methods_area = ""
    case ("grid_cell_average")
      cell_methods_area = "area: mean where land"
    case default
      call cable_abort("Unexpected reduction method '" // output_variable%reduction_method // &
        "' for variable '" // output_variable%field_name // "'", __FILE__, __LINE__)
    end select

    cell_methods = adjustl(trim(cell_methods_time) // " " // trim(cell_methods_area))

  end function

  subroutine define_variables(output_profile, restart)
    type(cable_output_profile_t), intent(inout) :: output_profile
    logical, intent(in), optional :: restart

    logical :: restart_local
    integer :: i, j

    character(MAX_LEN_DIM), allocatable :: required_dimensions(:), dim_names(:)
    character(len=128) :: variable_name

    type(cable_output_variable_t), allocatable :: all_output_variables(:)

    restart_local = .false.
    if (present(restart)) restart_local = restart

    all_output_variables = [ &
      output_profile%coordinate_variables, &
      output_profile%output_variables &
    ]

    do i = 1, size(all_output_variables)
      associate(output_var => all_output_variables(i))
        if (.not. allocated(output_var%data_shape)) cycle
        dim_names = infer_dim_names( &
          output_profile, &
          output_var, &
          time_axis=(.not. (restart_local .or. output_var%parameter)) &
        )
        if (.not. allocated(required_dimensions)) then
          required_dimensions = dim_names
        else
          required_dimensions = [ &
            required_dimensions, &
            pack(dim_names, [( &
              all(dim_names(j) /= required_dimensions), &
              j = 1, &
              size(dim_names) &
            )]) &
          ]
        end if
      end associate
    end do

    do i = 1, size(required_dimensions)
      select case (required_dimensions(i))
      case ("mp")
        call output_profile%output_file%def_dims(["mp"], [mp_global])
      case ("mland")
        call output_profile%output_file%def_dims(["mland"], [mland_global])
      case ("land")
        call output_profile%output_file%def_dims(["land"], [mland_global])
      case ("x")
        call output_profile%output_file%def_dims(["x"], [xdimsize])
      case ("y")
        call output_profile%output_file%def_dims(["y"], [ydimsize])
      case ("patch")
        call output_profile%output_file%def_dims(["patch"], [max_vegpatches])
      case ("soil")
        call output_profile%output_file%def_dims(["soil"], [ms])
      case ("snow")
        call output_profile%output_file%def_dims(["snow"], [msn])
      case ("rad")
        call output_profile%output_file%def_dims(["rad"], [nrb])
      case ("soil_carbon_pools")
        call output_profile%output_file%def_dims(["soil_carbon_pools"], [ncs])
      case ("soilcarbon")
        call output_profile%output_file%def_dims(["soilcarbon"], [ncs])
      case ("plant_carbon_pools")
        call output_profile%output_file%def_dims(["plant_carbon_pools"], [ncp])
      case ("plantcarbon")
        call output_profile%output_file%def_dims(["plantcarbon"], [ncp])
      case ("time")
        ! time dimension defined separately below
      case default
        call cable_abort("Unexpected dimension name '" // required_dimensions(i) // "'", __FILE__, __LINE__)
      end select
    end do

    if (output_profile%grid_type == "restart") then
      call output_profile%output_file%def_dims(["time"], [1])
    else
      call output_profile%output_file%def_dims(["time"], [CABLE_NETCDF_UNLIMITED])
    end if

    call output_profile%output_file%def_var("time", ["time"], CABLE_NETCDF_DOUBLE)
    call output_profile%output_file%put_att("time", "units", timeunits)
    call output_profile%output_file%put_att("time", "coordinate", time_coord)
    call output_profile%output_file%put_att("time", "calendar", calendar)

    if (output_profile%grid_type /= "restart") then
      call output_profile%output_file%def_dims(["nv"], [2])
      call output_profile%output_file%def_var("time_bnds", ["nv  ", "time"], CABLE_NETCDF_DOUBLE)
      call output_profile%output_file%put_att("time", "bounds", "time_bnds")
    end if

    do i = 1, size(output_profile%coordinate_variables)
      associate(coord_var => output_profile%coordinate_variables(i))
        variable_name = var_name(coord_var)
        call output_profile%output_file%def_var( &
          var_name=variable_name, &
          dim_names=infer_dim_names(output_profile, coord_var), &
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
        dim_names = infer_dim_names( &
          output_profile, &
          output_var, &
          time_axis=(.not. (restart_local .or. output_var%parameter)) &
        )
        call output_profile%output_file%def_var( &
          var_name=variable_name, &
          dim_names=dim_names, &
          type=netcdf_var_type(output_var) &
        )
        if (allocated(output_var%metadata)) then
          do j = 1, size(output_var%metadata)
            call output_profile%output_file%put_att(variable_name, output_var%metadata(j)%name, output_var%metadata(j)%value)
          end do
        end if
        select case (netcdf_var_type(output_var))
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
        call output_profile%output_file%put_att(variable_name, "cell_methods", infer_cell_methods(output_var))
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
