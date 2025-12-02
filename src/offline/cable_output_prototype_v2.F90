module cable_output_prototype_v2_mod
  use iso_fortran_env, only: int32, real32, real64

  use cable_def_types_mod, only: mp, mp_global
  use cable_def_types_mod, only: mland, mland_global
  use cable_def_types_mod, only: ms
  use cable_def_types_mod, only: msn
  use cable_def_types_mod, only: nrb
  use cable_def_types_mod, only: ncs
  use cable_def_types_mod, only: ncp
  use cable_def_types_mod, only: met_type

  use cable_abort_module, only: cable_abort

  use cable_io_vars_module, only: metGrid, patch_type, land_type, xdimsize, ydimsize, max_vegpatches
  use cable_io_vars_module, only: timeunits, calendar, time_coord
  use cable_io_vars_module, only: check, ON_TIMESTEP, ON_WRITE

  use cable_checks_module, only: check_range

  use cable_io_vars_module, only: output, patchout

  use cable_timing_utils_mod, only: time_step_matches

  use aggregator_mod, only: aggregator_mod_init
  use aggregator_mod, only: aggregator_mod_end
  use aggregator_mod, only: aggregator_t
  use aggregator_mod, only: aggregator_handle_t
  use aggregator_mod, only: aggregator_int32_1d_t
  use aggregator_mod, only: aggregator_int32_2d_t
  use aggregator_mod, only: aggregator_int32_3d_t
  use aggregator_mod, only: aggregator_real32_1d_t
  use aggregator_mod, only: aggregator_real32_2d_t
  use aggregator_mod, only: aggregator_real32_3d_t
  use aggregator_mod, only: aggregator_real64_1d_t
  use aggregator_mod, only: aggregator_real64_2d_t
  use aggregator_mod, only: aggregator_real64_3d_t
  use aggregator_mod, only: store_aggregator

  use cable_netcdf_mod, only: cable_netcdf_file_t
  use cable_netcdf_mod, only: cable_netcdf_decomp_t
  use cable_netcdf_mod, only: cable_netcdf_create_file
  use cable_netcdf_mod, only: CABLE_NETCDF_INT
  use cable_netcdf_mod, only: CABLE_NETCDF_FLOAT
  use cable_netcdf_mod, only: CABLE_NETCDF_DOUBLE
  use cable_netcdf_mod, only: CABLE_NETCDF_UNLIMITED
  use cable_netcdf_mod, only: MAX_LEN_VAR => CABLE_NETCDF_MAX_STR_LEN_VAR
  use cable_netcdf_mod, only: MAX_LEN_DIM => CABLE_NETCDF_MAX_STR_LEN_DIM

  use cable_grid_reductions_mod, only: grid_cell_average

  implicit none
  private

  public :: cable_output_mod_init
  public :: cable_output_mod_end
  public :: cable_output_add_variable
  public :: cable_output_commit
  public :: cable_output_update
  public :: cable_output_write_parameters
  public :: output
  public :: patchout
  public :: requires_x_y_output_grid
  public :: requires_land_output_grid

  integer(kind=int32), parameter :: FILL_VALUE_INT32  = -9999_int32
  real(kind=real32),   parameter :: FILL_VALUE_REAL32 = -1.0e+33_real32
  real(kind=real64),   parameter :: FILL_VALUE_REAL64 = -1.0e+33_real64

  character(len=*), parameter :: DEFAULT_ACCUMULATION_FREQUENCY = "all"

  type cable_output_variable_t
    character(len=MAX_LEN_VAR) :: name
    character(len=MAX_LEN_DIM), allocatable :: dims(:)
    integer :: var_type
    character(len=50) :: units
    character(len=100) :: long_name
    character(len=100) :: cell_methods
    character(len=10) :: accumulation_frequency
    logical :: active
    character(len=50) :: reduction_method
    real, dimension(2) :: range
    type(aggregator_handle_t) :: aggregator_handle
    class(cable_netcdf_decomp_t), pointer :: decomp => null()
    real(kind=real32), pointer :: temp_buffer_real32_1d(:)       => null()
    real(kind=real32), pointer :: temp_buffer_real32_2d(:, :)    => null()
    real(kind=real32), pointer :: temp_buffer_real32_3d(:, :, :) => null()
    real(kind=real64), pointer :: temp_buffer_real64_1d(:)       => null()
    real(kind=real64), pointer :: temp_buffer_real64_2d(:, :)    => null()
    real(kind=real64), pointer :: temp_buffer_real64_3d(:, :, :) => null()
  end type

  type cable_output_profile_t
    real :: previous_write_time = 0.0
    integer :: frame = 0
    class(cable_netcdf_file_t), allocatable :: output_file
    type(cable_output_variable_t), allocatable :: output_variables(:), output_parameters(:)

    type(aggregator_handle_t), allocatable :: aggregators_accumulate_time_step(:)
    type(aggregator_handle_t), allocatable :: aggregators_accumulate_daily(:)
  end type

  ! Temporary buffers for computing grid-cell averages for each variable class
  real(kind=real32), allocatable, target :: temp_buffer_land_real32(:)
  real(kind=real64), allocatable, target :: temp_buffer_land_real64(:)
  real(kind=real32), allocatable, target :: temp_buffer_land_soil_real32(:, :)
  real(kind=real64), allocatable, target :: temp_buffer_land_soil_real64(:, :)
  real(kind=real32), allocatable, target :: temp_buffer_land_snow_real32(:, :)
  real(kind=real64), allocatable, target :: temp_buffer_land_snow_real64(:, :)
  real(kind=real32), allocatable, target :: temp_buffer_land_rad_real32(:, :)
  real(kind=real64), allocatable, target :: temp_buffer_land_rad_real64(:, :)
  real(kind=real32), allocatable, target :: temp_buffer_land_plantcarbon_real32(:, :)
  real(kind=real64), allocatable, target :: temp_buffer_land_plantcarbon_real64(:, :)
  real(kind=real32), allocatable, target :: temp_buffer_land_soilcarbon_real32(:, :)
  real(kind=real64), allocatable, target :: temp_buffer_land_soilcarbon_real64(:, :)

  ! TODO(Sean): once cable_write.F90 is removed, move the output_inclusion_type
  ! from cable_io_vars_module to here (as this would no longer introduce a cyclic
  ! module dependency). Then uncomment declarations below:
  ! type(output_inclusion_t) :: output
  ! type(output_inclusion_t) :: patchout ! do we want patch-specific info

  type(cable_output_profile_t), allocatable :: global_profile

contains

  logical function requires_x_y_output_grid(output_grid, met_grid)
    character(len=*), intent(in) :: output_grid
    character(len=*), intent(in) :: met_grid
    requires_x_y_output_grid = (( &
      output_grid == 'default' .AND. met_grid == 'mask' &
    ) .OR. ( &
      output_grid == 'mask' .OR. output_grid == 'ALMA' &
    ))
  end function

  logical function requires_land_output_grid(output_grid, met_grid)
    character(len=*), intent(in) :: output_grid
    character(len=*), intent(in) :: met_grid
    requires_land_output_grid = ( &
      output_grid == 'land' .OR. (output_grid == 'default' .AND. met_grid == 'land') &
    )
  end function

  logical function check_invalid_frequency(sampling_frequency, accumulation_frequency)
    character(len=*), intent(in) :: sampling_frequency
    character(len=*), intent(in) :: accumulation_frequency

    check_invalid_frequency = .false.

    ! TODO(Sean): return true if sampling frequency is greater than accumulation frequency

  end function

  subroutine cable_output_mod_init()
    class(cable_netcdf_file_t), allocatable :: output_file

    ! Initialize temporary buffers for grid-cell averaging
    allocate(temp_buffer_land_real32(mland))
    allocate(temp_buffer_land_real64(mland))
    allocate(temp_buffer_land_soil_real32(mland, ms))
    allocate(temp_buffer_land_soil_real64(mland, ms))
    allocate(temp_buffer_land_snow_real32(mland, msn))
    allocate(temp_buffer_land_snow_real64(mland, msn))
    allocate(temp_buffer_land_rad_real32(mland, nrb))
    allocate(temp_buffer_land_rad_real64(mland, nrb))
    allocate(temp_buffer_land_plantcarbon_real32(mland, ncp))
    allocate(temp_buffer_land_plantcarbon_real64(mland, ncp))
    allocate(temp_buffer_land_soilcarbon_real32(mland, ncs))
    allocate(temp_buffer_land_soilcarbon_real64(mland, ncs))

    call aggregator_mod_init()

    allocate(cable_output_profile_t::global_profile)

  end subroutine

  subroutine cable_output_mod_end()

    if (allocated(global_profile%output_file)) call global_profile%output_file%close()

    deallocate(global_profile)

    call aggregator_mod_end()

    deallocate(temp_buffer_land_real32)
    deallocate(temp_buffer_land_real64)
    deallocate(temp_buffer_land_soil_real32)
    deallocate(temp_buffer_land_soil_real64)
    deallocate(temp_buffer_land_snow_real32)
    deallocate(temp_buffer_land_snow_real64)
    deallocate(temp_buffer_land_rad_real32)
    deallocate(temp_buffer_land_rad_real64)
    deallocate(temp_buffer_land_plantcarbon_real32)
    deallocate(temp_buffer_land_plantcarbon_real64)
    deallocate(temp_buffer_land_soilcarbon_real32)
    deallocate(temp_buffer_land_soilcarbon_real64)

  end subroutine

  subroutine cable_output_add_variable( &
    name, dims, var_type, units, long_name, active, reduction_method, &
    decomp, range, accumulation_frequency, aggregator, parameter &
  )
    character(len=*), intent(in) :: name
    character(len=*), dimension(:), intent(in) :: dims
    integer, intent(in) :: var_type
    character(len=*), intent(in) :: units
    character(len=*), intent(in) :: long_name
    logical, intent(in) :: active
    character(len=*), intent(in), optional :: reduction_method
    class(cable_netcdf_decomp_t), intent(in), target :: decomp
    real, dimension(2), intent(in) :: range
    character(len=*), intent(in), optional :: accumulation_frequency
    class(aggregator_t), intent(in) :: aggregator
    logical, intent(in), optional :: parameter

    logical :: is_parameter
    type(cable_output_variable_t) :: output_var

    is_parameter = .false.
    if (present(parameter)) is_parameter = parameter

    if (present(reduction_method)) then
      select case (reduction_method)
      case ("none")
        ! No additional checks needed for 'none' reduction
      case ("grid_cell_average")
        select type (aggregator)
        type is (aggregator_real32_1d_t)
          if (size(aggregator%source_data, 1) /= mp) call cable_abort("Incompatible source data size for grid reduction", __FILE__, __LINE__)
        type is (aggregator_real32_2d_t)
          if (size(aggregator%source_data, 1) /= mp) call cable_abort("Incompatible source data size for grid reduction", __FILE__, __LINE__)
        type is (aggregator_real32_3d_t)
          if (size(aggregator%source_data, 1) /= mp) call cable_abort("Incompatible source data size for grid reduction", __FILE__, __LINE__)
        type is (aggregator_real64_1d_t)
          if (size(aggregator%source_data, 1) /= mp) call cable_abort("Incompatible source data size for grid reduction", __FILE__, __LINE__)
        type is (aggregator_real64_2d_t)
          if (size(aggregator%source_data, 1) /= mp) call cable_abort("Incompatible source data size for grid reduction", __FILE__, __LINE__)
        type is (aggregator_real64_3d_t)
          if (size(aggregator%source_data, 1) /= mp) call cable_abort("Incompatible source data size for grid reduction", __FILE__, __LINE__)
        class default
          call cable_abort("Unexpected aggregator type", __FILE__, __LINE__)
        end select
      case default
        call cable_abort("Invalid reduction method", __FILE__, __LINE__)
      end select
    end if

    ! TODO(Sean): determine cell_methods based on reduction and aggregation method

    output_var%name = name
    output_var%dims = dims
    output_var%units = units
    output_var%long_name = long_name
    output_var%active = active
    output_var%range = range
    output_var%decomp => decomp
    output_var%var_type = var_type

    if (present(reduction_method)) then
      output_var%reduction_method = reduction_method
    else
      output_var%reduction_method = "none"
    end if

    if (present(accumulation_frequency)) then
      output_var%accumulation_frequency = accumulation_frequency
    else
      output_var%accumulation_frequency = DEFAULT_ACCUMULATION_FREQUENCY
    end if

    if (active) then

      if (check_invalid_frequency( &
        sampling_frequency=output%averaging, &
        accumulation_frequency=output_var%accumulation_frequency &
      )) then
        call cable_abort("Sampling frequency and accumulation frequency are incompatible", __FILE__, __LINE__)
      end if

      if (present(reduction_method)) then
        select type(aggregator)
        type is (aggregator_real32_1d_t)
          if (all(shape(aggregator%source_data) == [mp])) then
            output_var%temp_buffer_real32_1d => temp_buffer_land_real32
          else
            call cable_abort("Unexpected source data shape for grid reduction", __FILE__, __LINE__)
          end if
        type is (aggregator_real64_1d_t)
          if (all(shape(aggregator%source_data) == [mp])) then
            output_var%temp_buffer_real64_1d => temp_buffer_land_real64
          else
            call cable_abort("Unexpected source data shape for grid reduction", __FILE__, __LINE__)
          end if
        type is (aggregator_real32_2d_t)
          if (all(shape(aggregator%source_data) == [mp, ms])) then
            output_var%temp_buffer_real32_2d => temp_buffer_land_soil_real32
          else if (all(shape(aggregator%source_data) == [mp, nrb])) then
            output_var%temp_buffer_real32_2d => temp_buffer_land_rad_real32
          else if (all(shape(aggregator%source_data) == [mp, msn])) then
            output_var%temp_buffer_real32_2d => temp_buffer_land_snow_real32
          else if (all(shape(aggregator%source_data) == [mp, nrb])) then
            output_var%temp_buffer_real32_2d => temp_buffer_land_rad_real32
          else if (all(shape(aggregator%source_data) == [mp, ncp])) then
            output_var%temp_buffer_real32_2d => temp_buffer_land_plantcarbon_real32
          else if (all(shape(aggregator%source_data) == [mp, ncs])) then
            output_var%temp_buffer_real32_2d => temp_buffer_land_soilcarbon_real32
          else
            call cable_abort("Unexpected source data shape for grid reduction", __FILE__, __LINE__)
          end if
        type is (aggregator_real64_2d_t)
          if (all(shape(aggregator%source_data) == [mp, ms])) then
            output_var%temp_buffer_real64_2d => temp_buffer_land_soil_real64
          else if (all(shape(aggregator%source_data) == [mp, nrb])) then
            output_var%temp_buffer_real64_2d => temp_buffer_land_rad_real64
          else if (all(shape(aggregator%source_data) == [mp, msn])) then
            output_var%temp_buffer_real64_2d => temp_buffer_land_snow_real64
          else if (all(shape(aggregator%source_data) == [mp, nrb])) then
            output_var%temp_buffer_real64_2d => temp_buffer_land_rad_real64
          else if (all(shape(aggregator%source_data) == [mp, ncp])) then
            output_var%temp_buffer_real64_2d => temp_buffer_land_plantcarbon_real64
          else if (all(shape(aggregator%source_data) == [mp, ncs])) then
            output_var%temp_buffer_real64_2d => temp_buffer_land_soilcarbon_real64
          else
            call cable_abort("Unexpected source data shape for grid reduction", __FILE__, __LINE__)
          end if
        class default
          call cable_abort("Unexpected aggregator type", __FILE__, __LINE__)
        end select
      end if

      output_var%aggregator_handle = store_aggregator(aggregator)

      if (is_parameter) then
        call output_var%aggregator_handle%init()
        if (.not. allocated(global_profile%output_parameters)) then
          global_profile%output_parameters = [output_var]
        else
          global_profile%output_parameters = [global_profile%output_parameters, output_var]
        end if
      else
        select case(output_var%accumulation_frequency)
        case("all")
          if (.not. allocated(global_profile%aggregators_accumulate_time_step)) then
            global_profile%aggregators_accumulate_time_step = [output_var%aggregator_handle]
          else
            global_profile%aggregators_accumulate_time_step = [global_profile%aggregators_accumulate_time_step, output_var%aggregator_handle]
          end if
        case("daily")
          if (.not. allocated(global_profile%aggregators_accumulate_daily)) then
            global_profile%aggregators_accumulate_daily = [output_var%aggregator_handle]
          else
            global_profile%aggregators_accumulate_daily = [global_profile%aggregators_accumulate_daily, output_var%aggregator_handle]
          end if
        case default
          call cable_abort("Invalid accumulation frequency", __FILE__, __LINE__)
        end select

        if (.not. allocated(global_profile%output_variables)) then
          global_profile%output_variables = [output_var]
        else
          global_profile%output_variables = [global_profile%output_variables, output_var]
        end if

      end if

    end if

  end subroutine cable_output_add_variable

  subroutine cable_output_commit()
    class(cable_netcdf_file_t), allocatable :: output_file
    integer :: i

    output_file = cable_netcdf_create_file("test_output.nc") ! TODO(Sean): use filename from namelist

    call output_file%def_dims(["x", "y"], [xdimsize, ydimsize])
    call output_file%def_dims(["patch"], [max_vegpatches])
    call output_file%def_dims(["soil"], [ms])
    call output_file%def_dims(["rad"], [nrb])
    call output_file%def_dims(["soil_carbon_pools"], [ncs])
    call output_file%def_dims(["plant_carbon_pools"], [ncp])
    call output_file%def_dims(["time"], [CABLE_NETCDF_UNLIMITED])
    call output_file%def_dims(["nv"], [2])

    if (requires_x_y_output_grid(output%grid, metgrid)) then
      call output_file%def_dims(["z"], [1]) ! Atmospheric 'z' dim of size 1 to comply with ALMA grid type
    else if (requires_land_output_grid(output%grid, metgrid)) then
      call output_file%def_dims(["land"], [mland_global])
      call output_file%def_var("local_lat", ["land"], CABLE_NETCDF_FLOAT)
      call output_file%put_att("local_lat", "units", "degrees_north")
      call output_file%def_var("local_lon", ["land"], CABLE_NETCDF_FLOAT)
      call output_file%put_att("local_lon", "units", "degrees_east")
    else
      call cable_abort("Error: Unable to determine output grid type", __FILE__, __LINE__)
    end if

    call output_file%def_var("time", ["time"], CABLE_NETCDF_DOUBLE)
    call output_file%put_att("time", "units", timeunits)
    call output_file%put_att("time", "coordinate", time_coord)
    call output_file%put_att("time", "calendar", calendar)
    call output_file%put_att("time", "bounds", "time_bnds")
    call output_file%def_var("time_bnds", ["nv", "time"], CABLE_NETCDF_DOUBLE)

    ! Define latitude and longitude variable (ALMA):
    call output_file%def_var("latitude", ["x", "y"], CABLE_NETCDF_FLOAT)
    call output_file%put_att("latitude", "units", "degrees_north")
    call output_file%def_var("longitude", ["x", "y"], CABLE_NETCDF_FLOAT)
    call output_file%put_att("longitude", "units", "degrees_east")

    ! Write "cordinate variables" to enable reading by GrADS:
    call output_file%def_var("x", ["x"], CABLE_NETCDF_FLOAT)
    call output_file%put_att("x", "units", "degrees_east")
    call output_file%put_att("x", "comment", "x coordinate variable for GrADS compatibility")
    call output_file%def_var("y", ["y"], CABLE_NETCDF_FLOAT)
    call output_file%put_att("y", "units", "degrees_north")
    call output_file%put_att("y", "comment", "y coordinate variable for GrADS compatibility")

    ! TODO(Sean): define remaining coordinate variables

    ! TODO(Sean): add global attributes

    ! TODO(Sean): should we just have a single list of output variables instead
    ! of parameters and variables?

    do i = 1, size(global_profile%output_parameters)
      associate(output_var => global_profile%output_parameters(i))
        call output_file%def_var( &
          var_name=output_var%name, &
          dim_names=output_var%dims, &
          type=output_var%var_type &
        )
        call output_file%put_att(output_var%name, 'units', output_var%units)
        call output_file%put_att(output_var%name, 'long_name', output_var%long_name)
        select case (output_var%var_type)
        case (CABLE_NETCDF_INT)
          call output_file%put_att(output_var%name, '_FillValue', FILL_VALUE_INT32)
          call output_file%put_att(output_var%name, 'missing_value', FILL_VALUE_INT32)
        case (CABLE_NETCDF_FLOAT)
          call output_file%put_att(output_var%name, '_FillValue', FILL_VALUE_REAL32)
          call output_file%put_att(output_var%name, 'missing_value', FILL_VALUE_REAL32)
        case (CABLE_NETCDF_DOUBLE)
          call output_file%put_att(output_var%name, '_FillValue', FILL_VALUE_REAL64)
          call output_file%put_att(output_var%name, 'missing_value', FILL_VALUE_REAL64)
        end select
        ! TODO(Sean): set cell_methods attribute
      end associate
    end do

    do i = 1, size(global_profile%output_variables)
      associate(output_var => global_profile%output_variables(i))
        call output_file%def_var( &
          var_name=output_var%name, &
          dim_names=output_var%dims, &
          type=output_var%var_type &
        )
        call output_file%put_att(output_var%name, 'units', output_var%units)
        call output_file%put_att(output_var%name, 'long_name', output_var%long_name)
        select case (output_var%var_type)
        case (CABLE_NETCDF_INT)
          call output_file%put_att(output_var%name, '_FillValue', FILL_VALUE_INT32)
          call output_file%put_att(output_var%name, 'missing_value', FILL_VALUE_INT32)
        case (CABLE_NETCDF_FLOAT)
          call output_file%put_att(output_var%name, '_FillValue', FILL_VALUE_REAL32)
          call output_file%put_att(output_var%name, 'missing_value', FILL_VALUE_REAL32)
        case (CABLE_NETCDF_DOUBLE)
          call output_file%put_att(output_var%name, '_FillValue', FILL_VALUE_REAL64)
          call output_file%put_att(output_var%name, 'missing_value', FILL_VALUE_REAL64)
        end select
        ! TODO(Sean): set cell_methods attribute
      end associate
    end do

    global_profile%output_file = output_file

    ! Initialise all aggregators

    do i = 1, size(global_profile%aggregators_accumulate_time_step)
      associate(aggregator_handle => global_profile%aggregators_accumulate_time_step(i))
        call aggregator_handle%init()
      end associate
    end do

    do i = 1, size(global_profile%aggregators_accumulate_daily)
      associate(aggregator_handle => global_profile%aggregators_accumulate_daily(i))
        call aggregator_handle%init()
      end associate
    end do

  end subroutine

  subroutine cable_output_write_parameters(time_index, patch, landpt, met)
    integer, intent(in) :: time_index
    type(patch_type), intent(in) :: patch(:)
    type(land_type), intent(in) :: landpt(:)
    type(met_type), intent(in) :: met

    integer :: i

    do i = 1, size(global_profile%output_parameters)
      associate(output_variable => global_profile%output_parameters(i))
        call check_variable_range(output_variable, time_index, met)
        call output_variable%aggregator_handle%accumulate()
        select case (output_variable%reduction_method)
        case ("grid_cell_average")
          call write_variable_grid_cell_average(output_variable, global_profile%output_file, patch, landpt)
        case ("none")
          call write_variable(output_variable, global_profile%output_file)
        case default
          call cable_abort("Invalid reduction method", __FILE__, __LINE__)
        end select
        call output_variable%aggregator_handle%reset()
      end associate
    end do

  end subroutine cable_output_write_parameters

  subroutine cable_output_update(time_index, dels, leaps, start_year, patch, landpt, met)
    integer, intent(in) :: time_index
    real, intent(in) :: dels
    logical, intent(in) :: leaps
    integer, intent(in) :: start_year
    type(patch_type), intent(in) :: patch(:)
    type(land_type), intent(in) :: landpt(:)
    type(met_type), intent(in) :: met

    real :: current_time
    integer :: i

    if (check%ranges == ON_TIMESTEP) then
      do i = 1, size(global_profile%output_variables)
        call check_variable_range(global_profile%output_variables(i), time_index, met)
      end do
    end if

    do i = 1, size(global_profile%aggregators_accumulate_time_step)
      associate(aggregator_handle => global_profile%aggregators_accumulate_time_step(i))
        call aggregator_handle%accumulate()
      end associate
    end do

    if (time_step_matches(dels, time_index, "daily", leaps, start_year)) then
      do i = 1, size(global_profile%aggregators_accumulate_daily)
        associate(aggregator_handle => global_profile%aggregators_accumulate_daily(i))
          call aggregator_handle%accumulate()
        end associate
      end do
    end if

    if (time_step_matches(dels, time_index, output%averaging, leaps, start_year)) then

      do i = 1, size(global_profile%output_variables)
        associate(output_variable => global_profile%output_variables(i))
          if (check%ranges == ON_WRITE) call check_variable_range(output_variable, time_index, met)
          call output_variable%aggregator_handle%normalise()
          select case (output_variable%reduction_method)
          case ("grid_cell_average")
            call write_variable_grid_cell_average(output_variable, global_profile%output_file, patch, landpt, global_profile%frame + 1)
          case ("none")
            call write_variable(output_variable, global_profile%output_file, global_profile%frame + 1)
          case default
            call cable_abort("Invalid reduction method", __FILE__, __LINE__)
          end select
          call output_variable%aggregator_handle%reset()
        end associate
      end do

      current_time = time_index * dels
      call global_profile%output_file%put_var("time", (current_time + global_profile%previous_write_time) / 2.0, start=[global_profile%frame + 1])
      call global_profile%output_file%put_var("time_bnds", [global_profile%previous_write_time, current_time], start=[1, global_profile%frame + 1])
      global_profile%previous_write_time = current_time
      global_profile%frame = global_profile%frame + 1

    end if

  end subroutine cable_output_update

  subroutine check_variable_range(output_variable, time_index, met)
    type(cable_output_variable_t), intent(in) :: output_variable
    integer, intent(in) :: time_index
    type(met_type), intent(in) :: met

    select type (aggregator => output_variable%aggregator_handle%aggregator)
    type is (aggregator_int32_1d_t)
      ! TODO(Sean): implement range checking for integer types
    type is (aggregator_int32_2d_t)
      ! TODO(Sean): implement range checking for integer types
    type is (aggregator_int32_3d_t)
      ! TODO(Sean): implement range checking for integer types
    type is (aggregator_real32_1d_t)
      call check_range(output_variable%name, aggregator%source_data, output_variable%range, time_index, met)
    type is (aggregator_real32_2d_t)
      call check_range(output_variable%name, aggregator%source_data, output_variable%range, time_index, met)
    type is (aggregator_real32_3d_t)
      call check_range(output_variable%name, aggregator%source_data, output_variable%range, time_index, met)
    type is (aggregator_real64_1d_t)
      ! TODO(Sean): implement range checking for double precision types
    type is (aggregator_real64_2d_t)
      ! TODO(Sean): implement range checking for double precision types
    type is (aggregator_real64_3d_t)
      ! TODO(Sean): implement range checking for double precision types
    class default
      call cable_abort("Unexpected aggregator type", __FILE__, __LINE__)
    end select

  end subroutine check_variable_range

  subroutine write_variable(output_variable, output_file, time_index)
    type(cable_output_variable_t), intent(inout) :: output_variable
    class(cable_netcdf_file_t), intent(inout) :: output_file
    integer, intent(in), optional :: time_index

    select type (aggregator => output_variable%aggregator_handle%aggregator)
    type is (aggregator_int32_1d_t)
      call output_file%write_darray( &
            var_name=output_variable%name, &
            values=aggregator%aggregated_data, &
            decomp=output_variable%decomp, &
            frame=time_index)
    type is (aggregator_int32_2d_t)
      call output_file%write_darray( &
            var_name=output_variable%name, &
            values=aggregator%aggregated_data, &
            decomp=output_variable%decomp, &
            frame=time_index)
    type is (aggregator_int32_3d_t)
      call output_file%write_darray( &
            var_name=output_variable%name, &
            values=aggregator%aggregated_data, &
            decomp=output_variable%decomp, &
            frame=time_index)
    type is (aggregator_real32_1d_t)
      call output_file%write_darray( &
            var_name=output_variable%name, &
            values=aggregator%aggregated_data, &
            decomp=output_variable%decomp, &
            fill_value=FILL_VALUE_REAL32, &
            frame=time_index)
    type is (aggregator_real32_2d_t)
      call output_file%write_darray( &
            var_name=output_variable%name, &
            values=aggregator%aggregated_data, &
            decomp=output_variable%decomp, &
            fill_value=FILL_VALUE_REAL32, &
            frame=time_index)
    type is (aggregator_real32_3d_t)
      call output_file%write_darray( &
            var_name=output_variable%name, &
            values=aggregator%aggregated_data, &
            decomp=output_variable%decomp, &
            fill_value=FILL_VALUE_REAL32, &
            frame=time_index)
    type is (aggregator_real64_1d_t)
      call output_file%write_darray( &
            var_name=output_variable%name, &
            values=aggregator%aggregated_data, &
            decomp=output_variable%decomp, &
            fill_value=FILL_VALUE_REAL64, &
            frame=time_index)
    type is (aggregator_real64_2d_t)
      call output_file%write_darray( &
            var_name=output_variable%name, &
            values=aggregator%aggregated_data, &
            decomp=output_variable%decomp, &
            fill_value=FILL_VALUE_REAL64, &
            frame=time_index)
    type is (aggregator_real64_3d_t)
      call output_file%write_darray( &
            var_name=output_variable%name, &
            values=aggregator%aggregated_data, &
            decomp=output_variable%decomp, &
            fill_value=FILL_VALUE_REAL64, &
            frame=time_index)
    class default
      call cable_abort("Unexpected aggregator type", __FILE__, __LINE__)
    end select

  end subroutine write_variable

  subroutine write_variable_grid_cell_average(output_variable, output_file, patch, landpt, time_index)
    type(cable_output_variable_t), intent(inout) :: output_variable
    class(cable_netcdf_file_t), intent(inout) :: output_file
    type(patch_type), intent(in) :: patch(:)
    type(land_type), intent(in) :: landpt(:)
    integer, intent(in), optional :: time_index

    select type (aggregator => output_variable%aggregator_handle%aggregator)
    type is (aggregator_real32_1d_t)
      call grid_cell_average( &
            input_array=aggregator%aggregated_data, &
            output_array=output_variable%temp_buffer_real32_1d, &
            landpt=landpt, &
            patch=patch)
      call output_file%write_darray( &
            var_name=output_variable%name, &
            values=output_variable%temp_buffer_real32_1d, &
            decomp=output_variable%decomp, &
            fill_value=FILL_VALUE_REAL32, &
            frame=time_index)
    type is (aggregator_real32_2d_t)
      call grid_cell_average( &
            input_array=aggregator%aggregated_data, &
            output_array=output_variable%temp_buffer_real32_2d, &
            landpt=landpt, &
            patch=patch)
      call output_file%write_darray( &
            var_name=output_variable%name, &
            values=output_variable%temp_buffer_real32_2d, &
            decomp=output_variable%decomp, &
            fill_value=FILL_VALUE_REAL32, &
            frame=time_index)
    type is (aggregator_real32_3d_t)
      call grid_cell_average( &
            input_array=aggregator%aggregated_data, &
            output_array=output_variable%temp_buffer_real32_3d, &
            landpt=landpt, &
            patch=patch)
      call output_file%write_darray( &
            var_name=output_variable%name, &
            values=output_variable%temp_buffer_real32_3d, &
            decomp=output_variable%decomp, &
            fill_value=FILL_VALUE_REAL32, &
            frame=time_index)
    type is (aggregator_real64_1d_t)
      call grid_cell_average( &
            input_array=aggregator%aggregated_data, &
            output_array=output_variable%temp_buffer_real64_1d, &
            landpt=landpt, &
            patch=patch)
      call output_file%write_darray( &
            var_name=output_variable%name, &
            values=output_variable%temp_buffer_real64_1d, &
            decomp=output_variable%decomp, &
            fill_value=FILL_VALUE_REAL64, &
            frame=time_index)
    type is (aggregator_real64_2d_t)
      call grid_cell_average( &
            input_array=aggregator%aggregated_data, &
            output_array=output_variable%temp_buffer_real64_2d, &
            landpt=landpt, &
            patch=patch)
      call output_file%write_darray( &
            var_name=output_variable%name, &
            values=output_variable%temp_buffer_real64_2d, &
            decomp=output_variable%decomp, &
            fill_value=FILL_VALUE_REAL64, &
            frame=time_index)
    type is (aggregator_real64_3d_t)
      call grid_cell_average( &
            input_array=aggregator%aggregated_data, &
            output_array=output_variable%temp_buffer_real64_3d, &
            landpt=landpt, &
            patch=patch)
      call output_file%write_darray( &
            var_name=output_variable%name, &
            values=output_variable%temp_buffer_real64_3d, &
            decomp=output_variable%decomp, &
            fill_value=FILL_VALUE_REAL64, &
            frame=time_index)
    class default
      call cable_abort("Unexpected aggregator type", __FILE__, __LINE__)
    end select

  end subroutine write_variable_grid_cell_average

end module
