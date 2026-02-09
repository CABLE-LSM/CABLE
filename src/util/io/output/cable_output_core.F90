module cable_output_core_mod

  use iso_fortran_env, only: int32, real32, real64

  use cable_def_types_mod, only: met_type

  use cable_io_vars_module, only: patch_type
  use cable_io_vars_module, only: land_type
  use cable_io_vars_module, only: metgrid
  use cable_io_vars_module, only: output
  use cable_io_vars_module, only: check
  use cable_io_vars_module, only: ON_TIMESTEP
  use cable_io_vars_module, only: ON_WRITE

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

  use cable_netcdf_mod, only: cable_netcdf_decomp_t
  use cable_netcdf_mod, only: cable_netcdf_file_t
  use cable_netcdf_mod, only: cable_netcdf_create_file
  use cable_netcdf_mod, only: CABLE_NETCDF_INT
  use cable_netcdf_mod, only: CABLE_NETCDF_FLOAT
  use cable_netcdf_mod, only: CABLE_NETCDF_DOUBLE
  use cable_netcdf_mod, only: CABLE_NETCDF_IOTYPE_CLASSIC

  use cable_abort_module, only: cable_abort

  use cable_checks_module, only: check_range

  use cable_timing_utils_mod, only: time_step_matches

  use cable_grid_reductions_mod, only: grid_cell_average
  use cable_grid_reductions_mod, only: first_patch_in_grid_cell

  use cable_output_types_mod, only: cable_output_variable_t
  use cable_output_types_mod, only: cable_output_profile_t
  use cable_output_types_mod, only: FILL_VALUE_INT32
  use cable_output_types_mod, only: FILL_VALUE_REAL32
  use cable_output_types_mod, only: FILL_VALUE_REAL64

  use cable_output_reduction_buffers_mod, only: allocate_grid_reduction_buffers
  use cable_output_reduction_buffers_mod, only: deallocate_grid_reduction_buffers
  use cable_output_reduction_buffers_mod, only: associate_temp_buffer_int32
  use cable_output_reduction_buffers_mod, only: associate_temp_buffer_real32
  use cable_output_reduction_buffers_mod, only: associate_temp_buffer_real64

  use cable_output_decomp_mod, only: allocate_decompositions
  use cable_output_decomp_mod, only: deallocate_decompositions
  use cable_output_decomp_mod, only: associate_decomp_int32
  use cable_output_decomp_mod, only: associate_decomp_real32
  use cable_output_decomp_mod, only: associate_decomp_real64

  use cable_output_utils_mod, only: check_invalid_frequency
  use cable_output_utils_mod, only: dim_size
  use cable_output_utils_mod, only: define_variables
  use cable_output_utils_mod, only: set_global_attributes

  use cable_output_definitions_mod, only: coordinate_variables

  implicit none
  private

  public :: cable_output_mod_init
  public :: cable_output_mod_end
  public :: cable_output_register_output_variables
  public :: cable_output_profiles_init
  public :: cable_output_update
  public :: cable_output_write
  public :: cable_output_write_parameters
  public :: cable_output_write_restart

  type(cable_output_profile_t), allocatable :: global_profile

  type(cable_output_variable_t), allocatable :: registered_output_variables(:)

contains

  subroutine cable_output_mod_init()
    class(cable_netcdf_file_t), allocatable :: output_file

    call allocate_decompositions()
    call allocate_grid_reduction_buffers()

  end subroutine

  subroutine cable_output_mod_end()

    if (allocated(global_profile%output_file)) call global_profile%output_file%close()

    deallocate(global_profile)

    call deallocate_grid_reduction_buffers()
    call deallocate_decompositions()

  end subroutine

  subroutine cable_output_register_output_variables(output_variables)
    type(cable_output_variable_t), dimension(:), intent(in) :: output_variables
    integer :: i

    do i = 1, size(output_variables)
      associate(output_var => output_variables(i))
        if (all(output_var%reduction_method /= [character(32) :: "none", "grid_cell_average", "first_patch_in_grid_cell"])) then
          call cable_abort("Invalid reduction method for variable " // trim(output_var%name), __FILE__, __LINE__)
        end if
        if (all(output_var%aggregation_method /= [character(32) :: "point", "mean", "max", "min", "sum"])) then
          call cable_abort("Invalid aggregation method for variable " // trim(output_var%name), __FILE__, __LINE__)
        end if
        if (all(output_var%var_type /= [CABLE_NETCDF_INT, CABLE_NETCDF_FLOAT, CABLE_NETCDF_DOUBLE])) then
          call cable_abort("Invalid variable type for variable " // trim(output_var%name), __FILE__, __LINE__)
        end if
        if (count(output_var%name == output_variables(:)%name) > 1) then
          call cable_abort("Duplicate variable name found: " // trim(output_var%name), __FILE__, __LINE__)
        end if
        if (( &
          .not. allocated(output_var%data_shape) .and. output_var%aggregator%rank() /= 0 &
        ) .or. ( &
          allocated(output_var%data_shape) .and. any(dim_size(output_var%data_shape) /= output_var%aggregator%shape()) &
        )) then
          call cable_abort("Data shape does not match aggregator shape for variable " // trim(output_var%name), __FILE__, __LINE__)
        end if
        if (output_var%reduction_method /= "none" .and. .not. output_var%distributed) then
          call cable_abort("Grid cell reductions require distributed output for variable " // trim(output_var%name), __FILE__, __LINE__)
        end if
      end associate
    end do

    registered_output_variables = output_variables

  end subroutine cable_output_register_output_variables

  subroutine cable_output_profiles_init()
    class(cable_netcdf_file_t), allocatable :: output_file
    integer :: i

    character(32) :: grid_type

    if (output%grid == "land" .OR. (output%grid == "default" .AND. metgrid == "land")) then
      grid_type = "land"
    else if (( &
      output%grid == "default" .AND. metgrid == "mask" &
    ) .OR. ( &
      output%grid == "mask" .OR. output%grid == "ALMA" &
    )) then
      grid_type = "mask"
    else
      call cable_abort("Unable to determine output grid type.", __FILE__, __LINE__)
    end if

    global_profile = cable_output_profile_t( &
      sampling_frequency=output%averaging, &
      grid_type=grid_type, &
      file_name="test_output.nc", & ! TODO(Sean): use filename from namelist
      output_file=cable_netcdf_create_file("test_output.nc", iotype=CABLE_NETCDF_IOTYPE_CLASSIC), & ! TODO(Sean): use filename from namelist
      output_variables=[ &
        coordinate_variables(grid_type), &
        pack(registered_output_variables, registered_output_variables(:)%active) &
      ] &
    )

    do i = 1, size(global_profile%output_variables)
      associate(output_var => global_profile%output_variables(i))
        call check_invalid_frequency( &
          sampling_frequency=global_profile%sampling_frequency, &
          accumulation_frequency=output_var%accumulation_frequency, &
          var_name=output_var%name, &
          file_name=global_profile%file_name &
        )
      end associate
    end do

    call define_variables(global_profile)

    call set_global_attributes(global_profile)

    call global_profile%output_file%end_def()

    do i = 1, size(global_profile%output_variables)
      associate(output_variable => global_profile%output_variables(i))
        call output_variable%aggregator%init(method=output_variable%aggregation_method)
      end associate
    end do

  end subroutine

  subroutine cable_output_write_parameters(time_index, patch, landpt, met)
    integer, intent(in) :: time_index
    type(patch_type), intent(in) :: patch(:)
    type(land_type), intent(in) :: landpt(:)
    type(met_type), intent(in) :: met

    integer :: i

    do i = 1, size(global_profile%output_variables)
      associate(output_variable => global_profile%output_variables(i))
        if (.not. output_variable%parameter) cycle
        call check_variable_range(output_variable, time_index, met)
        call output_variable%aggregator%accumulate()
        call write_variable(global_profile, output_variable, patch, landpt)
        call output_variable%aggregator%reset()
      end associate
    end do

  end subroutine cable_output_write_parameters

  subroutine cable_output_update(time_index, dels, leaps, start_year, met)
    integer, intent(in) :: time_index
    real, intent(in) :: dels
    logical, intent(in) :: leaps
    integer, intent(in) :: start_year
    type(met_type), intent(in) :: met

    real :: current_time
    integer :: i

    if (check%ranges == ON_TIMESTEP) then
      do i = 1, size(global_profile%output_variables)
        call check_variable_range(global_profile%output_variables(i), time_index, met)
      end do
    end if

    do i = 1, size(global_profile%output_variables)
      associate(output_variable => global_profile%output_variables(i))
        if (time_step_matches(dels, time_index, output_variable%accumulation_frequency, leaps, start_year)) then
          call output_variable%aggregator%accumulate()
        end if
      end associate
    end do

  end subroutine cable_output_update

  subroutine cable_output_write(time_index, dels, leaps, start_year, met, patch, landpt)
    integer, intent(in) :: time_index
    real, intent(in) :: dels
    logical, intent(in) :: leaps
    integer, intent(in) :: start_year
    type(met_type), intent(in) :: met
    type(patch_type), intent(in) :: patch(:)
    type(land_type), intent(in) :: landpt(:)

    real :: current_time
    integer :: i

    if (time_step_matches(dels, time_index, global_profile%sampling_frequency, leaps, start_year)) then

      do i = 1, size(global_profile%output_variables)
        associate(output_variable => global_profile%output_variables(i))
          if (output_variable%parameter) cycle
          if (check%ranges == ON_WRITE) call check_variable_range(output_variable, time_index, met)
          call write_variable(global_profile, output_variable, patch, landpt, frame=global_profile%frame + 1)
          call output_variable%aggregator%reset()
        end associate
      end do

      current_time = time_index * dels

      if (global_profile%sampling_frequency == "all") then
        call global_profile%output_file%put_var("time", current_time, start=[global_profile%frame + 1])
      else
        call global_profile%output_file%put_var("time", (current_time + global_profile%previous_write_time) / 2.0, start=[global_profile%frame + 1])
      end if

      call global_profile%output_file%put_var("time_bnds", [global_profile%previous_write_time, current_time], start=[1, global_profile%frame + 1])

      global_profile%previous_write_time = current_time
      global_profile%frame = global_profile%frame + 1

    end if

  end subroutine cable_output_write

  subroutine cable_output_write_restart(current_time)
    real, intent(in) :: current_time !! Current simulation time

    type(cable_output_profile_t), allocatable :: restart_output_profile
    integer :: i

    restart_output_profile = cable_output_profile_t( &
      sampling_frequency="none", &
      grid_type="restart", &
      file_name="test_restart.nc", & ! TODO(Sean): use filename from namelist
      output_file=cable_netcdf_create_file("test_restart.nc", iotype=CABLE_NETCDF_IOTYPE_CLASSIC), & ! TODO(Sean): use filename from namelist
      output_variables=[ &
        coordinate_variables(grid_type="restart"), &
        pack(registered_output_variables, registered_output_variables(:)%restart) &
      ] &
    )

    call define_variables(restart_output_profile)

    call restart_output_profile%output_file%end_def()

    call restart_output_profile%output_file%put_var("time", [current_time])

    do i = 1, size(restart_output_profile%output_variables)
      call write_variable(restart_output_profile, restart_output_profile%output_variables(i), restart=.true.)
    end do

    call restart_output_profile%output_file%close()

  end subroutine cable_output_write_restart

  subroutine check_variable_range(output_variable, time_index, met)
    type(cable_output_variable_t), intent(in) :: output_variable
    integer, intent(in) :: time_index
    type(met_type), intent(in) :: met

    select type (aggregator => output_variable%aggregator)
    type is (aggregator_int32_0d_t)
      ! TODO(Sean): implement range checking for integer types
    type is (aggregator_int32_1d_t)
      ! TODO(Sean): implement range checking for integer types
    type is (aggregator_int32_2d_t)
      ! TODO(Sean): implement range checking for integer types
    type is (aggregator_int32_3d_t)
      ! TODO(Sean): implement range checking for integer types
    type is (aggregator_real32_0d_t)
      ! TODO(Sean): implement range checking for scalars
    type is (aggregator_real32_1d_t)
      call check_range(output_variable%name, aggregator%source_data, output_variable%range, time_index, met)
    type is (aggregator_real32_2d_t)
      call check_range(output_variable%name, aggregator%source_data, output_variable%range, time_index, met)
    type is (aggregator_real32_3d_t)
      call check_range(output_variable%name, aggregator%source_data, output_variable%range, time_index, met)
    type is (aggregator_real64_0d_t)
      ! TODO(Sean): implement range checking for double precision types
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

  subroutine write_variable(output_profile, output_variable, patch, landpt, frame, restart)
    type(cable_output_profile_t), intent(inout) :: output_profile
    type(cable_output_variable_t), intent(inout), target :: output_variable
    type(patch_type), intent(in), optional :: patch(:)
    type(land_type), intent(in), optional :: landpt(:)
    integer, intent(in), optional :: frame
    logical, intent(in), optional :: restart

    class(cable_netcdf_decomp_t), pointer :: decomp
    integer :: i, ndims
    logical :: restart_local

    integer(kind=int32), pointer :: write_buffer_int32_0d
    integer(kind=int32), pointer :: write_buffer_int32_1d(:)
    integer(kind=int32), pointer :: write_buffer_int32_2d(:, :)
    integer(kind=int32), pointer :: write_buffer_int32_3d(:, :, :)
    real(kind=real32),   pointer :: write_buffer_real32_0d
    real(kind=real32),   pointer :: write_buffer_real32_1d(:)
    real(kind=real32),   pointer :: write_buffer_real32_2d(:, :)
    real(kind=real32),   pointer :: write_buffer_real32_3d(:, :, :)
    real(kind=real64),   pointer :: write_buffer_real64_0d
    real(kind=real64),   pointer :: write_buffer_real64_1d(:)
    real(kind=real64),   pointer :: write_buffer_real64_2d(:, :)
    real(kind=real64),   pointer :: write_buffer_real64_3d(:, :, :)

    decomp => null()

    write_buffer_int32_0d  => null()
    write_buffer_int32_1d  => null()
    write_buffer_int32_2d  => null()
    write_buffer_int32_3d  => null()
    write_buffer_real32_0d => null()
    write_buffer_real32_1d => null()
    write_buffer_real32_2d => null()
    write_buffer_real32_3d => null()
    write_buffer_real64_0d => null()
    write_buffer_real64_1d => null()
    write_buffer_real64_2d => null()
    write_buffer_real64_3d => null()

    restart_local = .false.
    if (present(restart)) restart_local = restart

    if (.not. restart_local .and. output_variable%reduction_method /= "none") then
      if (.not. present(patch) .or. .not. present(landpt)) then
        call cable_abort("Optional arguments patch and landpt must be present for grid reductions", __FILE__, __LINE__)
      end if
    end if

    select type (aggregator => output_variable%aggregator)
    type is (aggregator_int32_0d_t)
      if (output_variable%reduction_method /= "none") then
        call cable_abort("Grid cell reductions are not supported for scalar variables", __FILE__, __LINE__)
      end if
      if (output_variable%distributed) then
        call cable_abort("Distributed writes are not supported for scalar variables", __FILE__, __LINE__)
      end if
      write_buffer_int32_0d => aggregator%aggregated_data
      if (restart_local) write_buffer_int32_0d => aggregator%source_data
      if (present(frame)) then
        call output_profile%output_file%inq_var_ndims(output_variable%name, ndims)
        call output_profile%output_file%put_var( &
              var_name=output_variable%name, &
              values=write_buffer_int32_0d, &
              start=[(1, i = 1, ndims - 1), frame])
      else
        call output_profile%output_file%put_var( &
              var_name=output_variable%name, &
              values=write_buffer_int32_0d)
      end if
    type is (aggregator_int32_1d_t)
      if (restart_local) then
        write_buffer_int32_1d => aggregator%source_data
      else if (output_variable%reduction_method == "none") then
        write_buffer_int32_1d => aggregator%aggregated_data
      else if (output_variable%reduction_method == "grid_cell_average") then
        call cable_abort("Reduction method grid_cell_average is not supported for integer variables", __FILE__, __LINE__)
      else if (output_variable%reduction_method == "first_patch_in_grid_cell") then
        call associate_temp_buffer_int32(output_variable, temp_buffer_int32_1d=write_buffer_int32_1d)
        call first_patch_in_grid_cell( &
              input_array=aggregator%aggregated_data, &
              output_array=write_buffer_int32_1d, &
              landpt=landpt)
      else
        call cable_abort("Invalid reduction method", __FILE__, __LINE__)
      end if
      if (output_variable%distributed) then
        call associate_decomp_int32(output_profile, output_variable, decomp)
        call output_profile%output_file%write_darray( &
              var_name=output_variable%name, &
              values=write_buffer_int32_1d, &
              decomp=decomp, &
              fill_value=FILL_VALUE_INT32, &
              frame=frame)
      else if (present(frame)) then
        call output_profile%output_file%inq_var_ndims(output_variable%name, ndims)
        call output_profile%output_file%put_var( &
              var_name=output_variable%name, &
              values=write_buffer_int32_1d, &
              start=[(1, i = 1, ndims - 1), frame])
      else
        call output_profile%output_file%put_var( &
              var_name=output_variable%name, &
              values=write_buffer_int32_1d)
      end if
    type is (aggregator_int32_2d_t)
      if (restart_local) then
        write_buffer_int32_2d => aggregator%source_data
      else if (output_variable%reduction_method == "none") then
        write_buffer_int32_2d => aggregator%aggregated_data
      else if (output_variable%reduction_method == "grid_cell_average") then
        call cable_abort("Reduction method grid_cell_average is not supported for integer variables", __FILE__, __LINE__)
      else if (output_variable%reduction_method == "first_patch_in_grid_cell") then
        call associate_temp_buffer_int32(output_variable, temp_buffer_int32_2d=write_buffer_int32_2d)
        call first_patch_in_grid_cell( &
              input_array=aggregator%aggregated_data, &
              output_array=write_buffer_int32_2d, &
              landpt=landpt)
      else
        call cable_abort("Invalid reduction method", __FILE__, __LINE__)
      end if
      if (output_variable%distributed) then
        call associate_decomp_int32(output_profile, output_variable, decomp)
        call output_profile%output_file%write_darray( &
              var_name=output_variable%name, &
              values=write_buffer_int32_2d, &
              decomp=decomp, &
              fill_value=FILL_VALUE_INT32, &
              frame=frame)
      else if (present(frame)) then
        call output_profile%output_file%inq_var_ndims(output_variable%name, ndims)
        call output_profile%output_file%put_var( &
              var_name=output_variable%name, &
              values=write_buffer_int32_2d, &
              start=[(1, i = 1, ndims - 1), frame])
      else
        call output_profile%output_file%put_var( &
              var_name=output_variable%name, &
              values=write_buffer_int32_2d)
      end if
    type is (aggregator_int32_3d_t)
      if (restart_local) then
        write_buffer_int32_3d => aggregator%source_data
      else if (output_variable%reduction_method == "none") then
        write_buffer_int32_3d => aggregator%aggregated_data
      else if (output_variable%reduction_method == "grid_cell_average") then
        call cable_abort("Reduction method grid_cell_average is not supported for integer variables", __FILE__, __LINE__)
      else if (output_variable%reduction_method == "first_patch_in_grid_cell") then
        call associate_temp_buffer_int32(output_variable, temp_buffer_int32_3d=write_buffer_int32_3d)
        call first_patch_in_grid_cell( &
              input_array=aggregator%aggregated_data, &
              output_array=write_buffer_int32_3d, &
              landpt=landpt)
      else
        call cable_abort("Invalid reduction method", __FILE__, __LINE__)
      end if
      if (output_variable%distributed) then
        call associate_decomp_int32(output_profile, output_variable, decomp)
        call output_profile%output_file%write_darray( &
              var_name=output_variable%name, &
              values=write_buffer_int32_3d, &
              decomp=decomp, &
              fill_value=FILL_VALUE_INT32, &
              frame=frame)
      else if (present(frame)) then
        call output_profile%output_file%inq_var_ndims(output_variable%name, ndims)
        call output_profile%output_file%put_var( &
              var_name=output_variable%name, &
              values=write_buffer_int32_3d, &
              start=[(1, i = 1, ndims - 1), frame])
      else
        call output_profile%output_file%put_var( &
              var_name=output_variable%name, &
              values=write_buffer_int32_3d)
      end if
    type is (aggregator_real32_0d_t)
      if (output_variable%reduction_method /= "none") then
        call cable_abort("Grid cell reductions are not supported for scalar variables", __FILE__, __LINE__)
      end if
      if (output_variable%distributed) then
        call cable_abort("Distributed writes are not supported for scalar variables", __FILE__, __LINE__)
      end if
      write_buffer_real32_0d => aggregator%aggregated_data
      if (restart_local) write_buffer_real32_0d => aggregator%source_data
      if (present(frame)) then
        call output_profile%output_file%inq_var_ndims(output_variable%name, ndims)
        call output_profile%output_file%put_var( &
              var_name=output_variable%name, &
              values=write_buffer_real32_0d, &
              start=[(1, i = 1, ndims - 1), frame])
      else
        call output_profile%output_file%put_var( &
              var_name=output_variable%name, &
              values=write_buffer_real32_0d)
      end if
    type is (aggregator_real32_1d_t)
      if (restart_local) then
        write_buffer_real32_1d => aggregator%source_data
      else if (output_variable%reduction_method == "none") then
        write_buffer_real32_1d => aggregator%aggregated_data
      else if (output_variable%reduction_method == "grid_cell_average") then
        call associate_temp_buffer_real32(output_variable, temp_buffer_real32_1d=write_buffer_real32_1d)
        call grid_cell_average( &
              input_array=aggregator%aggregated_data, &
              output_array=write_buffer_real32_1d, &
              landpt=landpt, &
              patch=patch)
      else if (output_variable%reduction_method == "first_patch_in_grid_cell") then
        call associate_temp_buffer_real32(output_variable, temp_buffer_real32_1d=write_buffer_real32_1d)
        call first_patch_in_grid_cell( &
              input_array=aggregator%aggregated_data, &
              output_array=write_buffer_real32_1d, &
              landpt=landpt)
      else
        call cable_abort("Invalid reduction method", __FILE__, __LINE__)
      end if
      if (output_variable%distributed) then
        call associate_decomp_real32(output_profile, output_variable, decomp)
        call output_profile%output_file%write_darray( &
              var_name=output_variable%name, &
              values=write_buffer_real32_1d, &
              decomp=decomp, &
              fill_value=FILL_VALUE_REAL32, &
              frame=frame)
      else if (present(frame)) then
        call output_profile%output_file%inq_var_ndims(output_variable%name, ndims)
        call output_profile%output_file%put_var( &
              var_name=output_variable%name, &
              values=write_buffer_real32_1d, &
              start=[(1, i = 1, ndims - 1), frame])
      else
        call output_profile%output_file%put_var( &
              var_name=output_variable%name, &
              values=write_buffer_real32_1d)
      end if
    type is (aggregator_real32_2d_t)
      if (restart_local) then
        write_buffer_real32_2d => aggregator%source_data
      else if (output_variable%reduction_method == "none") then
        write_buffer_real32_2d => aggregator%aggregated_data
      else if (output_variable%reduction_method == "grid_cell_average") then
        call associate_temp_buffer_real32(output_variable, temp_buffer_real32_2d=write_buffer_real32_2d)
        call grid_cell_average( &
              input_array=aggregator%aggregated_data, &
              output_array=write_buffer_real32_2d, &
              landpt=landpt, &
              patch=patch)
      else if (output_variable%reduction_method == "first_patch_in_grid_cell") then
        call associate_temp_buffer_real32(output_variable, temp_buffer_real32_2d=write_buffer_real32_2d)
        call first_patch_in_grid_cell( &
              input_array=aggregator%aggregated_data, &
              output_array=write_buffer_real32_2d, &
              landpt=landpt)
      else
        call cable_abort("Invalid reduction method", __FILE__, __LINE__)
      end if
      if (output_variable%distributed) then
        call associate_decomp_real32(output_profile, output_variable, decomp)
        call output_profile%output_file%write_darray( &
              var_name=output_variable%name, &
              values=write_buffer_real32_2d, &
              decomp=decomp, &
              fill_value=FILL_VALUE_REAL32, &
              frame=frame)
      else if (present(frame)) then
        call output_profile%output_file%inq_var_ndims(output_variable%name, ndims)
        call output_profile%output_file%put_var( &
              var_name=output_variable%name, &
              values=write_buffer_real32_2d, &
              start=[(1, i = 1, ndims - 1), frame])
      else
        call output_profile%output_file%put_var( &
              var_name=output_variable%name, &
              values=write_buffer_real32_2d)
      end if
    type is (aggregator_real32_3d_t)
      if (restart_local) then
        write_buffer_real32_3d => aggregator%source_data
      else if (output_variable%reduction_method == "none") then
        write_buffer_real32_3d => aggregator%aggregated_data
      else if (output_variable%reduction_method == "grid_cell_average") then
        call associate_temp_buffer_real32(output_variable, temp_buffer_real32_3d=write_buffer_real32_3d)
        call grid_cell_average( &
              input_array=aggregator%aggregated_data, &
              output_array=write_buffer_real32_3d, &
              landpt=landpt, &
              patch=patch)
      else if (output_variable%reduction_method == "first_patch_in_grid_cell") then
        call associate_temp_buffer_real32(output_variable, temp_buffer_real32_3d=write_buffer_real32_3d)
        call first_patch_in_grid_cell( &
              input_array=aggregator%aggregated_data, &
              output_array=write_buffer_real32_3d, &
              landpt=landpt)
      else
        call cable_abort("Invalid reduction method", __FILE__, __LINE__)
      end if
      if (output_variable%distributed) then
        call associate_decomp_real32(output_profile, output_variable, decomp)
        call output_profile%output_file%write_darray( &
              var_name=output_variable%name, &
              values=write_buffer_real32_3d, &
              decomp=decomp, &
              fill_value=FILL_VALUE_REAL32, &
              frame=frame)
      else if (present(frame)) then
        call output_profile%output_file%inq_var_ndims(output_variable%name, ndims)
        call output_profile%output_file%put_var( &
              var_name=output_variable%name, &
              values=write_buffer_real32_3d, &
              start=[(1, i = 1, ndims - 1), frame])
      else
        call output_profile%output_file%put_var( &
              var_name=output_variable%name, &
              values=write_buffer_real32_3d)
      end if
    type is (aggregator_real64_0d_t)
      if (output_variable%reduction_method /= "none") then
        call cable_abort("Grid cell reductions are not supported for scalar variables", __FILE__, __LINE__)
      end if
      if (output_variable%distributed) then
        call cable_abort("Distributed writes are not supported for scalar variables", __FILE__, __LINE__)
      end if
      write_buffer_real64_0d => aggregator%aggregated_data
      if (restart_local) write_buffer_real64_0d => aggregator%source_data
      if (present(frame)) then
        call output_profile%output_file%inq_var_ndims(output_variable%name, ndims)
        call output_profile%output_file%put_var( &
              var_name=output_variable%name, &
              values=write_buffer_real64_0d, &
              start=[(1, i = 1, ndims - 1), frame])
      else
        call output_profile%output_file%put_var( &
              var_name=output_variable%name, &
              values=write_buffer_real64_0d)
      end if
    type is (aggregator_real64_1d_t)
      if (restart_local) then
        write_buffer_real64_1d => aggregator%source_data
      else if (output_variable%reduction_method == "none") then
        write_buffer_real64_1d => aggregator%aggregated_data
      else if (output_variable%reduction_method == "grid_cell_average") then
        call associate_temp_buffer_real64(output_variable, temp_buffer_real64_1d=write_buffer_real64_1d)
        call grid_cell_average( &
              input_array=aggregator%aggregated_data, &
              output_array=write_buffer_real64_1d, &
              landpt=landpt, &
              patch=patch)
      else if (output_variable%reduction_method == "first_patch_in_grid_cell") then
        call associate_temp_buffer_real64(output_variable, temp_buffer_real64_1d=write_buffer_real64_1d)
        call first_patch_in_grid_cell( &
              input_array=aggregator%aggregated_data, &
              output_array=write_buffer_real64_1d, &
              landpt=landpt)
      else
        call cable_abort("Invalid reduction method", __FILE__, __LINE__)
      end if
      if (output_variable%distributed) then
        call associate_decomp_real64(output_profile, output_variable, decomp)
        call output_profile%output_file%write_darray( &
              var_name=output_variable%name, &
              values=write_buffer_real64_1d, &
              decomp=decomp, &
              fill_value=FILL_VALUE_REAL64, &
              frame=frame)
      else if (present(frame)) then
        call output_profile%output_file%inq_var_ndims(output_variable%name, ndims)
        call output_profile%output_file%put_var( &
              var_name=output_variable%name, &
              values=write_buffer_real64_1d, &
              start=[(1, i = 1, ndims - 1), frame])
      else
        call output_profile%output_file%put_var( &
              var_name=output_variable%name, &
              values=write_buffer_real64_1d)
      end if
    type is (aggregator_real64_2d_t)
      if (restart_local) then
        write_buffer_real64_2d => aggregator%source_data
      else if (output_variable%reduction_method == "none") then
        write_buffer_real64_2d => aggregator%aggregated_data
      else if (output_variable%reduction_method == "grid_cell_average") then
        call associate_temp_buffer_real64(output_variable, temp_buffer_real64_2d=write_buffer_real64_2d)
        call grid_cell_average( &
              input_array=aggregator%aggregated_data, &
              output_array=write_buffer_real64_2d, &
              landpt=landpt, &
              patch=patch)
      else if (output_variable%reduction_method == "first_patch_in_grid_cell") then
        call associate_temp_buffer_real64(output_variable, temp_buffer_real64_2d=write_buffer_real64_2d)
        call first_patch_in_grid_cell( &
              input_array=aggregator%aggregated_data, &
              output_array=write_buffer_real64_2d, &
              landpt=landpt)
      else
        call cable_abort("Invalid reduction method", __FILE__, __LINE__)
      end if
      if (output_variable%distributed) then
        call associate_decomp_real64(output_profile, output_variable, decomp)
        call output_profile%output_file%write_darray( &
              var_name=output_variable%name, &
              values=write_buffer_real64_2d, &
              decomp=decomp, &
              fill_value=FILL_VALUE_REAL64, &
              frame=frame)
      else if (present(frame)) then
        call output_profile%output_file%inq_var_ndims(output_variable%name, ndims)
        call output_profile%output_file%put_var( &
              var_name=output_variable%name, &
              values=write_buffer_real64_2d, &
              start=[(1, i = 1, ndims - 1), frame])
      else
        call output_profile%output_file%put_var( &
              var_name=output_variable%name, &
              values=write_buffer_real64_2d)
      end if
    type is (aggregator_real64_3d_t)
      if (restart_local) then
        write_buffer_real64_3d => aggregator%source_data
      else if (output_variable%reduction_method == "none") then
        write_buffer_real64_3d => aggregator%aggregated_data
      else if (output_variable%reduction_method == "grid_cell_average") then
        call associate_temp_buffer_real64(output_variable, temp_buffer_real64_3d=write_buffer_real64_3d)
        call grid_cell_average( &
              input_array=aggregator%aggregated_data, &
              output_array=write_buffer_real64_3d, &
              landpt=landpt, &
              patch=patch)
      else if (output_variable%reduction_method == "first_patch_in_grid_cell") then
        call associate_temp_buffer_real64(output_variable, temp_buffer_real64_3d=write_buffer_real64_3d)
        call first_patch_in_grid_cell( &
              input_array=aggregator%aggregated_data, &
              output_array=write_buffer_real64_3d, &
              landpt=landpt)
      else
        call cable_abort("Invalid reduction method", __FILE__, __LINE__)
      end if
      if (output_variable%distributed) then
        call associate_decomp_real64(output_profile, output_variable, decomp)
        call output_profile%output_file%write_darray( &
              var_name=output_variable%name, &
              values=write_buffer_real64_3d, &
              decomp=decomp, &
              fill_value=FILL_VALUE_REAL64, &
              frame=frame)
      else if (present(frame)) then
        call output_profile%output_file%inq_var_ndims(output_variable%name, ndims)
        call output_profile%output_file%put_var( &
              var_name=output_variable%name, &
              values=write_buffer_real64_3d, &
              start=[(1, i = 1, ndims - 1), frame])
      else
        call output_profile%output_file%put_var( &
              var_name=output_variable%name, &
              values=write_buffer_real64_3d)
      end if
    class default
      call cable_abort("Unexpected aggregator type", __FILE__, __LINE__)
    end select

  end subroutine write_variable

end module
