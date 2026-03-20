submodule (cable_output_mod:cable_output_internal) cable_output_implementation
  use cable_common_module, only: filename
  use cable_io_vars_module, only: metgrid
  use cable_io_vars_module, only: output
  use cable_io_vars_module, only: check
  use cable_io_vars_module, only: ON_TIMESTEP
  use cable_io_vars_module, only: ON_WRITE
  use cable_netcdf_mod, only: cable_netcdf_create_file
  use cable_netcdf_mod, only: CABLE_NETCDF_IOTYPE_CLASSIC
  use cable_timing_mod, only: frequency_matches => cable_timing_frequency_matches
  use cable_array_utils_mod, only: array_eq
  implicit none

  !> This flag forces time averaging to computed by summing the diagnostics at
  !! each accumulation step and then dividing by the number of samples at write
  !! time, rather than computing the average incrementally. This is required to
  !! demonstrate bitwise reproducibility with the previous output module.
  logical, parameter :: normalised_averaging = .true.

  !> Global output stream instance for the main cable output file.
  type(cable_output_stream_t), allocatable :: global_output_stream

  !> Registered output variables.
  type(cable_output_variable_t), allocatable :: registered_output_variables(:)

contains

  module subroutine cable_output_implementation_init()

    call cable_output_decomp_init()
    call cable_output_reduction_buffers_init()

  end subroutine

  module subroutine cable_output_implementation_end()

    if (allocated(global_output_stream%output_file)) call global_output_stream%output_file%close()

    deallocate(global_output_stream)

    call cable_output_reduction_buffers_free()
    call cable_output_decomp_free()

  end subroutine

  module subroutine cable_output_implementation_register_output_variables(output_variables)
    type(cable_output_variable_t), dimension(:), intent(in) :: output_variables
    integer :: i

    do i = 1, size(output_variables)
      associate(output_var => output_variables(i))
        if (count(output_var%field_name == output_variables(:)%field_name) > 1) then
          call cable_abort("Duplicate field_name found: " // output_var%field_name, __FILE__, __LINE__)
        end if
        if (all(output_var%reduction_method /= allowed_reduction_methods)) then
          call cable_abort("Invalid reduction method for variable " // output_var%field_name, __FILE__, __LINE__)
        end if
        if (all(output_var%aggregation_method /= allowed_aggregation_methods)) then
          call cable_abort("Invalid aggregation method for variable " // output_var%field_name, __FILE__, __LINE__)
        end if
        if (.not. allocated(output_var%aggregator)) then
          call cable_abort("Undefined aggregator for variable " // output_var%field_name, __FILE__, __LINE__)
        end if
        if (.not. allocated(output_var%data_shape) .and. output_var%aggregator%rank() /= 0) then
          call cable_abort("Data shape does not match aggregator shape for variable " // output_var%field_name, __FILE__, __LINE__)
        end if
        if (allocated(output_var%data_shape)) then
          if (.not. array_eq(output_var%data_shape(:)%size(), output_var%aggregator%shape())) then
            call cable_abort("Data shape does not match aggregator shape for variable " // output_var%field_name, __FILE__, __LINE__)
          end if
        end if
      end associate
    end do

    registered_output_variables = output_variables

  end subroutine cable_output_implementation_register_output_variables

  module subroutine cable_output_implementation_init_streams()
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

    global_output_stream = cable_output_stream_t( &
      sampling_frequency=output%averaging, &
      grid_type=grid_type, &
      file_name=filename%out, &
      output_file=cable_netcdf_create_file(filename%out, iotype=CABLE_NETCDF_IOTYPE_CLASSIC), &
      coordinate_variables=coordinate_variables_list(grid_type), &
      output_variables=pack(registered_output_variables, registered_output_variables(:)%active) &
    )

    do i = 1, size(global_output_stream%output_variables)
      associate(output_var => global_output_stream%output_variables(i))
        if (count(output_var%get_netcdf_name() == global_output_stream%output_variables(:)%get_netcdf_name()) > 1) then
          call cable_abort("Duplicate netCDF variable name in output stream: " // output_var%get_netcdf_name(), __FILE__, __LINE__)
        end if
        if (frequency_is_greater_than(global_output_stream%sampling_frequency, output_var%accumulation_frequency)) then
          call cable_abort( &
            "Output stream sampling frequency '" // global_output_stream%sampling_frequency // &
            "' is greater than accumulation frequency '" // output_var%accumulation_frequency // &
            "' for variable '" // output_var%field_name, __FILE__, __LINE__ &
          )
        end if
        if (output_var%patchout) output_var%reduction_method = "none"
        if (global_output_stream%sampling_frequency == "all") output_var%aggregation_method = "point"
      end associate
    end do

    call cable_output_define_variables(global_output_stream)

    call global_output_stream%output_file%end_def()

    do i = 1, size(global_output_stream%coordinate_variables)
      associate(coordinate_variable => global_output_stream%coordinate_variables(i))
        call coordinate_variable%aggregator%init(method="point")
        call coordinate_variable%aggregator%accumulate()
        call cable_output_write_variable(global_output_stream, coordinate_variable)
        call coordinate_variable%aggregator%reset()
      end associate
    end do

    do i = 1, size(global_output_stream%output_variables)
      associate(output_variable => global_output_stream%output_variables(i))
        if (normalised_averaging .and. output_variable%aggregation_method == "mean") then
          call output_variable%aggregator%init(method="sum")
        else
          call output_variable%aggregator%init(method=output_variable%aggregation_method)
        end if
      end associate
    end do

  end subroutine

  module subroutine cable_output_implementation_write_parameters(time_index, patch, landpt)
    integer, intent(in) :: time_index
    type(patch_type), intent(in) :: patch(:)
    type(land_type), intent(in) :: landpt(:)
    integer :: i

    do i = 1, size(global_output_stream%output_variables)
      associate(output_variable => global_output_stream%output_variables(i))
        if (.not. output_variable%parameter) cycle
        call check_variable_range(output_variable, time_index)
        call output_variable%aggregator%accumulate( &
          scale=output_variable%scale_by, &
          div=output_variable%divide_by, &
          offset=output_variable%offset_by &
        )
        call cable_output_write_variable(global_output_stream, output_variable, patch, landpt)
        call output_variable%aggregator%reset()
      end associate
    end do

  end subroutine

  module subroutine cable_output_implementation_update(time_index, dels, met)
    integer, intent(in) :: time_index
    real, intent(in) :: dels
    type(met_type), intent(in) :: met
    real :: current_time
    integer :: i

    if (check%ranges == ON_TIMESTEP) then
      do i = 1, size(global_output_stream%output_variables)
        call check_variable_range(global_output_stream%output_variables(i), time_index, met)
      end do
    end if

    do i = 1, size(global_output_stream%output_variables)
      associate(output_variable => global_output_stream%output_variables(i))
        if (frequency_matches(dels, time_index, output_variable%accumulation_frequency)) then
          call output_variable%aggregator%accumulate( &
            scale=output_variable%scale_by, &
            div=output_variable%divide_by, &
            offset=output_variable%offset_by &
          )
        end if
      end associate
    end do

  end subroutine

  module subroutine cable_output_implementation_write(time_index, dels, met, patch, landpt)
    integer, intent(in) :: time_index
    real, intent(in) :: dels
    type(met_type), intent(in) :: met
    type(patch_type), intent(in) :: patch(:)
    type(land_type), intent(in) :: landpt(:)

    real :: current_time
    integer :: i

    if (frequency_matches(dels, time_index, global_output_stream%sampling_frequency)) then

      do i = 1, size(global_output_stream%output_variables)
        associate(output_variable => global_output_stream%output_variables(i))
          if (output_variable%parameter) cycle
          if (check%ranges == ON_WRITE) call check_variable_range(output_variable, time_index, met)
          if (normalised_averaging .and. output_variable%aggregation_method == "mean") then
            call output_variable%aggregator%div(real(output_variable%aggregator%counter))
          end if
          call cable_output_write_variable(global_output_stream, output_variable, patch, landpt, frame=global_output_stream%frame + 1)
          call output_variable%aggregator%reset()
        end associate
      end do

      current_time = time_index * dels

      if (global_output_stream%sampling_frequency == "all") then
        call global_output_stream%output_file%put_var("time", current_time, start=[global_output_stream%frame + 1])
      else
        call global_output_stream%output_file%put_var("time", (current_time + global_output_stream%previous_write_time) / 2.0, start=[global_output_stream%frame + 1])
      end if

      global_output_stream%previous_write_time = current_time
      global_output_stream%frame = global_output_stream%frame + 1

    end if

  end subroutine cable_output_implementation_write

  module subroutine cable_output_implementation_write_restart(current_time)
    real, intent(in) :: current_time !! Current simulation time

    type(cable_output_stream_t), allocatable :: restart_output_stream
    integer :: i

    restart_output_stream = cable_output_stream_t( &
      sampling_frequency="none", &
      grid_type="restart", &
      file_name=filename%restart_out, &
      output_file=cable_netcdf_create_file(filename%restart_out, iotype=CABLE_NETCDF_IOTYPE_CLASSIC), &
      coordinate_variables=coordinate_variables_list(grid_type="restart"), &
      output_variables=pack(registered_output_variables, registered_output_variables(:)%restart) &
    )

    call cable_output_define_variables(restart_output_stream, restart=.true.)

    call restart_output_stream%output_file%end_def()

    call restart_output_stream%output_file%put_var("time", [current_time])

    do i = 1, size(restart_output_stream%coordinate_variables)
      call cable_output_write_variable(restart_output_stream, restart_output_stream%coordinate_variables(i), restart=.true.)
    end do

    do i = 1, size(restart_output_stream%output_variables)
      call cable_output_write_variable(restart_output_stream, restart_output_stream%output_variables(i), restart=.true.)
    end do

    call restart_output_stream%output_file%close()

  end subroutine cable_output_implementation_write_restart

end submodule cable_output_implementation
