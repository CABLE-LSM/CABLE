! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

submodule (cable_output_mod:cable_output_internal_smod) cable_output_impl_smod
  !! Implementation of the public interface procedures in [[cable_output_mod]].
  use cable_common_module, only: filename
  use cable_io_vars_module, only: metgrid
  use cable_io_vars_module, only: output
  use cable_io_vars_module, only: check
  use cable_io_vars_module, only: ON_TIMESTEP
  use cable_io_vars_module, only: ON_WRITE
  use cable_netcdf_mod, only: cable_netcdf_create_file
  use cable_netcdf_mod, only: CABLE_NETCDF_IOTYPE_CLASSIC
  use cable_timing_mod, only: frequency_matches => cable_timing_frequency_matches
  use cable_timing_mod, only: frequency_is_greater_than => cable_timing_frequency_is_greater_than
  use cable_array_utils_mod, only: array_eq
  implicit none

  !> This flag forces time averaging to computed by summing the diagnostics at
  !! each accumulation step and then dividing by the number of samples at write
  !! time, rather than computing the average incrementally. This is required to
  !! demonstrate bitwise reproducibility with the previous output module.
  logical, parameter :: normalised_averaging = .true.

  !> Global output stream instance for the main cable output file.
  type(cable_output_stream_t) :: global_output_stream

  !> Registered output variables.
  type(cable_output_variable_t), allocatable :: registered_output_variables(:)

contains

  module subroutine cable_output_impl_init()
    !* Module initialisation procedure for `cable_output_mod`.
    !
    ! This procedure must be called before any other procedures in
    ! `cable_output_mod`.

    call cable_output_decomp_init()
    call cable_output_reduction_buffers_init()

  end subroutine

  module subroutine cable_output_impl_end()
    !* Module finalization procedure for `cable_output_mod`.
    !
    ! This procedure should be called at the end of the simulation after all
    ! output has been written.

    if (allocated(global_output_stream%output_file)) call global_output_stream%output_file%close()

    call cable_output_reduction_buffers_free()
    call cable_output_decomp_free()

  end subroutine

  module subroutine cable_output_impl_register_output_variables(output_variables)
    !* Registers output variables with the output module. Note that
    ! registering an output variable does not necessarily mean that the variable
    ! will be written to an output stream - this can depend on whether the
    ! output variable is active, or if it is a restart variable. Output
    ! variables should be registered if their associated diagnostic working
    ! variables are initialised in the model as this can help provide the
    ! information on the diagnostics which are available.
    type(cable_output_variable_t), dimension(:), intent(in) :: output_variables
      !! An array of output variable definitions to be registered.
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
        if (allocated(output_var%range)) then
          if (output_var%range(1) >= output_var%range(2)) then
            call cable_abort("Invalid range specified for variable " // output_var%field_name, __FILE__, __LINE__)
          end if
        end if

      end associate
    end do

    registered_output_variables = output_variables

  end subroutine cable_output_impl_register_output_variables

  module subroutine cable_output_impl_init_streams(dels)
    !! Initialise output streams based on the current output configuration.
    real, intent(in) :: dels !! The current time step size in seconds.
    integer :: i
    character(32) :: grid_type

    if (.not. allocated(registered_output_variables)) then
      call cable_abort("Output variables must be registered before initialising output streams.", __FILE__, __LINE__)
    end if

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
        if (frequency_is_greater_than(global_output_stream%sampling_frequency, output_var%accumulation_frequency, dels)) then
          call cable_abort( &
            "Output stream sampling frequency '" // global_output_stream%sampling_frequency // &
            "' is greater than accumulation frequency '" // output_var%accumulation_frequency // &
            "' for variable '" // output_var%field_name, __FILE__, __LINE__ &
          )
        end if
        if (output_var%patchout) output_var%reduction_method = "none"
        if (global_output_stream%sampling_frequency == "all") output_var%aggregation_method = "point"
        if (allocated(output_var%range)) output_var%range_native = (output_var%range - output_var%offset_by) * output_var%divide_by / output_var%scale_by
      end associate
    end do

    call cable_output_define_stream(global_output_stream)

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

  module subroutine cable_output_impl_update(time_index, dels, met)
    !* Updates the time aggregation accumulation for any output variables that
    ! are active in an output stream with an accumulation frequency that matches
    ! the current time step.
    integer, intent(in) :: time_index !! The current time step index in the simulation.
    real, intent(in) :: dels !! The current time step size in seconds.
    type(met_type), intent(in) :: met
      !* Met variables at the current time step to provide informative error
      ! messages for CABLE range checks.
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

  module subroutine cable_output_impl_write(time_index, dels, met, patch, landpt)
    !* Writes output variables to disk for any output streams with a sampling
    ! frequency that matches the current time step.
    integer, intent(in) :: time_index !! The current time step index in the simulation.
    real, intent(in) :: dels !! The current time step size in seconds.
    type(met_type), intent(in) :: met
      !* Met variables at the current time step to provide informative error
      ! messages for CABLE range checks.
    type(patch_type), intent(in) :: patch(:)
      !! The patch type instance for performing grid reductions over the patch dimension if required.
    type(land_type), intent(in) :: landpt(:)
      !! The land type instance for performing grid reductions over the patch dimension if required.
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

  end subroutine cable_output_impl_write

  module subroutine cable_output_impl_write_parameters(time_index, patch, landpt)
    !* Writes non-time varying parameter output variables to disk. This is
    ! done on the first time step of the simulation after the output streams
    ! have been initialised.
    integer, intent(in) :: time_index !! The current time step index in the simulation.
    type(patch_type), intent(in) :: patch(:)
      !! The patch type instance for performing grid reductions over the patch dimension if required.
    type(land_type), intent(in) :: landpt(:)
      !! The land type instance for performing grid reductions over the patch dimension if required.
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

  module subroutine cable_output_impl_write_restart(current_time)
    !* Writes variables to the CABLE restart file. This is done at the end of
    ! the simulation.
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

    call cable_output_define_stream(restart_output_stream, restart=.true.)

    call restart_output_stream%output_file%end_def()

    call restart_output_stream%output_file%put_var("time", [current_time])

    do i = 1, size(restart_output_stream%coordinate_variables)
      call cable_output_write_variable(restart_output_stream, restart_output_stream%coordinate_variables(i), restart=.true.)
    end do

    do i = 1, size(restart_output_stream%output_variables)
      call cable_output_write_variable(restart_output_stream, restart_output_stream%output_variables(i), restart=.true.)
    end do

    call restart_output_stream%output_file%close()

  end subroutine cable_output_impl_write_restart

end submodule cable_output_impl_smod
