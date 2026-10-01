! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

submodule (cable_output_mod:cable_output_common_smod) cable_output_impl_smod
  !! Implementation of the public interface procedures in [[cable_output_mod]].
  use cable_common_module, only: filename
  use cable_io_vars_module, only: check
  use cable_io_vars_module, only: ON_TIMESTEP
  use cable_io_vars_module, only: ON_WRITE
  use cable_netcdf_mod, only: cable_netcdf_create_file
  use cable_netcdf_mod, only: CABLE_NETCDF_IOTYPE_CLASSIC
  use cable_netcdf_mod, only: CABLE_NETCDF_IOTYPE_NETCDF4C
  use cable_timing_mod, only: frequency_matches => cable_timing_frequency_matches
  use cable_array_utils_mod, only: array_eq
  implicit none

  !> This flag forces time averaging to computed by summing the diagnostics at
  !! each accumulation step and then dividing by the number of samples at write
  !! time, rather than computing the average incrementally. This is required to
  !! demonstrate bitwise reproducibility with the previous output module.
  logical, parameter :: normalised_averaging = .true.

  !> The output files being written, one element per file. The configuration
  !! decides how many there are.
  type(cable_output_stream_t), allocatable :: output_streams(:)

  !> Registered output variable definitions.
  type(cable_output_variable_definition_t), allocatable :: registered_definitions(:)

contains

  module subroutine cable_output_impl_init()
    !* Module initialisation procedure for `cable_output_mod`.
    !
    ! This procedure must be called before any other procedures in
    ! `cable_output_mod`.

    ! Both set up shared work areas used by every stream: the parallel I/O
    ! decompositions, and the temporary arrays that hold tile-reduced data.
    call cable_output_decomp_init()
    call cable_output_reduction_buffers_init()

  end subroutine

  module subroutine cable_output_impl_end()
    !* Module finalization procedure for `cable_output_mod`.
    !
    ! This procedure should be called at the end of the simulation after all
    ! output has been written.

    integer :: stream_index

    ! Close every file that was opened. Closing also flushes anything still
    ! buffered. A stream whose file was never created is skipped.
    if (allocated(output_streams)) then
      do stream_index = 1, size(output_streams)
        if (allocated(output_streams(stream_index)%output_file)) call output_streams(stream_index)%output_file%close()
      end do
    end if

    ! Release the shared work areas last, after the files that used them are closed.
    call cable_output_reduction_buffers_free()
    call cable_output_decomp_free()

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
    integer :: stream_index, i

    ! No streams means no output was requested, so there is nothing to do.
    if (.not. allocated(output_streams)) return

    do stream_index = 1, size(output_streams)
      associate(output_stream => output_streams(stream_index))
        ! Optional check that each value is physically sensible, done every
        ! step so that a bad value is caught as soon as it appears.
        if (check%ranges == ON_TIMESTEP) then
          do i = 1, size(output_stream%output_variables)
            call check_variable_range(output_stream%output_variables(i), time_index, met)
          end do
        end if

        ! Accumulation is decided per variable, not per file. A variable is
        ! sampled when this step matches how often the model itself updates
        ! it (every step for most, once a day for a few), which is independent
        ! of how often the file is written.
        do i = 1, size(output_stream%output_variables)
          associate(output_variable => output_stream%output_variables(i))
            if (frequency_matches(dels, time_index, output_variable%accumulation_frequency)) then
              ! The aggregator reads the model's working variable directly
              ! (through the pointer it holds), applies the unit conversion
              ! and folds the result into the running mean, sum, max, min
              ! or latest value, according to its aggregation method.
              call output_variable%aggregator%accumulate( &
                scale=output_variable%scale_by, &
                div=output_variable%divide_by, &
                offset=output_variable%offset_by &
              )
            end if
          end associate
        end do
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
    integer :: stream_index, i

    if (.not. allocated(output_streams)) return

    do stream_index = 1, size(output_streams)
      associate(output_stream => output_streams(stream_index))
        ! Each file has its own clock. `cycle` moves on to the next file unless
        ! this step is the end of this file's period (for example the last
        ! step of a month for a monthly file).
        if (.not. frequency_matches(dels, time_index, output_stream%sampling_frequency)) cycle

        do i = 1, size(output_stream%output_variables)
          associate(output_variable => output_stream%output_variables(i))
            ! Parameters have no time axis and were written once at start-up.
            if (output_variable%parameter) cycle
            if (check%ranges == ON_WRITE) call check_variable_range(output_variable, time_index, met)
            ! A mean is accumulated as a running sum (see start_stream) and
            ! turned into a mean only now, by dividing by how many samples
            ! were added. This order of operations is what makes results match
            ! the original output module bit for bit.
            if (normalised_averaging .and. output_variable%aggregation_method == "mean") then
              call output_variable%aggregator%div(real(output_variable%aggregator%counter))
            end if
            ! The write routine applies any tile reduction, then hands the data
            ! to the NetCDF layer. `frame` is the record along the time axis.
            call cable_output_write_variable(output_stream, output_variable, patch, landpt, frame=output_stream%frame + 1)
            ! Start the next period from a clean aggregator.
            call output_variable%aggregator%reset()
          end associate
        end do

        ! The time value recorded with this record. For every-step output it
        ! is simply the current time. For longer periods it is the midpoint of
        ! the period: the average of this write time and the previous one
        ! (which is 0 before the first write).
        current_time = time_index * dels

        if (output_stream%sampling_frequency == "timestep") then
          call output_stream%output_file%put_var("time", current_time, start=[output_stream%frame + 1])
        else
          call output_stream%output_file%put_var("time", (current_time + output_stream%previous_write_time) / 2.0, &
            start=[output_stream%frame + 1])
        end if

        ! Remember when this write happened, for the next midpoint, and move
        ! on to the next record.
        output_stream%previous_write_time = current_time
        output_stream%frame = output_stream%frame + 1
      end associate
    end do

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
    integer :: stream_index, i

    if (.not. allocated(output_streams)) return

    ! Parameters do not change with time, so they are sampled once and written
    ! immediately, instead of going through the per-step accumulate/write cycle.
    do stream_index = 1, size(output_streams)
      associate(output_stream => output_streams(stream_index))
        do i = 1, size(output_stream%output_variables)
          associate(output_variable => output_stream%output_variables(i))
            ! Only parameters here; the other variables wait for the time loop.
            if (.not. output_variable%parameter) cycle
            call check_variable_range(output_variable, time_index)
            ! One accumulate with the aggregator's "instant" behaviour copies
            ! the converted value into the aggregator so that the write
            ! routine can find it in the same place as for any other variable.
            call output_variable%aggregator%accumulate( &
              scale=output_variable%scale_by, &
              div=output_variable%divide_by, &
              offset=output_variable%offset_by &
            )
            ! No `frame` argument: a variable without a time axis is written whole.
            call cable_output_write_variable(output_stream, output_variable, patch, landpt)
            call output_variable%aggregator%reset()
          end associate
        end do
      end associate
    end do

  end subroutine

  module subroutine cable_output_impl_write_restart(current_time)
    !* Writes variables to the CABLE restart file. This is done at the end of
    ! the simulation.
    real, intent(in) :: current_time !! Current simulation time

    type(cable_output_stream_t), allocatable :: restart_output_stream
    type(cable_output_variable_t), allocatable :: restart_variables(:)
    integer :: i, restart_count, restart_index

    if (.not. allocated(registered_definitions)) then
      call cable_abort("Output variables must be registered before writing a restart file.", __FILE__, __LINE__)
    end if

    ! Restart variables are those definitions with a restart name. They are
    ! written as they are, whatever the output configuration selects.
    ! Count them first so the array can be allocated to the right size, then
    ! fill it in a second pass.
    restart_count = count(len_trim(registered_definitions(:)%restart_name) > 0)
    allocate(restart_variables(restart_count))
    restart_index = 0
    do i = 1, size(registered_definitions)
      if (len_trim(registered_definitions(i)%restart_name) == 0) cycle
      restart_index = restart_index + 1
      restart_variables(restart_index) = cable_output_restart_variable(registered_definitions(i))
    end do

    ! The restart file is a stream with no frequency: it is written once, at
    ! the end. It always uses the classic NetCDF format, and the "restart" grid
    ! type gives it the compact layout (patch and land dimensions as they are in
    ! memory) rather than a lat/lon grid.
    restart_output_stream = cable_output_stream_t( &
      sampling_frequency="none", &
      grid_type="restart", &
      file_name=filename%restart_out, &
      output_file=cable_netcdf_create_file(filename%restart_out, iotype=CABLE_NETCDF_IOTYPE_CLASSIC), &
      coordinate_variables=coordinate_variables_list(grid_type="restart"), &
      output_variables=restart_variables &
    )

    ! Define the file's contents, then leave define mode so data can be written.
    call cable_output_define_stream(restart_output_stream, restart=.true.)

    call restart_output_stream%output_file%end_def()

    ! A restart file has a single time value: when the run stopped.
    call restart_output_stream%output_file%put_var("time", [current_time])

    do i = 1, size(restart_output_stream%coordinate_variables)
      call cable_output_write_variable(restart_output_stream, restart_output_stream%coordinate_variables(i), restart=.true.)
    end do

    ! `restart=.true.` makes the write routine take the model's current value
    ! directly, skipping aggregation and tile reduction.
    do i = 1, size(restart_output_stream%output_variables)
      call cable_output_write_variable(restart_output_stream, restart_output_stream%output_variables(i), restart=.true.)
    end do

    call restart_output_stream%output_file%close()

  end subroutine cable_output_impl_write_restart

  ! ------------------------------------------------------------------ output configurations

  module subroutine cable_output_impl_register_definitions(definitions)
    !* Registers output variable definitions with the output module.
    type(cable_output_variable_definition_t), dimension(:), intent(in) :: definitions
      !! The output variable definitions to register.
    integer :: i, j

    ! These checks protect the rest of the output system from bad definitions.
    ! They are mistakes in CABLE's own catalogue or bindings, not in the user's
    ! configuration, so each one stops the run.
    do i = 1, size(definitions)
      associate(definition => definitions(i))
        if (len_trim(definition%field_name) == 0) then
          call cable_abort("Output variable definition without a field_name", __FILE__, __LINE__)
        end if
        ! Names must be unique. Comparing each definition only with those after
        ! it checks every pair once.
        do j = i + 1, size(definitions)
          if (definition%field_name == definitions(j)%field_name) then
            call cable_abort("Duplicate field_name found: " // trim(definition%field_name), __FILE__, __LINE__)
          end if
          if (len_trim(definition%restart_name) > 0 .and. definition%restart_name == definitions(j)%restart_name) then
            call cable_abort("Duplicate restart_name found: " // trim(definition%restart_name), __FILE__, __LINE__)
          end if
        end do
        ! A range with its bounds the wrong way round would flag every value.
        if (definition%range(1) >= definition%range(2)) then
          call cable_abort("Invalid range specified for variable " // trim(definition%field_name), __FILE__, __LINE__)
        end if
        ! The remaining checks look at the model variable itself, which only
        ! exists if the variable is available in this model configuration.
        if (definition%available) then
          if (.not. allocated(definition%aggregator)) then
            call cable_abort("Undefined aggregator for variable " // trim(definition%field_name), __FILE__, __LINE__)
          end if
          if (.not. allocated(definition%data_shape)) then
            call cable_abort("Undefined data shape for variable " // trim(definition%field_name), __FILE__, __LINE__)
          end if
          ! The shape the catalogue promises must equal the shape of the
          ! array the binding points at, or writes would read the wrong amount.
          if (.not. array_eq(definition%data_shape(:)%size(), definition%aggregator%shape())) then
            call cable_abort("Data shape does not match aggregator shape for variable " // &
              trim(definition%field_name), __FILE__, __LINE__)
          end if
        end if
      end associate
    end do

    ! Keep a copy for building streams and for the restart file.
    registered_definitions = definitions

  end subroutine cable_output_impl_register_definitions

  module function cable_output_impl_build_streams(config, definitions, grid_type) result(streams)
    !* Works out which files the configuration produces and which variables go
    ! into each, without creating any files.
    type(cable_output_config_t), intent(in) :: config
    type(cable_output_variable_definition_t), intent(in) :: definitions(:)
    character(len=*), intent(in) :: grid_type
    type(cable_output_stream_t), allocatable :: streams(:)
    type(cable_output_variable_t), allocatable :: stream_variables(:)
    integer, allocatable :: use_indexes(:)
    integer :: stream_index, use_index, file_count, file_index, variable_index

    ! Reject an unknown grid type now, rather than when the first file is defined.
    if (all(allowed_grid_types /= grid_type)) then
      call cable_abort("Invalid output grid type '" // trim(grid_type) // "'", __FILE__, __LINE__)
    end if

    ! Pass 1: count the files, so the result can be allocated once.
    ! `use_indexes` lists which variable uses belong to this stream: `pack`
    ! keeps the positions 1..n for which the use's stream number matches.
    file_count = 0
    do stream_index = 1, size(config%streams)
      associate(stream => config%streams(stream_index))
        use_indexes = pack([(use_index, use_index = 1, size(config%variables))], &
          config%variables(:)%stream_id == stream%stream_id)
        ! A stream that nothing was directed to makes no file.
        if (size(use_indexes) == 0) cycle
        ! One file for the whole stream, or one per variable if asked.
        if (stream%separate_file_per_variable) then
          file_count = file_count + size(use_indexes)
        else
          file_count = file_count + 1
        end if
      end associate
    end do

    ! Pass 2: build each file, in the same order as the count above.
    allocate(streams(file_count))
    file_index = 0
    do stream_index = 1, size(config%streams)
      associate(stream => config%streams(stream_index))
        use_indexes = pack([(use_index, use_index = 1, size(config%variables))], &
          config%variables(:)%stream_id == stream%stream_id)
        if (size(use_indexes) == 0) cycle
        if (stream%separate_file_per_variable) then
          ! Separate files: a file of its own for each variable use. The array
          ! `stream_variables` holds just that one variable each time round.
          do variable_index = 1, size(use_indexes)
            file_index = file_index + 1
            allocate(stream_variables(1))
            stream_variables(1) = variable_from_use(config%variables(use_indexes(variable_index)), definitions)
            ! Each variable gets a file named after it, next to where file_name would be.
            streams(file_index) = stream_from_config(stream, grid_type, stream_variables, &
              directory_of(stream%file_name) // trim(stream_variables(1)%get_netcdf_name()) // ".nc")
            deallocate(stream_variables)
          end do
        else
          ! Normal case: every variable of the stream goes into one file.
          allocate(stream_variables(size(use_indexes)))
          do variable_index = 1, size(use_indexes)
            stream_variables(variable_index) = variable_from_use(config%variables(use_indexes(variable_index)), definitions)
          end do
          file_index = file_index + 1
          streams(file_index) = stream_from_config(stream, grid_type, stream_variables, trim(stream%file_name))
          deallocate(stream_variables)
        end if
      end associate
    end do

  end function cable_output_impl_build_streams

  pure function directory_of(path) result(directory)
    !! The directory part of `path`, including the final `/`, or nothing if there is none.
    character(len=*), intent(in) :: path
    character(:), allocatable :: directory
    integer :: last_slash

    ! `back=.true.` finds the LAST "/". Keeping everything up to and including
    ! it gives the directory; with no "/" the result is empty, meaning the
    ! current directory.
    last_slash = index(trim(path), "/", back=.true.)
    directory = path(1:last_slash)
  end function directory_of

  function variable_from_use(variable_use, definitions) result(output_variable)
    !! The output variable for one use of a definition in the configuration.
    type(cable_output_variable_config_t), intent(in) :: variable_use
    type(cable_output_variable_definition_t), intent(in) :: definitions(:)
    type(cable_output_variable_t) :: output_variable
    integer :: definition_index

    ! Find the definition this use refers to. When the loop finishes without
    ! `exit`, the index has run one past the end, which is how "not found" is
    ! detected below. (The configuration was validated against the same
    ! definitions, so this is a safety net.)
    definition_index = 0
    do definition_index = 1, size(definitions)
      if (definitions(definition_index)%field_name == variable_use%field_name) exit
    end do
    if (definition_index > size(definitions)) then
      call cable_abort("Output configuration refers to unknown variable " // trim(variable_use%field_name), __FILE__, __LINE__)
    end if

    ! Combine the definition (what the variable is) with this use's choices
    ! (how it is written). The constructor copies the aggregator, so every use
    ! has its own running totals.
    output_variable = cable_output_variable_t( &
      definition=definitions(definition_index), &
      aggregation=trim(variable_use%aggregation), &
      reduction=trim(variable_use%reduction), &
      netcdf_name=trim(variable_use%netcdf_name), &
      metadata=variable_use%metadata &
    )
  end function variable_from_use

  function stream_from_config(stream, grid_type, stream_variables, file_name) result(output_stream)
    !! An output file with the settings of a configured stream.
    type(cable_output_stream_config_t), intent(in) :: stream
    character(len=*), intent(in) :: grid_type
    type(cable_output_variable_t), intent(in) :: stream_variables(:)
    character(len=*), intent(in) :: file_name
    type(cable_output_stream_t) :: output_stream

    ! Copy across the settings that belong to the file. The NetCDF file itself
    ! is not created here; that happens in start_stream.
    output_stream%sampling_frequency = stream%frequency
    output_stream%grid_type = grid_type
    output_stream%file_name = file_name
    output_stream%shuffle = stream%shuffle
    output_stream%compression_level = stream%compression_level
    output_stream%output_variables = stream_variables
    ! Global attributes are optional.
    if (allocated(stream%metadata)) output_stream%metadata = stream%metadata
  end function stream_from_config

  module subroutine cable_output_impl_init_streams_from_config(config, grid_type, dels)
    !* Creates the output files described by the configuration, defines their
    ! contents and writes their coordinate variables.
    type(cable_output_config_t), intent(in) :: config
    character(len=*), intent(in) :: grid_type
    real, intent(in) :: dels !! The current time step size in seconds.
    integer :: stream_index

    if (.not. allocated(registered_definitions)) then
      call cable_abort("Output variables must be registered before initialising output streams.", __FILE__, __LINE__)
    end if

    ! First work out what the files will contain (no files yet), then create
    ! and define each of them.
    output_streams = cable_output_impl_build_streams(config, registered_definitions, grid_type)
    do stream_index = 1, size(output_streams)
      call start_stream(output_streams(stream_index))
    end do

  end subroutine cable_output_impl_init_streams_from_config

  subroutine start_stream(output_stream)
    !! Creates an output file, defines its contents, and writes its coordinates.
    type(cable_output_stream_t), intent(inout) :: output_stream
    integer :: iotype, i

    ! Compressed files must be netCDF-4.
    iotype = CABLE_NETCDF_IOTYPE_CLASSIC
    if (output_stream%compression_level > 0) iotype = CABLE_NETCDF_IOTYPE_NETCDF4C
    output_stream%output_file = cable_netcdf_create_file(output_stream%file_name, iotype=iotype)
    ! Coordinates (latitude, longitude, ...) go into every file, and depend
    ! on the grid type.
    output_stream%coordinate_variables = coordinate_variables_list(output_stream%grid_type)

    ! NetCDF has two phases. In "define" mode the dimensions, variables and
    ! attributes are declared; `end_def` then switches to "data" mode, after
    ! which values can be written but the structure can no longer change.
    call cable_output_define_stream(output_stream)
    call output_stream%output_file%end_def()

    ! Coordinates never change, so each is sampled once with the "instant"
    ! behaviour and written straight away, then the aggregator is cleared.
    do i = 1, size(output_stream%coordinate_variables)
      associate(coordinate_variable => output_stream%coordinate_variables(i))
        call coordinate_variable%aggregator%init(method="instant")
        call coordinate_variable%aggregator%accumulate()
        call cable_output_write_variable(output_stream, coordinate_variable)
        call coordinate_variable%aggregator%reset()
      end associate
    end do

    ! Prepare each variable's aggregator for accumulating. `init` allocates the
    ! storage and selects the behaviour (mean, sum, max, min or instant). A
    ! mean is set up as a plain sum; the division by the sample count happens at
    ! write time (see cable_output_impl_write).
    do i = 1, size(output_stream%output_variables)
      associate(output_variable => output_stream%output_variables(i))
        if (normalised_averaging .and. output_variable%aggregation_method == "mean") then
          call output_variable%aggregator%init(method="sum")
        else
          call output_variable%aggregator%init(method=output_variable%aggregation_method)
        end if
      end associate
    end do
  end subroutine start_stream

end submodule cable_output_impl_smod
