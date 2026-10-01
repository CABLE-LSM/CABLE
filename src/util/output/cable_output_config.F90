! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

module cable_output_config_mod
  !! Reads and validates the output configuration file.
  !!
  !! The configuration is a YAML file with two central ideas: streams, each
  !! linked to one output file and one write frequency, and output variables,
  !! each directed to one stream. Groups and modules are shorthand for sets of
  !! variables. This module turns the file into a flat list of stream settings
  !! and a flat list of variable uses (one per variable, aggregation and stream),
  !! with all groups and modules expanded and all `{...}` templates filled in.
  !!
  !! Everything here is a function of its arguments. Problems are collected
  !! rather than stopping at the first, so that one run reports every mistake in
  !! the file.
  use cable_error_handler_mod, only: cable_abort
  use cable_yaml_mod, only: cable_yaml_node_t
  use cable_yaml_mod, only: cable_yaml_load_file
  use cable_yaml_mod, only: cable_yaml_parse_text
  use cable_output_mod, only: cable_output_variable_definition_t
  use cable_output_mod, only: cable_output_attribute_t
  use cable_output_mod, only: cable_output_stream_config_t
  use cable_output_mod, only: cable_output_variable_config_t
  use cable_output_mod, only: cable_output_config_t
  use cable_output_mod, only: CABLE_OUTPUT_CONFIG_LENGTH
  use cable_output_mod, only: allowed_aggregation_methods
  use cable_output_mod, only: allowed_reduction_methods
  use cable_timing_mod, only: cable_timing_frequencies
  use cable_timing_mod, only: cable_output_frequency_rank => cable_timing_frequency_rank
  use cable_netcdf_mod, only: CABLE_NETCDF_INT
  implicit none
  private

  public :: cable_output_config_load
  public :: cable_output_config_parse
  public :: cable_output_config_build
  public :: cable_output_substitute
  public :: cable_output_frequency_rank
  public :: cable_output_stream_config_t
  public :: cable_output_variable_config_t
  public :: cable_output_config_t
  public :: CABLE_OUTPUT_CONFIG_LENGTH

  !> Keys that may appear in the file, in a stream, and in a variable, group or module entry. These are named
  !! constants rather than array constructors written in the calls, because gfortran can pass the wrong string
  !! length for such actual arguments.
  character(len=32), parameter :: top_level_keys(4) = [character(len=32) :: &
    "streams", "variables", "groups", "modules" &
  ]
  character(len=32), parameter :: stream_keys(7) = [character(len=32) :: &
    "file_name", "frequency", "netcdf_name", "shuffle", "compression_level", "separate_file_per_variable", "metadata" &
  ]
  character(len=32), parameter :: specification_keys(6) = [character(len=32) :: &
    "name", "stream", "aggregation", "reduction", "netcdf_name", "metadata" &
  ]
  !> Substitutions available in the metadata of a stream.
  character(len=32), parameter :: stream_substitution_names(3) = [character(len=32) :: &
    "frequency", "start_date", "end_date" &
  ]

  !> Attributes that are always written from the variable definition and so cannot be set in the file.
  character(len=16), parameter :: reserved_attributes(4) = [character(len=16) :: &
    "standard_name", "long_name", "units", "cell_methods" &
  ]

contains

  ! ------------------------------------------------------------------ entry points

  function cable_output_config_load(file_name, definitions, start_date, end_date) result(config)
    !! Read the configuration file `file_name`. Aborts, listing every problem, if it is invalid.
    character(len=*), intent(in) :: file_name
    type(cable_output_variable_definition_t), intent(in) :: definitions(:)
    character(len=*), intent(in) :: start_date !! First simulation date, `YYYY-MM-DD`.
    character(len=*), intent(in) :: end_date !! Last simulation date, `YYYY-MM-DD`.
    type(cable_output_config_t) :: config
    type(cable_yaml_node_t) :: root
    character(:), allocatable :: problems

    ! Reading the file and building the configuration are separate steps: the
    ! build collects every problem it finds, and only here, at the boundary
    ! with the rest of the model, is a non-empty list turned into a stop.
    root = cable_yaml_load_file(file_name)
    config = cable_output_config_build(root, definitions, start_date, end_date, problems)
    if (len(problems) > 0) then
      call cable_abort("Invalid output configuration file " // trim(file_name) // ":" // problems, __FILE__, __LINE__)
    end if
  end function cable_output_config_load

  function cable_output_config_parse(text, definitions, start_date, end_date) result(config)
    !! As `cable_output_config_load`, for a configuration held in a string.
    character(len=*), intent(in) :: text
    type(cable_output_variable_definition_t), intent(in) :: definitions(:)
    character(len=*), intent(in) :: start_date
    character(len=*), intent(in) :: end_date
    type(cable_output_config_t) :: config
    type(cable_yaml_node_t) :: root
    character(:), allocatable :: problems

    ! Identical to the file version once the text has been parsed.
    root = cable_yaml_parse_text(text)
    config = cable_output_config_build(root, definitions, start_date, end_date, problems)
    if (len(problems) > 0) then
      call cable_abort("Invalid output configuration:" // problems, __FILE__, __LINE__)
    end if
  end function cable_output_config_parse

  function cable_output_config_build(root, definitions, start_date, end_date, problems) result(config)
    !! Build the configuration from a parsed YAML document.
    !!
    !! `problems` is empty if the configuration is valid; otherwise it holds one
    !! line per problem, each starting with a newline.
    type(cable_yaml_node_t), intent(in) :: root
    type(cable_output_variable_definition_t), intent(in) :: definitions(:)
    character(len=*), intent(in) :: start_date
    character(len=*), intent(in) :: end_date
    character(:), allocatable, intent(out) :: problems
    type(cable_output_config_t) :: config
    type(cable_output_variable_config_t), allocatable :: uses(:)
    type(cable_yaml_node_t) :: streams_node
    integer :: use_count

    ! `problems` starts empty and only ever grows: every check below appends a
    ! line and carries on, so one run reports every mistake in the file.
    problems = ""
    allocate(config%streams(0))
    ! `uses` is a growable buffer. It starts with room for 16 variable uses and
    ! doubles when full (see add_use); `use_count` is how many are filled.
    allocate(uses(16))
    use_count = 0

    ! Without a mapping at the top there is nothing to look inside, so this is
    ! the one problem that cannot be combined with others.
    if (.not. root%is_mapping()) then
      problems = problems // new_line('a') // "  the file must contain a mapping with streams, variables, groups and modules"
      allocate(config%variables(0))
      return
    end if
    ! A misspelt section name (for example `variable:`) would otherwise be
    ! silently ignored, leaving the user wondering why nothing is written.
    call check_keys(root, top_level_keys, "the file", problems)

    ! Streams are read first because every variable entry refers to one by number.
    if (root%has("streams")) then
      streams_node = root%get("streams")
      call read_streams(streams_node, start_date, end_date, config%streams, problems)
    else
      problems = problems // new_line('a') // "  the file defines no streams"
    end if

    ! Modules, then groups, then variables: this fixes the order of variables in the output files.
    call add_specifications(root, "modules", config%streams, definitions, start_date, end_date, uses, use_count, problems)
    call add_specifications(root, "groups", config%streams, definitions, start_date, end_date, uses, use_count, problems)
    call add_specifications(root, "variables", config%streams, definitions, start_date, end_date, uses, use_count, problems)

    ! Trim the buffer down to the uses actually filled.
    allocate(config%variables(use_count))
    if (use_count > 0) config%variables = uses(1:use_count)

    ! These two checks look across entries, so they run once everything has been
    ! expanded: the two rules cannot be seen from a single entry.
    call check_streams_unique(config%streams, problems)
    call check_rules(config, problems)
  end function cable_output_config_build

  ! ------------------------------------------------------------------ streams

  subroutine read_streams(streams_node, start_date, end_date, streams, problems)
    !! Read every stream. A subroutine because gfortran drops updates made to the
    !! deferred-length string `problems` inside a function called on the right of an assignment.
    type(cable_yaml_node_t), intent(in) :: streams_node
    character(len=*), intent(in) :: start_date, end_date
    type(cable_output_stream_config_t), allocatable, intent(out) :: streams(:)
    character(:), allocatable, intent(inout) :: problems
    type(cable_yaml_node_t) :: stream_node
    integer :: stream_index

    if (.not. streams_node%is_mapping()) then
      problems = problems // new_line('a') // "  streams must be a mapping from stream number to settings"
      allocate(streams(0))
      return
    end if
    ! One stream per entry of the mapping, in file order.
    allocate(streams(streams_node%size()))
    do stream_index = 1, streams_node%size()
      stream_node = streams_node%at(stream_index)
      call read_stream(stream_node, start_date, end_date, streams(stream_index), problems)
    end do
  end subroutine read_streams

  subroutine read_stream(stream_node, start_date, end_date, stream, problems)
    !! A subroutine rather than a function so that `problems`, a deferred-length
    !! string, is reliably updated for the caller by gfortran.
    type(cable_yaml_node_t), intent(in) :: stream_node
    character(len=*), intent(in) :: start_date, end_date
    type(cable_output_stream_config_t), intent(inout) :: stream
    character(:), allocatable, intent(inout) :: problems
    character(:), allocatable :: context, text
    character(len=32) :: substitution_values(3)
    integer :: read_status

    ! The key of each stream entry is its number. Keys are stored as text, so
    ! convert it here; a non-number (`two:`) is reported rather than guessed.
    read(stream_node%key, *, iostat=read_status) stream%stream_id
    ! `context` names the stream in every message that follows, so a problem
    ! can be found without counting lines in the file.
    context = "stream " // stream_node%key
    if (read_status /= 0) then
      problems = problems // new_line('a') // "  stream numbers must be integers, found '" // stream_node%key // "'"
      return
    end if
    if (.not. stream_node%is_mapping()) then
      problems = problems // new_line('a') // "  " // context // " must be a mapping"
      return
    end if
    call check_keys(stream_node, stream_keys, context, problems)

    ! file_name and frequency are the two required settings. require_scalar
    ! records the problem itself when one is missing or malformed, so the
    ! assignment is simply skipped in that case.
    if (require_scalar(stream_node, "file_name", context, problems)) stream%file_name = stream_node%get_string("file_name")
    if (require_scalar(stream_node, "frequency", context, problems)) then
      ! Only the frequencies the timing code understands are accepted; there
      ! are no aliases, so a typo cannot silently become something else.
      text = stream_node%get_string("frequency")
      if (any(cable_timing_frequencies == text)) then
        stream%frequency = text
      else
        problems = problems // new_line('a') // "  " // context // ": unknown frequency '" // text // "'"
      end if
    end if
    ! Optional settings keep the defaults set in the type unless the file
    ! gives a value.
    if (stream_node%has("netcdf_name")) stream%netcdf_name = stream_node%get_string("netcdf_name")
    if (stream_node%has("shuffle")) call read_logical(stream_node, "shuffle", context, stream%shuffle, problems)
    if (stream_node%has("separate_file_per_variable")) then
      call read_logical(stream_node, "separate_file_per_variable", context, stream%separate_file_per_variable, problems)
    end if
    if (stream_node%has("compression_level")) then
      call read_integer(stream_node, "compression_level", context, stream%compression_level, problems)
      ! NetCDF's deflate levels run 1 (fastest) to 9 (smallest); 0 means none.
      ! An out-of-range value is reset to 0 so later code never sees it.
      if (stream%compression_level < 0 .or. stream%compression_level > 9) then
        problems = problems // new_line('a') // "  " // context // ": compression_level must be from 0 to 9"
        stream%compression_level = 0
      end if
    end if
    ! Stream attributes can mention the frequency and the run dates, so those
    ! are the only substitutions offered here. The values sit in a named
    ! variable because gfortran can mis-pass an inline array constructor.
    substitution_values = [character(len=32) :: stream%frequency, start_date, end_date]
    call read_metadata(stream_node, context, stream_substitution_names, substitution_values, stream%metadata, problems)
  end subroutine read_stream

  subroutine check_streams_unique(streams, problems)
    !! Stream numbers must be unique, and no two streams may write to the same file.
    type(cable_output_stream_config_t), intent(in) :: streams(:)
    character(:), allocatable, intent(inout) :: problems
    character(len=16) :: number_text
    integer :: first, second

    ! Compare every pair of streams once (second always comes after first).
    do first = 1, size(streams)
      do second = first + 1, size(streams)
        if (streams(first)%stream_id == streams(second)%stream_id) then
          write(number_text, '(i0)') streams(first)%stream_id
          problems = problems // new_line('a') // "  stream " // trim(number_text) // " is defined twice"
        end if
        ! Two streams writing the same file would overwrite each other. Streams
        ! that write one file per variable are exempt: their file_name is only
        ! used for its directory, so sharing it is harmless.
        if (len_trim(streams(first)%file_name) > 0 .and. streams(first)%file_name == streams(second)%file_name &
            .and. .not. (streams(first)%separate_file_per_variable .or. streams(second)%separate_file_per_variable)) then
          problems = problems // new_line('a') // "  streams " // stream_label(streams(first)) // " and " // &
            stream_label(streams(second)) // " both write to file " // trim(streams(first)%file_name)
        end if
      end do
    end do
  end subroutine check_streams_unique

  ! ------------------------------------------------------------------ variables, groups and modules

  subroutine add_specifications(root, section, streams, definitions, start_date, end_date, uses, use_count, problems)
    !! Expand every entry of the `variables`, `groups` or `modules` section into variable uses.
    type(cable_yaml_node_t), intent(in) :: root
    character(len=*), intent(in) :: section
    type(cable_output_stream_config_t), intent(in) :: streams(:)
    type(cable_output_variable_definition_t), intent(in) :: definitions(:)
    character(len=*), intent(in) :: start_date, end_date
    type(cable_output_variable_config_t), allocatable, intent(inout) :: uses(:)
    integer, intent(inout) :: use_count
    character(:), allocatable, intent(inout) :: problems
    type(cable_yaml_node_t) :: section_node, entry
    integer :: entry_index

    ! A section that is absent is fine: a file may use only variables, or only modules.
    if (.not. root%has(section)) return
    section_node = root%get(section)
    ! Each section must be a list of entries. An empty value (`groups:` with
    ! nothing under it) is accepted as an empty list.
    if (.not. section_node%is_sequence()) then
      if (section_node%size() == 0) return
      problems = problems // new_line('a') // "  " // section // " must be a list"
      return
    end if
    do entry_index = 1, section_node%size()
      entry = section_node%at(entry_index)
      call add_specification(entry, section, entry_index, streams, definitions, start_date, end_date, &
        uses, use_count, problems)
    end do
  end subroutine add_specifications

  subroutine add_specification(entry, section, entry_index, streams, definitions, start_date, end_date, &
      uses, use_count, problems)
    type(cable_yaml_node_t), intent(in) :: entry
    character(len=*), intent(in) :: section
    integer, intent(in) :: entry_index
    type(cable_output_stream_config_t), intent(in) :: streams(:)
    type(cable_output_variable_definition_t), intent(in) :: definitions(:)
    character(len=*), intent(in) :: start_date, end_date
    type(cable_output_variable_config_t), allocatable, intent(inout) :: uses(:)
    integer, intent(inout) :: use_count
    character(:), allocatable, intent(inout) :: problems
    character(len=16) :: index_text
    character(:), allocatable :: context, name, aggregation, reduction, netcdf_name_template
    integer :: stream_id, stream_index, definition_index, member_count
    logical :: named_directly

    ! Until the entry's name is known, messages refer to it by position.
    write(index_text, '(i0)') entry_index
    context = section // " entry " // trim(index_text)
    if (.not. entry%is_mapping()) then
      problems = problems // new_line('a') // "  " // context // " must be a mapping"
      return
    end if
    call check_keys(entry, specification_keys, context, problems)
    ! name, stream and aggregation must all be present. Without them there is
    ! nothing more to check for this entry, so stop here; the `.and.` does not
    ! short-circuit, so every missing key is still reported.
    if (.not. (require_scalar(entry, "name", context, problems) .and. require_scalar(entry, "stream", context, problems) &
        .and. require_scalar(entry, "aggregation", context, problems))) return
    ! From here on messages can use the entry's name, which is more useful.
    name = entry%get_string("name")
    context = section // " entry '" // name // "'"

    ! The stream number must match one of the streams defined in the file.
    call read_integer(entry, "stream", context, stream_id, problems)
    stream_index = find_stream(streams, stream_id)
    if (stream_index == 0) then
      problems = problems // new_line('a') // "  " // context // ": stream " // entry%get_string("stream") // " is not defined"
      return
    end if

    ! Validate the two method names against the lists the engine supports.
    ! Returning early is safe: an unknown method makes the rest of the entry
    ! meaningless.
    aggregation = entry%get_string("aggregation")
    if (.not. any(allowed_aggregation_methods == aggregation)) then
      problems = problems // new_line('a') // "  " // context // ": unknown aggregation '" // aggregation // "'"
      return
    end if
    reduction = "none"
    if (entry%has("reduction")) reduction = entry%get_string("reduction")
    if (.not. any(allowed_reduction_methods == reduction)) then
      problems = problems // new_line('a') // "  " // context // ": unknown reduction '" // reduction // "'"
      return
    end if
    ! The NetCDF name template comes from the entry if it has one, otherwise
    ! from its stream (which defaults to "{field_name}").
    netcdf_name_template = trim(streams(stream_index)%netcdf_name)
    if (entry%has("netcdf_name")) netcdf_name_template = entry%get_string("netcdf_name")

    ! Expansion. Walk every registered definition and pick those the entry
    ! selects: by name for a variable, by group or module name otherwise. A
    ! group or module therefore expands to however many definitions belong to it.
    named_directly = section == "variables"
    member_count = 0
    do definition_index = 1, size(definitions)
      if (.not. is_member(definitions(definition_index), section, name)) cycle
      ! Restart-only entries are never output.
      if (.not. definitions(definition_index)%outputtable) cycle
      if (.not. definitions(definition_index)%available) then
        ! Groups and modules skip variables the model cannot provide; naming one directly is a mistake.
        if (named_directly) then
          problems = problems // new_line('a') // "  " // context // ": variable is not available in this model configuration"
        end if
        cycle
      end if
      ! A selected variable becomes one use in the entry's stream.
      member_count = member_count + 1
      call add_use(entry, context, streams(stream_index), definitions(definition_index), aggregation, reduction, &
        netcdf_name_template, start_date, end_date, uses, use_count, problems)
    end do

    ! Nothing selected is fine if the name exists but every member was
    ! unavailable (the group is just empty here). It is an error only if no
    ! definition has that name at all, which is almost always a typo.
    if (member_count == 0 .and. .not. any_definition_matches(definitions, section, name)) then
      problems = problems // new_line('a') // "  " // context // ": unknown " // singular(section) // " '" // name // "'"
    end if
  end subroutine add_specification

  subroutine add_use(entry, context, stream, definition, aggregation, reduction, netcdf_name_template, &
      start_date, end_date, uses, use_count, problems)
    !! Check one variable against one stream and append it to `uses`.
    type(cable_yaml_node_t), intent(in) :: entry
    character(len=*), intent(in) :: context
    type(cable_output_stream_config_t), intent(in) :: stream
    type(cable_output_variable_definition_t), intent(in) :: definition
    character(len=*), intent(in) :: aggregation, reduction, netcdf_name_template, start_date, end_date
    type(cable_output_variable_config_t), allocatable, intent(inout) :: uses(:)
    integer, intent(inout) :: use_count
    character(:), allocatable, intent(inout) :: problems
    type(cable_output_variable_config_t) :: variable_use
    type(cable_output_variable_config_t), allocatable :: grown(:)
    character(:), allocatable :: substitution_error, label
    character(len=16), parameter :: names(6) = [character(len=16) :: &
      "field_name", "frequency", "aggregation", "reduction", "start_date", "end_date"]
    character(len=64) :: values(6)

    ! `label` adds the variable's name to the entry's context, so a problem in
    ! a group points at the member that caused it.
    label = context // " (" // trim(definition%field_name) // ")"
    variable_use%field_name = definition%field_name
    variable_use%stream_id = stream%stream_id
    variable_use%aggregation = aggregation
    variable_use%reduction = reduction

    call check_use_against_definition(label, definition, stream, reduction, problems)

    ! Fill in the NetCDF name template. The values line up with `names`, so
    ! each {name} in the template is replaced by the value at the same position.
    values = [character(len=64) :: trim(definition%field_name), trim(stream%frequency), aggregation, reduction, &
      start_date, end_date]
    variable_use%netcdf_name = cable_output_substitute(netcdf_name_template, names, values, substitution_error)
    if (allocated(substitution_error)) then
      problems = problems // new_line('a') // "  " // label // ": netcdf_name: " // substitution_error
    end if
    ! Only check the finished name if the template itself was valid,
    ! otherwise the name is half-substituted and the check would add noise.
    ! NetCDF variable names cannot contain '/'.
    if (allocated(substitution_error) .eqv. .false.) then
      if (len_trim(variable_use%netcdf_name) == 0 .or. index(variable_use%netcdf_name, "/") > 0) then
        problems = problems // new_line('a') // "  " // label // ": netcdf_name '" // trim(variable_use%netcdf_name) // &
          "' must not be empty or contain '/'"
      end if
    end if
    call read_metadata(entry, label, names, values, variable_use%metadata, problems)

    ! Append to the growable buffer. When it is full, make one twice as big,
    ! copy the filled part across and swap it in, so appending stays cheap
    ! however many variables a group expands to.
    if (use_count == size(uses)) then
      allocate(grown(2 * size(uses)))
      grown(1:use_count) = uses(1:use_count)
      call move_alloc(grown, uses)
    end if
    use_count = use_count + 1
    uses(use_count) = variable_use
  end subroutine add_use

  subroutine check_use_against_definition(label, definition, stream, reduction, problems)
    !! Reductions and frequencies must make sense for the variable being written.
    character(len=*), intent(in) :: label
    type(cable_output_variable_definition_t), intent(in) :: definition
    type(cable_output_stream_config_t), intent(in) :: stream
    character(len=*), intent(in) :: reduction
    character(:), allocatable, intent(inout) :: problems
    character(len=64) :: first_dimension

    ! Reductions collapse the per-tile dimension to one value per grid cell, so
    ! they only make sense for a variable that has one.
    if (reduction /= "none") then
      if (.not. allocated(definition%data_shape)) then
        problems = problems // new_line('a') // "  " // label // ": a scalar cannot be reduced"
      else if (size(definition%data_shape) == 0) then
        problems = problems // new_line('a') // "  " // label // ": a scalar cannot be reduced"
      else
        ! The tile dimension is always the first one, and its name starts with
        ! "patch" (the code's older word for tile).
        first_dimension = definition%data_shape(1)%name()
        if (first_dimension(1:5) /= "patch") then
          problems = problems // new_line('a') // "  " // label // ": only variables defined per tile can be reduced"
        end if
      end if
      ! Averaging whole numbers such as vegetation type has no meaning.
      if (reduction == "grid_cell_average" .and. definition%var_type == CABLE_NETCDF_INT) then
        problems = problems // new_line('a') // "  " // label // ": grid_cell_average is not defined for integer variables"
      end if
    end if
    ! A file cannot be written more often than the model updates a variable.
    ! Frequencies are ranked from finest to coarsest, so a stream whose rank is
    ! lower than the variable's native rank is too fine. Parameters are
    ! written once whatever the frequency, so they are exempt.
    if (.not. definition%parameter .and. len_trim(stream%frequency) > 0) then
      if (cable_output_frequency_rank(stream%frequency) < cable_output_frequency_rank(definition%native_frequency)) then
        problems = problems // new_line('a') // "  " // label // ": the variable is only updated " // &
          trim(definition%native_frequency) // ", so it cannot be written every " // trim(stream%frequency)
      end if
    end if
  end subroutine check_use_against_definition

  ! ------------------------------------------------------------------ the two rules

  subroutine check_rules(config, problems)
    !! Rule 1: no two variables in a stream share a NetCDF name.
    !! Rule 2: no two variables in a stream are numerically identical.
    type(cable_output_config_t), intent(in) :: config
    character(:), allocatable, intent(inout) :: problems
    integer :: first, second

    ! Compare every pair of uses once. Only uses in the same stream can clash,
    ! because each stream is a separate file.
    do first = 1, size(config%variables)
      do second = first + 1, size(config%variables)
        if (config%variables(first)%stream_id /= config%variables(second)%stream_id) cycle
        associate (one => config%variables(first), other => config%variables(second))
          ! Rule 2 is tested first: the same field with the same aggregation
          ! and reduction is a duplicate whatever the NetCDF names are. The
          ! `else if` means a duplicate is reported once, as rule 2, instead
          ! of also as a name clash.
          if (one%field_name == other%field_name .and. one%aggregation == other%aggregation &
              .and. one%reduction == other%reduction) then
            problems = problems // new_line('a') // "  " // trim(one%field_name) // " is written twice to stream " // &
              stream_id_text(one%stream_id) // " with the same aggregation (" // trim(one%aggregation) // &
              ") and reduction (" // trim(one%reduction) // "), even if the NetCDF names differ"
          ! Rule 1: different variables (or the same one aggregated two ways)
          ! that ended up with the same NetCDF name.
          else if (one%netcdf_name == other%netcdf_name) then
            problems = problems // new_line('a') // "  stream " // stream_id_text(one%stream_id) // &
              " would contain two variables named '" // trim(one%netcdf_name) // "' (from " // &
              trim(one%field_name) // " and " // trim(other%field_name) // "); set netcdf_name to tell them apart"
          end if
        end associate
      end do
    end do
  end subroutine check_rules

  ! ------------------------------------------------------------------ substitution

  function cable_output_substitute(template, names, values, error) result(text)
    !! Replace each `{name}` in `template` with the matching entry of `values`.
    !!
    !! Unknown names and unmatched braces leave the text unchanged and are
    !! described in `error`, which is not allocated if the template is valid.
    character(len=*), intent(in) :: template
    character(len=*), intent(in) :: names(:)
    character(len=*), intent(in) :: values(:)
    character(:), allocatable, intent(out) :: error
    character(:), allocatable :: text
    integer :: position, closing, name_index

    ! Scan the template one character at a time, copying text through until a
    ! "{". Then find the matching "}" and replace everything between them by
    ! the value of that name. `position` is where scanning has reached.
    text = ""
    position = 1
    do while (position <= len(template))
      if (template(position:position) == "{") then
        ! `index` searched from `position`, so convert the result back to a
        ! position in the whole template.
        closing = index(template(position:), "}")
        ! No closing brace: keep the rest of the text as it is and stop.
        if (closing == 0) then
          call add_error("unmatched '{' in template '" // template // "'")
          text = text // template(position:)
          exit
        end if
        closing = position + closing - 1
        ! The name is the text between the braces. An unknown name is kept
        ! in the output (braces and all) and recorded as an error, so the caller
        ! can show both what was written and what was wrong.
        name_index = find_name(names, template(position + 1:closing - 1))
        if (name_index == 0) then
          call add_error("unknown substitution '{" // template(position + 1:closing - 1) // "}' in template '" // template // "'")
          text = text // template(position:closing)
        else
          text = text // trim(values(name_index))
        end if
        ! Carry on after the closing brace.
        position = closing + 1
      else
        ! Ordinary character: copy it across.
        text = text // template(position:position)
        position = position + 1
      end if
    end do

  contains

    subroutine add_error(message)
      character(len=*), intent(in) :: message
      ! Several problems in one template are joined with "; " into one message.
      if (allocated(error)) then
        error = error // "; " // message
      else
        error = message
      end if
    end subroutine add_error

  end function cable_output_substitute

  pure function find_name(names, name) result(name_index)
    character(len=*), intent(in) :: names(:)
    character(len=*), intent(in) :: name
    integer :: name_index
    integer :: candidate

    ! Exact match against each known name, ignoring only trailing padding.
    name_index = 0
    do candidate = 1, size(names)
      if (trim(names(candidate)) == name) then
        name_index = candidate
        return
      end if
    end do
  end function find_name

  ! ------------------------------------------------------------------ helpers

  subroutine read_metadata(node, context, names, values, attributes, problems)
    !! The `metadata` mapping of `node` as attributes, with templates filled in.
    type(cable_yaml_node_t), intent(in) :: node
    character(len=*), intent(in) :: context
    character(len=*), intent(in) :: names(:)
    character(len=*), intent(in) :: values(:)
    type(cable_output_attribute_t), allocatable, intent(out) :: attributes(:)
    character(:), allocatable, intent(inout) :: problems
    type(cable_yaml_node_t) :: metadata, attribute_node
    character(:), allocatable :: substitution_error, value
    integer :: attribute_index

    ! Metadata is optional; no key means no extra attributes.
    if (.not. node%has("metadata")) then
      allocate(attributes(0))
      return
    end if
    metadata = node%get("metadata")
    ! `metadata:` with nothing under it reads as an empty value, which is
    ! treated as no attributes instead of an error.
    if (metadata%is_scalar()) then
      if (len(metadata%as_string()) == 0) then
        allocate(attributes(0))
        return
      end if
    end if
    if (.not. metadata%is_mapping()) then
      problems = problems // new_line('a') // "  " // context // ": metadata must be a mapping"
      allocate(attributes(0))
      return
    end if
    ! One attribute per entry, in file order. Each is checked, then its value
    ! has any {name} templates filled in.
    allocate(attributes(metadata%size()))
    do attribute_index = 1, metadata%size()
      attribute_node = metadata%at(attribute_index)
      attributes(attribute_index)%name = attribute_node%key
      ! units, long_name, standard_name and cell_methods always come from
      ! the variable definition, so letting the file set them would make two
      ! sources disagree.
      if (any(reserved_attributes == attribute_node%key)) then
        problems = problems // new_line('a') // "  " // context // ": metadata cannot set '" // attribute_node%key // &
          "'; it is always taken from the variable definition"
      end if
      ! An attribute is a single piece of text; a nested mapping cannot be
      ! written to a NetCDF attribute. `cycle` skips substitution for it.
      if (.not. attribute_node%is_scalar()) then
        problems = problems // new_line('a') // "  " // context // ": metadata '" // attribute_node%key // "' must be a single value"
        cycle
      end if
      value = cable_output_substitute(attribute_node%as_string(), names, values, substitution_error)
      if (allocated(substitution_error)) then
        problems = problems // new_line('a') // "  " // context // ": metadata '" // attribute_node%key // "': " // &
          substitution_error
      end if
      attributes(attribute_index)%value = value
    end do
  end subroutine read_metadata

  subroutine check_keys(node, allowed, context, problems)
    !! Report every key of the mapping `node` that is not in `allowed`.
    type(cable_yaml_node_t), intent(in) :: node
    character(len=*), intent(in) :: allowed(:)
    character(len=*), intent(in) :: context
    character(:), allocatable, intent(inout) :: problems
    type(cable_yaml_node_t) :: child
    integer :: child_index

    ! Every key present must be one of the allowed names; this catches typos
    ! that would otherwise be ignored without a word.
    do child_index = 1, node%size()
      child = node%at(child_index)
      if (.not. any(allowed == child%key)) then
        problems = problems // new_line('a') // "  " // context // ": unknown key '" // child%key // "'"
      end if
    end do
  end subroutine check_keys

  function require_scalar(node, key, context, problems) result(present_and_valid)
    !! True if `node` has `key` holding a single value; otherwise records a problem.
    type(cable_yaml_node_t), intent(in) :: node
    character(len=*), intent(in) :: key
    character(len=*), intent(in) :: context
    character(:), allocatable, intent(inout) :: problems
    logical :: present_and_valid
    type(cable_yaml_node_t) :: child

    ! Start from "not valid" so every early return reports failure.
    present_and_valid = .false.
    if (.not. node%has(key)) then
      problems = problems // new_line('a') // "  " // context // ": missing required key '" // key // "'"
      return
    end if
    child = node%get(key)
    if (.not. child%is_scalar()) then
      problems = problems // new_line('a') // "  " // context // ": '" // key // "' must be a single value"
      return
    end if
    present_and_valid = .true.
  end function require_scalar

  subroutine read_integer(node, key, context, value, problems)
    type(cable_yaml_node_t), intent(in) :: node
    character(len=*), intent(in) :: key
    character(len=*), intent(in) :: context
    integer, intent(out) :: value
    character(:), allocatable, intent(inout) :: problems
    character(:), allocatable :: text
    integer :: read_status

    ! Take the text of the value and try to read a whole number from it. A
    ! non-zero status means it was not one; report it and fall back to 0 so the
    ! caller never receives an undefined number.
    value = 0
    text = node%get_string(key)
    read(text, *, iostat=read_status) value
    if (read_status /= 0) then
      problems = problems // new_line('a') // "  " // context // ": '" // key // "' must be an integer, found '" // text // "'"
      value = 0
    end if
  end subroutine read_integer

  subroutine read_logical(node, key, context, value, problems)
    type(cable_yaml_node_t), intent(in) :: node
    character(len=*), intent(in) :: key
    character(len=*), intent(in) :: context
    logical, intent(out) :: value
    character(:), allocatable, intent(inout) :: problems
    character(:), allocatable :: text

    ! Same three spellings as cable_yaml_mod accepts. Anything else is reported,
    ! and the value stays false.
    value = .false.
    text = node%get_string(key)
    select case (text)
    case ("true", "True", "TRUE")
      value = .true.
    case ("false", "False", "FALSE")
      value = .false.
    case default
      problems = problems // new_line('a') // "  " // context // ": '" // key // "' must be true or false, found '" // text // "'"
    end select
  end subroutine read_logical

  pure function find_stream(streams, stream_id) result(stream_index)
    type(cable_output_stream_config_t), intent(in) :: streams(:)
    integer, intent(in) :: stream_id
    integer :: stream_index
    integer :: candidate

    ! Streams are identified by their number, which need not be 1, 2, 3...
    ! so look the number up rather than using it as a position.
    stream_index = 0
    do candidate = 1, size(streams)
      if (streams(candidate)%stream_id == stream_id) then
        stream_index = candidate
        return
      end if
    end do
  end function find_stream

  pure function is_member(definition, section, name) result(member)
    !! Whether `definition` is selected by the entry `name` of the given section.
    type(cable_output_variable_definition_t), intent(in) :: definition
    character(len=*), intent(in) :: section
    character(len=*), intent(in) :: name
    logical :: member

    ! What an entry's `name` refers to depends on the section it is in: a
    ! variable's own name, or the group or module the variable belongs to.
    select case (section)
    case ("variables")
      member = trim(definition%field_name) == name
    case ("groups")
      member = trim(definition%group_name) == name
    case ("modules")
      member = trim(definition%module_name) == name
    case default
      member = .false.
    end select
  end function is_member

  pure function any_definition_matches(definitions, section, name) result(found)
    type(cable_output_variable_definition_t), intent(in) :: definitions(:)
    character(len=*), intent(in) :: section
    character(len=*), intent(in) :: name
    logical :: found
    integer :: candidate

    ! Used to tell "this group exists but has nothing available" apart from
    ! "no such group": only the second is a mistake in the file.
    found = .false.
    do candidate = 1, size(definitions)
      if (definitions(candidate)%outputtable .and. is_member(definitions(candidate), section, name)) then
        found = .true.
        return
      end if
    end do
  end function any_definition_matches

  pure function singular(section) result(word)
    character(len=*), intent(in) :: section
    character(:), allocatable :: word
    ! Drop the trailing "s" (groups -> group, modules -> module) for messages.
    ! The explicit "variables" line below gives the same answer as that rule; it
    ! is redundant and only spells out the most common case.
    word = section(1:len(section) - 1)
    if (section == "variables") word = "variable"
  end function singular

  function stream_id_text(stream_id) result(text)
    integer, intent(in) :: stream_id
    character(:), allocatable :: text
    character(len=16) :: buffer

    write(buffer, '(i0)') stream_id
    text = trim(buffer)
  end function stream_id_text

  function stream_label(stream) result(text)
    type(cable_output_stream_config_t), intent(in) :: stream
    character(:), allocatable :: text
    text = stream_id_text(stream%stream_id)
  end function stream_label

end module cable_output_config_mod
