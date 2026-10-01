! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

module cable_output_catalogue_mod
  !! Builds output variable definitions from the output catalogue.
  !!
  !! The catalogue (`output_catalogue.yaml`, embedded into the executable at
  !! build time) says what each output variable is. The bindings (see
  !! `cable_output_bindings_mod`) say where its data comes from. This module
  !! joins the two by name.
  use cable_error_handler_mod, only: cable_abort
  use cable_netcdf_mod, only: CABLE_NETCDF_INT
  use cable_netcdf_mod, only: CABLE_NETCDF_FLOAT
  use cable_yaml_mod, only: cable_yaml_node_t
  use cable_yaml_mod, only: cable_yaml_parse_text
  use cable_output_catalogue_text_mod, only: catalogue_text => cable_output_catalogue_text_mod_get
  use cable_output_mod, only: cable_output_variable_definition_t
  use cable_output_mod, only: cable_output_binding_t
  use cable_output_mod, only: cable_output_attribute_t
  use cable_output_mod, only: cable_output_get_dimension
  implicit none
  private

  public :: cable_output_catalogue_definitions
  public :: cable_output_catalogue_definitions_from_text

contains

  function cable_output_catalogue_definitions(bindings) result(definitions)
    !! Definitions of every variable in the built-in output catalogue.
    type(cable_output_binding_t), intent(in) :: bindings(:)
      !! Bindings from all model components. Must match the catalogue one to one.
    type(cable_output_variable_definition_t), allocatable :: definitions(:)

    ! The catalogue text was embedded into the executable at build time, so
    ! there is no file to locate at run time.
    definitions = cable_output_catalogue_definitions_from_text(catalogue_text(), bindings)
  end function cable_output_catalogue_definitions

  function cable_output_catalogue_definitions_from_text(text, bindings) result(definitions)
    !! Definitions of every variable in the catalogue given as YAML `text`.
    !!
    !! Aborts, listing every problem found, if a catalogue entry has no binding
    !! or a binding has no catalogue entry.
    character(len=*), intent(in) :: text
    type(cable_output_binding_t), intent(in) :: bindings(:)
    type(cable_output_variable_definition_t), allocatable :: definitions(:)
    type(cable_yaml_node_t) :: root, entries, entry
    integer :: entry_index, binding_index
    character(:), allocatable :: problems
    logical :: binding_used(size(bindings))

    root = cable_yaml_parse_text(text)
    entries = root%get("variables")
    ! One definition per catalogue entry, in catalogue order.
    allocate(definitions(entries%size()))
    ! `binding_used` remembers which bindings found a partner, so leftovers
    ! can be reported after the loop.
    binding_used = .false.
    problems = ""

    ! Direction 1: every catalogue entry needs a binding with the same name.
    do entry_index = 1, entries%size()
      entry = entries%at(entry_index)
      binding_index = find_binding(bindings, entry%get_string("name"))
      ! Keep going after a miss so that all missing names are listed together.
      if (binding_index == 0) then
        problems = problems // new_line('a') // "  catalogue entry has no binding: " // entry%get_string("name")
        cycle
      end if
      binding_used(binding_index) = .true.
      definitions(entry_index) = build_definition(entry, bindings(binding_index))
    end do

    ! Direction 2: every binding must have been claimed by a catalogue entry.
    ! Together the two directions make the catalogue and bindings match one to one.
    do binding_index = 1, size(bindings)
      if (.not. binding_used(binding_index)) then
        problems = problems // new_line('a') // "  binding has no catalogue entry: " // trim(bindings(binding_index)%name)
      end if
    end do

    ! A mismatch is a mistake in CABLE's own files rather than in the user's
    ! configuration, so it stops the run.
    if (len(problems) > 0) then
      call cable_abort("Output catalogue and bindings do not match:" // problems, __FILE__, __LINE__)
    end if
  end function cable_output_catalogue_definitions_from_text

  pure function find_binding(bindings, name) result(binding_index)
    !! Position of the binding called `name`, or 0 if there is none.
    type(cable_output_binding_t), intent(in) :: bindings(:)
    character(len=*), intent(in) :: name
    integer :: binding_index
    integer :: candidate

    ! Linear search: there are a few hundred bindings and this runs once.
    binding_index = 0
    do candidate = 1, size(bindings)
      if (trim(bindings(candidate)%name) == name) then
        binding_index = candidate
        return
      end if
    end do
  end function find_binding

  function build_definition(entry, binding) result(definition)
    !! Combine one catalogue entry with its binding.
    type(cable_yaml_node_t), intent(in) :: entry
    type(cable_output_binding_t), intent(in) :: binding
    type(cable_output_variable_definition_t) :: definition
    type(cable_yaml_node_t) :: dimension_names, dimension_node, metadata, attribute_node
    integer :: index

    ! Part 1: descriptive settings from the catalogue entry. Every key except
    ! name and type is optional, and a missing key leaves the default set in the
    ! definition type (output true, not a parameter, distributed, updated every
    ! time step, no group or module, no restart name).
    definition%field_name = entry%get_string("name")
    if (entry%has("output")) definition%outputtable = entry%get_logical("output")
    if (entry%has("restart_name")) definition%restart_name = entry%get_string("restart_name")
    if (entry%has("parameter")) definition%parameter = entry%get_logical("parameter")
    if (entry%has("distributed")) definition%distributed = entry%get_logical("distributed")
    if (entry%has("native_frequency")) definition%native_frequency = entry%get_string("native_frequency")
    if (entry%has("group")) definition%group_name = entry%get_string("group")
    if (entry%has("module")) definition%module_name = entry%get_string("module")
    definition%var_type = netcdf_type(entry%get_string("type"), definition%field_name)

    ! Dimensions are stored as names (`patch`, `soil`). Turn each into the
    ! output module's dimension object, which also carries its current size.
    dimension_names = entry%get("dimensions")
    allocate(definition%data_shape(dimension_names%size()))
    do index = 1, dimension_names%size()
      dimension_node = dimension_names%at(index)
      definition%data_shape(index) = cable_output_get_dimension(dimension_node%as_string())
    end do

    ! Variable attributes (units, long_name, ...) are copied in file order.
    if (entry%has("metadata")) then
      metadata = entry%get("metadata")
      allocate(definition%metadata(metadata%size()))
      do index = 1, metadata%size()
        attribute_node = metadata%at(index)
        definition%metadata(index)%name = attribute_node%key
        definition%metadata(index)%value = attribute_node%as_string()
      end do
    end if

    ! Part 2: what only Fortran can supply, taken from the binding. A binding
    ! for a variable the model cannot provide has no aggregator, because there
    ! is no working variable to point at.
    if (allocated(binding%aggregator)) allocate(definition%aggregator, source=binding%aggregator)
    definition%available = binding%available
    definition%scale_by = binding%scale_by
    definition%divide_by = binding%divide_by
    definition%offset_by = binding%offset_by
    if (allocated(binding%range)) then
      ! The binding states the range in output units (what the user sees). The
      ! range check compares against the model's own working variable, so undo
      ! the conversion the aggregator applies on output:
      !   output = scale * native / divide + offset
      ! which rearranges to
      !   native = (output - offset) * divide / scale
      definition%range = (binding%range - definition%offset_by) * definition%divide_by / definition%scale_by
    end if
  end function build_definition

  function netcdf_type(type_name, field_name) result(type_code)
    !! NetCDF type code for the catalogue type name `float` or `int`.
    character(len=*), intent(in) :: type_name
    character(len=*), intent(in) :: field_name
    integer :: type_code

    ! Only two types exist in the catalogue. Anything else is a typo in the
    ! catalogue, reported with the name of the entry that has it.
    select case (type_name)
    case ("float")
      type_code = CABLE_NETCDF_FLOAT
    case ("int")
      type_code = CABLE_NETCDF_INT
    case default
      call cable_abort("Unknown type '" // type_name // "' for catalogue entry " // field_name, __FILE__, __LINE__)
    end select
  end function netcdf_type

end module cable_output_catalogue_mod
