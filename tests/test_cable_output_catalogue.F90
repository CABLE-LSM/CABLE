! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

module test_cable_output_catalogue
  !! Tests for the output catalogue and the definitions built from it.
  use fortuno_interface_mod, only: test_list_t
  use fortuno_interface_mod, only: test_case
  use fortuno_interface_mod, only: check
  use cable_yaml_mod, only: cable_yaml_node_t
  use cable_yaml_mod, only: cable_yaml_parse_text
  use cable_output_catalogue_text_mod, only: catalogue_text => cable_output_catalogue_text_mod_get
  use cable_output_catalogue_mod, only: cable_output_catalogue_definitions_from_text
  use cable_output_mod, only: cable_output_binding_t
  use cable_output_mod, only: cable_output_variable_definition_t
  use cable_netcdf_mod, only: CABLE_NETCDF_FLOAT
  use cable_netcdf_mod, only: CABLE_NETCDF_INT
  use aggregator_mod, only: new_aggregator
  implicit none
  private

  public :: cable_output_catalogue_test_list

  character(len=32), parameter :: allowed_dimensions(7) = [character(len=32) :: &
    "patch", "soil", "rad", "snow", "plant_carbon_pools", "soil_carbon_pools", "land_global" &
  ]
  character(len=32), parameter :: allowed_groups(10) = [character(len=32) :: &
    "met", "flux", "radiation", "carbon", "soil", "snow", "veg", "params", "balances", "casa" &
  ]

  character(len=*), parameter :: small_catalogue_text = &
    "variables:" // new_line('a') // &
    "  - name: ""Rainf""" // new_line('a') // &
    "    dimensions: [patch]" // new_line('a') // &
    "    type: float" // new_line('a') // &
    "    group: met" // new_line('a') // &
    "    metadata:" // new_line('a') // &
    "      units: ""kg/m^2/s""" // new_line('a') // &
    "      long_name: ""Rainfall+snowfall""" // new_line('a') // &
    "  - name: ""iveg""" // new_line('a') // &
    "    dimensions: [patch]" // new_line('a') // &
    "    type: int" // new_line('a') // &
    "    parameter: true" // new_line('a') // &
    "    restart_name: ""iveg""" // new_line('a') // &
    "  - name: ""mvtype""" // new_line('a') // &
    "    output: false" // new_line('a') // &
    "    dimensions: []" // new_line('a') // &
    "    type: int" // new_line('a') // &
    "    distributed: false" // new_line('a') // &
    "    restart_name: ""mvtype""" // new_line('a')

contains

  function cable_output_catalogue_test_list() result(test_list)
    type(test_list_t) :: test_list

    test_list = test_list_t([ &
      test_case("test_catalogue_entries_are_well_formed", test_catalogue_entries_are_well_formed), &
      test_case("test_catalogue_names_are_unique", test_catalogue_names_are_unique), &
      test_case("test_catalogue_entry_counts", test_catalogue_entry_counts), &
      test_case("test_definitions_take_metadata_from_catalogue", test_definitions_take_metadata_from_catalogue), &
      test_case("test_definitions_take_data_from_bindings", test_definitions_take_data_from_bindings), &
      test_case("test_definitions_flags_and_restart_names", test_definitions_flags_and_restart_names) &
    ])
  end function cable_output_catalogue_test_list

  function contains_name(names, name) result(found)
    character(len=*), intent(in) :: names(:)
    character(len=*), intent(in) :: name
    logical :: found
    found = any(names == name)
  end function contains_name

  subroutine test_catalogue_entries_are_well_formed()
    type(cable_yaml_node_t) :: root, entries, entry, dimension_names, dimension_node
    integer :: entry_index, dimension_index

    root = cable_yaml_parse_text(catalogue_text())
    entries = root%get("variables")
    call check(entries%is_sequence(), msg="variables should be a sequence")
    do entry_index = 1, entries%size()
      entry = entries%at(entry_index)
      call check(entry%has("name"), msg="entry needs a name")
      call check(entry%has("dimensions"), msg="entry needs dimensions: " // entry%get_string("name"))
      call check(entry%has("type"), msg="entry needs a type: " // entry%get_string("name"))
      call check(any([character(len=5) :: "float", "int"] == entry%get_string("type")), &
        msg="type should be float or int: " // entry%get_string("name"))
      dimension_names = entry%get("dimensions")
      do dimension_index = 1, dimension_names%size()
        dimension_node = dimension_names%at(dimension_index)
        call check(contains_name(allowed_dimensions, dimension_node%as_string()), &
          msg="unknown dimension in " // entry%get_string("name"))
      end do
      if (entry%has("group")) then
        call check(contains_name(allowed_groups, entry%get_string("group")), &
          msg="unknown group in " // entry%get_string("name"))
        call check(entry%has("module"), msg="entry with a group needs a module: " // entry%get_string("name"))
      end if
      if (entry%has("output")) then
        call check(.not. entry%get_logical("output"), msg="output is only ever set to false")
        call check(entry%has("restart_name"), msg="restart-only entry needs restart_name: " // entry%get_string("name"))
      end if
    end do
  end subroutine test_catalogue_entries_are_well_formed

  subroutine test_catalogue_names_are_unique()
    type(cable_yaml_node_t) :: root, entries, entry, other
    integer :: first, second

    root = cable_yaml_parse_text(catalogue_text())
    entries = root%get("variables")
    do first = 1, entries%size()
      entry = entries%at(first)
      do second = first + 1, entries%size()
        other = entries%at(second)
        call check(entry%get_string("name") /= other%get_string("name"), &
          msg="duplicate catalogue name " // entry%get_string("name"))
        if (entry%has("restart_name") .and. other%has("restart_name")) then
          call check(entry%get_string("restart_name") /= other%get_string("restart_name"), &
            msg="duplicate restart name " // entry%get_string("restart_name"))
        end if
      end do
    end do
  end subroutine test_catalogue_names_are_unique

  subroutine test_catalogue_entry_counts()
    !! Guards against entries being lost when the catalogue is edited: these are
    !! the numbers of variables defined by the legacy diagnostics modules.
    type(cable_yaml_node_t) :: root, entries, entry
    integer :: entry_index, restart_count, restart_only_count

    root = cable_yaml_parse_text(catalogue_text())
    entries = root%get("variables")
    restart_count = 0
    restart_only_count = 0
    do entry_index = 1, entries%size()
      entry = entries%at(entry_index)
      if (entry%has("restart_name")) restart_count = restart_count + 1
      if (entry%has("output")) restart_only_count = restart_only_count + 1
    end do
    call check(entries%size() == 160, msg="expected 160 catalogue entries")
    call check(restart_count == 42, msg="expected 42 restart variables")
    call check(restart_only_count == 26, msg="expected 26 restart-only variables")
  end subroutine test_catalogue_entry_counts

  function small_bindings(rainfall, vegetation_type, vegetation_types) result(bindings)
    real, intent(inout) :: rainfall(:)
    integer, intent(inout) :: vegetation_type(:)
    integer, intent(inout) :: vegetation_types
    type(cable_output_binding_t), allocatable :: bindings(:)

    bindings = [ &
      cable_output_binding_t(name="Rainf", aggregator=new_aggregator(rainfall), divide_by=2.0, &
        offset_by=1.0, range=[0.0, 100.0], available=.false.), &
      cable_output_binding_t(name="iveg", aggregator=new_aggregator(vegetation_type)), &
      cable_output_binding_t(name="mvtype", aggregator=new_aggregator(vegetation_types)) &
    ]
  end function small_bindings

  subroutine test_definitions_take_metadata_from_catalogue()
    real :: rainfall(3)
    integer :: vegetation_type(3), vegetation_types
    type(cable_output_variable_definition_t), allocatable :: definitions(:)

    definitions = cable_output_catalogue_definitions_from_text( &
      small_catalogue_text, small_bindings(rainfall, vegetation_type, vegetation_types))
    call check(size(definitions) == 3, msg="expected three definitions")
    call check(trim(definitions(1)%field_name) == "Rainf", msg="field name")
    call check(definitions(1)%var_type == CABLE_NETCDF_FLOAT, msg="float type")
    call check(definitions(2)%var_type == CABLE_NETCDF_INT, msg="int type")
    call check(size(definitions(1)%data_shape) == 1, msg="one dimension")
    call check(size(definitions(3)%data_shape) == 0, msg="scalar has no dimensions")
    call check(trim(definitions(1)%group_name) == "met", msg="group name")
    call check(size(definitions(1)%metadata) == 2, msg="two attributes")
    call check(trim(definitions(1)%metadata(1)%name) == "units", msg="first attribute name")
    call check(trim(definitions(1)%metadata(1)%value) == "kg/m^2/s", msg="first attribute value")
    call check(trim(definitions(1)%metadata(2)%name) == "long_name", msg="second attribute name")
    call check(.not. allocated(definitions(2)%metadata), msg="no metadata when the catalogue has none")
  end subroutine test_definitions_take_metadata_from_catalogue

  subroutine test_definitions_take_data_from_bindings()
    real :: rainfall(3)
    integer :: vegetation_type(3), vegetation_types
    type(cable_output_variable_definition_t), allocatable :: definitions(:)

    definitions = cable_output_catalogue_definitions_from_text( &
      small_catalogue_text, small_bindings(rainfall, vegetation_type, vegetation_types))
    call check(definitions(1)%divide_by == 2.0, msg="divide_by")
    call check(definitions(1)%offset_by == 1.0, msg="offset_by")
    call check(definitions(1)%scale_by == 1.0, msg="scale_by defaults to 1")
    ! Range converts to working-variable units: (range - offset) * divide_by / scale_by.
    call check(definitions(1)%range(1) == -2.0, msg="lower range in native units")
    call check(definitions(1)%range(2) == 198.0, msg="upper range in native units")
    call check(.not. definitions(1)%available, msg="availability comes from the binding")
    call check(definitions(2)%available, msg="available by default")
    call check(definitions(2)%range(1) == -huge(0.0), msg="unbounded range by default")
    call check(allocated(definitions(1)%aggregator), msg="aggregator is copied from the binding")
  end subroutine test_definitions_take_data_from_bindings

  subroutine test_definitions_flags_and_restart_names()
    real :: rainfall(3)
    integer :: vegetation_type(3), vegetation_types
    type(cable_output_variable_definition_t), allocatable :: definitions(:)

    definitions = cable_output_catalogue_definitions_from_text( &
      small_catalogue_text, small_bindings(rainfall, vegetation_type, vegetation_types))
    call check(definitions(1)%outputtable, msg="outputtable by default")
    call check(len_trim(definitions(1)%restart_name) == 0, msg="no restart name unless stated")
    call check(definitions(2)%parameter, msg="parameter flag")
    call check(definitions(2)%distributed, msg="distributed by default")
    call check(trim(definitions(2)%restart_name) == "iveg", msg="restart name")
    call check(.not. definitions(3)%outputtable, msg="restart-only entry is not outputtable")
    call check(.not. definitions(3)%distributed, msg="distributed flag can be cleared")
    call check(trim(definitions(1)%native_frequency) == "timestep", msg="native frequency default")
  end subroutine test_definitions_flags_and_restart_names

end module test_cable_output_catalogue
