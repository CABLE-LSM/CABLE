! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

module test_cable_output_config
  !! Tests for cable_output_config_mod, using the examples from the output
  !! configuration design.
  use fortuno_interface_mod, only: test_list_t
  use fortuno_interface_mod, only: test_case
  use fortuno_interface_mod, only: check
  use cable_yaml_mod, only: cable_yaml_node_t
  use cable_yaml_mod, only: cable_yaml_parse_text
  use cable_output_mod, only: cable_output_variable_definition_t
  use cable_output_mod, only: cable_output_get_dimension
  use cable_netcdf_mod, only: CABLE_NETCDF_FLOAT
  use cable_netcdf_mod, only: CABLE_NETCDF_INT
  use cable_output_config_mod
  implicit none
  private

  public :: cable_output_config_test_list

  character(len=*), parameter :: nl = new_line('a')

contains

  function cable_output_config_test_list() result(test_list)
    type(test_list_t) :: test_list

    test_list = test_list_t([ &
      test_case("test_config_minimum_specification", test_minimum_specification), &
      test_case("test_config_maximum_specification", test_maximum_specification), &
      test_case("test_config_group_and_module_expansion", test_group_and_module_expansion), &
      test_case("test_config_multiple_streams", test_multiple_streams), &
      test_case("test_config_rule1_duplicate_names", test_rule1_duplicate_names), &
      test_case("test_config_rule1_resolved_by_netcdf_name", test_rule1_resolved_by_netcdf_name), &
      test_case("test_config_rule2_identical_variables", test_rule2_identical_variables), &
      test_case("test_config_diurnal_example_is_valid", test_diurnal_example_is_valid), &
      test_case("test_config_unavailable_in_group_is_skipped", test_unavailable_in_group_is_skipped), &
      test_case("test_config_unavailable_by_name_is_an_error", test_unavailable_by_name_is_an_error), &
      test_case("test_config_reports_every_problem", test_reports_every_problem), &
      test_case("test_config_reduction_checks", test_reduction_checks), &
      test_case("test_config_frequency_checks", test_frequency_checks), &
      test_case("test_config_separate_files_use_dates", test_separate_files_use_dates), &
      test_case("test_config_stream_checks", test_stream_checks), &
      test_case("test_config_substitution", test_substitution), &
      test_case("test_config_value_and_metadata_problems", test_value_and_metadata_problems) &
    ])
  end function cable_output_config_test_list

  ! ------------------------------------------------------------------ helpers

  function make_definition(name, group, module_name, dimension_name, is_integer, is_parameter, &
      native_frequency, available, outputtable) result(definition)
    character(len=*), intent(in) :: name
    character(len=*), intent(in), optional :: group, module_name, dimension_name, native_frequency
    logical, intent(in), optional :: is_integer, is_parameter, available, outputtable
    type(cable_output_variable_definition_t) :: definition

    definition%field_name = name
    definition%var_type = CABLE_NETCDF_FLOAT
    if (present(group)) definition%group_name = group
    if (present(module_name)) definition%module_name = module_name
    if (present(is_integer)) then
      if (is_integer) definition%var_type = CABLE_NETCDF_INT
    end if
    if (present(is_parameter)) definition%parameter = is_parameter
    if (present(native_frequency)) definition%native_frequency = native_frequency
    if (present(available)) definition%available = available
    if (present(outputtable)) definition%outputtable = outputtable
    if (present(dimension_name)) then
      allocate(definition%data_shape(1))
      definition%data_shape(1) = cable_output_get_dimension(dimension_name)
    else
      allocate(definition%data_shape(0))
    end if
  end function make_definition

  function test_definitions() result(definitions)
    type(cable_output_variable_definition_t), allocatable :: definitions(:)

    definitions = [ &
      make_definition("GPP", "carbon_cycle", "biogeochemistry", "patch"), &
      make_definition("labile_carbon", "carbon_pools", "biogeochemistry", "patch"), &
      make_definition("wood_carbon", "carbon_pools", "biogeochemistry", "patch"), &
      make_definition("evaporation", "water_cycle", "biogeophysics", "patch"), &
      make_definition("transpiration", "water_cycle", "biogeophysics", "patch"), &
      make_definition("wtd", "water_cycle", "biogeophysics", "patch", available=.false.), &
      make_definition("Qle", "energy_cycle", "biogeophysics", "patch"), &
      make_definition("Qh", "energy_cycle", "biogeophysics", "patch"), &
      make_definition("iveg", "params", "parameters", "patch", is_integer=.true., is_parameter=.true.), &
      make_definition("Tmx", "veg", "biogeophysics", "patch", native_frequency="daily"), &
      make_definition("Area", "params", "parameters", "land_global", is_parameter=.true.), &
      make_definition("nap", outputtable=.false.) &
    ]
  end function test_definitions

  function build(text, problems) result(config)
    character(len=*), intent(in) :: text
    character(:), allocatable, intent(out) :: problems
    type(cable_output_config_t) :: config
    type(cable_yaml_node_t) :: root

    root = cable_yaml_parse_text(text)
    config = cable_output_config_build(root, test_definitions(), "2000-01-01", "2000-12-31", problems)
  end function build

  function mentions(problems, text) result(found)
    character(len=*), intent(in) :: problems, text
    logical :: found
    found = index(problems, text) > 0
  end function mentions

  function count_lines(problems) result(count)
    character(len=*), intent(in) :: problems
    integer :: count, position
    count = 0
    do position = 1, len(problems)
      if (problems(position:position) == nl) count = count + 1
    end do
  end function count_lines

  function names_of(config) result(text)
    !! "name,name,..." of the config's variables, in order.
    type(cable_output_config_t), intent(in) :: config
    character(:), allocatable :: text
    integer :: i

    text = ""
    do i = 1, size(config%variables)
      if (i > 1) text = text // ","
      text = text // trim(config%variables(i)%netcdf_name)
    end do
  end function names_of

  ! ------------------------------------------------------------------ tests

  subroutine test_minimum_specification()
    type(cable_output_config_t) :: config
    character(:), allocatable :: problems

    config = build( &
      "streams:" // nl // "    1:" // nl // "        file_name: cable_output.nc" // nl // "        frequency: daily" // nl // &
      "variables:" // nl // "    - name: GPP" // nl // "      stream: 1" // nl // "      aggregation: mean" // nl, problems)
    call check(len(problems) == 0, msg="minimum specification should be valid: " // problems)
    call check(size(config%streams) == 1, msg="one stream")
    call check(trim(config%streams(1)%file_name) == "cable_output.nc", msg="file name")
    call check(trim(config%streams(1)%frequency) == "daily", msg="frequency")
    ! Defaults from the design.
    call check(trim(config%streams(1)%netcdf_name) == "{field_name}", msg="default netcdf_name template")
    call check(config%streams(1)%shuffle, msg="shuffle defaults to true")
    call check(config%streams(1)%compression_level == 1, msg="compression_level defaults to 1")
    call check(.not. config%streams(1)%separate_file_per_variable, msg="separate files default to false")
    call check(size(config%variables) == 1, msg="one variable")
    call check(trim(config%variables(1)%field_name) == "GPP", msg="variable name")
    call check(trim(config%variables(1)%netcdf_name) == "GPP", msg="default NetCDF name is the field name")
    call check(trim(config%variables(1)%reduction) == "none", msg="reduction defaults to none")
    call check(config%variables(1)%stream_id == 1, msg="variable stream")
  end subroutine test_minimum_specification

  subroutine test_maximum_specification()
    type(cable_output_config_t) :: config
    character(:), allocatable :: problems

    config = build( &
      "streams:" // nl // "    1:" // nl // "        file_name: cable_output.nc" // nl // "        frequency: monthly" // nl // &
      "        netcdf_name: ""{field_name}_{aggregation}""" // nl // "        shuffle: false" // nl // &
      "        compression_level: 4" // nl // "        separate_file_per_variable: false" // nl // &
      "        metadata:" // nl // "            model: CABLE" // nl // &
      "            experiment: example_experiment_id" // nl // "            period: ""{start_date} to {end_date}""" // nl // &
      "variables:" // nl // "    - name: GPP" // nl // "      stream: 1" // nl // "      aggregation: mean" // nl // &
      "      reduction: grid_cell_average" // nl // "      netcdf_name: GPP_mean" // nl // "      metadata:" // nl // &
      "        description: Computed using X algorithm." // nl // "        note: ""{aggregation} of {field_name}""" // nl, &
      problems)
    call check(len(problems) == 0, msg="maximum specification should be valid: " // problems)
    call check(.not. config%streams(1)%shuffle, msg="shuffle can be turned off")
    call check(config%streams(1)%compression_level == 4, msg="compression level")
    call check(size(config%streams(1)%metadata) == 3, msg="three global attributes")
    call check(trim(config%streams(1)%metadata(1)%name) == "model", msg="first global attribute name")
    call check(trim(config%streams(1)%metadata(1)%value) == "CABLE", msg="first global attribute value")
    call check(trim(config%streams(1)%metadata(3)%value) == "2000-01-01 to 2000-12-31", msg="dates are substituted")
    call check(trim(config%variables(1)%netcdf_name) == "GPP_mean", msg="variable netcdf_name")
    call check(trim(config%variables(1)%reduction) == "grid_cell_average", msg="reduction")
    call check(size(config%variables(1)%metadata) == 2, msg="two variable attributes")
    call check(trim(config%variables(1)%metadata(2)%value) == "mean of GPP", msg="variable attribute substitution")
  end subroutine test_maximum_specification

  subroutine test_group_and_module_expansion()
    type(cable_output_config_t) :: config
    character(:), allocatable :: problems

    config = build( &
      "streams:" // nl // "    1:" // nl // "        file_name: a.nc" // nl // "        frequency: daily" // nl // &
      "groups:" // nl // "    - name: carbon_pools" // nl // "      stream: 1" // nl // "      aggregation: mean" // nl // &
      "      netcdf_name: ""{field_name}_{aggregation}""" // nl // &
      "modules:" // nl // "    - name: parameters" // nl // "      stream: 1" // nl // "      aggregation: instant" // nl, problems)
    call check(len(problems) == 0, msg="expansion should be valid: " // problems)
    ! Modules are expanded before groups, each in definition order.
    call check(names_of(config) == "iveg,Area,labile_carbon_mean,wood_carbon_mean", msg="expanded names: " // names_of(config))
    call check(trim(config%variables(3)%aggregation) == "mean", msg="group aggregation applies to members")
    call check(trim(config%variables(1)%aggregation) == "instant", msg="module aggregation applies to members")
  end subroutine test_group_and_module_expansion

  subroutine test_multiple_streams()
    type(cable_output_config_t) :: config
    character(:), allocatable :: problems

    config = build( &
      "streams:" // nl // "    1:" // nl // "        file_name: biogeophysics.nc" // nl // "        frequency: monthly" // nl // &
      "    2:" // nl // "        file_name: biogeochemistry.nc" // nl // "        frequency: monthly" // nl // &
      "modules:" // nl // "    - name: biogeophysics" // nl // "      stream: 1" // nl // "      aggregation: mean" // nl // &
      "    - name: biogeochemistry" // nl // "      stream: 2" // nl // "      aggregation: mean" // nl, problems)
    call check(len(problems) == 0, msg="multiple streams should be valid: " // problems)
    call check(size(config%streams) == 2, msg="two streams")
    call check(config%streams(2)%stream_id == 2, msg="second stream number")
    ! Biogeophysics: evaporation, transpiration, Qle, Qh, Tmx (wtd is unavailable). Biogeochemistry: 3.
    call check(size(config%variables) == 8, msg="expected 8 variables")
    call check(all(config%variables(1:5)%stream_id == 1), msg="biogeophysics goes to stream 1")
    call check(all(config%variables(6:8)%stream_id == 2), msg="biogeochemistry goes to stream 2")
  end subroutine test_multiple_streams

  subroutine test_rule1_duplicate_names()
    !! A group and a variable in it, with different aggregation, both using the same NetCDF name.
    type(cable_output_config_t) :: config
    character(:), allocatable :: problems

    config = build( &
      "streams:" // nl // "    1:" // nl // "        file_name: a.nc" // nl // "        frequency: daily" // nl // &
      "groups:" // nl // "    - name: carbon_pools" // nl // "      stream: 1" // nl // "      aggregation: mean" // nl // &
      "variables:" // nl // "    - name: labile_carbon" // nl // "      stream: 1" // nl // "      aggregation: instant" // nl, &
      problems)
    call check(mentions(problems, "two variables named 'labile_carbon'"), msg="rule 1 should be reported: " // problems)
    call check(count_lines(problems) == 1, msg="only rule 1 is broken")
  end subroutine test_rule1_duplicate_names

  subroutine test_rule1_resolved_by_netcdf_name()
    type(cable_output_config_t) :: config
    character(:), allocatable :: problems

    config = build( &
      "streams:" // nl // "    1:" // nl // "        file_name: a.nc" // nl // "        frequency: daily" // nl // &
      "groups:" // nl // "    - name: carbon_pools" // nl // "      stream: 1" // nl // "      aggregation: mean" // nl // &
      "      netcdf_name: ""{field_name}_{aggregation}""" // nl // &
      "variables:" // nl // "    - name: labile_carbon" // nl // "      stream: 1" // nl // "      aggregation: instant" // nl, &
      problems)
    call check(len(problems) == 0, msg="distinct names resolve rule 1: " // problems)
    call check(names_of(config) == "labile_carbon_mean,wood_carbon_mean,labile_carbon", msg=names_of(config))
  end subroutine test_rule1_resolved_by_netcdf_name

  subroutine test_rule2_identical_variables()
    !! Same field, aggregation and reduction in one stream is a duplicate even with different names.
    type(cable_output_config_t) :: config
    character(:), allocatable :: problems

    config = build( &
      "streams:" // nl // "    1:" // nl // "        file_name: a.nc" // nl // "        frequency: daily" // nl // &
      "groups:" // nl // "    - name: carbon_pools" // nl // "      stream: 1" // nl // "      aggregation: mean" // nl // &
      "variables:" // nl // "    - name: labile_carbon" // nl // "      stream: 1" // nl // "      aggregation: mean" // nl // &
      "      netcdf_name: labile_carbon_mean" // nl, problems)
    call check(mentions(problems, "labile_carbon is written twice to stream 1"), msg="rule 2 should be reported: " // problems)
    call check(count_lines(problems) == 1, msg="only rule 2 is broken")
  end subroutine test_rule2_identical_variables

  subroutine test_diurnal_example_is_valid()
    type(cable_output_config_t) :: config
    character(:), allocatable :: problems

    config = build( &
      "streams:" // nl // "    1:" // nl // "        file_name: cable_water_diurnal.nc" // nl // "        frequency: 3hrly" // nl // &
      "groups:" // nl // "    - name: water_cycle" // nl // "      stream: 1" // nl // "      aggregation: sum" // nl // &
      "      netcdf_name: ""{field_name}_{aggregation}""" // nl // &
      "    - name: energy_cycle" // nl // "      stream: 1" // nl // "      aggregation: sum" // nl // &
      "      netcdf_name: ""{field_name}_{aggregation}""" // nl // &
      "variables:" // nl // "    - name: evaporation" // nl // "      stream: 1" // nl // "      aggregation: instant" // nl // &
      "      netcdf_name: evaporation_instant" // nl // "    - name: Qle" // nl // "      stream: 1" // nl // &
      "      aggregation: instant" // nl // "      netcdf_name: Qle_instant" // nl, problems)
    call check(len(problems) == 0, msg="diurnal example should be valid: " // problems)
    ! water_cycle (wtd unavailable): 2 variables, energy_cycle: 2, plus 2 instantaneous.
    call check(size(config%variables) == 6, msg="expected 6 variables")
    call check(trim(config%variables(1)%netcdf_name) == "evaporation_sum", msg="first name")
  end subroutine test_diurnal_example_is_valid

  subroutine test_unavailable_in_group_is_skipped()
    type(cable_output_config_t) :: config
    character(:), allocatable :: problems
    integer :: i

    config = build( &
      "streams:" // nl // "    1:" // nl // "        file_name: a.nc" // nl // "        frequency: daily" // nl // &
      "    2:" // nl // "        file_name: b.nc" // nl // "        frequency: daily" // nl // &
      "groups:" // nl // "    - name: water_cycle" // nl // "      stream: 1" // nl // "      aggregation: mean" // nl // &
      "modules:" // nl // "    - name: biogeophysics" // nl // "      stream: 2" // nl // "      aggregation: max" // nl, problems)
    call check(len(problems) == 0, msg="an unavailable member of a group is not an error: " // problems)
    do i = 1, size(config%variables)
      call check(trim(config%variables(i)%field_name) /= "wtd", msg="unavailable variable must be excluded")
    end do
    call check(size(config%variables) == 2 + 5, msg="water_cycle gives 2, biogeophysics gives 5")
  end subroutine test_unavailable_in_group_is_skipped

  subroutine test_unavailable_by_name_is_an_error()
    type(cable_output_config_t) :: config
    character(:), allocatable :: problems

    config = build( &
      "streams:" // nl // "    1:" // nl // "        file_name: a.nc" // nl // "        frequency: daily" // nl // &
      "variables:" // nl // "    - name: wtd" // nl // "      stream: 1" // nl // "      aggregation: mean" // nl, problems)
    call check(mentions(problems, "not available in this model configuration"), msg="expected an error: " // problems)
    call check(mentions(problems, "wtd"), msg="error should name the variable: " // problems)
  end subroutine test_unavailable_by_name_is_an_error

  subroutine test_reports_every_problem()
    type(cable_output_config_t) :: config
    character(:), allocatable :: problems

    config = build( &
      "streams:" // nl // "    1:" // nl // "        file_name: a.nc" // nl // "        frequency: hourly" // nl // &
      "        colour: red" // nl // &
      "variables:" // nl // "    - name: no_such_variable" // nl // "      stream: 1" // nl // "      aggregation: mean" // nl // &
      "    - name: GPP" // nl // "      stream: 7" // nl // "      aggregation: mean" // nl // &
      "    - name: GPP" // nl // "      stream: 1" // nl // "      aggregation: median" // nl // &
      "    - name: GPP" // nl // "      stream: 1" // nl // &
      "groups:" // nl // "    - name: no_such_group" // nl // "      stream: 1" // nl // "      aggregation: mean" // nl, problems)
    call check(mentions(problems, "unknown frequency 'hourly'"), msg="frequency: " // problems)
    call check(mentions(problems, "unknown key 'colour'"), msg="unknown key: " // problems)
    call check(mentions(problems, "unknown variable 'no_such_variable'"), msg="variable: " // problems)
    call check(mentions(problems, "stream 7 is not defined"), msg="stream: " // problems)
    call check(mentions(problems, "unknown aggregation 'median'"), msg="aggregation: " // problems)
    call check(mentions(problems, "missing required key 'aggregation'"), msg="missing key: " // problems)
    call check(mentions(problems, "unknown group 'no_such_group'"), msg="group: " // problems)
    call check(count_lines(problems) == 7, msg="all seven problems reported together")
  end subroutine test_reports_every_problem

  subroutine test_reduction_checks()
    type(cable_output_config_t) :: config
    character(:), allocatable :: problems
    character(len=*), parameter :: streams = &
      "streams:" // nl // "    1:" // nl // "        file_name: a.nc" // nl // "        frequency: daily" // nl

    config = build(streams // "variables:" // nl // "    - name: GPP" // nl // "      stream: 1" // nl // &
      "      aggregation: mean" // nl // "      reduction: dominant_tile" // nl, problems)
    call check(len(problems) == 0, msg="dominant_tile is valid for a per-tile float: " // problems)

    config = build(streams // "variables:" // nl // "    - name: iveg" // nl // "      stream: 1" // nl // &
      "      aggregation: instant" // nl // "      reduction: grid_cell_average" // nl, problems)
    call check(mentions(problems, "not defined for integer"), msg="integer average: " // problems)

    config = build(streams // "variables:" // nl // "    - name: Area" // nl // "      stream: 1" // nl // &
      "      aggregation: instant" // nl // "      reduction: first_tile_on_cell" // nl, problems)
    call check(mentions(problems, "only variables defined per tile"), msg="non-tile variable: " // problems)

    config = build(streams // "variables:" // nl // "    - name: GPP" // nl // "      stream: 1" // nl // &
      "      aggregation: mean" // nl // "      reduction: first_patch_in_grid_cell" // nl, problems)
    call check(mentions(problems, "unknown reduction 'first_patch_in_grid_cell'"), msg="old name is rejected: " // problems)
  end subroutine test_reduction_checks

  subroutine test_frequency_checks()
    type(cable_output_config_t) :: config
    character(:), allocatable :: problems

    ! Tmx is only updated daily, so it cannot be written every 3 hours...
    config = build( &
      "streams:" // nl // "    1:" // nl // "        file_name: a.nc" // nl // "        frequency: 3hrly" // nl // &
      "variables:" // nl // "    - name: Tmx" // nl // "      stream: 1" // nl // "      aggregation: mean" // nl, problems)
    call check(mentions(problems, "only updated daily"), msg="too fine a frequency: " // problems)

    ! ...but a monthly mean of it is fine.
    config = build( &
      "streams:" // nl // "    1:" // nl // "        file_name: a.nc" // nl // "        frequency: monthly" // nl // &
      "variables:" // nl // "    - name: Tmx" // nl // "      stream: 1" // nl // "      aggregation: mean" // nl, problems)
    call check(len(problems) == 0, msg="monthly mean of a daily variable: " // problems)

    ! Parameters are written once, whatever the stream frequency.
    config = build( &
      "streams:" // nl // "    1:" // nl // "        file_name: a.nc" // nl // "        frequency: timestep" // nl // &
      "variables:" // nl // "    - name: iveg" // nl // "      stream: 1" // nl // "      aggregation: instant" // nl, problems)
    call check(len(problems) == 0, msg="parameters ignore frequency: " // problems)
    call check(cable_output_frequency_rank("timestep") < cable_output_frequency_rank("yearly"), msg="rank order")
    call check(cable_output_frequency_rank("fortnightly") == 0, msg="unknown frequency has rank 0")
  end subroutine test_frequency_checks

  subroutine test_separate_files_use_dates()
    type(cable_output_config_t) :: config
    character(:), allocatable :: problems

    config = build( &
      "streams:" // nl // "    1:" // nl // "        file_name: unused_file_name.nc" // nl // "        frequency: daily" // nl // &
      "        separate_file_per_variable: true" // nl // &
      "        netcdf_name: ""{field_name}_{aggregation}_{start_date}-{end_date}""" // nl // &
      "groups:" // nl // "    - name: water_cycle" // nl // "      stream: 1" // nl // "      aggregation: mean" // nl, problems)
    call check(len(problems) == 0, msg="separate files should be valid: " // problems)
    call check(config%streams(1)%separate_file_per_variable, msg="separate_file_per_variable")
    call check(trim(config%variables(1)%netcdf_name) == "evaporation_mean_2000-01-01-2000-12-31", &
      msg="file name from dates: " // trim(config%variables(1)%netcdf_name))
  end subroutine test_separate_files_use_dates

  subroutine test_stream_checks()
    type(cable_output_config_t) :: config
    character(:), allocatable :: problems

    config = build( &
      "streams:" // nl // "    1:" // nl // "        file_name: same.nc" // nl // "        frequency: daily" // nl // &
      "    2:" // nl // "        file_name: same.nc" // nl // "        frequency: monthly" // nl // &
      "    3:" // nl // "        file_name: c.nc" // nl // "        frequency: daily" // nl // "        compression_level: 12" // nl, &
      problems)
    call check(mentions(problems, "both write to file same.nc"), msg="file clash: " // problems)
    call check(mentions(problems, "compression_level must be from 0 to 9"), msg="compression range: " // problems)

    config = build("variables:" // nl // "    - name: GPP" // nl // "      stream: 1" // nl // "      aggregation: mean" // nl, &
      problems)
    call check(mentions(problems, "defines no streams"), msg="no streams: " // problems)
  end subroutine test_stream_checks

  subroutine test_value_and_metadata_problems()
    !! Every kind of problem in a value or in metadata must reach the caller.
    type(cable_output_config_t) :: config
    character(:), allocatable :: problems

    config = build( &
      "streams:" // nl // "    1:" // nl // "        file_name: a.nc" // nl // "        frequency: daily" // nl // &
      "        shuffle: maybe" // nl // "        separate_file_per_variable: sometimes" // nl // &
      "        compression_level: high" // nl // "        metadata:" // nl // "            title: ""{colour}""" // nl // &
      "    two:" // nl // "        file_name: b.nc" // nl // "        frequency: daily" // nl // &
      "variables:" // nl // "    - name: GPP" // nl // "      stream: first" // nl // "      aggregation: mean" // nl // &
      "    - name: GPP" // nl // "      stream: 1" // nl // "      aggregation: mean" // nl // &
      "      netcdf_name: ""{oops}""" // nl // "      metadata:" // nl // "        units: mm" // nl // &
      "        comment: ""{oops}""" // nl // "        nested:" // nl // "          a: b" // nl, problems)
    call check(mentions(problems, "'shuffle' must be true or false, found 'maybe'"), msg="logical: " // problems)
    call check(mentions(problems, "'separate_file_per_variable' must be true or false"), msg="second logical: " // problems)
    call check(mentions(problems, "'compression_level' must be an integer, found 'high'"), msg="integer: " // problems)
    call check(mentions(problems, "metadata 'title': unknown substitution '{colour}'"), msg="stream metadata: " // problems)
    call check(mentions(problems, "stream numbers must be integers, found 'two'"), msg="stream number: " // problems)
    call check(mentions(problems, "'stream' must be an integer, found 'first'"), msg="variable stream: " // problems)
    call check(mentions(problems, "netcdf_name: unknown substitution '{oops}'"), msg="netcdf_name: " // problems)
    call check(mentions(problems, "metadata cannot set 'units'"), msg="reserved attribute: " // problems)
    call check(mentions(problems, "metadata 'comment': unknown substitution '{oops}'"), msg="variable metadata: " // problems)
    call check(mentions(problems, "metadata 'nested' must be a single value"), msg="nested metadata: " // problems)

    config = build( &
      "streams:" // nl // "    1:" // nl // "        file_name: a.nc" // nl // "        frequency: daily" // nl // &
      "variables:" // nl // "    - name: GPP" // nl // "      stream: 1" // nl // "      aggregation: mean" // nl // &
      "      netcdf_name: ""dir/GPP""" // nl, problems)
    call check(mentions(problems, "must not be empty or contain '/'"), msg="slash in a NetCDF name: " // problems)
  end subroutine test_value_and_metadata_problems

  subroutine test_substitution()
    character(len=16), parameter :: names(2) = [character(len=16) :: "field_name", "frequency"]
    character(len=16), parameter :: values(2) = [character(len=16) :: "GPP", "daily"]
    character(:), allocatable :: error

    call check(cable_output_substitute("{field_name}_{frequency}", names, values, error) == "GPP_daily", msg="both names")
    call check(.not. allocated(error), msg="no error for a valid template")
    call check(cable_output_substitute("plain", names, values, error) == "plain", msg="no substitutions")
    call check(cable_output_substitute("{field_name}{field_name}", names, values, error) == "GPPGPP", msg="repeated name")
    call check(cable_output_substitute("", names, values, error) == "", msg="empty template")

    call check(cable_output_substitute("{aggregation_method}", names, values, error) == "{aggregation_method}", &
      msg="unknown names are left alone")
    call check(allocated(error), msg="unknown name is an error")
    call check(index(error, "unknown substitution '{aggregation_method}'") > 0, msg="error names the substitution")

    call check(cable_output_substitute("{field_name", names, values, error) == "{field_name", msg="unmatched brace text")
    call check(allocated(error), msg="unmatched brace is an error")
  end subroutine test_substitution

end module test_cable_output_config
