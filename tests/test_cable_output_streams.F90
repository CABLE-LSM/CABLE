! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

module test_cable_output_streams
  !! Tests for the output files that an output configuration produces.
  use fortuno_interface_mod, only: test_list_t
  use fortuno_interface_mod, only: test_case
  use fortuno_interface_mod, only: check
  use cable_yaml_mod, only: cable_yaml_parse_text
  use cable_output_mod, only: cable_output_variable_definition_t
  use cable_output_mod, only: cable_output_attribute_t
  use cable_output_mod, only: cable_output_stream_t
  use cable_output_mod, only: cable_output_variable_t
  use cable_output_mod, only: cable_output_get_dimension
  use cable_output_mod, only: cable_output_build_streams
  use cable_output_config_mod, only: cable_output_config_parse
  use cable_output_config_mod, only: cable_output_config_t
  use cable_netcdf_mod, only: CABLE_NETCDF_FLOAT
  use aggregator_mod, only: new_aggregator
  implicit none
  private

  public :: cable_output_streams_test_list

  character(len=*), parameter :: nl = new_line('a')
  real, target :: working_values(4)

contains

  function cable_output_streams_test_list() result(test_list)
    type(test_list_t) :: test_list

    test_list = test_list_t([ &
      test_case("test_streams_one_file_per_stream", test_one_file_per_stream), &
      test_case("test_streams_variable_settings", test_variable_settings), &
      test_case("test_streams_metadata_order_and_cell_methods", test_metadata_order_and_cell_methods), &
      test_case("test_streams_compression_settings_are_copied", test_compression_settings_are_copied), &
      test_case("test_streams_separate_file_per_variable", test_separate_file_per_variable), &
      test_case("test_streams_without_variables_make_no_file", test_streams_without_variables_make_no_file) &
    ])
  end function cable_output_streams_test_list

  function make_definition(name, group, is_parameter, native_frequency, units) result(definition)
    character(len=*), intent(in) :: name, group
    logical, intent(in), optional :: is_parameter
    character(len=*), intent(in), optional :: native_frequency, units
    type(cable_output_variable_definition_t) :: definition

    definition%field_name = name
    definition%group_name = group
    definition%var_type = CABLE_NETCDF_FLOAT
    definition%divide_by = 2.0
    allocate(definition%data_shape(1))
    definition%data_shape(1) = cable_output_get_dimension("patch")
    allocate(definition%aggregator, source=new_aggregator(working_values))
    allocate(definition%metadata(2))
    definition%metadata(1) = cable_output_attribute_t("units", "W/m^2")
    if (present(units)) definition%metadata(1)%value = units
    definition%metadata(2) = cable_output_attribute_t("long_name", "Long name of " // name)
    if (present(is_parameter)) definition%parameter = is_parameter
    if (present(native_frequency)) definition%native_frequency = native_frequency
  end function make_definition

  function test_definitions() result(definitions)
    type(cable_output_variable_definition_t), allocatable :: definitions(:)

    definitions = [ &
      make_definition("GPP", "carbon_cycle"), &
      make_definition("Qle", "energy_cycle", units="W/m^2"), &
      make_definition("Qh", "energy_cycle"), &
      make_definition("iveg", "params", is_parameter=.true.), &
      make_definition("Tmx", "veg", native_frequency="daily") &
    ]
  end function test_definitions

  function configure(text) result(config)
    character(len=*), intent(in) :: text
    type(cable_output_config_t) :: config
    config = cable_output_config_parse(text, test_definitions(), "2000-01-01", "2000-12-31")
  end function configure

  subroutine test_one_file_per_stream()
    type(cable_output_config_t) :: config
    type(cable_output_stream_t), allocatable :: streams(:)

    config = configure( &
      "streams:" // nl // "    1:" // nl // "        file_name: monthly.nc" // nl // "        frequency: monthly" // nl // &
      "    2:" // nl // "        file_name: diurnal.nc" // nl // "        frequency: 3hrly" // nl // &
      "variables:" // nl // "    - name: GPP" // nl // "      stream: 1" // nl // "      aggregation: mean" // nl // &
      "    - name: Qle" // nl // "      stream: 1" // nl // "      aggregation: max" // nl // &
      "    - name: GPP" // nl // "      stream: 2" // nl // "      aggregation: instant" // nl // &
      "      netcdf_name: GPP_instant" // nl)
    streams = cable_output_build_streams(config, test_definitions(), "land")
    call check(size(streams) == 2, msg="one file for each stream")
    call check(trim(streams(1)%file_name) == "monthly.nc", msg="first file name")
    call check(trim(streams(1)%sampling_frequency) == "monthly", msg="first frequency")
    call check(size(streams(1)%output_variables) == 2, msg="two variables in the first file")
    call check(trim(streams(2)%file_name) == "diurnal.nc", msg="second file name")
    call check(trim(streams(2)%sampling_frequency) == "3hrly", msg="second frequency")
    call check(size(streams(2)%output_variables) == 1, msg="one variable in the second file")
    call check(trim(streams(1)%grid_type) == "land", msg="grid type is passed on")
    ! The same definition is used in both files, each use with its own aggregator.
    call check(allocated(streams(1)%output_variables(1)%aggregator), msg="aggregator in first use")
    call check(allocated(streams(2)%output_variables(1)%aggregator), msg="aggregator in second use")
  end subroutine test_one_file_per_stream

  subroutine test_variable_settings()
    type(cable_output_config_t) :: config
    type(cable_output_stream_t), allocatable :: streams(:)

    config = configure( &
      "streams:" // nl // "    1:" // nl // "        file_name: a.nc" // nl // "        frequency: monthly" // nl // &
      "variables:" // nl // "    - name: GPP" // nl // "      stream: 1" // nl // "      aggregation: mean" // nl // &
      "      reduction: grid_cell_average" // nl // "      netcdf_name: GPP_monthly_mean" // nl // &
      "    - name: Tmx" // nl // "      stream: 1" // nl // "      aggregation: max" // nl // &
      "    - name: iveg" // nl // "      stream: 1" // nl // "      aggregation: instant" // nl)
    streams = cable_output_build_streams(config, test_definitions(), "mask")
    associate(gpp => streams(1)%output_variables(1), tmx => streams(1)%output_variables(2), iveg => streams(1)%output_variables(3))
      call check(trim(gpp%field_name) == "GPP", msg="field name")
      call check(trim(gpp%netcdf_name) == "GPP_monthly_mean", msg="NetCDF name from the configuration")
      call check(trim(gpp%aggregation_method) == "mean", msg="aggregation")
      call check(trim(gpp%reduction_method) == "grid_cell_average", msg="reduction")
      call check(gpp%divide_by == 2.0, msg="conversion is taken from the definition")
      call check(size(gpp%data_shape) == 1, msg="shape is taken from the definition")
      call check(trim(gpp%accumulation_frequency) == "timestep", msg="accumulated every timestep by default")
      call check(trim(tmx%accumulation_frequency) == "daily", msg="native frequency becomes the accumulation frequency")
      call check(trim(tmx%reduction_method) == "none", msg="no reduction by default")
      call check(iveg%parameter, msg="parameter flag is taken from the definition")
    end associate
  end subroutine test_variable_settings

  subroutine test_metadata_order_and_cell_methods()
    type(cable_output_config_t) :: config
    type(cable_output_stream_t), allocatable :: streams(:)

    config = configure( &
      "streams:" // nl // "    1:" // nl // "        file_name: a.nc" // nl // "        frequency: monthly" // nl // &
      "variables:" // nl // "    - name: GPP" // nl // "      stream: 1" // nl // "      aggregation: mean" // nl // &
      "      reduction: grid_cell_average" // nl // "      metadata:" // nl // "        description: ""Computed using X.""" // nl // &
      "    - name: Qle" // nl // "      stream: 1" // nl // "      aggregation: sum" // nl // &
      "    - name: Qh" // nl // "      stream: 1" // nl // "      aggregation: min" // nl // &
      "    - name: Tmx" // nl // "      stream: 1" // nl // "      aggregation: instant" // nl // &
      "    - name: iveg" // nl // "      stream: 1" // nl // "      aggregation: instant" // nl)
    streams = cable_output_build_streams(config, test_definitions(), "land")
    associate(gpp => streams(1)%output_variables(1)%metadata)
      call check(size(gpp) == 4, msg="units, long_name, cell_methods and the user attribute")
      call check(trim(gpp(1)%name) == "units" .and. trim(gpp(1)%value) == "W/m^2", msg="units come first")
      call check(trim(gpp(2)%name) == "long_name", msg="then long_name")
      call check(trim(gpp(3)%name) == "cell_methods", msg="then cell_methods")
      call check(trim(gpp(3)%value) == "area: mean time: mean", msg="area and time methods: " // trim(gpp(3)%value))
      call check(trim(gpp(4)%name) == "description" .and. trim(gpp(4)%value) == "Computed using X.", msg="user attribute last")
    end associate
    call check(trim(streams(1)%output_variables(2)%metadata(3)%value) == "time: sum", msg="sum")
    call check(trim(streams(1)%output_variables(3)%metadata(3)%value) == "time: minimum", msg="minimum")
    call check(trim(streams(1)%output_variables(4)%metadata(3)%value) == "time: point", msg="instant")
    call check(size(streams(1)%output_variables(5)%metadata) == 2, msg="parameters have no cell_methods")
  end subroutine test_metadata_order_and_cell_methods

  subroutine test_compression_settings_are_copied()
    type(cable_output_config_t) :: config
    type(cable_output_stream_t), allocatable :: streams(:)

    config = configure( &
      "streams:" // nl // "    1:" // nl // "        file_name: a.nc" // nl // "        frequency: daily" // nl // &
      "        shuffle: false" // nl // "        compression_level: 5" // nl // &
      "        metadata:" // nl // "            model: CABLE" // nl // &
      "    2:" // nl // "        file_name: b.nc" // nl // "        frequency: daily" // nl // &
      "variables:" // nl // "    - name: GPP" // nl // "      stream: 1" // nl // "      aggregation: mean" // nl // &
      "    - name: GPP" // nl // "      stream: 2" // nl // "      aggregation: mean" // nl)
    streams = cable_output_build_streams(config, test_definitions(), "land")
    call check(.not. streams(1)%shuffle, msg="shuffle off")
    call check(streams(1)%compression_level == 5, msg="compression level")
    call check(streams(2)%shuffle, msg="default shuffle")
    call check(streams(2)%compression_level == 1, msg="default compression level")
    call check(size(streams(1)%metadata) == 1, msg="global attribute")
    call check(trim(streams(1)%metadata(1)%name) == "model" .and. trim(streams(1)%metadata(1)%value) == "CABLE", &
      msg="global attribute value")
    call check(size(streams(2)%metadata) == 0, msg="no global attributes by default")
  end subroutine test_compression_settings_are_copied

  subroutine test_separate_file_per_variable()
    type(cable_output_config_t) :: config
    type(cable_output_stream_t), allocatable :: streams(:)

    config = configure( &
      "streams:" // nl // "    1:" // nl // "        file_name: some/dir/unused.nc" // nl // "        frequency: daily" // nl // &
      "        separate_file_per_variable: true" // nl // &
      "        netcdf_name: ""{field_name}_{aggregation}_{start_date}-{end_date}""" // nl // &
      "        compression_level: 3" // nl // &
      "groups:" // nl // "    - name: energy_cycle" // nl // "      stream: 1" // nl // "      aggregation: mean" // nl)
    streams = cable_output_build_streams(config, test_definitions(), "land")
    call check(size(streams) == 2, msg="one file for each of Qle and Qh")
    call check(trim(streams(1)%file_name) == "some/dir/Qle_mean_2000-01-01-2000-12-31.nc", msg="file name in the directory of file_name: " // trim(streams(1)%file_name))
    call check(trim(streams(2)%file_name) == "some/dir/Qh_mean_2000-01-01-2000-12-31.nc", msg="second file name")
    call check(size(streams(1)%output_variables) == 1, msg="one variable in each file")
    call check(streams(2)%compression_level == 3, msg="stream settings apply to every file")
  end subroutine test_separate_file_per_variable

  subroutine test_streams_without_variables_make_no_file()
    type(cable_output_config_t) :: config
    type(cable_output_stream_t), allocatable :: streams(:)

    config = configure( &
      "streams:" // nl // "    1:" // nl // "        file_name: a.nc" // nl // "        frequency: daily" // nl // &
      "    2:" // nl // "        file_name: unused.nc" // nl // "        frequency: monthly" // nl // &
      "variables:" // nl // "    - name: GPP" // nl // "      stream: 1" // nl // "      aggregation: mean" // nl)
    streams = cable_output_build_streams(config, test_definitions(), "land")
    call check(size(streams) == 1, msg="the empty stream makes no file")
    call check(trim(streams(1)%file_name) == "a.nc", msg="only the used stream remains")
  end subroutine test_streams_without_variables_make_no_file

end module test_cable_output_streams
