! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

module test_cable_yaml
  !! Tests for cable_yaml_mod.
  use fortuno_interface_mod, only: test_list_t
  use fortuno_interface_mod, only: test_case
  use fortuno_interface_mod, only: check
  use cable_yaml_mod
  implicit none
  private

  public :: cable_yaml_test_list

  !> The stream/variable structure used in the output configuration design.
  character(len=*), parameter :: output_config_text = &
    "streams:" // new_line('a') // &
    "    1:" // new_line('a') // &
    "        file_name: cable_output.nc" // new_line('a') // &
    "        frequency: monthly" // new_line('a') // &
    "        netcdf_name: ""{field_name}_{aggregation}""" // new_line('a') // &
    "        shuffle: true" // new_line('a') // &
    "        compression_level: 4" // new_line('a') // &
    "        metadata:" // new_line('a') // &
    "            model: CABLE" // new_line('a') // &
    "    2:" // new_line('a') // &
    "        file_name: cable_evaporation_3hrly.nc" // new_line('a') // &
    "        frequency: 3hrly" // new_line('a') // &
    "variables:" // new_line('a') // &
    "    - name: GPP" // new_line('a') // &
    "      stream: 1" // new_line('a') // &
    "      aggregation: mean" // new_line('a') // &
    "    - name: Qle" // new_line('a') // &
    "      stream: 2" // new_line('a') // &
    "      aggregation: instant" // new_line('a')

contains

  function cable_yaml_test_list() result(test_list)
    type(test_list_t) :: test_list

    test_list = test_list_t([ &
      test_case("test_cable_yaml_mapping_with_integer_keys", test_mapping_with_integer_keys), &
      test_case("test_cable_yaml_list_of_mappings", test_list_of_mappings), &
      test_case("test_cable_yaml_scalar_conversions", test_scalar_conversions), &
      test_case("test_cable_yaml_quoted_braces_preserved", test_quoted_braces_preserved), &
      test_case("test_cable_yaml_nested_mapping", test_nested_mapping), &
      test_case("test_cable_yaml_key_order_preserved", test_key_order_preserved), &
      test_case("test_cable_yaml_has_reports_missing_keys", test_has_reports_missing_keys) &
    ])
  end function cable_yaml_test_list

  subroutine test_mapping_with_integer_keys()
    type(cable_yaml_node_t) :: root, streams, first_stream, second_stream

    root = cable_yaml_parse_text(output_config_text)
    streams = root%get("streams")
    first_stream = streams%at(1)
    second_stream = streams%at(2)
    call check(streams%is_mapping(), msg="streams should be a mapping")
    call check(streams%size() == 2, msg="expected two streams")
    call check(first_stream%key == "1", msg="first stream key should be 1")
    call check(second_stream%key == "2", msg="second stream key should be 2")
    call check(second_stream%get_string("file_name") == "cable_evaporation_3hrly.nc", msg="second stream file name")
  end subroutine test_mapping_with_integer_keys

  subroutine test_list_of_mappings()
    type(cable_yaml_node_t) :: root, variables, first_variable, second_variable

    root = cable_yaml_parse_text(output_config_text)
    variables = root%get("variables")
    first_variable = variables%at(1)
    second_variable = variables%at(2)
    call check(variables%is_sequence(), msg="variables should be a sequence")
    call check(variables%size() == 2, msg="expected two variables")
    call check(first_variable%get_string("name") == "GPP", msg="first variable name")
    call check(first_variable%get_integer("stream") == 1, msg="first variable stream")
    call check(first_variable%get_string("aggregation") == "mean", msg="first variable aggregation")
    ! Items must not be merged: the second item keeps its own values.
    call check(second_variable%get_string("name") == "Qle", msg="second variable name")
    call check(second_variable%get_integer("stream") == 2, msg="second variable stream")
    call check(second_variable%get_string("aggregation") == "instant", msg="second variable aggregation")
  end subroutine test_list_of_mappings

  subroutine test_scalar_conversions()
    type(cable_yaml_node_t) :: root, streams, stream, frequency

    root = cable_yaml_parse_text(output_config_text)
    streams = root%get("streams")
    stream = streams%at(1)
    frequency = stream%get("frequency")
    call check(stream%get_logical("shuffle"), msg="shuffle should be true")
    call check(stream%get_integer("compression_level") == 4, msg="compression level")
    call check(stream%get_string("compression_level") == "4", msg="integer as string")
    call check(stream%get_real("compression_level") == 4.0, msg="integer as real")
    call check(frequency%is_scalar(), msg="frequency should be a scalar")
  end subroutine test_scalar_conversions

  subroutine test_quoted_braces_preserved()
    type(cable_yaml_node_t) :: root, streams, stream

    root = cable_yaml_parse_text(output_config_text)
    streams = root%get("streams")
    stream = streams%at(1)
    call check(stream%get_string("netcdf_name") == "{field_name}_{aggregation}", &
      msg="quotes should be removed and braces kept")
  end subroutine test_quoted_braces_preserved

  subroutine test_nested_mapping()
    type(cable_yaml_node_t) :: root, streams, stream, metadata

    root = cable_yaml_parse_text(output_config_text)
    streams = root%get("streams")
    stream = streams%at(1)
    metadata = stream%get("metadata")
    call check(metadata%is_mapping(), msg="metadata should be a mapping")
    call check(metadata%get_string("model") == "CABLE", msg="metadata model")
  end subroutine test_nested_mapping

  subroutine test_key_order_preserved()
    type(cable_yaml_node_t) :: root, streams, stream, first_entry, second_entry, fifth_entry

    root = cable_yaml_parse_text(output_config_text)
    streams = root%get("streams")
    stream = streams%at(1)
    first_entry = stream%at(1)
    second_entry = stream%at(2)
    fifth_entry = stream%at(5)
    call check(first_entry%key == "file_name", msg="first key")
    call check(second_entry%key == "frequency", msg="second key")
    call check(fifth_entry%key == "compression_level", msg="fifth key")
  end subroutine test_key_order_preserved

  subroutine test_has_reports_missing_keys()
    type(cable_yaml_node_t) :: root, streams, second_stream, variables

    root = cable_yaml_parse_text(output_config_text)
    streams = root%get("streams")
    second_stream = streams%at(2)
    variables = root%get("variables")
    call check(root%has("streams"), msg="streams should be present")
    call check(.not. root%has("groups"), msg="groups should be absent")
    call check(.not. second_stream%has("netcdf_name"), msg="netcdf_name should be absent")
    call check(.not. variables%has("name"), msg="a sequence has no keys")
  end subroutine test_has_reports_missing_keys

end module test_cable_yaml
