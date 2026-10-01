! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

!> Entry point for running Fortuno tests.
program testapp
  use fortuno_interface_mod, only: execute_cmd_app
  use fortuno_interface_mod, only: test_list => test_list_t
  use test_cable_netcdf, only: cable_netcdf_test_list
  use test_cable_yaml, only: cable_yaml_test_list
  use test_cable_netcdf_compression, only: cable_netcdf_compression_test_list
  use test_cable_timing, only: cable_timing_test_list
  use test_cable_grid_reductions, only: cable_grid_reductions_test_list
  use test_cable_output_catalogue, only: cable_output_catalogue_test_list
  use test_cable_output_config, only: cable_output_config_test_list
  use test_cable_output_streams, only: cable_output_streams_test_list
  implicit none

  call execute_cmd_app(test_list([ &
      cable_netcdf_test_list(), &
      cable_yaml_test_list(), &
      cable_netcdf_compression_test_list(), &
      cable_timing_test_list(), &
      cable_grid_reductions_test_list(), &
      cable_output_catalogue_test_list(), &
      cable_output_config_test_list(), &
      cable_output_streams_test_list() &
  ]))

end program testapp
