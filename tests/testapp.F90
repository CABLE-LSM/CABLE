! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

!> Entry point for running Fortuno tests.
program testapp
  use fortuno_interface_mod, only: execute_cmd_app
  use fortuno_interface_mod, only: test_list => test_list_t
  use test_cable_netcdf, only: cable_netcdf_test_list
  implicit none

  call execute_cmd_app(test_list([ &
      cable_netcdf_test_list() &
  ]))

end program testapp
