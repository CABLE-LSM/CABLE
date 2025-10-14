module cable_tests
  use fortuno_interface_m, only: test_list
  use test_cable_netcdf, only: cable_netcdf_test_list
  implicit none

contains

  function tests()
    type(test_list) :: tests
    tests = test_list([&
      cable_netcdf_test_list()&
    ])
  end function tests

end module cable_tests

program test_cable
  
  use fortuno_interface_m, only: execute_cmd_app
  use cable_tests, only: tests

  implicit none

  call execute_cmd_app(tests())

end program test_cable
