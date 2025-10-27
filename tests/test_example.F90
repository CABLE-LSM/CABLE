! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

module test_example_mod
  use fortuno_interface_mod, only: test_list_t
  use fortuno_interface_mod, only: test_case
  use fortuno_interface_mod, only: check
  implicit none
  private

  public :: test_example

contains

  function test_example() result(example_tests)
    type(test_list_t) :: example_tests
    example_tests = test_list_t([ &
      test_case("example", example) &
    ])
  end function test_example

  subroutine example()
    !! Example test.
    call check(1 /= 2, msg="1 should not equal 2")
  end subroutine example

end module
