!! Copyright (C) 2024 Alex Buccheri
!!
!! This Source Code Form is subject to the terms of the Mozilla Public
!! License, v. 2.0. If a copy of the MPL was not distributed with this
!! file, You can obtain one at https://mozilla.org/MPL/2.0/.
!!

!> @brief Expose Fortuno serail data types and routines through common aliases.
!!
!!      Alias                Serial routine                      Description
!! ------------------------------------------------------------------------------------------------------
!! execute_cmd_app        execute_serial_cmd_app     Accepts an array of test_item, and runs them.
!! test_case_t            serial_case_base           Base type for representing a test case.
!! test_case              serial_case_item           Returns a test case instace as a generic test item.
!! suite                  serial_suite_item          Returns a suite instance wrapped as test_item.
!! check                  serial_check               Perform a logical check (assertion) on a condition.
!! check_failed           serial_check_failed        Returns .true. if the previous check failed, .false. otherwise.
!!

module fortuno_interface_m

  use fortuno_serial, only : &
    execute_cmd_app => execute_serial_cmd_app, &
    test_case_t => serial_case_base, &
    test_case => serial_case_item, &
    suite => serial_suite_item, &
    check => serial_check, &
    check_failed => serial_check_failed, &
    is_equal, &
    all_equal, &
    all_close, &
    test_item, &
    test_list

  implicit none

  ! Scope is deliberately public

contains

  function global_comm()
    !! Serial stub for global_comm function.
    use cable_mpi_stub_types_mod, only: MPI_COMM, MPI_COMM_NULL
    type(MPI_COMM) global_comm
    global_comm = MPI_COMM_NULL
  end function

  integer function num_ranks()
    !! Serial stub for num_ranks function.
    num_ranks = 1
  end function

  integer function this_rank()
    !! Serial stub for this_rank function.
    this_rank = 0
  end function

end module fortuno_interface_m
