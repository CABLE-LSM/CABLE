! Copyright (C) 2024 Alex Buccheri
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at https://mozilla.org/MPL/2.0/.

module fortuno_interface_mod
  !! Expose Fortuno serial data types and routines through common aliases.

  use fortuno_serial, only: test_case_t => serial_case_base
  use fortuno_serial, only: test_item_t => test_item
  use fortuno_serial, only: test_list_t => test_list

  use fortuno_serial, only: execute_cmd_app => execute_serial_cmd_app
  use fortuno_serial, only: test_suite => serial_suite_item
  use fortuno_serial, only: test_case => serial_case_item

  use fortuno_serial, only: check => serial_check
  use fortuno_serial, only: check_failed => serial_check_failed
  use fortuno_serial, only: skip => serial_skip
  use fortuno_serial, only: is_equal
  use fortuno_serial, only: all_equal
  use fortuno_serial, only: all_close

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

end module fortuno_interface_mod
