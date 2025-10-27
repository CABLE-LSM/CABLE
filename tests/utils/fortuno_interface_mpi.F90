! Copyright (C) 2024 Alex Buccheri
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at https://mozilla.org/MPL/2.0/.

module fortuno_interface_mod
  !! Expose Fortuno MPI data types and routines through common aliases.

  use fortuno_mpi, only: test_case_t => mpi_case_base
  use fortuno_mpi, only: test_item_t => test_item
  use fortuno_mpi, only: test_list_t => test_list

  use fortuno_mpi, only: execute_cmd_app => execute_mpi_cmd_app
  use fortuno_mpi, only: test_suite => mpi_suite_item
  use fortuno_mpi, only: test_case => mpi_case_item

  use fortuno_mpi, only: check => mpi_check
  use fortuno_mpi, only: check_failed => mpi_check_failed
  use fortuno_mpi, only: skip => mpi_skip
  use fortuno_mpi, only: is_equal
  use fortuno_mpi, only: all_equal
  use fortuno_mpi, only: all_close

  use fortuno_mpi, only: global_comm
  use fortuno_mpi, only: num_ranks
  use fortuno_mpi, only: this_rank

  implicit none

  ! Scope is deliberately public

end module fortuno_interface_mod
