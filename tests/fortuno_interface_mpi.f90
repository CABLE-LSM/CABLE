!! Copyright (C) 2024 Alex Buccheri
!!
!! This Source Code Form is subject to the terms of the Mozilla Public
!! License, v. 2.0. If a copy of the MPL was not distributed with this
!! file, You can obtain one at https://mozilla.org/MPL/2.0/.
!!

!> @brief Expose Fortuno MPI data types and routines through common aliases.
!!
!!      Alias                   MPI routine                       Description
!! ------------------------------------------------------------------------------------------------------
!! execute_cmd_app            execute_mpi_cmd_app     Accepts an array of test_item, and runs them.
!! test_case_t                mpi_case_base           Base type for representing a test case.
!! test_case                  mpi_case_item           Returns a test case instace as a generic test item.
!! suite                      mpi_suite_item          Returns a suite instance wrapped as test_item.
!! check                      mpi_check               Perform a logical check (assertion) on a condition.
!! check_failed               mpi_check_failed        Returns .true. if the previous check failed, .false. otherwise.
!!

module fortuno_interface_m

  use fortuno_mpi, only : &
    execute_cmd_app => execute_mpi_cmd_app, &
    test_case_t => mpi_case_base, &
    test_case => mpi_case_item, &
    suite => mpi_suite_item, &
    check => mpi_check, &
    check_failed => mpi_check_failed, &
    is_equal, &
    all_equal, &
    all_close, &
    test_item, &
    test_list, &
    global_comm, &
    num_ranks, &
    this_rank

  implicit none

  ! Scope is deliberately public

end module fortuno_interface_m
