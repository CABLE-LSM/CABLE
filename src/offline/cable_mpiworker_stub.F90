! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

MODULE cable_mpiworker
  !! Stub for the worker driver when MPI is not available.
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: mpidrv_worker

CONTAINS

  SUBROUTINE mpidrv_worker(comm)
    !! Stub for when MPI is not available
    INTEGER, INTENT(IN) :: comm

    ! This should never be called!
    STOP

  END SUBROUTINE mpidrv_worker

END MODULE cable_mpiworker
