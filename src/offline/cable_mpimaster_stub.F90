! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

MODULE cable_mpimaster
  !! Stub for the master driver when MPI is not available.
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: mpidrv_master

CONTAINS

  SUBROUTINE mpidrv_master(comm)
    !! Stub for when MPI is not available
    INTEGER, INTENT(IN) :: comm

    ! This should never be called!
    STOP

  END SUBROUTINE mpidrv_master

END MODULE cable_mpimaster
