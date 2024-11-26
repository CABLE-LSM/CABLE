! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

MODULE cable_mpimaster
  !! Stub for the master driver when MPI is not available.
  USE CABLE_PLUME_MIP, ONLY : PLUME_MIP_TYPE
  USE CABLE_CRU, ONLY : CRU_TYPE
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: mpidrv_master

CONTAINS

  SUBROUTINE mpidrv_master(comm, trunk_sumbal, dels, koffset, kend, PLUME, CRU)
    !! Stub for when MPI is not available
    INTEGER, INTENT(IN) :: comm
    DOUBLE PRECISION, INTENT(IN) :: trunk_sumbal
    REAL, INTENT(INOUT) :: dels
    INTEGER, INTENT(INOUT) :: koffset
    INTEGER, INTENT(INOUT) :: kend
    TYPE(PLUME_MIP_TYPE), INTENT(IN) :: PLUME
    TYPE(CRU_TYPE), INTENT(IN) :: CRU

    ! This should never be called!
    STOP

  END SUBROUTINE mpidrv_master

END MODULE cable_mpimaster
