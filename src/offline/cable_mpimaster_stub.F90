! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.

MODULE cable_mpimaster
  !! Stub for the master driver when MPI is not available.
  USE cable_mpi_mod, ONLY : mpi_grp_t
  USE CABLE_PLUME_MIP, ONLY : PLUME_MIP_TYPE
  USE CABLE_CRU, ONLY : CRU_TYPE
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: mpidrv_master

CONTAINS

  SUBROUTINE mpidrv_master(comm, dels, koffset, kend, PLUME, CRU, mpi_grp_master)
    !! Stub for when MPI is not available
    INTEGER, INTENT(IN) :: comm
    REAL, INTENT(INOUT) :: dels
    INTEGER, INTENT(INOUT) :: koffset
    INTEGER, INTENT(INOUT) :: kend
    TYPE(PLUME_MIP_TYPE), INTENT(IN) :: PLUME
    TYPE(CRU_TYPE), INTENT(IN) :: CRU
    TYPE(mpi_grp_t), INTENT(INOUT) :: mpi_grp_master

    ! This should never be called!
    STOP

  END SUBROUTINE mpidrv_master

END MODULE cable_mpimaster
