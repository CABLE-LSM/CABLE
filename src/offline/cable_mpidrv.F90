!==============================================================================
! This source code is part of the
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
!
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) is located
! in each directory containing CABLE code.
!
! ==============================================================================
! Purpose: Bare bones MPI driver for CABLE
!
! Contact: Bernard.Pak@csiro.au
!
! History: MPI wrapper developed by Maciej Golebiewski (2012)
!
! ==============================================================================
!
PROGRAM mpi_driver

  USE mpi
  USE cable_driver_init_mod

  USE cable_mpicommon, ONLY : comm, rank
  USE cable_mpimaster
  USE cable_mpiworker

  IMPLICIT NONE

  INTEGER :: ierr
  REAL    :: etime ! Declare the type of etime()

  CALL cable_driver_init()

  IF (rank == 0) THEN
     CALL mpidrv_master (comm)
  ELSE
     CALL mpidrv_worker (comm)
  END IF

  CALL MPI_Finalize (ierr)

  CALL CPU_TIME(etime)
  PRINT *, 'Finished. ', etime, ' seconds needed for '

END PROGRAM mpi_driver
