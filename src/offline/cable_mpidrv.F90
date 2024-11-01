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
  USE cable_mpi_mod, ONLY : mpi_grp_t, mpi_mod_init, mpi_mod_end, mpi_check_error
  USE cable_driver_init_mod

  USE cable_mpimaster
  USE cable_mpiworker

  IMPLICIT NONE

  REAL    :: etime ! Declare the type of etime()
  TYPE(mpi_grp_t) :: mpi_grp

  call mpi_mod_init()
  mpi_grp = mpi_grp_t()

  CALL cable_driver_init(mpi_grp)

  IF (mpi_grp%size < 2) THEN
    WRITE(*,*) 'This program needs at least 2 processes to run!'
    CALL mpi_grp%abort()
  END IF

  IF (mpi_grp%rank == 0) THEN
     CALL mpidrv_master (mpi_grp%comm)
  ELSE
     CALL mpidrv_worker (mpi_grp%comm)
  END IF

  CALL mpi_mod_end()

  CALL CPU_TIME(etime)
  PRINT *, 'Finished. ', etime, ' seconds needed for '

END PROGRAM mpi_driver
