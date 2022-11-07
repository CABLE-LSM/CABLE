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

  USE cable_mpicommon
  USE cable_mpimaster
  USE cable_mpiworker
  USE cable_namelist_util, ONLY: get_namelist_file_name,&
       CABLE_NAMELIST,arg_not_namelist

  IMPLICIT NONE

  INTEGER :: comm, np, rank, ierr
  REAL    :: etime ! Declare the type of etime()



  CALL MPI_Init (ierr)
  CALL MPI_Comm_dup (MPI_COMM_WORLD, comm, ierr)
  CALL MPI_Comm_size (comm, np, ierr)

  !check to see if first argument passed to cable is
  !the name of the namelist file
  !if not use cable.nml
  CALL get_namelist_file_name()

  IF (np < 2) THEN
     WRITE (*,*) 'This program needs at least 2 processes to run!'
     CALL MPI_Abort (comm, 0, ierr)
  END IF

  CALL MPI_Comm_rank (comm, rank, ierr)

  IF (rank == 0) THEN
     CALL mpidrv_master (comm)
  ELSE
     CALL mpidrv_worker (comm)
  END IF

  CALL MPI_Finalize (ierr)

  CALL CPU_TIME(etime)
  PRINT *, 'Finished. ', etime, ' seconds needed for '

END PROGRAM mpi_driver
