MODULE cable_driver_init_mod
  !! Module for CABLE offline driver initialisation.
  USE cable_namelist_util, ONLY : get_namelist_file_name
#ifdef __MPI__
  USE mpi
  USE cable_mpicommon, ONLY : comm, rank
#endif
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: cable_driver_init

CONTAINS

  SUBROUTINE cable_driver_init()
    !! Model initialisation routine for the CABLE offline driver.
#ifdef __MPI__
    INTEGER :: np, ierr
#endif

    !check to see if first argument passed to cable is
    !the name of the namelist file
    !if not use cable.nml
    CALL get_namelist_file_name()

#ifdef __MPI__
    CALL MPI_Init(ierr)
    CALL MPI_Comm_dup(MPI_COMM_WORLD, comm, ierr)
    CALL MPI_Comm_size(comm, np, ierr)

    IF (np < 2) THEN
      WRITE (*,*) 'This program needs at least 2 processes to run!'
      CALL MPI_Abort(comm, 0, ierr)
    END IF

    CALL MPI_Comm_rank(comm, rank, ierr)
#endif

  END SUBROUTINE cable_driver_init

END MODULE cable_driver_init_mod
