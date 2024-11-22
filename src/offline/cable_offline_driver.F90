! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.
PROGRAM cable_offline_driver
  USE iso_fortran_env, ONLY : error_unit
  USE cable_mpi_mod, ONLY : mpi_grp_t, mpi_mod_init, mpi_mod_end
  USE cable_driver_init_mod
  USE cable_serial
  USE cable_mpimaster
  USE cable_mpiworker
  USE cable_common_module, ONLY : cable_user

  IMPLICIT NONE

  REAL    :: etime ! Declare the type of etime()
  TYPE(mpi_grp_t) :: mpi_grp
  DOUBLE PRECISION :: trunk_sumbal
    !! Reference value for quasi-bitwise reproducibility checks.
  INTEGER :: NRRRR !! Number of repeated spin-up cycles

  call mpi_mod_init()
  mpi_grp = mpi_grp_t()

  CALL cable_driver_init(mpi_grp, trunk_sumbal, NRRRR)

  SELECT CASE(TRIM(cable_user%MetType))
  CASE('gswp', 'gswp3')
    CALL cable_driver_init_gswp(mpi_grp)
  CASE('plum')
  CASE('cru')
  CASE('site')
  CASE('')
  CASE DEFAULT
    WRITE(error_unit,*) "Error: unknown value for cable_user%MetType (", TRIM(cable_user%MetType), ")."
    STOP
  END SELECT

  IF (mpi_grp%size == 1) THEN
    CALL serialdrv(trunk_sumbal, NRRRR)
  ELSE
    IF (mpi_grp%rank == 0) THEN
      CALL mpidrv_master(mpi_grp%comm, trunk_sumbal)
    ELSE
      CALL mpidrv_worker(mpi_grp%comm)
    END IF
  END IF

  CALL mpi_mod_end()

  CALL CPU_TIME(etime)
  PRINT *, 'Finished. ', etime, ' seconds needed for '

END PROGRAM cable_offline_driver
