! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation
! (CSIRO) ABN 41 687 119 230.
PROGRAM cable_offline_driver
  USE iso_fortran_env, ONLY : error_unit
  USE cable_mpi_mod, ONLY : mpi_grp_t, mpi_mod_init, mpi_mod_end
  USE cable_driver_common_mod, ONLY: &
      cable_driver_init,             &
      cable_driver_init_gswp,        &
      cable_driver_init_plume,       &
      cable_driver_init_cru,         &
      cable_driver_init_site,        &
      cable_driver_init_default
  USE cable_serial
  USE cable_mpimaster
  USE cable_mpiworker
  USE cable_common_module, ONLY : cable_user
  USE CABLE_PLUME_MIP, ONLY : PLUME_MIP_TYPE
  USE CABLE_CRU, ONLY : CRU_TYPE
  USE CABLE_site, ONLY : site_TYPE

  IMPLICIT NONE

  REAL    :: etime ! Declare the type of etime()
  TYPE(mpi_grp_t) :: mpi_grp
  DOUBLE PRECISION :: trunk_sumbal
    !! Reference value for quasi-bitwise reproducibility checks.
  INTEGER :: NRRRR !! Number of repeated spin-up cycles
  REAL :: dels !! Time step size in seconds
  INTEGER :: koffset = 0 !! Timestep to start at
  INTEGER :: kend !! No. of time steps in run
  INTEGER, ALLOCATABLE :: GSWP_MID(:,:) !! NetCDF file IDs for GSWP met forcing
  TYPE(PLUME_MIP_TYPE) :: PLUME
  TYPE(CRU_TYPE) :: CRU
  TYPE (site_TYPE) :: site

  call mpi_mod_init()
  mpi_grp = mpi_grp_t()

  CALL cable_driver_init(mpi_grp, trunk_sumbal, NRRRR)

  SELECT CASE(TRIM(cable_user%MetType))
  CASE('gswp')
    CALL cable_driver_init_gswp(mpi_grp, GSWP_MID, NRRRR)
  CASE('gswp3')
    CALL cable_driver_init_gswp(mpi_grp)
  CASE('prin')
    CALL cable_driver_init_gswp(mpi_grp, GSWP_MID, NRRRR)
  CASE('plum')
    CALL cable_driver_init_plume(dels, koffset, PLUME)
  CASE('cru')
    CALL cable_driver_init_cru(dels, koffset, CRU)
  CASE('site')
    CALL cable_driver_init_site(site)
    CALL cable_driver_init_default(dels, koffset, kend)
  CASE('')
    CALL cable_driver_init_default(dels, koffset, kend)
  CASE DEFAULT
    WRITE(error_unit,*) "Error: unknown value for cable_user%MetType (", TRIM(cable_user%MetType), ")."
    STOP
  END SELECT

  IF (mpi_grp%size == 1) THEN
    CALL serialdrv(trunk_sumbal, NRRRR, dels, koffset, kend, GSWP_MID, PLUME, CRU, site)
  ELSE
    IF (mpi_grp%rank == 0) THEN
      CALL mpidrv_master(mpi_grp%comm, trunk_sumbal, dels, koffset, kend, PLUME, CRU)
    ELSE
      CALL mpidrv_worker(mpi_grp%comm)
    END IF
  END IF

  CALL mpi_mod_end()

  CALL CPU_TIME(etime)
  PRINT *, 'Finished. ', etime, ' seconds needed for '

END PROGRAM cable_offline_driver
