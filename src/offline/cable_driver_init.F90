MODULE cable_driver_init_mod
  !! Module for CABLE offline driver initialisation.
  USE cable_common_module, ONLY : &
    filename,                     &
    calcsoilalbedo,               &
    redistrb,                     &
    wiltParam,                    &
    satuParam,                    &
    cable_user,                   &
    gw_params
  USE cable_IO_vars_module, ONLY : &
    soilparmnew,                   &
    output,                        &
    patchout,                      &
    check,                         &
    verbose,                       &
    leaps,                         &
    logn,                          &
    fixedCO2,                      &
    ncciy,                         &
    gswpfile,                      &
    globalMetfile
  USE casadimension, ONLY : icycle
  USE casavariable, ONLY : casafile
  USE cable_namelist_util, ONLY : &
    get_namelist_file_name,       &
    CABLE_NAMELIST
#ifdef __MPI__
  USE mpi
  USE cable_mpicommon, ONLY : comm, rank
#endif
  IMPLICIT NONE
  PRIVATE

  LOGICAL, SAVE, PUBLIC :: vegparmnew    = .FALSE. ! using new format input file (BP dec 2007)
  LOGICAL, SAVE, PUBLIC :: spinup        = .FALSE. ! model spinup to soil state equilibrium?
  LOGICAL, SAVE, PUBLIC :: spincasa      = .FALSE. ! TRUE: CASA-CNP Will spin mloop times, FALSE: no spin up
  LOGICAL, SAVE, PUBLIC :: l_casacnp     = .FALSE. ! using CASA-CNP with CABLE
  LOGICAL, SAVE, PUBLIC :: l_landuse     = .FALSE. ! using CASA-CNP with CABLE
  LOGICAL, SAVE, PUBLIC :: l_laiFeedbk   = .FALSE. ! using prognostic LAI
  LOGICAL, SAVE, PUBLIC :: l_vcmaxFeedbk = .FALSE. ! using prognostic Vcmax

  REAL, SAVE, PUBLIC :: delsoilM ! allowed variation in soil moisture for spin up
  REAL, SAVE, PUBLIC :: delsoilT ! allowed variation in soil temperature for spin up
  REAL, SAVE, PUBLIC :: delgwM = 1e-4

  NAMELIST /CABLE/  &
    filename,       & ! TYPE, containing input filenames
    vegparmnew,     & ! use new soil param. method
    soilparmnew,    & ! use new soil param. method
    calcsoilalbedo, & ! albedo considers soil color Ticket #27
    spinup,         & ! spinup model (soil) to steady state
    delsoilM,       &
    delsoilT,       &
    delgwM,         &
    output,         &
    patchout,       &
    check,          &
    verbose,        &
    leaps,          &
    logn,           &
    fixedCO2,       &
    spincasa,       &
    l_casacnp,      &
    l_landuse,      &
    l_laiFeedbk,    &
    l_vcmaxFeedbk,  &
    icycle,         &
    casafile,       &
    ncciy,          &
    gswpfile,       &
    globalMetfile,  &
    redistrb,       &
    wiltParam,      &
    satuParam,      &
    cable_user,     & ! additional USER switches
    gw_params

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

#ifndef __MPI__
    WRITE(*,*) "THE NAME LIST IS ", CABLE_NAMELIST
#endif
    ! Open, read and close the namelist file.
    OPEN(10, FILE=CABLE_NAMELIST, STATUS="OLD", ACTION="READ")
    READ(10, NML=CABLE)
    CLOSE(10)

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
