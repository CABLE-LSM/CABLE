! CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
! Copyright (c) 2015, Commonwealth Scientific and Industrial Research Organisation 
! (CSIRO) ABN 41 687 119 230.

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
  USE cable_mpi_mod, ONLY : mpi_grp_t
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

  SUBROUTINE cable_driver_init(mpi_grp)
    TYPE(mpi_grp_t), INTENT(IN) :: mpi_grp !! MPI group to use

    !! Model initialisation routine for the CABLE offline driver.

    !check to see if first argument passed to cable is
    !the name of the namelist file
    !if not use cable.nml
    CALL get_namelist_file_name()

    IF (mpi_grp%rank == 0) THEN
      WRITE(*,*) "THE NAME LIST IS ", CABLE_NAMELIST
    END IF

    ! Open, read and close the namelist file.
    OPEN(10, FILE=CABLE_NAMELIST, STATUS="OLD", ACTION="READ")
    READ(10, NML=CABLE)
    CLOSE(10)

    cable_runtime%offline = .TRUE.

  END SUBROUTINE cable_driver_init

END MODULE cable_driver_init_mod
