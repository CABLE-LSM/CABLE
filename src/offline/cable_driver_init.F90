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
    gw_params,                    &
    cable_runtime
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
    globalMetfile,                 &
    set_group_output_values
  USE casadimension, ONLY : icycle
  USE casavariable, ONLY : casafile
  USE cable_namelist_util, ONLY : &
    get_namelist_file_name,       &
    CABLE_NAMELIST,               &
    arg_not_namelist
  USE cable_mpi_mod, ONLY : mpi_grp_t
  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER :: CASAONLY_ICYCLE_MIN = 10

  LOGICAL, SAVE, PUBLIC :: vegparmnew    = .FALSE. ! using new format input file (BP dec 2007)
  LOGICAL, SAVE, PUBLIC :: spinup        = .FALSE. ! model spinup to soil state equilibrium?
  LOGICAL, SAVE, PUBLIC :: spincasa      = .FALSE. ! TRUE: CASA-CNP Will spin mloop times, FALSE: no spin up
  LOGICAL, SAVE, PUBLIC :: CASAONLY      = .FALSE. ! ONLY Run CASA-CNP
  LOGICAL, SAVE, PUBLIC :: l_casacnp     = .FALSE. ! using CASA-CNP with CABLE
  LOGICAL, SAVE, PUBLIC :: l_landuse     = .FALSE. ! using CASA-CNP with CABLE
  LOGICAL, SAVE, PUBLIC :: l_laiFeedbk   = .FALSE. ! using prognostic LAI
  LOGICAL, SAVE, PUBLIC :: l_vcmaxFeedbk = .FALSE. ! using prognostic Vcmax

  REAL, SAVE, PUBLIC :: delsoilM ! allowed variation in soil moisture for spin up
  REAL, SAVE, PUBLIC :: delsoilT ! allowed variation in soil temperature for spin up
  REAL, SAVE, PUBLIC :: delgwM = 1e-4

  INTEGER, SAVE, PUBLIC :: LALLOC = 0 ! allocation coefficient for passing to spincasa

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
    CASAONLY,       &
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

  SUBROUTINE cable_driver_init(mpi_grp, trunk_sumbal, NRRRR)
    !! Model initialisation routine for the CABLE offline driver.
    TYPE(mpi_grp_t), INTENT(IN) :: mpi_grp !! MPI group to use
    DOUBLE PRECISION, INTENT(OUT) :: trunk_sumbal
      !! Reference value for quasi-bitwise reproducibility checks.
    INTEGER, INTENT(OUT) :: NRRRR !! Number of repeated spin-up cycles

    INTEGER :: ioerror, unit
    CHARACTER(len=4) :: cRank ! for worker-logfiles

    !check to see if first argument passed to cable is
    !the name of the namelist file
    !if not use cable.nml
    CALL get_namelist_file_name()

    IF (mpi_grp%rank == 0) THEN
      WRITE(*,*) "THE NAME LIST IS ", CABLE_NAMELIST
    END IF

    ! Open, read and close the namelist file.
    OPEN(NEWUNIT=unit, FILE=CABLE_NAMELIST, STATUS="OLD", ACTION="READ")
    READ(unit, NML=CABLE)
    CLOSE(unit)

    cable_runtime%offline = .TRUE.

    ! Open, read and close the consistency check file.
    ! Check triggered by cable_user%consistency_check = .TRUE. in cable.nml
    IF (mpi_grp%rank == 0 .AND. cable_user%consistency_check) THEN
      OPEN(NEWUNIT=unit, FILE=filename%trunk_sumbal, STATUS='old', ACTION='READ', IOSTAT=ioerror)
      IF(ioerror == 0) THEN
        READ(unit, *) trunk_sumbal  ! written by previous trunk version
      END IF
      CLOSE(unit)
    END IF

    ! Open log file:
    IF (mpi_grp%rank == 0) THEN
      OPEN(logn, FILE=filename%log)
    ELSE IF (cable_user%logworker) THEN
      WRITE(cRank, FMT='(I4.4)') mpi_grp%rank
      OPEN(NEWUNIT=logn, FILE="cable_log_"//cRank, STATUS="REPLACE")
    ELSE
      OPEN(NEWUNIT=logn, FILE="/dev/null")
    END IF

    IF (IARGC() > 0 .AND. arg_not_namelist) THEN
      CALL GETARG(1, filename%met)
      CALL GETARG(2, casafile%cnpipool)
    END IF

    ! Initialise flags to output individual variables according to group
    ! options from the namelist file
    IF (mpi_grp%rank == 0) THEN
      CALL set_group_output_values()
    END IF

    IF (TRIM(cable_user%POPLUC_RunType) == 'static') THEN
      cable_user%POPLUC= .FALSE.
    END IF

    ! TODO(Sean): we should not be setting namelist parameters in the following if
    ! block - all options are all configurable via the namelist file and is
    ! unclear that these options are being overwritten. A better approach would be
    ! to error for bad combinations of namelist parameters.
    IF (icycle > CASAONLY_ICYCLE_MIN) THEN
      icycle                     = icycle - CASAONLY_ICYCLE_MIN
      CASAONLY                   = .TRUE.
      CABLE_USER%CASA_DUMP_READ  = .TRUE.
      CABLE_USER%CASA_DUMP_WRITE = .FALSE.
    ELSE IF (icycle == 0) THEN
      CABLE_USER%CASA_DUMP_READ  = .FALSE.
      spincasa                   = .FALSE.
      CABLE_USER%CALL_POP        = .FALSE.
    END IF

    ! TODO(Sean): overwriting l_casacnp defeats the purpose of it being a namelist
    ! option - we should either deprecate the l_casacnp option or not overwrite
    ! its value.
    l_casacnp = icycle > 0

    IF (l_casacnp .AND. (icycle == 0 .OR. icycle > 3)) THEN
      STOP 'icycle must be 1 to 3 when using casaCNP'
    END IF
    IF ((l_laiFeedbk .OR. l_vcmaxFeedbk) .AND. (.NOT. l_casacnp)) THEN
      STOP 'casaCNP required to get prognostic LAI or Vcmax'
    END IF
    IF (l_vcmaxFeedbk .AND. (icycle < 2 .OR. icycle > 3)) THEN
      STOP 'icycle must be 2 to 3 to get prognostic Vcmax'
    END IF
    IF (icycle > 0 .AND. (.NOT. soilparmnew)) THEN
      STOP 'casaCNP must use new soil parameters'
    END IF

    ! vh_js ! suggest LALLOC should ulitmately be a switch in the .nml file
    IF (cable_user%CALL_POP) THEN
      LALLOC = 3 ! for use with POP: makes use of pipe model to partition between stem and leaf
    END IF

    NRRRR = MERGE(MAX(cable_user%CASA_NREP,1), 1, CASAONLY)

  END SUBROUTINE cable_driver_init

END MODULE cable_driver_init_mod
