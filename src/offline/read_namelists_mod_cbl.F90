MODULE read_namelists_mod_cbl

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='READ_NAMELISTS_MOD_CBL'

CONTAINS

SUBROUTINE     read_cable_namelist( unitnumber, cable_namelist, vegparmnew,     &
                                             spinup    ,      &
                                             spincasa  ,      &
                                             CASAONLY  ,      &
                                             l_casacnp ,      &
                                             l_landuse ,      &
                                             l_laiFeedbk ,    &
                                             l_vcmaxFeedbk    )
! Description:
!  Read the cable namelist

USE cable_common_module, ONLY: &
  filename,                    &
  calcsoilalbedo,              &
  redistrb,                    &
  wiltParam,                   &
  satuParam,                   &
  snmin,                       &
  cable_user,                  &
  gw_params
  
USE cable_IO_vars_module, ONLY: &
  soilparmnew,                  &
  output,                       &
  patchout,                     &
  check,                        &
  verbose,                      &
  leaps,                        &
  logn,                         &
  fixedCO2,                     &
  ncciy,                        &
  gswpfile,                     &
  globalMetfile
  
USE casadimension, ONLY: icycle
USE casavariable,  ONLY: casafile
   
IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: unitnumber
CHARACTER(LEN=*), INTENT(IN) :: cable_namelist

! additional declarations to satisfy namelist dec
CHARACTER(LEN=*), PARAMETER :: RoutineName='read_cable_namelist'
LOGICAL, INTENT(OUT) :: vegparmnew     ! using new format input file (BP dec 2007)
LOGICAL, INTENT(OUT) :: spinup         ! model spinup to soil state equilibrium?
LOGICAL, INTENT(OUT) :: spincasa       ! TRUE: CASA-CNP Will spin mloop times, FALSE: no spin up
LOGICAL, INTENT(OUT) :: CASAONLY       ! ONLY Run CASA-CNP
LOGICAL, INTENT(OUT) :: l_casacnp      ! using CASA-CNP with CABLE
LOGICAL, INTENT(OUT) :: l_landuse      ! using CASA-CNP with CABLE
LOGICAL, INTENT(OUT) :: l_laiFeedbk    ! using prognostic LAI
LOGICAL, INTENT(OUT) :: l_vcmaxFeedbk  ! using prognostic Vcmax

REAL :: delsoilM ! allowed variation in soil moisture for spin up
REAL :: delsoilT ! allowed variation in soil temperature for spin up
REAL :: delgwM = 1e-4

INTEGER :: LALLOC = 0            ! alloc coeff passed to spincasa

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
snmin,          &
cable_user,     & ! additional USER switches
gw_params

! Open, read and close the namelist file.
OPEN( UNIT=unitnumber, FILE=CABLE_NAMELIST, STATUS="OLD", ACTION="READ" )
  READ( unitnumber, NML=CABLE )
CLOSE( unitnumber )

RETURN
END SUBROUTINE read_cable_namelist

!SUBROUTINE read_cable_model_environment(unitnumber)
!
!! Description:
!! Read the cable namelist
!
!USE cable_model_env_mod, ONLY: read_nml_cable_model_env
!USE cable_model_env_mod, ONLY: set_derived_variables_cable_model_env
!
!IMPLICIT NONE
!
!! Subroutine arguments
!INTEGER, INTENT(IN) :: unitnumber
!
!CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_CABLE_MODEL_ENVIRONMENT'
!
!CALL read_nml_cable_model_env(unitnumber)
!
!CALL set_derived_variables_cable_model_env()
!
!RETURN
!END SUBROUTINE read_cable_model_environment

END MODULE read_namelists_mod_cbl
