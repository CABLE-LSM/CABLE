MODULE read_namelists_mod_cbl
  
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='READ_NAMELISTS_MOD_CBL'

CONTAINS

SUBROUTINE read_cable_namelists( unitnumber )

! Description:
!  Organize reading ALL cable namelists

USE cable_common_module, ONLY: cable_runtime
USE cable_namelist_util, ONLY: cable_namelist
    
IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: unitnumber

CHARACTER(LEN=*), PARAMETER :: RoutineName='read_cable_namelists'

CALL read_cable_namelist( unitnumber, cable_namelist )

IF( cable_runtime%offline ) THEN
  CALL read_offline_namelist( "offline.nml" )
ENDIF

RETURN
END SUBROUTINE read_cable_namelists



SUBROUTINE read_cable_namelist( unitnumber, cable_namelist )
! Description:
!  Read the cable namelist

!!USE cable_common_module_temp, ONLY: calcsoilalbedo, redistrb, wiltParam,       &
!!                                    satuParam, snmin, cable_user, gw_params
  
USE cable_common_module, ONLY: &
  calcsoilalbedo,              &
  redistrb,                    &
  wiltParam,                   &
  satuParam,                   &
  snmin,                       &
  cable_user,                  &
  gw_params
  
USE temp_module,   ONLY: l_casacnp,   &
                       l_landuse,     &
                       l_laiFeedbk,   &
                       l_vcmaxFeedbk
   
USE casadimension, ONLY: icycle
USE casavariable,  ONLY: casafile
   
IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: unitnumber
CHARACTER(LEN=*), INTENT(IN) :: cable_namelist

CHARACTER(LEN=*), PARAMETER :: RoutineName='read_cable_namelist'

NAMELIST /CABLE/  &
calcsoilalbedo, & ! albedo considers soil color Ticket #27
l_casacnp,      &
l_landuse,      &
l_laiFeedbk,    &
l_vcmaxFeedbk,  &
icycle,         &
redistrb,       &
wiltParam,      &
satuParam,      &
snmin,          &
cable_user,     & 
gw_params

! Open, read and close the namelist file.
OPEN( UNIT=unitnumber, FILE=CABLE_NAMELIST, STATUS="OLD", ACTION="READ" )
  READ( unitnumber, NML=CABLE )
CLOSE( unitnumber )

RETURN
END SUBROUTINE read_cable_namelist



SUBROUTINE read_offline_namelist( offline_namelist )
! Description:
!  Read the cable namelist for offline apps only

USE cable_common_module, ONLY: filename
USE casavariable,        ONLY: casafile
  
USE cable_IO_vars_module, ONLY: &
  soilparmnew,                  &
  fixedCO2,                     &
  output,                       &
  patchout,                     &
  check,                        &
  verbose,                      &
  leaps,                        &
  logn,                         &
  ncciy,                        &
  gswpfile,                     &
  globalMetfile
  
USE temp_module,   ONLY: vegparmnew,    &
                       spinup,        &
                       spincasa,      &
                       CASAONLY,      &
                       delsoilM,      &
                       delsoilT,      &
                       delgwM,        &
                       LALLOC
IMPLICIT NONE

! Subroutine arguments
CHARACTER(LEN=*), INTENT(IN) :: offline_namelist 

INTEGER                     :: unitnumber
CHARACTER(LEN=*), PARAMETER :: RoutineName='read_offline_cable_namelist'

NAMELIST /offline/  &
filename,       & ! TYPE, containing input filenames
vegparmnew,     & ! use new soil param. method
soilparmnew,    & ! use new soil param. method
spinup,         & ! spinup model (soil) to steady state
delsoilM,       &
delsoilT,       &
delgwM,         &
fixedCO2,       &
output,         &
patchout,       &
check,          &
verbose,        &
leaps,          &
logn,           &
spincasa,       &
CASAONLY,       &
casafile,       &
ncciy,          &
gswpfile,       &
globalMetfile

! Open, read and close the namelist file.
OPEN( NEWUNIT=unitnumber, FILE=offline_namelist, STATUS="OLD", ACTION="READ" )
  READ( unitnumber, NML=CABLE )
CLOSE( unitnumber )

RETURN
END SUBROUTINE read_offline_namelist



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
